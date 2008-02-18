/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_static_drt.H"
#include "../io/io_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "stru_resulttest.H"


/*----------------------------------------------------------------------*
 |                                                         maf 05/07    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
 |                                                         maf 05/07    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                           maf 05/07
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                         maf 05/07    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                         maf 05/07    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*
 |                                                         maf 05/07    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;

/*----------------------------------------------------------------------*
  | structural nonlinear static solution routine             maf 05/07  |
 *----------------------------------------------------------------------*/
void stru_static_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // get a communicator and myrank
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank = Comm.MyPID();

  //----------------------------------------------------- get error file
  FILE* errfile = allfiles.out_err;

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[0];

  //-----------------------------------------------------create a solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();

  // -------------------------------------------------------------------
  // create empty stiffness matrix
  // -------------------------------------------------------------------
  // `81' is an initial guess for the bandwidth of the matrices
  // A better guess will be determined later.
  RefCountPtr<LINALG::SparseMatrix> stiff_mat = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81));

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // a zero vector of full length
  RefCountPtr<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  RefCountPtr<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  RefCountPtr<Epetra_Vector> invtoggle = LINALG::CreateVector(*dofrowmap,false);
  // displacements D_{n} at last time
  RefCountPtr<Epetra_Vector> dis = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  RefCountPtr<Epetra_Vector> disn = LINALG::CreateVector(*dofrowmap,true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  RefCountPtr<Epetra_Vector> disi = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  RefCountPtr<Epetra_Vector> fint = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_ext at last times
  RefCountPtr<Epetra_Vector> fext = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_{n+1} at new time
  RefCountPtr<Epetra_Vector> fextn = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known at out-of-balance-force
  RefCountPtr<Epetra_Vector> fresm = LINALG::CreateVector(*dofrowmap,false);

  // nodal normal stresses
  RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*dofrowmap,true);
  // nodal shear stresses
  RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*dofrowmap,true);

  if (statvar->nr_controltyp != control_load) dserror("Only load control implemented");

  /*
  ** solution control parameters are inherited from dynamic routine:
  ** dt     = stepsize
  ** istep  = load step index
  ** time   = redundant, equals istep*dt
  */
  //------------------------------------------ time integration parameters
  const double dt = statvar->stepsize;
  int istep = 0;
  double time = 0.0;  // we should add an input parameter
  double timen;

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(actdis);
  if (genprob.restart){
    int restartstep = genprob.restart;
    RefCountPtr<DRT::Discretization> rcpdiscret(actdis);
    rcpdiscret.release();
    IO::DiscretizationReader reader(rcpdiscret,restartstep);
    double rtime  = reader.ReadDouble("time");
    int    rstep = reader.ReadInt("step");
    if (rstep != restartstep) dserror("Time step on file not equal to given step");

    reader.ReadVector(dis, "displacement");
    //reader.ReadVector(fext, "fexternal");
    //reader.ReadMesh(restartstep);

    // override current time and step with values from file
    time = rtime;
    istep = rstep;
  }

  //---------------------------------------------- do "stress" calculation
  int mod_stress = istep % statvar->resevry_stress;
  if (!mod_stress && ioflags.struct_stress==1)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",zeros);
    actdis->SetState("displacement",dis);
    actdis->Evaluate(p,null,null,null,null,null);
    actdis->ClearState();
  }


  // write mesh always at beginning of calc or restart
  output.WriteMesh(istep,time);


  //-------------------------------------- parameters for volume constraint
  // required volumes (can change over time with a given loadcurve
  //vector<double> referencevolumes;
  Epetra_SerialDenseVector referencevolumes;
  Epetra_SerialDenseVector startvolumes;
  Epetra_SerialDenseVector actvol;
  //vector<double> actvol;
  Epetra_SerialDenseVector volerr;
  //min, max Condition IDs, and total number if volconstrID
  int minConstrID;
  int maxConstrID;
  int numConstrID;
  //"Loadfactor" for restricted volume
  double fact;
  //lagrange multiplier and increments (-> Uzawa)
  RefCountPtr<Epetra_SerialDenseVector> lagrMultVec;
  RefCountPtr<Epetra_SerialDenseVector> lagrMultInc;
  //every volume gets is on dofrowmap
  map<int, RCP<Epetra_Map> > voldofrowmaps;
  //vector contains ALL first "contrain derivative", i.e. additional column in stiff_mat
  //We will deal with non overlapping volumes, so "constrain derivative" for a specific volume
  //can be determined using voldofrowmaps
  RefCountPtr<Epetra_Vector> constrVec = LINALG::CreateVector(*dofrowmap,true);


  //-------------------------------- calculate external force distribution
  //---- which is scaled by the load factor lambda (itself stays constant)
  {
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_eleload");
    // choose what to assemble
    params.set("assemble matrix 1",false);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);

    //other parameters needed by the elements
    params.set("total time",time);
    params.set("delta time",dt);

    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("displacement",dis);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    actdis->EvaluateDirichlet(params,*dis,*dirichtoggle);
    actdis->ClearState();
    actdis->SetState("displacement",dis);
    // predicted rhs
    actdis->EvaluateNeumann(params,*fext); // *fext holds external force vector
    params.set("action","calc_struct_constrvol");
    actdis->EvaluateCondition(params,"VolumeConstraint_3D");
    synchronize_VolConstraint(Comm, params,startvolumes,minConstrID,maxConstrID,numConstrID);
    //inititalize lambda vector and increments

    if (numConstrID!=0)
    {
    	if (!myrank)
    		cout<<"Reference Volume: "<<startvolumes[0]<<endl;
    	lagrMultVec=rcp(new Epetra_SerialDenseVector(numConstrID));
    	lagrMultInc=rcp(new Epetra_SerialDenseVector(numConstrID));
    	lagrMultVec->Scale(0.0);
    	lagrMultInc->Scale(0.0);
    }

    setupVolDofrowmaps(actdis,voldofrowmaps,numConstrID,"VolumeConstraint_3D");
    actdis->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle->PutScalar(1.0);
  invtoggle->Update(-1.0,*dirichtoggle,1.0);

  //------------------------------------------------- output initial state
  output.NewStep(istep, time);
  output.WriteVector("displacement", dis);
  //---------------------------------------------- output stresses
  if (!mod_stress && ioflags.struct_stress==1)
  {
    output.WriteElementData();

    if (0)
    {
      // create the parameters for the discretization
      ParameterList p;
      p.set("action","calc_struct_stress_nodal");
      RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*(actdis->DofRowMap()),true);
      RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*(actdis->DofRowMap()),true);
      actdis->Evaluate(p,null,null,normal_stresses,shear_stresses,null);
      actdis->ClearState();
      output.WriteStressVector("nodal_stresses_xyz", normal_stresses, shear_stresses);
    }
  }
  //---------------------------------------------end of output initial state

  //========================================== start of time/loadstep loop
  while ( istep < statvar->nstep)
  {
    //------------------------------------------------------- current time
    // we are at t_{n} == time; the new time is t_{n+1} == time+dt
    timen = time + dt;

    //--------------------------------------------------- predicting state
    // constant predictor : displacement in domain
    disn->Update(1.0, *dis, 0.0);

    // eval fint and stiffness matrix at current istep
    // and apply new displacements at DBCs
    {
      // destroy and create new matrix
      stiff_mat->Zero();
      // create the parameters for the discretization
      ParameterList params;
      // action for elements
      params.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      params.set("assemble matrix 1",true);
      params.set("assemble matrix 2",false);
      params.set("assemble vector 1",true);
      params.set("assemble vector 2",false);
      params.set("assemble vector 3",false);
      // other parameters needed by the elements
      params.set("total time",timen);  // load factor (pseudo time)
      params.set("delta time",dt);  // load factor increment (pseudo time increment)
      // set vector values needed by elements
      actdis->ClearState();
      actdis->SetState("residual displacement",disi);
      // predicted dirichlet values
      // disn then also holds prescribed new dirichlet displacements
      actdis->EvaluateDirichlet(params,*disn,*dirichtoggle);
      actdis->SetState("displacement",disn);
      fint->PutScalar(0.0);  // initialise internal force vector
      actdis->Evaluate(params,stiff_mat,null,fint,null,null);
      // predicted rhs
      fextn->PutScalar(0.0);  // initialize external force vector (load vect)
      actdis->EvaluateNeumann(params,*fextn); // *fext holds external force vector at current step
      params.set("action","calc_struct_volconstrstiff");
      params.set("assemble vector 2",true);
      params.set("MinID",minConstrID);
      params.set("MaxID",maxConstrID);
      params.set("NumberofID",numConstrID);
      params.set("LagrMultVector",lagrMultVec);
      constrVec->PutScalar(0.0);
      actdis->EvaluateCondition(params,stiff_mat,fint,constrVec,"VolumeConstraint_3D");
      synchronize_VolConstraint(Comm, params,actvol,minConstrID,maxConstrID,numConstrID);
      volerr.Size(numConstrID);
      fact=params.get("LoadCurveFactor",1.0);
      referencevolumes=startvolumes;
      if (numConstrID!=0)
      {
    	  	referencevolumes.Scale(fact);
       		if (!myrank)
       		{
       			cout<<"New Reference Volume: "<<referencevolumes[0]<<endl;
       		}
      }
      for (int iter = 0; iter < numConstrID; ++iter)
      {
    	  volerr[iter]=referencevolumes[iter]-actvol[iter];
      }
      actdis->ClearState();
      }
    // complete stiffness matrix
    stiff_mat->Complete();

    double stiffnorm;
    stiffnorm = stiff_mat->NormFrobenius();

    // evaluate residual at current istep
    // R{istep,numiter=0} = F_int{istep,numiter=0} - F_ext{istep}
    fresm->Update(1.0,*fint,-1.0,*fextn,0.0);

    // blank residual at DOFs on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm);
      fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
    }

    //------------------------------------------------ build residual norm
    double norm;
    fresm->Norm2(&norm);
    if (!myrank) cout << " Predictor residual forces " << norm << endl; fflush(stdout);

    double norm_vol;
    if (numConstrID!=0)
    {
    	norm_vol=volerr.Norm2();
    }
    else
    {
    	norm_vol=0.0;
    }
    //=================================================== equilibrium loop
    int numiter=0;
    while ((norm > statvar->toldisp||norm_vol>statvar->toldisp) && numiter < statvar->maxiter)
    {
      //----------------------- apply dirichlet BCs to system of equations
      fresm->Scale(-1.0);     // rhs = -R = -fresm
      disi->PutScalar(0.0);   // Useful? depends on solver and more
      LINALG::ApplyDirichlettoSystem(stiff_mat,disi,fresm,zeros,dirichtoggle);

      if (numConstrID==0)
      {
	    	//Do usual newton step
	    	// solve for disi
	      	// Solve K . IncD = -R  ===>  IncD_{n+1}
	      	if (numiter==0)
	      	{
                  solver.Solve(stiff_mat->Matrix(),disi,fresm,true,true);
	      	}
	      	else
	      	{
                  solver.Solve(stiff_mat->Matrix(),disi,fresm,true,false);
	      	}
      }
      else
      {
    	  	//===================================================uzawa loop
	        // For every iteration step an uzawa algorithm is used to solve the linear system.
	        //Preparation of uzawa method to solve the linear system.
	        double norm_uzawa;
	        double norm_vol_uzawa;
	        int numiter_uzawa=0;
	        lagrMultInc->Scale(0.0);
	        Epetra_Vector constrVecWeight(*constrVec);
	        Epetra_SerialDenseVector dotprod(numConstrID);

	        // Compute residual of the uzawa algorithm
	        for (int foo = 0; foo < numConstrID; ++foo)
      	  	{
	        	Epetra_Vector onlyvol_foo(*(voldofrowmaps[foo]));
      		  	Epetra_Vector onlydis_foo(*(voldofrowmaps[foo]));
      		  	Epetra_Import importer1(*(voldofrowmaps[foo]),constrVec->Map());
      		  	Epetra_Import importer3(*(voldofrowmaps[foo]),disi->Map());
      		  	onlyvol_foo.Import(*constrVec,importer1,Insert);
      		  	onlydis_foo.Import(*disi,importer3,Insert);
      		  	onlydis_foo.Dot(onlyvol_foo,&(dotprod[foo]));
      		  	onlyvol_foo.Scale(-(*lagrMultInc)[foo]);
      		  	Epetra_Import importer2(constrVec->Map(),*(voldofrowmaps[foo]));
      		  	constrVecWeight.Import(onlyvol_foo,importer2,Insert);
      	  	}

	        RCP<Epetra_Vector> fresmcopy=rcp(new Epetra_Vector(*fresm));
	      	fresmcopy->Update(1.0,constrVecWeight,1.0);
	      	Epetra_Vector uzawa_res(*fresmcopy);
	      	stiff_mat->Multiply(false,*disi,uzawa_res);
	      	uzawa_res.Update(1.0,*fresmcopy,-1.0);
	      	// blank residual DOFs which are on Dirichlet BC
	      	{
	      	    Epetra_Vector rescopy(uzawa_res);
	      	    uzawa_res.Multiply(1.0,*invtoggle,rescopy,0.0);
	      	}
	      	uzawa_res.Norm2(&norm_uzawa);

	      	Epetra_SerialDenseVector vol_res(numConstrID);
	      	for (int foo = 0; foo < numConstrID; ++foo)
	      	{
	      	  vol_res[foo]=dotprod[foo]+volerr[foo];
	      	}
	      	norm_vol_uzawa=vol_res.Norm2();

	        //Solve one iteration step with augmented lagrange
	      	//Since we calculate displacement norm as well, at least one step has to be taken
	        while (((norm_uzawa > statvar->toldisp||norm_vol_uzawa>statvar->toldisp)
	        		&& numiter_uzawa < statvar->maxiter)||numiter_uzawa<1)
	        {

	      	  LINALG::ApplyDirichlettoSystem(stiff_mat,disi,fresmcopy,zeros,dirichtoggle);
	      	  // solve for disi
	      	  // Solve K . IncD = -R  ===>  IncD_{n+1}
	      	  if (numiter_uzawa==0&&numiter==0)
	      	  {
                    solver.Solve(stiff_mat->Matrix(),disi,fresmcopy,true,true);
	      	  }
	      	  else
	      	  {
                    solver.Solve(stiff_mat->Matrix(),disi,fresmcopy,true,false);
	      	  }

	      	  //compute lagrange multiplier increments
	      	  const double alpha=10;


	      	  for (int foo = 0; foo < numConstrID; ++foo)
	      	  {
	      		  Epetra_Vector onlyvol_foo(*(voldofrowmaps[foo]));
	      		  Epetra_Vector onlydis_foo(*(voldofrowmaps[foo]));
	      		  Epetra_Import importer1(*(voldofrowmaps[foo]),constrVec->Map());
	      		  onlyvol_foo.Import(*constrVec,importer1,Insert);
	      		  onlydis_foo.Import(*disi,importer1,Insert);
	      		  onlydis_foo.Dot(onlyvol_foo,&(dotprod[foo]));
	      		  (*lagrMultInc)[foo]+=alpha*(dotprod[foo]+volerr[foo]);
	      	  }
	      	  //Compute residual of the uzawa algorithm
	      	  constrVecWeight.PutScalar(0.0);
	      	  for (int foo = 0; foo < numConstrID; ++foo)
	      	  {
	      		  Epetra_Vector onlyvol_foo(*(voldofrowmaps[foo]));
	      		  Epetra_Vector onlydis_foo(*(voldofrowmaps[foo]));
	  		      Epetra_Import importer1(*(voldofrowmaps[foo]),constrVec->Map());
	  		      Epetra_Import importer3(*(voldofrowmaps[foo]),disi->Map());
	  		      onlyvol_foo.Import(*constrVec,importer1,Insert);
	  		      onlyvol_foo.Scale(-(*lagrMultInc)[foo]);
	  		      onlydis_foo.Import(*disi,importer3,Insert);
	  		      Epetra_Import importer2(constrVec->Map(),*(voldofrowmaps[foo]));
	  		      constrVecWeight.Import(onlyvol_foo,importer2,Insert);
	      	  }
	      	  fresmcopy->Update(1.0,constrVecWeight,1.0,*fresm,0.0);
	      	  Epetra_Vector uzawa_res(*fresmcopy);
	      	  stiff_mat->Multiply(false,*disi,uzawa_res);
	      	  uzawa_res.Update(1.0,*fresmcopy,-1.0);
	      	  // blank residual DOFs which are on Dirichlet BC
	      	  {
	      	      Epetra_Vector rescopy(uzawa_res);
	      		  uzawa_res.Multiply(1.0,*invtoggle,rescopy,0.0);
	      	  }
	      	  uzawa_res.Norm2(&norm_uzawa);
	      	  Epetra_SerialDenseVector vol_res(numConstrID);
	      	  for (int foo = 0; foo < numConstrID; ++foo)
	      	  {
	      		  vol_res[foo]=dotprod[foo]+volerr[foo];
	      	  }
	      	  norm_vol_uzawa=vol_res.Norm2();
	      	  numiter_uzawa++;
	        }//Uzawa loop

	        if (numiter_uzawa==statvar->maxiter)
	        {
	              dserror("Uzawa unconverged in %d iterations",numiter_uzawa);
	        }

	        //update lagrange multiplier
	        for (int foo = 0; foo < numConstrID; ++foo)
	        {
	        	  (*lagrMultVec)[foo]=(*lagrMultVec)[foo]+(*lagrMultInc)[foo];
	        }
      }

      // update displacements
      // D_{istep,numiter+1} := D_{istep,numiter} + IncD_{numiter}
      disn->Update(1.0, *disi, 1.0);

      // compute internal forces and stiffness at current iterate numiter
      {
        // zero out stiffness
        stiff_mat->Zero();
        // create the parameters for the discretization
        ParameterList params;
        // action for elements
        params.set("action","calc_struct_nlnstiff");
        // choose what to assemble
        params.set("assemble matrix 1",true);
        params.set("assemble matrix 2",false);
        params.set("assemble vector 1",true);
        params.set("assemble vector 2",false);
        params.set("assemble vector 3",false);
        // other parameters needed by the elements
        params.set("total time",timen);  // load factor (pseudo time)
        params.set("delta time",dt);  // load factor increment (pseudo time increment)
        // set vector values needed by elements
        actdis->ClearState();
        actdis->SetState("residual displacement",disi);
        actdis->SetState("displacement",disn);
        fint->PutScalar(0.0);  // initialise internal force vector
        actdis->Evaluate(params,stiff_mat,null,fint,null,null);

        params.set("action","calc_struct_volconstrstiff");
        params.set("assemble vector 2",true);
        params.set("MinID",minConstrID);
        params.set("MaxID",maxConstrID);
        params.set("NumberofID",numConstrID);
        params.set("LagrMultVector",lagrMultVec);
        constrVec->PutScalar(0.0);
        actdis->EvaluateCondition(params,stiff_mat,fint,constrVec,"VolumeConstraint_3D");
        synchronize_VolConstraint(Comm, params,actvol,minConstrID,maxConstrID,numConstrID);
        if (!myrank && numConstrID!=0)
        {
        	cout<<"Computed Volume after Newton Step: "<<actvol[0]<<endl;
        }

        for (int iter = 0; iter < numConstrID; ++iter)
        {
        	volerr[iter]=referencevolumes[iter]-actvol[iter];
        }

        actdis->ClearState();
      }
      // complete stiffness matrix
      stiff_mat->Complete();

      // evaluate new residual fresm at current iterate numiter
      // R{istep,numiter} = F_int{istep,numiter} - F_ext{istep}
      fresm->Update(1.0,*fint,-1.0,*fextn,0.0);

      // blank residual DOFs which are on Dirichlet BC
      {
        Epetra_Vector fresmcopy(*fresm);
        fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
      }

      //---------------------------------------------- build residual norm
      double disinorm;
      disi->Norm2(&disinorm);

      norm_vol=volerr.Norm2();

      fresm->Norm2(&norm);
      // a short message
      if (!myrank)
      {
        printf("numiter %d res-norm %e dis-norm %e \n",numiter+1, norm, disinorm);
        fprintf(errfile,"numiter %d res-norm %e dis-norm %e\n",numiter+1, norm, disinorm);
        fflush(stdout);
        fflush(errfile);
      }
      // criteria to stop Newton iteration
      norm = disinorm;

      for (int iter = 0; iter < numConstrID; ++iter)
      {
           	(*lagrMultInc)[iter]=0.0;
      }

      //--------------------------------- increment equilibrium loop index
      ++numiter;
    } //
    //============================================= end equilibrium loop

    //-------------------------------- test whether max iterations was hit
    if (statvar->maxiter == 1 && statvar->nstep == 1)
      printf("computed 1 step with 1 iteration: STATIC LINEAR SOLUTION\n");
    else if (numiter==statvar->maxiter)
      dserror("Newton unconverged in %d iterations",numiter);

    //---------------------------- determine new end-quantities and update
    // new displacements at t_{n+1} -> t_n
    // D_{n} := D_{n+1}
    dis->Update(1.0, *disn, 0.0);

    //----- update anything that needs to be updated at the element level
    {
      // create the parameters for the discretization
      ParameterList params;
      // action for elements
      params.set("action","calc_struct_update_istep");
      // choose what to assemble
      params.set("assemble matrix 1",false);
      params.set("assemble matrix 2",false);
      params.set("assemble vector 1",false);
      params.set("assemble vector 2",false);
      params.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      params.set("total time",timen);
      params.set("delta time",dt);
      actdis->Evaluate(params,null,null,null,null,null);
    }

    //---------------------------------------------- do stress calculation
    int mod_stress = istep % statvar->resevry_stress;
    if (!mod_stress && ioflags.struct_stress==1)
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_stress");
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      actdis->ClearState();
      actdis->SetState("residual displacement",zeros);
      actdis->SetState("displacement",dis);
      actdis->Evaluate(p,null,null,null,null,null);
      actdis->ClearState();
    }


    //------------------------------------------ increment time/load step
    ++istep;      // load step n := n + 1
    time += dt;   // load factor / pseudo time  t_n := t_{n+1} = t_n + Delta t

    //------------------------------------------------- write restart step
    bool isdatawritten = false;
    if (istep % statvar->resevery_restart==0)
    {
      output.WriteMesh(istep,time);
      output.NewStep(istep, time);
      output.WriteVector("displacement",dis);
      //output.WriteVector("fexternal", fext);
      isdatawritten = true;

      if (!myrank)
      {
        cout << "====== Restart written in step " << istep << endl;
        fflush(stdout);
        fprintf(errfile,"====== Restart written in step %d\n",istep);
        fflush(errfile);
      }
    }

    //----------------------------------------------------- output results
    int mod_disp   = istep % statvar->resevry_disp;
    if (!mod_disp && ioflags.struct_disp==1 && !isdatawritten)
    {
      output.NewStep(istep, time);
      output.WriteVector("displacement", dis);
      isdatawritten = true;
    }

    //---------------------------------------------------- output stresses
    if (!mod_stress && ioflags.struct_stress==1)
    {
      if (!isdatawritten) output.NewStep(istep, timen);
      isdatawritten = true;
      output.WriteElementData();

      if (0)
      {
        // create the parameters for the discretization
        ParameterList p;
        // action for elements
        p.set("action","calc_struct_stress_nodal");
        RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*(actdis->DofRowMap()),true);
        RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*(actdis->DofRowMap()),true);
        actdis->Evaluate(p,null,null,normal_stresses,shear_stresses,null);
        actdis->ClearState();
        if (!isdatawritten) output.NewStep(istep, timen);
        output.WriteStressVector("nodal_stresses_xyz", normal_stresses, shear_stresses);
        isdatawritten = true;
      }
    }

//
//
//    ofstream f_system("stresses.gmsh");
//    stringstream gmshfilecontent;
//    gmshfilecontent << "View \" Solid Elements stresses \" {" << endl;
//    // plot elements
//    for (int i=0; i < actdis->NumMyColElements(); ++i)
//    {
//      if (actdis->lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
//      DRT::ELEMENTS::So_sh8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8*>(actdis->lColElement(i));
//      if (!actele) dserror("cast to So_sh8* failed");
//      // plot elements
//      gmshfilecontent << IO::GMSH::elementToString(actele->thickdir_, actele) << endl;
//      Epetra_SerialDenseMatrix* mystresses(8,8);
//      mystresses = actele->data_.Get<Epetra_SerialDenseMatrix>("Stresses");
//      cout << *mystresses;
//    }
//    gmshfilecontent << "};" << endl;

    //---------------------------------------------------------- print out
    if (!myrank)
    {
      printf("step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
             istep,statvar->nstep,timen,dt,numiter);
      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
              istep,statvar->nstep,timen,dt,numiter);
      printf("----------------------------------------------------------------------------------\n");
      fprintf(errfile,"----------------------------------------------------------------------------------\n");
      fflush(stdout);
      fflush(errfile);
    }

  }  //=============================================end time/loadstep loop

  // Structure Resulttests
  DRT::ResultTestManager testmanager(actdis->Comm());
  testmanager.AddFieldTest(rcp(new StruResultTest(actdis,dis)));
  testmanager.TestAll();

  //----------------------------- this is the end my lonely friend the end
  return;
} // end of stru_static_drt()

/*----------------------------------------------------------------------*
 |                                                          tk 10/07    |
 |small subroutine to synchronize processors after evaluating the      |
 |volume constraint volume                                              |
 *----------------------------------------------------------------------*/
void synchronize_VolConstraint(const Epetra_Comm& Comm,
								ParameterList& params,
								Epetra_SerialDenseVector& vect,
								int& minID,
								int& maxID,
								int& numID)
{
	Comm.MaxAll(&(params.get("MaxID",0)),&maxID,1);
	Comm.MinAll(&(params.get("MinID",100000)),&minID,1);
	if (minID==100000)
	{
		minID=1;
		numID=0;
	}
	else
	{
		numID=1+maxID-minID;
		vect.Size(1+maxID-minID);
		for (int i = 0; i <= maxID-minID; ++i)
		{
			char volname[30];
			sprintf(volname,"computed volume %d",i+1);
			double currvol;
			Comm.SumAll(&(params.get(volname,0.0)),&(currvol),1);
			vect[i]=currvol;
		}
	}
	return;
}

/*----------------------------------------------------------------------*
 |                                                          tk 10/07    |
 |subroutine to setup separate dofrowmap for any boundary condition ID  |
 *----------------------------------------------------------------------*/
void setupVolDofrowmaps(RefCountPtr<DRT::Discretization> actdis,
		  map<int, RCP<Epetra_Map> >& voldofrowmaps,
		  const int numconstID,
		  const string& condstring)
{
    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    map<int, vector<int> > dofdata;

    for (int iter = 0; iter < numconstID; ++iter)
    {
		(dofdata[iter]).reserve(actdis->NumMyRowNodes()*3);
	}

    //loop over all my nodes
    for (int i=0; i<actdis->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = actdis->lRowNode(i);
    //loop through the conditions and build dofdata vectors
      DRT::Condition* cond = node->GetCondition(condstring);
      if (cond)
      {
    	  const vector<int>*    CondIDVec  = cond->Get<vector<int> >("ConditionID");
  		  int CondID=(*CondIDVec)[0];

  		  vector<int> dof = actdis->Dof(node);
  		  int numdofs = dof.size();
  		  for (int j=0; j<numdofs; ++j)
  		  {
  			  // add this velocity dof to the velmapdata vector
  		  	  (dofdata[CondID-1]).push_back(dof[j]);
  		  }//dof loop

      }

    }
    for (int voliter = 0; voliter < numconstID; ++voliter)
    {
    	voldofrowmaps[voliter]= rcp(new Epetra_Map(-1,
                                    dofdata[voliter].size(),&dofdata[voliter][0],0,
                                    actdis->Comm()));
    }
	return;
}
#endif  // #ifdef CCADISCRET
