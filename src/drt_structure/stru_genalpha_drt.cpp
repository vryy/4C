/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
  | structural nonlinear dynamics (gen-alpha)              m.gee 12/06  |
 *----------------------------------------------------------------------*/
void stru_genalpha_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // get a communicator and myrank
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank = Comm.MyPID();

  //----------------------------------------------------- get error file
  FILE* errfile = DRT::Problem::Instance()->ErrorFile()->Handle();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& solveparams = DRT::Problem::Instance()->StructSolverParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  STRUCT_DYN_CALC dynvar;
  memset(&dynvar, 0, sizeof(STRUCT_DYN_CALC));
  double          time = 0.0;  // we should add an input parameter
  double          timen;  // new time t_{n+1}

  //-----------------------------------------------------create a solver
  LINALG::Solver solver(solveparams,actdis->Comm(),errfile);
  actdis->ComputeNullSpaceIfNecessary(solver.Params());

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();

  // -------------------------------------------------------------------
  // create empty matrices
  // -------------------------------------------------------------------
  // `81' is an initial guess for the bandwidth of the matrices
  // A better guess will be determined later.
  RefCountPtr<LINALG::SparseMatrix> stiff_mat = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81));
  RefCountPtr<LINALG::SparseMatrix> mass_mat  = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81));
  RefCountPtr<LINALG::SparseMatrix> damp_mat  = null;
  bool damping = false;
  if (sdyn.get<int>("DAMPING")==1)
  {
    damping = true;
    damp_mat = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81));
  }

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
  // velocities V_{n} at last time
  RefCountPtr<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n} at last time
  RefCountPtr<Epetra_Vector> acc = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  RefCountPtr<Epetra_Vector> disn = LINALG::CreateVector(*dofrowmap,true);
  // velocities V_{n+1} at new time
  RefCountPtr<Epetra_Vector> veln = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n+1} at new time
  RefCountPtr<Epetra_Vector> accn = LINALG::CreateVector(*dofrowmap,true);

  // mid-displacements D_{n+1-alpha_f}
  RefCountPtr<Epetra_Vector> dism = LINALG::CreateVector(*dofrowmap,true);
  // mid-velocities V_{n+1-alpha_f}
  RefCountPtr<Epetra_Vector> velm = LINALG::CreateVector(*dofrowmap,true);
  // mid-accelerations A_{n+1-alpha_m}
  RefCountPtr<Epetra_Vector> accm = LINALG::CreateVector(*dofrowmap,true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  RefCountPtr<Epetra_Vector> disi = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  RefCountPtr<Epetra_Vector> fint = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_ext at last times
  RefCountPtr<Epetra_Vector> fext = LINALG::CreateVector(*dofrowmap,true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  RefCountPtr<Epetra_Vector> fextm = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_{n+1} at new time
  RefCountPtr<Epetra_Vector> fextn = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known at out-of-balance-force
  RefCountPtr<Epetra_Vector> fresm = LINALG::CreateVector(*dofrowmap,false);

  //-------------------------------------------- calculate external forces
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
    // other parameters needed by the elements
    params.set("total time",time);
    params.set("delta time",sdyn.get<double>("TIMESTEP"));
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("displacement",dis);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    actdis->EvaluateDirichlet(params,dis,null,null,dirichtoggle);
    actdis->ClearState();
    actdis->SetState("displacement",dis);
    // predicted rhs
    actdis->EvaluateNeumann(params,*fext);
    actdis->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle->PutScalar(1.0);
  invtoggle->Update(-1.0,*dirichtoggle,1.0);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_nlnstiffmass");
    // choose what to assemble
    params.set("assemble matrix 1",true);
    params.set("assemble matrix 2",true);
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",time);
    params.set("delta time",sdyn.get<double>("TIMESTEP"));
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",zeros);
    actdis->SetState("displacement",dis);
    //actdis->SetState("velocity",vel); // not used at the moment
    actdis->Evaluate(params,stiff_mat,mass_mat,fint,null,null);
    actdis->ClearState();
  }

  // complete stiffness and mass matrix
  stiff_mat->Complete();
  mass_mat->Complete();

  // build damping matrix if neccessary
  if (damping)
  {
    damp_mat->Add(*stiff_mat,false,sdyn.get<double>("K_DAMP"),0.0);
    damp_mat->Add(*mass_mat ,false,sdyn.get<double>("M_DAMP"),1.0);
    damp_mat->Complete();
  }

  //--------------------------- calculate consistent initial accelerations
  {
    RefCountPtr<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
    if (damping) damp_mat->Multiply(false,*vel,*rhs);
    rhs->Update(-1.0,*fint,1.0,*fext,-1.0);
    Epetra_Vector rhscopy(*rhs);
    rhs->Multiply(1.0,*invtoggle,rhscopy,0.0);
    solver.Solve(mass_mat->EpetraMatrix(),acc,rhs,true,true);
  }

  //------------------------------------------ time integration parameters
  const double beta   = sdyn.get<double>("BETA");
  const double gamma  = sdyn.get<double>("GAMMA");
  const double alpham = sdyn.get<double>("ALPHA_M");
  const double alphaf = sdyn.get<double>("ALPHA_F");
  const double dt     = sdyn.get<double>("TIMESTEP");
  int istep = 0;  // time step index

  //------------------------------------------------- output initial state
  output.NewStep(istep, time);
  output.WriteVector("displacement", dis);
  output.WriteVector("velocity", vel);
  output.WriteVector("acceleration", acc);

  //==================================================== start of timeloop
  while (time<=sdyn.get<double>("MAXTIME") && istep<=sdyn.get<int>("NUMSTEP"))
  {
    //------------------------------------------------------- current time
    // we are at t_{n} == time; the new time is t_{n+1} == time+dt
    timen = time + dt;

    //--------------------------------------------------- predicting state
    // constant predictor : displacement in domain
    disn->Update(1.0, *dis, 0.0);

    // apply new displacements at DBCs
    // and get new external force vector
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
      // other parameters needed by the elements
      params.set("total time",timen);
      params.set("delta time",dt);
      // set vector values needed by elements
      actdis->ClearState();
      actdis->SetState("displacement",disn);
      // predicted dirichlet values
      // disn then also holds prescribed new dirichlet displacements
      actdis->EvaluateDirichlet(params,disn,null,null,dirichtoggle);
      actdis->ClearState();
      actdis->SetState("displacement",disn);
      fextn->PutScalar(0.0);  // initialize external force vector (load vect)
      actdis->EvaluateNeumann(params,*fextn);
      actdis->ClearState();
    }

    // consistent predictor
    // predicting velocity V_{n+1} (veln)
    // V_{n+1} := gamma/(beta*dt) * (D_{n+1} - D_n)
    //          + (beta-gamma)/beta * V_n
    //          + (2.*beta-gamma)/(2.*beta) * A_n
    //veln->Update(1.0,*disn,-1.0,*dis,0.0);
    //veln->Update((beta-gamma)/beta,*vel,
    //            (2.*beta-gamma)/(2.*beta),*acc,
    //             gamma/(beta*dt));
    // predicting accelerations A_{n+1} (accn)
    // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
    //          - 1./(beta*dt) * V_n
    //          + (2.*beta-1.)/(2.*beta) * A_n
    //accn->Update(1.0,*disn,-1.0,*dis,0.0);
    //accn->Update(-1./(beta*dt),*vel,
    //            (2.*beta-1.)/(2.*beta),*acc,
    //            1./(beta*dt*dt));

    // constant predictor
    veln->Update(1.0,*vel,0.0);
    accn->Update(1.0,*acc,0.0);

    //------------------------------ compute interpolated dis, vel and acc
    // consistent predictor
    // mid-displacements D_{n+1-alpha_f} (dism)
    //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
    //dism->Update(1.-alphaf,*disn,alphaf,*dis,0.0);
    // mid-velocities V_{n+1-alpha_f} (velm)
    //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
    //velm->Update(1.-alphaf,*veln,alphaf,*vel,0.0);
    // mid-accelerations A_{n+1-alpha_m} (accm)
    //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
    //accm->Update(1.-alpham,*accn,alpham,*acc,0.0);

    // constant predictor
    dism->Update(1.0,*dis,0.0);
    velm->Update(1.0,*vel,0.0);
    accm->Update(1.0,*accm,0.0);

    //------------------------------- compute interpolated external forces
    // external mid-forces F_{ext;n+1-alpha_f} (fextm)
    //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
    //                         + alpha_f * F_{ext;n}
    fextm->Update(1.-alphaf,*fextn,alphaf,*fext,0.0);

    //------------- eval fint at interpolated state, eval stiffness matrix
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
      // other parameters that might be needed by the elements
      params.set("total time",timen);
      params.set("delta time",dt);
      // set vector values needed by elements
      actdis->ClearState();
      actdis->SetState("residual displacement",disi);
      actdis->SetState("displacement",dism);
      //actdis->SetState("velocity",velm); // not used at the moment
      fint->PutScalar(0.0);  // initialise internal force vector
      actdis->Evaluate(params,stiff_mat,null,fint,null,null);
      actdis->ClearState();
      // do NOT finalize the stiffness matrix to mass and damping to it later
    }

    //-------------------------------------------- compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add mid-inertial force
    mass_mat->Multiply(false,*accm,*fresm);
    // add mid-viscous damping force
    if (damping)
    {
    	RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    	damp_mat->Multiply(false,*velm,*fviscm);
    	fresm->Update(1.0,*fviscm,1.0);
    }

    // add static mid-balance
    fresm->Update(1.0,*fint,-1.0,*fextm,1.0);

    // blank residual at DOFs on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm);
      fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
    }

    //------------------------------------------------ build residual norm
    double disinorm = 1.0e6;
    double fresmnorm = 1.0e6;
    fresm->Norm2(&fresmnorm);
    if (!myrank) cout << "Predictor residual forces " << fresmnorm << endl; fflush(stdout);

    //=================================================== equilibrium loop
    int numiter=0;
    while (disinorm>sdyn.get<double>("TOLDISP") && numiter<=sdyn.get<int>("MAXITER"))
    {
      //------------------------------------------- effective rhs is fresm
      //---------------------------------------------- build effective lhs
      // (using matrix stiff_mat as effective matrix)
      stiff_mat->Add(*mass_mat,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
      if (damping)
        stiff_mat->Add(*damp_mat,false,(1.-alphaf)*gamma/(beta*dt),1.0);
      stiff_mat->Complete();

      //----------------------- apply dirichlet BCs to system of equations
      fresm->Scale(-1.0);  // delete this by building fresm with other sign
      disi->PutScalar(0.0);  // Useful? depends on solver and more
      LINALG::ApplyDirichlettoSystem(stiff_mat,disi,fresm,zeros,dirichtoggle);

      //--------------------------------------------------- solve for disi
      // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
      if (numiter==0)
      {
        solver.Solve(stiff_mat->EpetraMatrix(),disi,fresm,true,true);
      }
      else
      {
        solver.Solve(stiff_mat->EpetraMatrix(),disi,fresm,true,false);
      }

      //---------------------------------- update mid configuration values
      // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
      dism->Update(1.-alphaf, *disi, 1.0);
      // V_{n+1-alpha_f} := V_{n+1-alpha_f}
      //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
      velm->Update((1.-alphaf)*gamma/(beta*dt), *disi, 1.0);
      // A_{n+1-alpha_m} := A_{n+1-alpha_m}
      //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
      accm->Update((1.-alpham)/(beta*dt*dt), *disi, 1.0);

      //---------------------------- compute internal forces and stiffness
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
        // other parameters that might be needed by the elements
        params.set("total time",timen);
        params.set("delta time",dt);
        // set vector values needed by elements
        actdis->ClearState();
        actdis->SetState("residual displacement",disi);
        actdis->SetState("displacement",dism);
        //actdis->SetState("velocity",velm); // not used at the moment
        fint->PutScalar(0.0);  // initialise internal force vector
        actdis->Evaluate(params,stiff_mat,null,fint,null,null);
        actdis->ClearState();
        // do NOT finalize the stiffness matrix to add masses to it later
      }

      //------------------------------------------ compute residual forces
      // Res = M . A_{n+1-alpha_m}
      //     + C . V_{n+1-alpha_f}
      //     + F_int(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}
      // add inertia mid-forces
      mass_mat->Multiply(false,*accm,*fresm);
      // add viscous mid-forces
      if (damping)
      {
        RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
        damp_mat->Multiply(false,*velm,*fviscm);
        fresm->Update(1.0,*fviscm,1.0);
      }
      // add static mid-balance
      fresm->Update(1.0,*fint,-1.0,*fextm,1.0);
      // blank residual DOFs with are on Dirichlet BC
      {
        Epetra_Vector fresmcopy(*fresm);
        fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
      }

      //---------------------------------------------- build residual norm
      disi->Norm2(&disinorm);

      fresm->Norm2(&fresmnorm);
      // a short message
      if (!myrank)
      {
        printf("numiter %d res-norm %e dis-norm %e\n",numiter+1, fresmnorm, disinorm);
        fprintf(errfile,"numiter %d res-norm %e dis-norm %e\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
        fflush(errfile);
      }

      //--------------------------------- increment equilibrium loop index
      ++numiter;
    } //============================================= end equilibrium loop

    //-------------------------------- test whether max iterations was hit
    if (numiter==sdyn.get<int>("MAXITER")) dserror("Newton unconverged in %d iterations",numiter);

    //---------------------------- determine new end-quantities and update
    // new displacements at t_{n+1} -> t_n
    //    D_{n} := D_{n+1} = 1./(1.-alphaf) * D_{n+1-alpha_f}
    //                     - alphaf/(1.-alphaf) * D_n
    dis->Update(1./(1.-alphaf), *dism, -alphaf/(1.-alphaf));
    // new velocities at t_{n+1} -> t_n
    //    V_{n} := V_{n+1} = 1./(1.-alphaf) * V_{n+1-alpha_f}
    //                     - alphaf/(1.-alphaf) * V_n
    vel->Update(1./(1.-alphaf), *velm, -alphaf/(1.-alphaf));
    // new accelerations at t_{n+1} -> t_n
    //    A_{n} := A_{n+1} = 1./(1.-alpham) * A_{n+1-alpha_m}
    //                     - alpham/(1.-alpham) * A_n
    acc->Update(1./(1.-alpham), *accm, -alpham/(1.-alpham));
    // update new external force
    //    F_{ext;n} := F_{ext;n+1}
    fext->Update(1.0, *fextn, 0.0);

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

    //------------------------------------------------ increment time step
    ++istep;  // step n := n + 1
    time += dt;  // time t_n := t_{n+1} = t_n + Delta t

    //----------------------------------------------------- output results
    int mod_disp   = istep % sdyn.get<int>("RESEVRYDISP");
    if (!mod_disp && ioflags.struct_disp==1)
    {
      output.NewStep(istep, time);
      output.WriteVector("displacement", dis);
      output.WriteVector("velocity", vel);
      output.WriteVector("acceleration", acc);
    }

    //---------------------------------------------- do stress calculation
    int mod_stress = istep % sdyn.get<int>("RESEVRYSTRS");
    if (!mod_stress && ioflags.struct_stress==1)
    {
      // create the parameters for the discretization
      ParameterList params;
      // action for elements
      params.set("action","calc_struct_stress");
      // choose what to assemble
      params.set("assemble matrix 1",false);
      params.set("assemble matrix 2",false);
      params.set("assemble vector 1",false);
      params.set("assemble vector 2",false);
      params.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      params.set("total time",timen);
      params.set("delta time",dt);
      // set vector values needed by elements
      actdis->ClearState();
      actdis->SetState("residual displacement",zeros);
      actdis->SetState("displacement",dis);
      actdis->Evaluate(params,null,null,null,null,null);
      actdis->ClearState();
    }

    //---------------------------------------------------------- print out
    if (!myrank)
    {
      printf("step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
             istep,sdyn.get<int>("NUMSTEP"),timen,dt,numiter);
      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
              istep,sdyn.get<int>("NUMSTEP"),timen,dt,numiter);
      printf("----------------------------------------------------------------------------------\n");
      fprintf(errfile,"----------------------------------------------------------------------------------\n");
      fflush(stdout);
      fflush(errfile);
    }

  }  //===================================================== end time loop



  //----------------------------- this is the end my lonely friend the end
  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef CCADISCRET
