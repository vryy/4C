/*!----------------------------------------------------------------------
\file microstatic.cpp
\brief Static control for  microstructural problems in case of multiscale
analyses

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Epetra_LinearProblem.h>
#include <Amesos_Klu.h>

#include "microstatic.H"

#include <vector>

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../io/io_drt_micro.H"

#include "../drt_mfsi/mfsi_debug.H"

using namespace IO;

/*----------------------------------------------------------------------*
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |  ctor (public)|
 *----------------------------------------------------------------------*/
MicroStatic::MicroStatic(RCP<ParameterList> params,
                         RCP<DRT::Discretization> dis,
                         RCP<LINALG::Solver> solver):
params_(params),
discret_(dis),
solver_(solver)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time    = params_->get<double>("total time"      ,0.0);
  double dt      = params_->get<double>("delta time"      ,0.01);
//   int istep      = params_->get<int>   ("step"            ,0);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  if (!discret_->Filled()) discret_->FillComplete();
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  myrank_ = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // create empty matrices
  // -------------------------------------------------------------------
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,true,false));

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap,true);
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  invtoggle_ = LINALG::CreateVector(*dofrowmap,false);

  // displacements D_{n} at last time
  dis_ = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  disn_ = LINALG::CreateVector(*dofrowmap,true);

  // mid-displacements D_{n+1-alpha_f}
  dism_ = LINALG::CreateVector(*dofrowmap,true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  fint_ = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap,true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known at out-of-balance-force
  fresm_ = LINALG::CreateVector(*dofrowmap,false);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // holding also boundary forces due to Dirichlet/Microboundary
  // conditions
  fresm_dirich_ = LINALG::CreateVector(*dofrowmap,false);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",zeros_);
    discret_->SetState("displacement",dis_);
    discret_->Evaluate(p,stiff_,null,fint_,null,null);
    discret_->ClearState();
  }

  //------------------------------------------------------ time step index
//   istep = 0;
//   params_->set<int>("step",istep);

  // Determine dirichtoggle_ and its inverse since boundary conditions for
  // microscale simulations are due to the MicroBoundary condition
  // (and not Dirichlet BC)

  MicroStatic::DetermineToggle();
  MicroStatic::SetUpHomogenization();

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

  // -------------------------- Calculate initial volume of microstructure
  ParameterList p;
  // action for elements
  p.set("action","calc_init_vol");
  p.set("V0", 0.0);
//   discret_->Evaluate(p,null,null,null,null,null);
  discret_->EvaluateCondition(p, null, null, null, null, null, "MicroBoundary");
  V0_ = p.get<double>("V0", -1.0);
  if (V0_ == -1.0)
    dserror("Calculation of initial volume failed");

  // ------------------------- Calculate initial density of microstructure
  // the macroscopic density has to be averaged over the entire
  // microstructural reference volume

  // create the parameters for the discretization
  ParameterList par;
  // action for elements
  par.set("action","calc_homog_stressdens");
  // set density to zero
  par.set("homogdens", 0.0);
  // set flag that only density has to be calculated
  par.set("onlydens", true);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement",zeros_);

  discret_->SetState("displacement",dism_);

  discret_->Evaluate(par,null,null,null,null,null);
  discret_->ClearState();

  density_ = 1/V0_*par.get<double>("homogdens", 0.0);
  if (density_ == 0.0)
    dserror("Density determined from homogenization procedure equals zero!");

  // Check for surface stress conditions due to interfacial phenomena
  vector<DRT::Condition*> surfstresscond(0);
  discret_->GetCondition("SurfaceStress",surfstresscond);
  if (surfstresscond.size())
  {
    surf_stress_man_=rcp(new DRT::SurfStressManager(*discret_));
  }

  return;
} // MicroStatic::MicroStatic


MFSI::Debug dbg;


/*----------------------------------------------------------------------*
 |  do predictor step (public)                               mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStatic::Predictor(const Epetra_SerialDenseMatrix* defgrd)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_->get<double>("total time"     ,0.0);
  double dt          = params_->get<double>("delta time"     ,0.01);
  string convcheck   = params_->get<string>("convcheck"      ,"AbsRes_Or_AbsDis");
  int istep          = params_->get<int>   ("step"            ,0);
  //bool   printscreen = params_->get<bool>  ("print to screen",false);

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (istep != 0 && (convcheck != "AbsRes_And_AbsDis" || convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  // apply new displacements at DBCs -> this has to be done with the
  // mid-displacements since the given macroscopic deformation
  // gradient is evaluated at the mid-point!
  {
    // dism then also holds prescribed new dirichlet displacements
    EvaluateMicroBC(defgrd);
    discret_->ClearState();
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
  }

  //------------------------------- compute interpolated external forces
  fextm_->Scale(0.0);  // we do not have any external forces in the microproblem!

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    disi_->Scale(0.0);
    discret_->SetState("residual displacement",disi_);
    discret_->SetState("displacement",dism_);
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_->Evaluate(p,stiff_,null,fint_,null,null);
    discret_->ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // complete stiffness matrix
    stiff_->Complete();
  }

  //-------------------------------------------- compute residual forces
  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);

    // save this vector for homogenization
    fresm_dirich_->Update(1.0, fresmcopy, 0.0);

    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //dbg.DumpVector("fresm", *discret_, *fresm_);
  //dbg.DumpVector("dism", *discret_, *dism_);


  //------------------------------------------------ build residual norm
//   double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (istep == 0 && (convcheck != "AbsRes_And_AbsDis" || convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

//   if (printscreen)
//     fresm_->Norm2(&fresmnorm);
//   if (!myrank_ && printscreen)
//   {
//     PrintPredictor(convcheck, fresmnorm);
//   }

  return;
} // MicroStatic::Predictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStatic::FullNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_->get<double>("total time"             ,0.0);
  double dt        = params_->get<double>("delta time"             ,0.01);
  int    maxiter   = params_->get<int>   ("max iterations"         ,10);
  string convcheck = params_->get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_->get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_->get<double>("tolerance residual"     ,1.0e-07);
  //bool printscreen = params_->get<bool>  ("print to screen",true);
  //double toldisp = 1.0e-08;

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_->Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<=maxiter)
  {

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_->Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,true);
    else
      solver_->Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,false);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.0,*disi_,1.0);

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",true);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",time);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("residual displacement",disi_);
      discret_->SetState("displacement",dism_);
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_->Evaluate(p,stiff_,null,fint_,null,null);
      discret_->ClearState();

      if (surf_stress_man_!=null)
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
      }
    }

    // complete stiffness matrix
    stiff_->Complete();

    //------------------------------------------ compute residual forces
    // add static mid-balance
    //fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.);
    // blank residual DOFs which are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);

      // save this vector for homogenization
      fresm_dirich_->Update(1.0, fresmcopy, 0.0);

      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    // a short message
//     if (!myrank_ && printscreen)
//     {
//       PrintNewton(printscreen,print_unconv,timer,numiter,maxiter,
//                   fresmnorm,disinorm,convcheck);
//     }

    //--------------------------------- increment equilibrium loop index
    ++numiter;
  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
//      if (!myrank_ && printscreen)
//      {
//        PrintNewton(printscreen,print_unconv,timer,numiter,maxiter,
//                    fresmnorm,disinorm,convcheck);
//      }
  }

  params_->set<int>("num iterations",numiter);

  return;
} // MicroStatic::FullNewton()


/*----------------------------------------------------------------------*
 |  write output (public)                                       lw 02/08|
 *----------------------------------------------------------------------*/
void MicroStatic::Output(RefCountPtr<MicroDiscretizationWriter> output,
                         const double time,
                         const int istep,
                         const double dt)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------

  bool   iodisp        = params_->get<bool>  ("io structural disp"     ,true);
  int    updevrydisp   = params_->get<int>   ("io disp every nstep"    ,1);
  bool   iostress      = params_->get<bool>  ("io structural stress"   ,false);
  int    updevrystress = params_->get<int>   ("io stress every nstep"  ,10);
  bool   iostrain      = params_.get<bool>  ("io structural strain"   ,false);
  int    writeresevry  = params_->get<int>   ("write restart every"    ,0);

  bool isdatawritten = false;

  //------------------------------------------------- write restart step
  if (writeresevry and istep%writeresevry==0)
  {
    output->WriteMesh(istep,time);
    output->NewStep(istep, time);
    output->WriteVector("displacement",dis_);
    isdatawritten = true;

    if (surf_stress_man_!=null)
    {
      RCP<Epetra_Map> surfrowmap=surf_stress_man_->GetSurfRowmap();
      RCP<Epetra_Vector> A_old=rcp(new Epetra_Vector(*surfrowmap, true));
      RCP<Epetra_Vector> con_quot=rcp(new Epetra_Vector(*surfrowmap, true));
      surf_stress_man_->GetHistory(A_old, con_quot);
      output->WriteVector("Aold", A_old);
      output->WriteVector("conquot", con_quot);
    }
  }

  //----------------------------------------------------- output results
  if (iodisp && updevrydisp && istep%updevrydisp==0 && !isdatawritten)
  {
    output->NewStep(istep, time);
    output->WriteVector("displacement",dis_);
    isdatawritten = true;
  }

  //------------------------------------- do stress calculation and output
  if (updevrystress and !(istep%updevrystress) and iostress)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
    Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stress);
    p.set("strain", strain);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",zeros_);
    discret_->SetState("displacement",dis_);
    discret_->Evaluate(p,null,null,null,null,null);
    discret_->ClearState();
    if (!isdatawritten) output->NewStep(istep, time);
    isdatawritten = true;
    output->WriteVector("gauss_stresses_xyz",*stress,*discret_->ElementColMap());
    if (iostrain)
    {
      output_.WriteVector("gauss_strains_xyz",*strain,*discret_.ElementColMap());
    }
  }

  return;
} // MicroStatic::Output()


/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStatic::SetDefaults(ParameterList& params)
{
  params.set<bool>  ("print to screen"        ,false);
  params.set<bool>  ("print to err"           ,false);
  params.set<FILE*> ("err file"               ,NULL);
  params.set<bool>  ("damping"                ,false);
  params.set<double>("damping factor K"       ,0.00001);
  params.set<double>("damping factor M"       ,0.00001);
  params.set<double>("beta"                   ,0.292);
  params.set<double>("gamma"                  ,0.581);
  params.set<double>("alpha m"                ,0.378);
  params.set<double>("alpha f"                ,0.459);
  params.set<double>("total time"             ,0.0);
  params.set<double>("delta time"             ,0.01);
  params.set<int>   ("step"                   ,0);
  params.set<int>   ("nstep"                  ,5);
  params.set<int>   ("max iterations"         ,10);
  params.set<int>   ("num iterations"         ,-1);
  params.set<double>("tolerance displacements",1.0e-07);
  params.set<bool>  ("io structural disp"     ,true);
  params.set<int>   ("io disp every nstep"    ,1);
  params.set<bool>  ("io structural stress"   ,false);
  params.set<int>   ("restart"                ,0);
  params.set<int>   ("write restart every"    ,0);
  // takes values "constant" consistent"
  params.set<string>("predictor"              ,"constant");
  // takes values "full newton" , "modified newton" , "nonlinear cg"
  params.set<string>("equilibrium iteration"  ,"full newton");
  return;
}

/*----------------------------------------------------------------------*
 |  read restart (public)                                       lw 03/08|
 *----------------------------------------------------------------------*/
void MicroStatic::ReadRestart(int step, RCP<Epetra_Vector> dis, string name)
{
  IO::MicroDiscretizationReader reader(discret_, step, name);
  double time  = reader.ReadDouble("time");
  int    rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  reader.ReadVector(dis, "displacement");
  // It does not make any sense to read the mesh and corresponding
  // element based data because we only have one instance of the
  // microscale discretization and surely have different element based
  // data at every Gauss point -> when implementing EAS, we have to
  // think about how to handle element based EAS history data
  // reader.ReadMesh(step);

  // override current time and step with values from file
  params_->set<double>("total time",time);
  params_->set<int>   ("step",rstep);

  if (surf_stress_man_!=null)
  {
    RCP<Epetra_Map> surfmap=surf_stress_man_->GetSurfRowmap();
    RCP<Epetra_Vector> A_old = LINALG::CreateVector(*surfmap,true);
    RCP<Epetra_Vector> con_quot = LINALG::CreateVector(*surfmap,true);
    reader.ReadVector(A_old, "Aold");
    reader.ReadVector(con_quot, "conquot");
    surf_stress_man_->SetHistory(A_old, con_quot);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
MicroStatic::~MicroStatic()
{
  return;
}


void MicroStatic::DetermineToggle()
{
  int np = 0;   // number of prescribed (=boundary) dofs needed for the
                // creation of vectors and matrices for homogenization
                // procedure

  RefCountPtr<DRT::Discretization> dis =
    DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  vector<DRT::Condition*> conds;
  dis->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!dis->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = dis->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      vector<int> dofs = dis->Dof(actnode);
      const unsigned numdf = dofs.size();

      for (unsigned j=0; j<numdf; ++j)
      {
        const int gid = dofs[j];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);

        if ((*dirichtoggle_)[lid] != 1.0)  // be careful not to count dofs more
                                           // than once since nodes belong to
                                           // several surfaces simultaneously
          ++np;

        (*dirichtoggle_)[lid] = 1.0;
      }
    }
  }

  np_ = np;
}

void MicroStatic::EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd)
{
  vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("MicroBoundary condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int j=0; j<nnode; ++j)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[j])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[j]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[j]);

      // nodal coordinates
      const double* x = actnode->X();

      // boundary displacements are prescribed via the macroscopic
      // deformation gradient
      double dism_prescribed[3];
      Epetra_SerialDenseMatrix Du(*defgrd);
      Epetra_SerialDenseMatrix I(3,3);
      I(0,0)=-1.0;
      I(1,1)=-1.0;
      I(2,2)=-1.0;
      Du+=I;

      for (int k=0; k<3;k++)
      {
        double dis = 0.;

        for (int l=0;l<3;l++)
        {
          dis += Du(k, l) * x[l];
        }

        dism_prescribed[k] = dis;
      }

      vector<int> dofs = discret_->Dof(actnode);
      //cout << "dofs:\n" << dofs[0] << "\n" << dofs[1] << "\n" << dofs[2] << endl;

      for (int l=0; l<3; ++l)
      {
        const int gid = dofs[l];

        const int lid = dism_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*dism_)[lid] = dism_prescribed[l];
      }
    }
  }
}

void MicroStatic::SetOldState(RefCountPtr<Epetra_Vector> dis,
                              RefCountPtr<Epetra_Vector> dism,
                              RefCountPtr<DRT::SurfStressManager> surfman)
{
  dis_ = dis;
  dism_ = dism;
  surf_stress_man_ = surfman;
  fext_->PutScalar(0.);     // we do not have any external loads on
                            // the microscale, so assign all components
                            // to zero
}

void MicroStatic::UpdateNewTimeStep(RefCountPtr<Epetra_Vector> dis,
                                    RefCountPtr<Epetra_Vector> dism)
{
  double alphaf = params_->get<double>("alpha f",0.459);
  dis->Update(1.0/(1.0-alphaf), *dism, -alphaf/(1.0-alphaf));
  dism->Update(1.0, *dis, 0.0);

  if (surf_stress_man_!=null)
  {
    surf_stress_man_->Update();
  }
}

void MicroStatic::SetTime(double timen, int istep)
{
  params_->set<double>("total time", timen);
  params_->set<int>   ("step", istep);
}

RefCountPtr<Epetra_Vector> MicroStatic::ReturnNewDism() { return rcp(new Epetra_Vector(*dism_)); }

void MicroStatic::ClearState()
{
  dis_ = null;
  dism_ = null;
}

void MicroStatic::SetUpHomogenization()
{
  int indp = 0;
  int indf = 0;

  ndof_ = discret_->DofRowMap()->NumMyElements();

  std::vector <int>   pdof(np_);
  std::vector <int>   fdof(ndof_-np_);        // changed this, previously this
                                              // has been just fdof(np_),
                                              // but how should that
                                              // work for ndof_-np_>np_???

  for (int it=0; it<ndof_; ++it)
  {
    if ((*dirichtoggle_)[it] == 1.0)
    {
      pdof[indp]=discret_->DofRowMap()->GID(it);
      ++indp;
    }
    else
    {
      fdof[indf]=discret_->DofRowMap()->GID(it);
      ++indf;
    }
  }

  // create map based on the determined dofs of prescribed and free nodes
  pdof_ = rcp(new Epetra_Map(-1, np_, &pdof[0], 0, discret_->Comm()));
  fdof_ = rcp(new Epetra_Map(-1, ndof_-np_, &fdof[0], 0, discret_->Comm()));

  // create importer
  importp_ = rcp(new Epetra_Import(*pdof_, *(discret_->DofRowMap())));
  importf_ = rcp(new Epetra_Import(*fdof_, *(discret_->DofRowMap())));

  // create vector containing material coordinates of prescribed nodes
  Epetra_Vector Xp_temp(*pdof_);

  vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("MicroBoundary condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);

      // nodal coordinates
      const double* x = actnode->X();

      vector<int> dofs = discret_->Dof(actnode);

      for (int k=0; k<3; ++k)
      {
        const int gid = dofs[k];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);

        for (int l=0;l<np_;++l)
        {
          if (pdof[l]==gid)
            Xp_temp[l]=x[k];
        }
      }
    }
  }

  Xp_ = LINALG::CreateVector(*pdof_,true);
  *Xp_ = Xp_temp;

  // now create D and its transpose DT (following Miehe et al., 2002)

  Epetra_Map Dmap(9, 0, Epetra_SerialComm());
  D_ = rcp(new Epetra_MultiVector(Dmap, np_));

  for (int n=0;n<np_/3;++n)
  {
    Epetra_Vector* temp1 = (*D_)(3*n);
    (*temp1)[0] = (*Xp_)[3*n];
    (*temp1)[3] = (*Xp_)[3*n+1];
    (*temp1)[6] = (*Xp_)[3*n+2];
    Epetra_Vector* temp2 = (*D_)(3*n+1);
    (*temp2)[1] = (*Xp_)[3*n+1];
    (*temp2)[4] = (*Xp_)[3*n+2];
    (*temp2)[7] = (*Xp_)[3*n];
    Epetra_Vector* temp3 = (*D_)(3*n+2);
    (*temp3)[2] = (*Xp_)[3*n+2];
    (*temp3)[5] = (*Xp_)[3*n];
    (*temp3)[8] = (*Xp_)[3*n+1];
  }

  Epetra_MultiVector DT(*pdof_, 9);

  for (int n=0;n<np_/3;++n)
  {
    (*(DT(0)))[3*n]   = (*Xp_)[3*n];
    (*(DT(1)))[3*n+1] = (*Xp_)[3*n+1];
    (*(DT(2)))[3*n+2] = (*Xp_)[3*n+2];
    (*(DT(3)))[3*n]   = (*Xp_)[3*n+1];
    (*(DT(4)))[3*n+1] = (*Xp_)[3*n+2];
    (*(DT(5)))[3*n+2] = (*Xp_)[3*n];
    (*(DT(6)))[3*n]   = (*Xp_)[3*n+2];
    (*(DT(7)))[3*n+1] = (*Xp_)[3*n];
    (*(DT(8)))[3*n+2] = (*Xp_)[3*n+1];
  }

  rhs_ = rcp(new Epetra_MultiVector(*(discret_->DofRowMap()), 9));

  for (int i=0;i<9;++i)
  {
    ((*rhs_)(i))->Export(*(DT(i)), *importp_, Insert);
  }
}


/*----------------------------------------------------------------------*
 |  check convergence of Newton iteration (public)              lw 12/07|
 *----------------------------------------------------------------------*/
bool MicroStatic::Converged(const string type, const double disinorm,
                            const double resnorm, const double toldisp,
                            const double tolres)
{
  if (type == "AbsRes_Or_AbsDis")
  {
    return (disinorm<toldisp or resnorm<tolres);
  }
  else if (type == "AbsRes_And_AbsDis")
  {
    return (disinorm<toldisp and resnorm<tolres);
  }
  else if (type == "RelRes_Or_AbsDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    return (disinorm<toldisp or (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_And_AbsDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    return (disinorm<toldisp and (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_Or_RelDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    if (ref_disnorm_ == 0.) ref_disnorm_ = 1.0;
    return ((disinorm/ref_disnorm_)<toldisp or (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_And_RelDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    if (ref_disnorm_ == 0.) ref_disnorm_ = 1.0;

    return ((disinorm/ref_disnorm_)<toldisp and (resnorm/ref_fnorm_)<tolres);
  }
  else
  {
    dserror("Requested convergence check not (yet) implemented");
    return true;
  }
}

/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks   lw 12/07|
 *----------------------------------------------------------------------*/
void MicroStatic::CalcRefNorms()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  dis_->Norm2(&ref_disnorm_);


  double fintnorm, fextnorm;
  fint_->Norm2(&fintnorm);
  fextm_->Norm2(&fextnorm);

  ref_fnorm_=max(fintnorm, fextnorm);
}

/*----------------------------------------------------------------------*
 |  print to screen and/or error file                           lw 12/07|
 *----------------------------------------------------------------------*/
void MicroStatic::PrintNewton(bool printscreen, bool print_unconv,
                              Epetra_Time timer, int numiter,
                              int maxiter, double fresmnorm, double disinorm,
                              string convcheck)
{
  bool relres        = (convcheck == "RelRes_And_AbsDis" || convcheck == "RelRes_Or_AbsDis");
  bool relres_reldis = (convcheck == "RelRes_And_RelDis" || convcheck == "RelRes_Or_RelDis");

  if (relres)
  {
    fresmnorm /= ref_fnorm_;
  }
  if (relres_reldis)
  {
    fresmnorm /= ref_fnorm_;
    disinorm  /= ref_disnorm_;
  }

  if (print_unconv)
  {
    if (printscreen)
    {
      if (relres)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
      else
        {
        printf("      MICROSCALE numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
    }
  }
  else
  {
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e absolute dis-norm %e time %10.5f\n\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e scaled dis-norm %e time %10.5f\n\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d absolute res-norm %e absolute dis-norm %e time %10.5f\n\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
  }
}

/*----------------------------------------------------------------------*
 |  print to screen                                             lw 12/07|
 *----------------------------------------------------------------------*/
void MicroStatic::PrintPredictor(string convcheck, double fresmnorm)
{
  if (convcheck != "AbsRes_And_AbsDis" && convcheck != "AbsRes_Or_AbsDis")
  {
    fresmnorm /= ref_fnorm_;
    cout << "      MICROSCALE Predictor scaled res-norm " << fresmnorm << endl;
  }
  else
  {
    cout << "      MICROSCALE Predictor absolute res-norm " << fresmnorm << endl;
  }
  fflush(stdout);
}




// after having finished all the testing, remove cmat from the input
// parameters, since no constitutive tensor is calculated here.

void MicroStatic::Homogenization(Epetra_SerialDenseVector* stress,
                                       Epetra_SerialDenseMatrix* cmat,
                                       double *density,
                                       const Epetra_SerialDenseMatrix* defgrd,
                                       const string action)
{
  // determine macroscopic parameters via averaging (homogenization) of
  // microscopic features
  // this was implemented against the background of serial usage
  // -> if a parallel version of microscale simulations is EVER wanted,
  // carefully check if/what/where things have to change

  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  p.set("action","calc_homog_stressdens");
  // choose what to assemble
  p.set("assemble matrix 1",false);
  p.set("assemble matrix 2",false);
  p.set("assemble vector 1",false);
  p.set("assemble vector 2",false);
  p.set("assemble vector 3",false);
  // set stresses and densities to zero
//   p.set("homogdens", 0.0);

  p.set("homogP11", 0.0);
  p.set("homogP12", 0.0);
  p.set("homogP13", 0.0);
  p.set("homogP21", 0.0);
  p.set("homogP22", 0.0);
  p.set("homogP23", 0.0);
  p.set("homogP31", 0.0);
  p.set("homogP32", 0.0);
  p.set("homogP33", 0.0);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement",zeros_);

  // we have to distinguish here whether we are in a homogenization
  // procedure during the nonlinear solution or in the post-processing
  // (macroscopic stress calculation) phase -> in the former case, the
  // displacement has to be chosen at the generalized mid-point, in
  // the latter one we need to evaluate everything at the end of the
  // time step
  if (action == "stress_calc")
    discret_->SetState("displacement",dis_);
  else
    discret_->SetState("displacement",dism_);
  discret_->Evaluate(p,null,null,null,null,null);
  discret_->ClearState();

  // homogenized density was already determined in the constructor
  *density = density_;

  Epetra_SerialDenseMatrix P(3, 3);
  P(0, 0) = 1/V0_*p.get<double>("homogP11", 0.0);
  P(0, 1) = 1/V0_*p.get<double>("homogP12", 0.0);
  P(0, 2) = 1/V0_*p.get<double>("homogP13", 0.0);
  P(1, 0) = 1/V0_*p.get<double>("homogP21", 0.0);
  P(1, 1) = 1/V0_*p.get<double>("homogP22", 0.0);
  P(1, 2) = 1/V0_*p.get<double>("homogP23", 0.0);
  P(2, 0) = 1/V0_*p.get<double>("homogP31", 0.0);
  P(2, 1) = 1/V0_*p.get<double>("homogP32", 0.0);
  P(2, 2) = 1/V0_*p.get<double>("homogP33", 0.0);

  // determine inverse of deformation gradient

  Epetra_SerialDenseMatrix F_inv(3,3);

  double detF= (*defgrd)(0,0) * (*defgrd)(1,1) * (*defgrd)(2,2)
             + (*defgrd)(0,1) * (*defgrd)(1,2) * (*defgrd)(2,0)
             + (*defgrd)(0,2) * (*defgrd)(1,0) * (*defgrd)(2,1)
             - (*defgrd)(0,0) * (*defgrd)(1,2) * (*defgrd)(2,1)
             - (*defgrd)(0,1) * (*defgrd)(1,0) * (*defgrd)(2,2)
             - (*defgrd)(0,2) * (*defgrd)(1,1) * (*defgrd)(2,0);

  F_inv(0,0) = ((*defgrd)(1,1)*(*defgrd)(2,2)-(*defgrd)(1,2)*(*defgrd)(2,1))/detF;
  F_inv(0,1) = ((*defgrd)(0,2)*(*defgrd)(2,1)-(*defgrd)(2,2)*(*defgrd)(0,1))/detF;
  F_inv(0,2) = ((*defgrd)(0,1)*(*defgrd)(1,2)-(*defgrd)(1,1)*(*defgrd)(0,2))/detF;
  F_inv(1,0) = ((*defgrd)(1,2)*(*defgrd)(2,0)-(*defgrd)(2,2)*(*defgrd)(1,0))/detF;
  F_inv(1,1) = ((*defgrd)(0,0)*(*defgrd)(2,2)-(*defgrd)(2,0)*(*defgrd)(0,2))/detF;
  F_inv(1,2) = ((*defgrd)(0,2)*(*defgrd)(1,0)-(*defgrd)(1,2)*(*defgrd)(0,0))/detF;
  F_inv(2,0) = ((*defgrd)(1,0)*(*defgrd)(2,1)-(*defgrd)(2,0)*(*defgrd)(1,1))/detF;
  F_inv(2,1) = ((*defgrd)(0,1)*(*defgrd)(2,0)-(*defgrd)(2,1)*(*defgrd)(0,0))/detF;
  F_inv(2,2) = ((*defgrd)(0,0)*(*defgrd)(1,1)-(*defgrd)(1,0)*(*defgrd)(0,1))/detF;

  for (int i=0; i<3; ++i)
  {
    (*stress)[0] += F_inv(0, i)*P(i,0);                     // S11
    (*stress)[1] += F_inv(1, i)*P(i,1);                     // S22
    (*stress)[2] += F_inv(2, i)*P(i,2);                     // S33
    (*stress)[3] += F_inv(0, i)*P(i,1);                     // S12
    (*stress)[4] += F_inv(1, i)*P(i,2);                     // S23
    (*stress)[5] += F_inv(0, i)*P(i,2);                     // S13
  }

  // for testing reasons only!!!!!!!!!! begin testing region
  const double Emod  = 100.0;
  const double nu  = 0.2;

  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

  // Right Cauchy-Green tensor = F^T * F
  Epetra_SerialDenseMatrix cauchygreen(3,3);
  cauchygreen.Multiply('T','N',1.0,*defgrd,*defgrd,1.0);

  // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  Epetra_SerialDenseVector glstrain(6);
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);

  Epetra_SerialDenseVector ref_stress(6);
  (*cmat).Multiply('N',glstrain,ref_stress);
  /// end testing region!!!!!!!!!
}

void MicroStatic::StaticHomogenization(Epetra_SerialDenseVector* stress,
                                             Epetra_SerialDenseMatrix* cmat,
                                             double *density,
                                             const Epetra_SerialDenseMatrix* defgrd)
{
  // determine macroscopic parameters via averaging (homogenization) of
  // microscopic features accoring to Kouznetsova, Miehe etc.
  // this was implemented against the background of serial usage
  // -> if a parallel version of microscale simulations is EVER wanted,
  // carefully check if/what/where things have to change

  // split microscale stiffness and residual forces into parts
  // corresponding to prescribed and free dofs -> see thesis
  // of Kouznetsova (Computational homogenization for the multi-scale
  // analysis of multi-phase materials, Eindhoven, 2002)

  // split residual forces -> we want to extract fp

  RefCountPtr<DRT::Discretization> dis =
    DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  Epetra_Vector fp(*pdof_);

  int err = fp.Import(*fresm_dirich_, *importp_, Insert);
  if (err)
    dserror("Importing external forces of prescribed dofs using importer returned err=%d",err);

  // Now we have all forces in the material description acting on the
  // boundary nodes together in one vector
  // -> for calculating the stresses, we need to choose the
  // right three components corresponding to a single node and
  // take the inner product with the material coordinates of this
  // node. The sum over all boundary nodes delivers the first
  // Piola-Kirchhoff macroscopic stress which has to be transformed
  // into the second Piola-Kirchhoff counterpart.
  // All these complicated conversions are necessary since only for
  // the energy-conjugated pair of first Piola-Kirchhoff and
  // deformation gradient the averaging integrals can be transformed
  // into integrals over the boundaries only in case of negligible
  // inertial forces (which simplifies matters significantly) whereas
  // the calling macroscopic material routine demands a second
  // Piola-Kirchhoff stress tensor.

  // IMPORTANT: the RVE has to be centered around (0,0,0), otherwise
  // this approach does not work. This was also confirmed by
  // Kouznetsova in a discussion during USNCCM 9.

  fp.Scale(-1.0);

  Epetra_SerialDenseMatrix P(3,3);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      for (int n=0; n<np_/3; ++n)
      {
        P(i,j) += fp[n*3+i]*(*Xp_)[n*3+j];
      }
      P(i,j) /= V0_;
    }
  }

  // determine inverse of deformation gradient

  Epetra_SerialDenseMatrix F_inv(3,3);

  double detF= (*defgrd)(0,0) * (*defgrd)(1,1) * (*defgrd)(2,2)
             + (*defgrd)(0,1) * (*defgrd)(1,2) * (*defgrd)(2,0)
             + (*defgrd)(0,2) * (*defgrd)(1,0) * (*defgrd)(2,1)
             - (*defgrd)(0,0) * (*defgrd)(1,2) * (*defgrd)(2,1)
             - (*defgrd)(0,1) * (*defgrd)(1,0) * (*defgrd)(2,2)
             - (*defgrd)(0,2) * (*defgrd)(1,1) * (*defgrd)(2,0);

  F_inv(0,0) = ((*defgrd)(1,1)*(*defgrd)(2,2)-(*defgrd)(1,2)*(*defgrd)(2,1))/detF;
  F_inv(0,1) = ((*defgrd)(0,2)*(*defgrd)(2,1)-(*defgrd)(2,2)*(*defgrd)(0,1))/detF;
  F_inv(0,2) = ((*defgrd)(0,1)*(*defgrd)(1,2)-(*defgrd)(1,1)*(*defgrd)(0,2))/detF;
  F_inv(1,0) = ((*defgrd)(1,2)*(*defgrd)(2,0)-(*defgrd)(2,2)*(*defgrd)(1,0))/detF;
  F_inv(1,1) = ((*defgrd)(0,0)*(*defgrd)(2,2)-(*defgrd)(2,0)*(*defgrd)(0,2))/detF;
  F_inv(1,2) = ((*defgrd)(0,2)*(*defgrd)(1,0)-(*defgrd)(1,2)*(*defgrd)(0,0))/detF;
  F_inv(2,0) = ((*defgrd)(1,0)*(*defgrd)(2,1)-(*defgrd)(2,0)*(*defgrd)(1,1))/detF;
  F_inv(2,1) = ((*defgrd)(0,1)*(*defgrd)(2,0)-(*defgrd)(2,1)*(*defgrd)(0,0))/detF;
  F_inv(2,2) = ((*defgrd)(0,0)*(*defgrd)(1,1)-(*defgrd)(1,0)*(*defgrd)(0,1))/detF;

  // convert to second Piola-Kirchhoff stresses and store them in
  // vector format
  // assembly of stresses (cf Solid3 Hex8): S11,S22,S33,S12,S23,S13

  Epetra_SerialDenseVector S(6);

  for (int i=0; i<3; ++i)
  {
    S[0] += F_inv(0, i)*P(i,0);                     // S11
    S[1] += F_inv(1, i)*P(i,1);                     // S22
    S[2] += F_inv(2, i)*P(i,2);                     // S33
    S[3] += F_inv(0, i)*P(i,1);                     // S12
    S[4] += F_inv(1, i)*P(i,2);                     // S23
    S[5] += F_inv(0, i)*P(i,2);                     // S13
  }

  for (int i=0; i<6; ++i)
  {
    (*stress)[i]=S[i];
  }

  // The calculation of the consistent macroscopic constitutive tensor
  // follows
  //
  // C. Miehe, Computational micro-to-macro transitions for
  // discretized micro-structures of heterogeneous materials at finite
  // strains based on a minimization of averaged incremental energy.
  // Computer Methods in Applied Mechanics and Engineering 192: 559-591, 2003.

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Epetra_MultiVector cmatpf(D_->Map(), 9);

  Epetra_Vector x(*dofrowmap);
  Epetra_Vector y(*dofrowmap);

  // make a copy
  stiff_dirich_ = Teuchos::rcp(new LINALG::SparseMatrix(*stiff_));

  stiff_->ApplyDirichlet(dirichtoggle_);

  Epetra_LinearProblem linprob(&(*stiff_->EpetraMatrix()), &x, &y);
  int error=linprob.CheckInput();
  if (error)
    dserror("Input for linear problem inconsistent");
  Amesos_Klu solver(linprob);
  err = solver.NumericFactorization();   // LU decomposition of stiff_ only once
  if (err)
    dserror("Numeric factorization of stiff_ for homogenization failed");

  for (int i=0;i<9;++i)
  {
    x.PutScalar(0.0);
    y.Update(1.0, *((*rhs_)(i)), 0.0);
    solver.Solve();

    Epetra_Vector f(*dofrowmap);
    stiff_dirich_->Multiply(false, x, f);
    Epetra_Vector fexp(*pdof_);
    int err = fexp.Import(f, *importp_, Insert);
    if (err)
      dserror("Export of boundary 'forces' failed with err=%d", err);

    (cmatpf(i))->Multiply('N', 'N', 1.0/V0_, *D_, fexp, 0.0);
  }

  // We now have to transform the calculated constitutive tensor
  // relating first Piola-Kirchhoff stresses to the deformation
  // gradient into a constitutive tensor relating second
  // Piola-Kirchhoff stresses to Green-Lagrange strains.

  ConvertMat(cmatpf, F_inv, *stress, *cmat);

  // homogenized density was already determined in the constructor

  *density = density_;

  //cout << "cmat:\n" << *cmat << "\nstress:\n" << *stress << "\n";
}


#endif
