/*!----------------------------------------------------------------------
\file microstrugenalpha.cpp
\brief Generalized Alpha time integration for microstructural problems
in case of multiscale analyses

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "microstrugenalpha.H"

#include <vector>

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../io/io_drt_micro.H"

using namespace IO;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
MicroStruGenAlpha::MicroStruGenAlpha(RefCountPtr<ParameterList> params,
                                     RefCountPtr<DRT::Discretization> dis,
                                     RefCountPtr<LINALG::Solver> solver) :
params_(params),
discret_(dis),
solver_(solver)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time    = params_->get<double>("total time"      ,0.0);
  double dt      = params_->get<double>("delta time"      ,0.01);
  bool damping   = params_->get<bool>  ("damping"         ,false);
  double kdamp   = params_->get<double>("damping factor K",0.0);
  double mdamp   = params_->get<double>("damping factor M",0.0);
  int istep      = params_->get<int>   ("step"            ,0);

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
  stiff_ = LINALG::CreateMatrix(*dofrowmap,81);
  mass_  = LINALG::CreateMatrix(*dofrowmap,81);
  if (damping) damp_ = LINALG::CreateMatrix(*dofrowmap,81);

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
  // velocities V_{n} at last time
  vel_ = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n} at last time
  acc_ = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  disn_ = LINALG::CreateVector(*dofrowmap,true);
  // velocities V_{n+1} at new time
  veln_ = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n+1} at new time
  accn_ = LINALG::CreateVector(*dofrowmap,true);

  // mid-displacements D_{n+1-alpha_f}
  dism_ = LINALG::CreateVector(*dofrowmap,true);
  // mid-velocities V_{n+1-alpha_f}
  velm_ = LINALG::CreateVector(*dofrowmap,true);
  // mid-accelerations A_{n+1-alpha_m}
  accm_ = LINALG::CreateVector(*dofrowmap,true);

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

  //-------------------------------------------- calculate external forces
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("displacement",dis_);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    discret_->EvaluateDirichlet(p,*dis_,*dirichtoggle_);
    discret_->ClearState();
    discret_->SetState("displacement",dis_);
    // predicted rhs
    discret_->EvaluateNeumann(p,*fext_);
    discret_->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiffmass");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",true);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",zeros_);
    discret_->SetState("displacement",dis_);
    //discret_.SetState("velocity",vel_); // not used at the moment
    discret_->Evaluate(p,stiff_,mass_,fint_,null,null);
    discret_->ClearState();
  }
  // build damping matrix if desired
  LINALG::Complete(*mass_);
  maxentriesperrow_ = mass_->MaxNumEntries();
  if (damping)
  {
    LINALG::Complete(*stiff_);
    LINALG::Add(*stiff_,false,kdamp,*damp_,0.0);
    stiff_ = null;
    LINALG::Add(*mass_,false,mdamp,*damp_,1.0);
    LINALG::Complete(*damp_);
  }

  //--------------------------- calculate consistent initial accelerations
  {
    RefCountPtr<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
    if (damping) damp_->Multiply(false,*vel_,*rhs);
    rhs->Update(-1.0,*fint_,1.0,*fext_,-1.0);
    Epetra_Vector rhscopy(*rhs);
    rhs->Multiply(1.0,*invtoggle_,rhscopy,0.0);
    solver->Solve(mass_,acc_,rhs,true,true);
  }

  //------------------------------------------------------ time step index
  istep = 0;
  params_->set<int>("step",istep);



  // Determine dirichtoggle_ and its inverse since boundary conditions for
  // microscale simulations are due to the MicroBoundary condition and are
  // therefore not taken into account in the Constructor of StruGenAlpha

  MicroStruGenAlpha::DetermineToggle();

  return;
} // MicroStruGenAlpha::MicroStruGenAlpha


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStruGenAlpha::ConstantPredictor()
{
  //cout << "ConstantPredictor\n";

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_->get<double>("total time"     ,0.0);
  double dt          = params_->get<double>("delta time"     ,0.01);
  int    istep       = params_->get<int>   ("step"           ,0);
  bool   damping     = params_->get<bool>  ("damping"        ,false);
  double alphaf      = params_->get<double>("alpha f"        ,0.459);
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  cout << "total time: " << time << " step: " << istep << endl;

  // increment time and step
  double timen = time += dt;
  istep++;
  params_->set<double>("total time",timen);
  params_->set<int>   ("step"      ,istep);

  //--------------------------------------------------- predicting state
  // constant predictor : displacement in domain
  disn_->Update(1.0,*dis_,0.0);

  // apply new displacements at DBCs
  // and get new external force vector
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_->EvaluateDirichlet(p,*disn_,*dirichtoggle_);
    discret_->ClearState();
    discret_->SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_->EvaluateNeumann(p,*fextn_);
    discret_->ClearState();
  }

  // constant predictor
  veln_->Update(1.0,*vel_,0.0);
  accn_->Update(1.0,*acc_,0.0);

  //------------------------------ compute interpolated dis, vel and acc
  // constant predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  velm_->Update(1.0,*vel_,0.0);
  accm_->Update(1.0,*acc_,0.0);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
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
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",disi_);
    discret_->SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_->Evaluate(p,stiff_,null,fint_,null,null);
    discret_->ClearState();
    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*fresm_);
  // add mid-viscous damping force
  if (damping)
  {
      RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
      damp_->Multiply(false,*velm_,*fviscm);
      fresm_->Update(1.0,*fviscm,1.0);
  }

  // add static mid-balance
  fresm_->Update(1.0,*fint_,-1.0,*fextm_,1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  fresm_->Norm2(&norm_);

  return;
} // MicroStruGenAlpha::ConstantPredictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStruGenAlpha::FullNewton()
{
  //cout << "FullNewton\n";

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_->get<double>("total time"             ,0.0);
  double dt        = params_->get<double>("delta time"             ,0.01);
  int    maxiter   = params_->get<int>   ("max iterations"         ,10);
  bool   damping   = params_->get<bool>  ("damping"                ,false);
  double beta      = params_->get<double>("beta"                   ,0.292);
  double gamma     = params_->get<double>("gamma"                  ,0.581);
  double alpham    = params_->get<double>("alpha m"                ,0.378);
  double alphaf    = params_->get<double>("alpha f"                ,0.459);
  double toldisp   = params_->get<double>("tolerance displacements",1.0e-07);
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  while (norm_>toldisp && numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    LINALG::Add(*mass_,false,(1.-alpham)/(beta*dt*dt),*stiff_,1.-alphaf);
    if (damping)
      LINALG::Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),*stiff_,1.0);
    LINALG::Complete(*stiff_);

    //----------------------- apply dirichlet BCs to system of equations
    fresm_->Scale(-1.0);   // delete this by building fresm with other sign
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_->Solve(stiff_,disi_,fresm_,true,true);
    else
      solver_->Solve(stiff_,disi_,fresm_,true,false);
    stiff_ = null;

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    //velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
    // accelerations
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    //accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
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
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_->Evaluate(p,stiff_,null,fint_,null,null);
      discret_->ClearState();
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*fresm_);
    // add viscous mid-forces
    if (damping)
    {
      RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fviscm);
      fresm_->Update(1.0,*fviscm,1.0);
    }
    // add static mid-balance
    fresm_->Update(1.0,*fint_,-1.0,*fextm_,1.0);
    // blank residual DOFs with are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    double disinorm;
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&norm_);

    norm_ = disinorm;

    //--------------------------------- increment equilibrium loop index
    ++numiter;
  } // while (norm_>toldisp && numiter<=maxiter)
  //============================================= end equilibrium loop

  //-------------------------------- test whether max iterations was hit
  if (numiter==maxiter) dserror("Newton unconverged in %d iterations",numiter);
  params_->set<int>("num iterations",numiter);

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // MicroStruGenAlpha::FullNewton()


/*----------------------------------------------------------------------*
 |  do Update (public)                                       mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStruGenAlpha::Update()
{
  //cout << "Update\n";

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time          = params_->get<double>("total time"             ,0.0);
  double dt            = params_->get<double>("delta time"             ,0.01);

  double alpham        = params_->get<double>("alpha m"                ,0.378);
  double alphaf        = params_->get<double>("alpha f"                ,0.459);

  //---------------------------- determine new end-quantities and update
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1} = 1./(1.-alphaf) * D_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * D_n
  dis_->Update(1./(1.-alphaf),*dism_,-alphaf/(1.-alphaf));
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1} = 1./(1.-alphaf) * V_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * V_n
  vel_->Update(1./(1.-alphaf),*velm_,-alphaf/(1.-alphaf));
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1} = 1./(1.-alpham) * A_{n+1-alpha_m}
  //                     - alpham/(1.-alpham) * A_n
  acc_->Update(1./(1.-alpham),*accm_,-alpham/(1.-alpham));
  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0,*fextn_,0.0);

  //----- update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_update_istep");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",false);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    discret_->Evaluate(p,null,null,null,null,null);
  }
  return;
} // MicroStruGenAlpha::Update()


/*----------------------------------------------------------------------*
 |  write output (public)                                    mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStruGenAlpha::Output(RefCountPtr<MicroDiscretizationWriter> output,
                               const double time,
                               const int istep)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------

  bool   iodisp        = params_->get<bool>  ("io structural disp"     ,true);
  int    updevrydisp   = params_->get<int>   ("io disp every nstep"    ,1);

  //----------------------------------------------------- output results
  if (iodisp && updevrydisp && istep%updevrydisp==0)
  {
    cout << "Output for step: " << istep << " at time: " << time << endl;

    output->NewStep(istep, time);
    output->WriteVector("displacement",dis_);
    //output->WriteVector("velocity",vel_);
    //output->WriteVector("acceleration",acc_);
  }

  return;
} // MicroStruGenAlpha::Output()


/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void MicroStruGenAlpha::SetDefaults(ParameterList& params)
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
 |  dtor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
MicroStruGenAlpha::~MicroStruGenAlpha()
{
  return;
}


void MicroStruGenAlpha::DetermineToggle()
{
  //cout << "DetermineToggle\n";

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
        (*dirichtoggle_)[lid] = 1.0;

        //if (lid == 4)
        //  cout << "dirichtoggle of dof 4: 1 \n";
      }
    }
  }

  //cout << "dirichtoggle: " << *dirichtoggle_ << "\n";
}

void MicroStruGenAlpha::EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd)
{
  //cout << "EvaluateMicroBC\n";

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

      // nodal coordinates
      const double* x = actnode->X();

      // boundary displacements are prescribed via the macroscopic
      // deformation gradient
      double disn_prescribed[3];
      for (int i=0; i<3;i++)
      {
        double dis = 0.;

        for (int j=0;j<3;j++)
        {
          dis += (*defgrd)(i, j) * x[j];
        }

        disn_prescribed[i] = dis - x[i];
      }

      vector<int> dofs = dis->Dof(actnode);

      for (int k=0; k<3; ++k)
      {
        const int gid = dofs[k];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*disn_)[lid] = disn_prescribed[k];

        //if (lid == 4)
        //{
          //cout << "deformation gradient: " << *defgrd << endl;
          //cout << "prescribed dof " << lid << ": " << (*disn_)[lid] << endl;
        //}
      }
    }
  }

}

void MicroStruGenAlpha::SetOldState(RefCountPtr<Epetra_Vector> disp,
                                    RefCountPtr<Epetra_Vector> vel,
                                    RefCountPtr<Epetra_Vector> acc)
{
  //cout << "SetOldState\n";

  dis_ = disp;
  vel_ = vel;
  acc_ = acc;
}

void MicroStruGenAlpha::SetTime(double timen, int istep)
{
  params_->set<double>("total time", timen);
  params_->set<int>   ("step", istep);
}

RefCountPtr<Epetra_Vector> MicroStruGenAlpha::ReturnNewDisp() { return rcp(new Epetra_Vector(*disn_)); }

RefCountPtr<Epetra_Vector> MicroStruGenAlpha::ReturnNewVel() { return rcp(new Epetra_Vector(*veln_)); }

RefCountPtr<Epetra_Vector> MicroStruGenAlpha::ReturnNewAcc() { return rcp(new Epetra_Vector(*accn_)); }

void MicroStruGenAlpha::ClearState()
{
  //cout << "ClearState\n";

  dis_ = null;
  vel_ = null;
  acc_ = null;
}

#endif
#endif
