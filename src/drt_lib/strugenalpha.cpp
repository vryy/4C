/*!----------------------------------------------------------------------
\file strugenalpha.cpp
\brief Generalized Alpha time integration for structural problems

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "strugenalpha.H"
#include "iostream"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
StruGenAlpha::StruGenAlpha(ParameterList& params,
                           DRT::Discretization& dis,
                           LINALG::Solver& solver,
                           IO::DiscretizationWriter& output) :
params_(params),
discret_(dis),
solver_(solver),
output_(output)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time    = params_.get<double>("total time"      ,0.0);
  double dt      = params_.get<double>("delta time"      ,0.01);
  bool damping   = params_.get<bool>  ("damping"         ,false);
  double kdamp   = params_.get<double>("damping factor K",0.0);
  double mdamp   = params_.get<double>("damping factor M",0.0);
  int istep      = params_.get<int>   ("step"            ,0);
  bool outerr    = params_.get<bool>  ("print to err"    ,false);
  FILE* errfile  = params_.get<FILE*> ("err file"        ,NULL);
  if (!errfile) outerr = false;

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  if (!discret_.Filled()) discret_.FillComplete();
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
  myrank_ = discret_.Comm().MyPID();

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
    discret_.ClearState();
    discret_.SetState("displacement",dis_);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,*dis_,*dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",dis_);
    // predicted rhs
    discret_.EvaluateNeumann(p,*fext_);
    discret_.ClearState();
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
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    //discret_.SetState("velocity",vel_); // not used at the moment
    discret_.Evaluate(p,stiff_,mass_,fint_,null,null);
    discret_.ClearState();
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
    solver.Solve(mass_,acc_,rhs,true,true);
  }

  //------------------------------------------------------ time step index
  istep = 0;
  params_.set<int>("step",istep);


  return;
} // StruGenAlpha::StruGenAlpha


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ConstantPredictor()
{

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    istep       = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // increment time and step
  double timen = time += dt;
  istep++;
  params_.set<double>("total time",timen);
  params_.set<int>   ("step"      ,istep);

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
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,*disn_,*dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
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
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();
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
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  fresm_->Norm2(&norm_);
  if (!myrank_ && printscreen)
  {
    cout << "Predictor residual forces " << norm_ << endl;
    fflush(stdout);
  }

  return;
} // StruGenAlpha::ConstantPredictor()


/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ConsistentPredictor()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    istep       = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  double alpham      = params_.get<double>("alpha m"        ,0.378);
  double beta        = params_.get<double>("beta"           ,0.292);
  double gamma       = params_.get<double>("gamma"          ,0.581);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // increment time and step
  double timen = time += dt;
  istep++;
  params_.set<double>("total time",timen);
  params_.set<int>   ("step"      ,istep);

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
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,*disn_,*dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // consistent predictor
  // predicting velocity V_{n+1} (veln)
  // V_{n+1} := gamma/(beta*dt) * (D_{n+1} - D_n)
  //          + (beta-gamma)/beta * V_n
  //          + (2.*beta-gamma)/(2.*beta) * A_n
  veln_->Update(1.0,*disn_,-1.0,*dis_,0.0);
  veln_->Update((beta-gamma)/beta,*vel_,(2.*beta-gamma)*dt/(2.*beta),*acc_,gamma/(beta*dt));
  // predicting accelerations A_{n+1} (accn)
  // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
  //          - 1./(beta*dt) * V_n
  //          + (2.*beta-1.)/(2.*beta) * A_n
  accn_->Update(1.0,*disn_,-1.0,*dis_,0.0);
  accn_->Update(-1./(beta*dt),*vel_,(2.*beta-1.)/(2.*beta),*acc_,1./(beta*dt*dt));

  //------------------------------ compute interpolated dis, vel and acc
  // consistent predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
  velm_->Update(1.-alphaf,*veln_,alphaf,*vel_,0.0);
  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
  accm_->Update(1.-alpham,*accn_,alpham,*acc_,0.0);

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
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();
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
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  fresm_->Norm2(&norm_);
  if (!myrank_ && printscreen)
  {
    cout << "Predictor residual forces " << norm_ << endl;
    fflush(stdout);
  }

  return;
} // StruGenAlpha::ConsistentPredictor()

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::FullNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm;
  double disinorm;
  fresm_->Norm2(&fresmnorm);
  while (norm_>toldisp && fresmnorm>toldisp && numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    LINALG::Add(*mass_,false,(1.-alpham)/(beta*dt*dt),*stiff_,1.-alphaf);
    if (damping)
      LINALG::Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),*stiff_,1.0);
    LINALG::Complete(*stiff_);

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_,disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_,disi_,fresm_,true,false);
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
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();
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
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);
    // a short message
    if (!myrank_)
    {
      if (printscreen)
      {
        printf("numiter %d res-norm %e dis-norm %e\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
      if (printerr)
      {
        fprintf(errfile,"numiter %d res-norm %e dis-norm %e\n",numiter+1, fresmnorm, disinorm);
        fflush(errfile);
      }
    }
    // criteria to stop Newton iteration
    norm_ = disinorm;

    //--------------------------------- increment equilibrium loop index
    ++numiter;

  } // while (norm_>toldisp && fresmnorm>toldisp && numiter<=maxiter)
  //=================================================================== end equilibrium loop

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter) dserror("Newton unconverged in %d iterations",numiter);
  params_.set<int>("num iterations",numiter);

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // StruGenAlpha::FullNewton()


/*----------------------------------------------------------------------*
 |  do modified Newton iteration (public)                    mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ModifiedNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  int numiter=0;
  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  LINALG::Add(*mass_,false,(1.-alpham)/(beta*dt*dt),*stiff_,1.-alphaf);
  if (damping)
    LINALG::Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),*stiff_,1.0);
  LINALG::Complete(*stiff_);
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

  //=================================================== equilibrium loop
  double fresmnorm;
  fresm_->Norm2(&fresmnorm);
  while (norm_>toldisp && fresmnorm>toldisp  && numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_,disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_,disi_,fresm_,false,false);

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

    //----------------------------------------- compute internal forces
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",false);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",time);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,null,null,fint_,null,null);
      discret_.ClearState();
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
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs with are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    double disinorm;
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&norm_);
    // a short message
    if (!myrank_)
    {
      if (printscreen)
      {
        printf("numiter %d res-norm %e dis-norm %e\n",numiter+1, norm_, disinorm);
        fflush(stdout);
      }
      if (printerr)
      {
        fprintf(errfile,"numiter %d res-norm %e dis-norm %e\n",numiter+1, norm_, disinorm);
        fflush(errfile);
      }
    }
    // criteria to stop Newton iteration
    norm_ = disinorm;

    //--------------------------------- increment equilibrium loop index
    ++numiter;
  } // while (norm_>toldisp && fresmnorm>toldisp && numiter<=maxiter)
  //============================================= end equilibrium loop

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter) dserror("Newton unconverged in %d iterations",numiter);
  params_.set<int>("num iterations",numiter);

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // StruGenAlpha::ModifiedNewton()


/*----------------------------------------------------------------------*
 |  do nonlinear cg iteration (public)                       mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::NonlinearCG()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  //double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  //bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;

  // get ml stuff from input file
  bool hasml = solver_.Params().isSublist("ML Parameters");
  if (!hasml) dserror("Solver does not have ML configured from input file");
  ParameterList& linearmllist = solver_.Params().sublist("ML Parameters");
  int outlevel = linearmllist.get<int>("output",0);
  int maxlevel = linearmllist.get<int>("max levels",3);
  int maxcsize = linearmllist.get<int>("coarse: max size",1);

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_==null)     dserror("stiffness matrix = null");
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //---------------------------------------------- set initial guess
  disi_->PutScalar(0.0);

  // create the nox parameter list if it does not exist
  ParameterList& noxparams = params_.sublist("nox parameters");
  RefCountPtr<Teuchos::ParameterList> rcpparams = rcp(&noxparams);
  rcpparams.release();

  ParameterList& printParams = noxparams.sublist("Printing");
  printParams.set("MyPID",myrank_);
  printParams.set("Output Precision", 9);
  printParams.set("Output Processor", 0);
  if (outlevel)
    printParams.set("Output Information",
                    NOX::Utils::OuterIteration +
                    //NOX::Utils::OuterIterationStatusTest +
                    //NOX::Utils::InnerIteration +
                    //NOX::Utils::Parameters +
                    //NOX::Utils::Details +
                    NOX::Utils::Warning
                    );
  else
    printParams.set("Output Information",0);

  // Set the nonlinear solver method as line search
  noxparams.set("Nonlinear Solver","Line Search Based");

  // get sublist for type of linesearch
  ParameterList& searchParams = noxparams.sublist("Line Search");
  searchParams.set("Method","NonlinearCG");

  // Sublist for direction
  ParameterList& dirParams = noxparams.sublist("Direction");
  dirParams.set("Method", "NonlinearCG");

  ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
  nlcgParams.set("Precondition","On");
  nlcgParams.set("Orthogonalize","Polak-Ribiere");
  //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
  nlcgParams.set("Restart Frequency", 25);

  ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 100);
  lsParams.set("Tolerance", 1e-11);
  lsParams.set("Output Frequency", 10);
  //lsParams.set("Preconditioning", "None");
  //lsParams.set("Preconditioner","None");
  lsParams.set("Preconditioning", "User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

  // create the nonlinear ml parameter list
  ParameterList& mlparams = params_.sublist("ml parameters");
  mlparams.set("nlnML output",                                      outlevel   ); // ML-output-level (0-10)
  mlparams.set("nlnML max levels",                                  maxlevel   ); // max. # levels (minimum = 2 !)
  mlparams.set("nlnML coarse: max size",                            maxcsize   ); // the size ML stops generating coarser levels
  mlparams.set("nlnML is linear preconditioner",                    true      );
  mlparams.set("nlnML is matrixfree",                               false      );
  mlparams.set("nlnML apply constraints",                           false      );
  mlparams.set("nlnML Jacobian fix diagonal",                       false      );
  mlparams.set("nlnML finite difference fine level",                false      );
  mlparams.set("nlnML finite difference alpha",                     1.0e-08    );
  mlparams.set("nlnML finite difference beta",                      1.0e-07    );
  mlparams.set("nlnML finite difference centered",                  false      );

  mlparams.set("nlnML absolute residual tolerance",                 toldisp    );
  mlparams.set("nlnML max cycles",                                  maxiter    );
  mlparams.set("nlnML adaptive recompute",                          0.0        ); // recompute if residual is larger then this value
  mlparams.set("nlnML offset recompute",                            0          ); // every offset this preconditioner is recomputed
  mlparams.set("nlnML additional adaptive nullspace",               0          ); // compute adaptive nullspace (additional kernel vectors)
  mlparams.set("nlnML PDE equations",                               6          ); // dof per node
  mlparams.set("nlnML null space: dimension",                       6          ); // dimension of nullspace
  mlparams.set("nlnML spatial dimension",                           3          );
  mlparams.set("nlnML coarse: type",                                "Uncoupled"); // Uncoupled METIS VBMETIS
  mlparams.set("nlnML nodes per aggregate",                         9          ); // # nodes per agg for coarsening METIS and VBMETIS

  mlparams.set("nlnML use nlncg on fine level",                     true); // use nlnCG or mod. Newton's method
  mlparams.set("nlnML use nlncg on medium level",                   true);
  mlparams.set("nlnML use nlncg on coarsest level",                 true);

  mlparams.set("nlnML max iterations newton-krylov fine level",     100); // # iterations of lin. CG in mod. Newton's method
  mlparams.set("nlnML max iterations newton-krylov medium level" ,  50);
  mlparams.set("nlnML max iterations newton-krylov coarsest level", 150);

  mlparams.set("nlnML linear smoother type fine level",             "MLS"); // MLS SGS BSGS Jacobi MLS Bcheby AmesosKLU
  mlparams.set("nlnML linear smoother type medium level",           "MLS");
  mlparams.set("nlnML linear smoother type coarsest level",         "AmesosKLU");
  mlparams.set("nlnML linear smoother sweeps fine level",           6);
  mlparams.set("nlnML linear smoother sweeps medium level",         6);
  mlparams.set("nlnML linear smoother sweeps coarsest level",       1);

  mlparams.set("nlnML nonlinear presmoothing sweeps fine level",    0);
  mlparams.set("nlnML nonlinear presmoothing sweeps medium level",  0);
  mlparams.set("nlnML nonlinear smoothing sweeps coarse level",     15);
  mlparams.set("nlnML nonlinear postsmoothing sweeps medium level", 3);
  mlparams.set("nlnML nonlinear postsmoothing sweeps fine level",   10);

  // create the fine level interface if it does not exist
  if (fineinterface_==null)
  {
    int printlevel = mlparams.get("nlnML output",6);
    fineinterface_ = rcp(new NoxInterface(*this,printlevel));
  }

  // create the nonlinear ml preconditioner if it does not exist
  if (prec_==null)
    prec_ = rcp(new NLNML::NLNML_Preconditioner(fineinterface_,mlparams,discret_.Comm()));
  // tell preconditioner to recompute from scratch
  prec_->setinit(false);

#if 0 // use the nonlinear preconditioner as a solver without the outer nox loop
  {
    Epetra_Time timer(discret_.Comm());
    double t0 = timer.ElapsedTime();
    //prec_->solve_variant();
    prec_->solve();
    double t1 = timer.ElapsedTime();

    // get status and print output message
    if (outlevel && myrank_==0)
    {
      printf("NOX/ML :============solve time incl. setup : %15.4f sec\n",t1-t0);
      double appltime = fineinterface_->getsumtime();
      printf("NOX/ML :===========of which time in ccarat : %15.4f sec\n",appltime);
      cout << "NOX/ML :======number calls to computeF in this solve : "
           << fineinterface_->getnumcallscomputeF() << "\n\n\n";
      fflush(stdout);
    }
    fineinterface_->resetsumtime();
    fineinterface_->setnumcallscomputeF(0);

  //---------------------------------- update mid configuration values
  // displacements
  // incremental (disi is now D_{n+1}-D_{n})
  dism_->Update((1.-alphaf),*disi_,1.0,*dis_,0.0);

  // velocities
  // incremental (required for constant predictor)
  velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                gamma/(beta*dt));

  // accelerations
  // incremental (required for constant predictor)
  accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                (2.*beta-1.+alpham)/(2.*beta),*acc_,
                (1.-alpham)/((1.-alphaf)*beta*dt*dt));

    stiff_ = null;
    return;
  }
#endif


  // create a matrix free operator if it does not exist
  if (matfreeoperator_==null)
    matfreeoperator_ = rcp(new NOX::Epetra::MatrixFree(printParams,fineinterface_,*disi_,false));

  // create a linear system if it does not exist
  if (rcpazlinsys_==null)
  {
    NOX::Epetra::Vector initialGuess(*disi_);
    const RefCountPtr<NOX::Epetra::Interface::Jacobian>       ijac  = matfreeoperator_;
    const RefCountPtr<Epetra_Operator>                        jac   = matfreeoperator_;
    const RefCountPtr<NOX::Epetra::Interface::Preconditioner> iprec = prec_;
    const RefCountPtr<Epetra_Operator>                        prec  = prec_;
    rcpazlinsys_ = rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,lsParams,
                                                            ijac,
                                                            jac,
                                                            iprec,
                                                            prec,
                                                            initialGuess));

    //const RefCountPtr<NOX::Epetra::Interface::Required>       ireq  = fineinterface_;
    //const RefCountPtr<NOX::Epetra::Interface::Preconditioner> iprec = prec_;
    //const RefCountPtr<Epetra_Operator>                        prec  = prec_;
    //rcpazlinsys_ = rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,lsParams,
    //                                                        ireq,
    //                                                        iprec,
    //                                                        prec,
    //                                                        initialGuess));
  }

  // create a group if it does not exist
  NOX::Epetra::Vector initialGuess(*disi_);
  RefCountPtr<NOX::Epetra::Group> rcpgrp = rcp(new NOX::Epetra::Group(printParams,
                                                                      fineinterface_,
                                                                      initialGuess,
                                                                      rcpazlinsys_));

  // create a convergence test if it does not exist
  if (combo_==null)
  {
    RefCountPtr<NOX::StatusTest::NormF> absresid =
      rcp( new NOX::StatusTest::NormF(toldisp));
    RefCountPtr<NOX::StatusTest::NormUpdate> nupdate =
      rcp(new NOX::StatusTest::NormUpdate(toldisp));
    RefCountPtr<NOX::StatusTest::Combo> converged =
      rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(nupdate);
    RefCountPtr<NOX::StatusTest::FiniteValue> fv =
      rcp(new NOX::StatusTest::FiniteValue());
    int maxcycle = mlparams.get("nlnML max cycles",200);
    RefCountPtr<NOX::StatusTest::MaxIters> maxiters =
      rcp(new NOX::StatusTest::MaxIters(maxcycle));
    combo_ = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo_->addStatusTest(maxiters);
    combo_->addStatusTest(converged);
    combo_->addStatusTest(fv);
  }

  // create nox solver manager if it does not exist
  RCP<NOX::Solver::Generic> noxsolver = NOX::Solver::buildSolver(rcpgrp,combo_,rcpparams);
  prec_->SetNoxSolver(noxsolver);

  // solve nonlinear problem
  Epetra_Time timer(discret_.Comm());
  double t0 = timer.ElapsedTime();
  NOX::StatusTest::StatusType status = noxsolver->solve();
  double t1 = timer.ElapsedTime();

  bool converged = true;
  if (status != NOX::StatusTest::Converged) converged = false;

  // get status and print output message
  if (outlevel && myrank_==0)
  {
    printf("NOX/ML :============solve time incl. setup : %15.4f sec\n",t1-t0);
    double appltime = fineinterface_->getsumtime();
    printf("NOX/ML :===========of which time in ccarat : %15.4f sec\n",appltime);
    cout << "NOX/ML :======number calls to computeF in this solve : "
         << fineinterface_->getnumcallscomputeF() << "\n\n\n";
    if (!converged)
      cout << "***WRN***: NOX not converged!\n";
    fflush(stdout);
  }
  fineinterface_->resetsumtime();
  fineinterface_->setnumcallscomputeF(0);

  //---------------------------------- update mid configuration values
  // displacements
  // incremental (disi is now D_{n+1}-D_{n})
  dism_->Update((1.-alphaf),*disi_,1.0,*dis_,0.0);

  // velocities
  // incremental (required for constant predictor)
  velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                gamma/(beta*dt));

  // accelerations
  // incremental (required for constant predictor)
  accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                (2.*beta-1.+alpham)/(2.*beta),*acc_,
                (1.-alpham)/((1.-alphaf)*beta*dt*dt));

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // StruGenAlpha::NonlinearCG()



/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PTC()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");


  // hard wired ptc parameters
  double ptcdt = 1.0e-02;
  double nc;
  fresm_->NormInf(&nc);
  double dti = 1/ptcdt;
  double dti0 = dti;
  RCP<Epetra_Vector> x0 = rcp(new Epetra_Vector(*disi_));
  double np = nc;

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm;
  fresm_->Norm2(&fresmnorm);
  while (norm_>toldisp && fresmnorm>toldisp && numiter<=maxiter)
  {
    double dtim = dti0;
    dti0 = dti;
    RCP<Epetra_Vector> xm = rcp(new Epetra_Vector(*x0));
    x0->Update(1.0,*disi_,0.0);
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    LINALG::Add(*mass_,false,(1.-alpham)/(beta*dt*dt),*stiff_,1.-alphaf);
    if (damping)
      LINALG::Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),*stiff_,1.0);
    LINALG::Complete(*stiff_);
    
    //------------------------------- do ptc modification to effective LHS
    {
      RCP<Epetra_Vector> tmp = LINALG::CreateVector(stiff_->RowMap(),false);
      tmp->PutScalar(dti);
      RCP<Epetra_Vector> diag = LINALG::CreateVector(stiff_->RowMap(),false);
      stiff_->ExtractDiagonalCopy(*diag);
      diag->Update(1.0,*tmp,1.0);
      stiff_->ReplaceDiagonalValues(*diag);
    }

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_,disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_,disi_,fresm_,true,false);
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
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();
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
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // compute inf norm of residual
    double np;
    fresm_->NormInf(&np);

    //---------------------------------------------- build residual norm
    double disinorm;
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);
    // a short message
    if (!myrank_)
    {
      if (printscreen)
      {
        printf("numiter %d res-norm %e dis-norm %e dti %e\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(stdout);
      }
      if (printerr)
      {
        fprintf(errfile,"numiter %d res-norm %e dis-norm %e dti %e\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(errfile);
      }
    }
    // criteria to stop Newton iteration
    norm_ = disinorm;

    //--------------------------------- increment equilibrium loop index
    ++numiter;

    //------------------------------------ PTC update of artificial time 
#if 1
    // SER step size control
    dti *= (np/nc);
    dti = max(dti,0.0);
    nc = np;
#endif    

#if 0
    {
      double ttau=0.75;	
      RCP<Epetra_Vector> d1 = LINALG::CreateVector(stiff_->RowMap(),false);
      d1->Update(1.0,*disi_,-1.0,*x0,0.0);
      d1->Scale(dti0);
      RCP<Epetra_Vector> d0 = LINALG::CreateVector(stiff_->RowMap(),false);
      d0->Update(1.0,*x0,-1.0,*xm,0.0);
      d0->Scale(dtim);
      double dt0 = 1/dti0;
      double dtm = 1/dtim;
      RCP<Epetra_Vector> xpp = LINALG::CreateVector(stiff_->RowMap(),false);
      xpp->Update(2.0/(dt0+dtm),*d1,-2.0/(dt0+dtm),*d0,0.0);
      RCP<Epetra_Vector> xtt = LINALG::CreateVector(stiff_->RowMap(),false);
      for (int i=0; i<xtt->MyLength(); ++i) (*xtt)[i] = abs((*xpp)[i])/(1.0+abs((*disi_)[i]));
      double ett;
      xtt->MaxValue(&ett);
      ett = ett / (2.*ttau);
      dti = sqrt(ett);
      nc = np;
    }
#endif    
  } // while (norm_>toldisp && fresmnorm>toldisp && numiter<=maxiter)
  //============================================= end equilibrium loop

  //-------------------------------- test whether max iterations was hit
  if (numiter==maxiter) dserror("Ptc unconverged in %d iterations",numiter);
  params_.set<int>("num iterations",numiter);

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // StruGenAlpha::PTC()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::UpdateandOutput()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time          = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  int    istep         = params_.get<int>   ("step"                   ,0);
  int    nstep         = params_.get<int>   ("nstep"                  ,5);
  int    numiter       = params_.get<int>   ("num iterations"         ,-1);

  double alpham        = params_.get<double>("alpha m"                ,0.378);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);

  bool   iodisp        = params_.get<bool>  ("io structural disp"     ,true);
  int    updevrydisp   = params_.get<int>   ("io disp every nstep"    ,10);
  bool   iostress      = params_.get<bool>  ("io structural stress"   ,false);
  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);

  int    writeresevry  = params_.get<int>   ("write restart every"    ,0);

  bool   printscreen   = params_.get<bool>  ("print to screen"        ,true);
  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
  if (!errfile) printerr = false;

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
    discret_.Evaluate(p,null,null,null,null,null);
  }

  //---------------------------------------------- do stress calculation
  if (updevrystress && !istep%updevrystress && iostress)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",false);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
  }

  //------------------------------------------------- write restart step
  bool isdatawritten = false;
  if (writeresevry && istep%writeresevry==0)
  {
    output_.NewStep(istep, time);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    output_.WriteMesh(istep,time);
    isdatawritten = true;

    if (discret_.Comm().MyPID()==0 && printscreen)
    {
      cout << "====== Restart written in step " << istep << endl;
      fflush(stdout);
    }
    if (errfile && printerr)
    {
      fprintf(errfile,"====== Restart written in step %d\n",istep);
      fflush(errfile);
    }
  }

  //----------------------------------------------------- output results
  if (iodisp && updevrydisp && istep%updevrydisp==0 && !isdatawritten)
  {
    output_.NewStep(istep, time);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    isdatawritten = true;
  }

  //---------------------------------------------------------- print out
  if (!myrank_)
  {
    if (printscreen)
    {
      printf("step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
             istep,nstep,time,dt,numiter);
      printf("----------------------------------------------------------------------------------\n");
      fflush(stdout);
    }
    if (printerr)
    {
      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
              istep,nstep,time,dt,numiter);
      fprintf(errfile,"----------------------------------------------------------------------------------\n");
      fflush(errfile);
    }
  }

  return;
} // StruGenAlpha::UpdateandOutput()




/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::Integrate()
{
  int    istep = params_.get<int>   ("step" ,0);
  int    nstep = params_.get<int>   ("nstep",5);

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  if (equil=="full newton")
  {
    for (int i=istep; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      FullNewton();
      UpdateandOutput();
    }
  }
  else if (equil=="modified newton")
  {
    for (int i=istep; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      ModifiedNewton();
      UpdateandOutput();
    }
  }
  else if (equil=="nonlinear cg")
  {
    for (int i=istep; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      NonlinearCG();
      UpdateandOutput();
    }
  }
  else if (equil=="ptc")
  {
    for (int i=istep; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      PTC();
      UpdateandOutput();
    }
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
}



/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::SetDefaults(ParameterList& params)
{
  params.set<bool>  ("print to screen"        ,true);
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
  params.set<bool>  ("io structural disp"     ,false);
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<bool>  ("io structural stress"   ,false);
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<int>   ("restart"                ,0);
  params.set<int>   ("write restart every"    ,0);
  // takes values "constant" consistent"
  params.set<string>("predictor"              ,"constant");
  // takes values "full newton" , "modified newton" , "nonlinear cg"
  params.set<string>("equilibrium iteration"  ,"full newton");
  return;
}

/*----------------------------------------------------------------------*
 |  read restart (public)                                    mwgee 06/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ReadRestart(int step)
{
  RefCountPtr<DRT::Discretization> rcpdiscret = rcp(&discret_);
  rcpdiscret.release();
  IO::DiscretizationReader reader(rcpdiscret,step);
  double time  = reader.ReadDouble("time");
  int    rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  reader.ReadVector(dis_, "displacement");
  reader.ReadVector(vel_, "velocity");
  reader.ReadVector(acc_, "acceleration");
  reader.ReadVector(fext_,"fexternal");
  reader.ReadMesh(step);

  // iverride current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
StruGenAlpha::~StruGenAlpha()
{
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
