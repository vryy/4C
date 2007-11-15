/*!----------------------------------------------------------------------
\file contactstrugenalpha.cpp
\brief Generalized Alpha time integration for structural problems with contact

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "contactstrugenalpha.H"
#include "iostream"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/07|
 *----------------------------------------------------------------------*/
ContactStruGenAlpha::ContactStruGenAlpha(ParameterList& params,
                                         DRT::Discretization& dis,
                                         LINALG::Solver& solver,
                                         IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output,false)
{
  havecontact_ = true; // set base class variable

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time    = params_.get<double>("total time"      ,0.0);
  double dt      = params_.get<double>("delta time"      ,0.01);
  bool damping   = params_.get<bool>  ("damping"         ,false);
  double kdamp   = params_.get<double>("damping factor K",0.0);
  double mdamp   = params_.get<double>("damping factor M",0.0);
  int step       = params_.get<int>   ("step"            ,0);
  bool outerr    = params_.get<bool>  ("print to err"    ,false);
  FILE* errfile  = params_.get<FILE*> ("err file"        ,NULL);
  if (!errfile) outerr = false;

  // -------------------------------------------------------------------
  // see whether we have contact boundary conditions
  // and create contact manager if so
  // -------------------------------------------------------------------
  {
    vector<DRT::Condition*> contactconditions(0);
    discret_.GetCondition("Contact",contactconditions);
    if (contactconditions.size()) havecontact_ = true;
    else dserror("No contact boundary conditions present");
    contactmanager_ = rcp(new CONTACT::Manager(discret_));
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  if (!discret_.Filled()) discret_.FillComplete();
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

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
  // also known as out-of-balance-force
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
  // close mass matrix
  LINALG::Complete(*mass_);
  maxentriesperrow_ = mass_->MaxNumEntries();

  // build damping matrix if desired  
  if (damping)
  {
    LINALG::Complete(*stiff_);
    LINALG::Add(*stiff_,false,kdamp,*damp_,0.0);
    LINALG::Add(*mass_,false,mdamp,*damp_,1.0);
    LINALG::Complete(*damp_);
    stiff_ = null;
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
  step = 0;
  params_.set<int>("step",step);


  return;
} // StruGenAlpha::StruGenAlpha




/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::ConsistentPredictor()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  //int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  double alpham      = params_.get<double>("alpha m"        ,0.378);
  double beta        = params_.get<double>("beta"           ,0.292);
  double gamma       = params_.get<double>("gamma"          ,0.581);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // increment time and step
  double timen = time + dt;  // t_{n+1}
  //int istep = step + 1;  // n+1

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
  veln_->Update((beta-gamma)/beta,*vel_,
                (2.*beta-gamma)*dt/(2.*beta),*acc_,gamma/(beta*dt));
  // predicting accelerations A_{n+1} (accn)
  // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
  //          - 1./(beta*dt) * V_n
  //          + (2.*beta-1.)/(2.*beta) * A_n
  accn_->Update(1.0,*disn_,-1.0,*dis_,0.0);
  accn_->Update(-1./(beta*dt),*vel_,
                (2.*beta-1.)/(2.*beta),*acc_,1./(beta*dt*dt));

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
  if (printscreen)
    fresm_->Norm2(&norm_);
  if (!myrank_ && printscreen)
  {
    cout << "Predictor residual forces " << norm_ << endl;
    fflush(stdout);
  }

  return;
} // ContactStruGenAlpha::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::ConstantPredictor()
{

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  //int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // increment time and step
  double timen = time += dt;
  //int istep = step + 1;

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
  if (printscreen)
    fresm_->Norm2(&norm_);
  if (!myrank_ && printscreen)
  {
    cout << "Predictor residual forces " << norm_ << endl;
    fflush(stdout);
  }

  return;
} // StruGenAlpha::ConstantPredictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::FullNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
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
    {
      LINALG::Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),*stiff_,1.0);
    }
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
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

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
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
     if ( (myrank_ ==  0) && (printscreen) )
     {
        printf("Newton iteration converged: numiter %d res-norm %e dis-norm %e\n", 
               numiter+1, fresmnorm, disinorm);
        fflush(stdout);
     }  
  }
  params_.set<int>("num iterations",numiter);

  //-------------------------------------- don't need this at the moment
  stiff_ = null;

  return;
} // ContactStruGenAlpha::FullNewton()



/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::UpdateandOutput()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time          = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  double timen         = time + dt;  // t_{n+1}
  int    step          = params_.get<int>   ("step"                   ,0);
  int    istep         = step + 1;  // n+1
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

  //----------------------------------------------- update time and step
  params_.set<double>("total time", timen);
  params_.set<int>("step", istep);

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
    p.set("total time",timen);
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
    p.set("total time",timen);
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
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    output_.WriteMesh(istep,timen);
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
    output_.NewStep(istep, timen);
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
             istep,nstep,timen,dt,numiter);
      printf("----------------------------------------------------------------------------------\n");
      fflush(stdout);
    }
    if (printerr)
    {
      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
              istep,nstep,timen,dt,numiter);
      fprintf(errfile,"----------------------------------------------------------------------------------\n");
      fflush(errfile);
    }
  }

  return;
} // ContactStruGenAlpha::UpdateandOutput()




/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::Integrate()
{
  int    step  = params_.get<int>   ("step" ,0);
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
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      FullNewton();
      UpdateandOutput();
    }
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void ContactStruGenAlpha::Integrate()



/*----------------------------------------------------------------------*
 |  read restart (public)                                    mwgee 06/07|
 *----------------------------------------------------------------------*/
void ContactStruGenAlpha::ReadRestart(int step)
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


#endif  // #ifdef CCADISCRET
