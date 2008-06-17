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
#include "contactdefines.H"

#include "../drt_lib/drt_globalproblem.H"

#include <iostream>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/07|
 *----------------------------------------------------------------------*/
CONTACT::ContactStruGenAlpha::ContactStruGenAlpha(ParameterList& params,
                                                  DRT::Discretization& dis,
                                                  LINALG::Solver& solver,
                                                  IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output)
{
  // -------------------------------------------------------------------
  // see whether we have contact boundary conditions
  // and create contact manager if so
  // -------------------------------------------------------------------
  {
    vector<DRT::Condition*> contactconditions(0);
    discret_.GetCondition("Contact",contactconditions);
    if (!contactconditions.size()) dserror("No contact boundary conditions present");
    double alphaf = params_.get<double>("alpha f",0.459);
    contactmanager_ = rcp(new CONTACT::Manager(discret_,alphaf));
  }
  return;
} // ContactStruGenAlpha::ContactStruGenAlpha


/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::ConsistentPredictor()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  double mdamp       = params_.get<double>("damping factor M",0.0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  double alpham      = params_.get<double>("alpha m"        ,0.378);
  double beta        = params_.get<double>("beta"           ,0.292);
#ifdef STRUGENALPHA_BE
  double delta       = params_.get<double>("delta"          ,beta);
#endif
  double gamma       = params_.get<double>("gamma"          ,0.581);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

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
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    p.set("damping factor M",mdamp);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  //cout << *disn_ << endl;

  // consistent predictor
  // predicting velocity V_{n+1} (veln)
  // V_{n+1} := gamma/(beta*dt) * (D_{n+1} - D_n)
  //          + (beta-gamma)/beta * V_n
  //          + (2.*beta-gamma)/(2.*beta) * A_n
  veln_->Update(1.0,*disn_,-1.0,*dis_,0.0);
#ifdef STRUGENALPHA_BE
  veln_->Update((delta-gamma)/delta,*vel_,
                (-gamma-2.*delta*gamma+2.*beta*gamma+2.*delta)*dt/(2.*delta),*acc_,gamma/(delta*dt));
#else
  veln_->Update((beta-gamma)/beta,*vel_,
                (2.*beta-gamma)*dt/(2.*beta),*acc_,gamma/(beta*dt));
#endif


#ifdef STRUGENALPHA_STRONGDBC
  // apply new velocities at DBCs
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    //p.set("time derivative degree",1);  // we want velocities
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("velocity",veln_);
    // predicted dirichlet values
    // veln_ then also holds prescribed new Dirichlet velocities
    discret_.EvaluateDirichlet(p,null,veln_,null,dirichtoggle_);
    discret_.ClearState();
  }
#endif

  //cout << *veln_ << endl;

  // predicting accelerations A_{n+1} (accn)
  // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
  //          - 1./(beta*dt) * V_n
  //          + (2.*beta-1.)/(2.*beta) * A_n
  accn_->Update(1.0,*disn_,-1.0,*dis_,0.0);
#ifdef STRUGENALPHA_BE
  accn_->Update(-1./(delta*dt),*vel_,
                (2.*beta-1.)/(2.*delta),*acc_,1./(delta*dt*dt));
#else
  accn_->Update(-1./(beta*dt),*vel_,
                (2.*beta-1.)/(2.*beta),*acc_,1./(beta*dt*dt));
#endif

#ifdef STRUGENALPHA_STRONGDBC
  // apply new accelerations at DBCs
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    //p.set("time derivative degree",2);  // we want accelerations
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("acceleration",accn_);
    // predicted dirichlet values
    // accn_ then also holds prescribed new Dirichlet accelerations
    discret_.EvaluateDirichlet(p,null,null,accn_,dirichtoggle_);
    discret_.ClearState();
  }
#endif

  //cout << *accn_ << endl;


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
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    // set vector values needed by elements
    discret_.ClearState();
    disi_->PutScalar(0.0);
    discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
    discret_.SetState("displacement",disn_);
#else
    discret_.SetState("displacement",dism_);
#endif
    //discret_.SetState("velocity",velm_); // not used at the moment
#ifdef STRUGENALPHA_FINTLIKETR
    fintn_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
    discret_.ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     + F_c(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  
  // FIXME: Strictly speaking we have to include the contact forces
  // here as well, but it does not matter for the following calculations,
  // so maybe we could just remove it and don't care about the wrong
  // predictor res-norm...?
  
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  
  // add mid-viscous damping force
  if (damping)
  {
    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    damp_->Multiply(false,*velm_,*fvisc_);
    fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
  fresm_->Update(1.0,*fextm_,-1.0);
  fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------------------ contact
  RCP<Epetra_Vector> z = contactmanager_->LagrMult();
  RCP<Epetra_Vector> zold = contactmanager_->LagrMultOld();
  
  // update current LM (predictor = last converged state)
  // this resetting is necessary due to multiple active set steps
  z->Update(1.0,*zold,0.0);
  contactmanager_->StoreNodalQuantities("lmcurrent");
  
  // evaluate Mortar coupling matrices for contact forces
  contactmanager_->Initialize(0);
  contactmanager_->SetState("displacement",disn_);
  contactmanager_->EvaluateMortar();
    
  // add contact forces
  contactmanager_->ContactForces(fresm_);
  RCP<Epetra_Vector> fc = contactmanager_->GetContactForces();
  Epetra_Vector fccopy(*fc);
  fc->Multiply(1.0,*invtoggle_,fccopy,0.0);
  if (fc!=null) fresm_->Update(-1.0,*fc,1.0);
    
  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step (possibly after a
  // restart)
  if (firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
    firststep_=false;
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }

  // remove contact forces from equilibrium again
  if (fc!=null) fresm_->Update(1.0,*fc,1.0);
    
  return;
} // ContactStruGenAlpha::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::ConstantPredictor()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  double mdamp       = params_.get<double>("damping factor M",0.0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  // increment time and step
  double timen = time + dt;
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
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    p.set("damping factor M",mdamp);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
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
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    // set vector values needed by elements
    discret_.ClearState();
    disi_->PutScalar(0.0);
    discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
    discret_.SetState("displacement",disn_);
#else
    discret_.SetState("displacement",dism_);
#endif
    //discret_.SetState("velocity",velm_); // not used at the moment
#ifdef STRUGENALPHA_FINTLIKETR
    fintn_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
    discret_.ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     + F_c(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  
  // FIXME: Strictly speaking we have to include the contact forces
  // here as well, but it does not matter for the following calculations,
  // so maybe we could just remove it and don't care about the wrong
  // predictor res-norm...?
  
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  // add mid-viscous damping force
  if (damping)
  {
    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    damp_->Multiply(false,*velm_,*fvisc_);
    fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
  fresm_->Update(1.0,*fextm_,-1.0);
  fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------------------ contact
  RCP<Epetra_Vector> z = contactmanager_->LagrMult();
  RCP<Epetra_Vector> zold = contactmanager_->LagrMultOld();
  
  // update current LM (predictor = last converged state)
  // this resetting is necessary due to multiple active set steps
  z->Update(1.0,*zold,0.0);
  contactmanager_->StoreNodalQuantities("lmcurrent");
  
  // evaluate Mortar coupling matrices for contact forces
  contactmanager_->Initialize(0);
  contactmanager_->SetState("displacement",disn_);
  contactmanager_->EvaluateMortar();
    
  // add contact forces
  contactmanager_->ContactForces(fresm_);
  RCP<Epetra_Vector> fc = contactmanager_->GetContactForces();
  Epetra_Vector fccopy(*fc);
  fc->Multiply(1.0,*invtoggle_,fccopy,0.0);
  if (fc!=null) fresm_->Update(-1.0,*fc,1.0);
    
  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step (possibly after a
  // restart)
  if (firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
    firststep_=false;
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }

  // remove contact forces from equilibrium again
  if (fc!=null) fresm_->Update(1.0,*fc,1.0);
    
  return;
} // ContactStruGenAlpha::ConstantPredictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::FullNewton()
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
#ifdef STRUGENALPHA_BE
  double delta     = params_.get<double>("delta"                  ,beta);
#endif
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  bool structrobin = params_.get<bool>  ("structrobin"            ,false);
  if (!errfile) printerr = false;
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
#ifdef STRUGENALPHA_BE
    stiff_->Add(*mass_,false,(1.-alpham)/(delta*dt*dt),1.-alphaf);
#else
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
#endif
    if (damping)
    {
#ifdef STRUGENALPHA_BE
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(delta*dt),1.0);
#else
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
#endif
    }
    stiff_->Complete();

    //-------------------------make contact modifications to lhs and rhs
    {
      contactmanager_->Initialize(numiter);
      contactmanager_->SetState("displacement",disn_);

      // (almost) all contact stuff is done here!
      contactmanager_->Evaluate(stiff_,fresm_,numiter);
    }
        
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    //------------------------------------ -- recover disi and Lag. Mult.
    {
      contactmanager_->Recover(disi_);
    }
        
    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(1.0,*disi_,1.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
    disn_->Update(1.0,*disi_,1.0);
    dism_->Update(1.-alphaf,*disi_,1.0);
#endif
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
#ifdef STRUGENALPHA_BE
    velm_->Update((delta-(1.0-alphaf)*gamma)/delta,*vel_,
                  (1.0-alphaf)*(-gamma-2.*delta*gamma+2.*beta*gamma+2.*delta)*dt/(2.*delta),*acc_,
                  gamma/(delta*dt));
#else
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
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
#ifdef STRUGENALPHA_BE
    accm_->Update(-(1.-alpham)/(delta*dt),*vel_,
                  (2.*beta-1.+alpham-2.*alpham*beta+2.*alpham*delta)/(2.*delta),*acc_,
                  (1.-alpham)/((1.-alphaf)*delta*dt*dt));
#else
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif
#endif

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      p.set("alpha f",alphaf);
      // set vector values needed by elements
      discret_.ClearState();
#ifdef STRUGENALPHA_FINTLIKETR
#else
      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      disi_->Scale(1.-alphaf);
#endif
      discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
#else
      discret_.SetState("displacement",dism_);
#endif
      //discret_.SetState("velocity",velm_); // not used at the moment
#ifdef STRUGENALPHA_FINTLIKETR
      fintn_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
      discret_.ClearState();

      if (surf_stress_man_!=null)
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
      }

      if (constrMan_->HaveConstraint())
      {
        constrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
      }

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddCondVector(alphas/dt,
                                   fsisurface_->ExtractCondVector(dism_),
                                   fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->CondMap();
        int numrdofs = robinmap.NumMyElements();
        int* rdofs = robinmap.MyGlobalElements();
        for (int lid=0; lid<numrdofs; ++lid)
        {
          int gid = rdofs[lid];
          // Note: This assemble might fail if we have a block matrix here.
          stiff_->Assemble(scale,gid,gid);
        }
      }
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //      + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // FIXME:
    // Strictly speaking, we would have to evaluate all contact stuff
    // here, to build the corect residuum.
    // this could be done by rearranging the code such that we move
    // the whole contact initialize - set state - evaluate to here!
      
    // add contact forces
    contactmanager_->ContactForces(fresm_);
    RCP<Epetra_Vector> fc = contactmanager_->GetContactForces();
    Epetra_Vector fccopy(*fc);
    fc->Multiply(1.0,*invtoggle_,fccopy,0.0);
    if (fc!=null) fresm_->Update(-1.0,*fc,1.0);
        
    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    //remove contact forces from equilibrium again
    if (fc!=null) fresm_->Update(1.0,*fc,1.0);
        
    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;
    
    //-------------------- update active set for semi-smooth Newton case
#ifdef CONTACTSEMISMOOTH
    contactmanager_->UpdateActiveSetSemiSmooth(disn_);
#endif // #ifdef CONTACTSEMISMOOTH
    
  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
    if (constrMan_->HaveMonitor())
    {
      constrMan_->ComputeMonitorValues(dism_);
    }
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // ContactStruGenAlpha::FullNewton()


/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::UpdateandOutput()
{
  Update();
  Output();
  return;
} // ContactStruGenAlpha::UpdateandOutput()

/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::Update()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time          = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  double timen         = time + dt;  // t_{n+1}
  int    step          = params_.get<int>   ("step"                   ,0);
  int    istep         = step + 1;  // n+1

  double alpham        = params_.get<double>("alpha m"                ,0.378);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);

  string iostress      = params_.get<string>("io structural stress"   ,"none");
  string iostrain      = params_.get<string>("io structural strain"   ,"none");

  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
  if (!errfile) printerr = false;

  //----------------------------------------------- update time and step
  params_.set<double>("total time", timen);
  params_.set<int>("step", istep);

  //---------------------- determine new end-point quantities and update
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
#ifdef STRUGENALPHA_FINTLIKETR
  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0,*fintn_,0.0);
#endif

  //---------------------------------------------- print contact to screen
  contactmanager_->PrintActiveSet();
    
  //---------------------------------- store Lagrange multipliers, D and M
  RCP<Epetra_Vector> z = contactmanager_->LagrMult();
  RCP<Epetra_Vector> zold = contactmanager_->LagrMultOld();
  zold->Update(1.0,*z,0.0);
  contactmanager_->StoreNodalQuantities("lmold");
  contactmanager_->StoreDM("old");

#ifdef PRESTRESS
  //----------- save the current green-lagrange strains in the material
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_prestress_update_green_lagrange");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    discret_.SetState("displacement",dis_);
    discret_.SetState("residual displacement",zeros_);
    discret_.Evaluate(p,null,null,null,null,null);
  }

  //----------------------------- reset the current disp/vel/acc to zero
  // (the structure does not move while prestraining it )
  // (prestraining with DBCs != 0 not allowed!)
  //dis_->Update(1.0,disold,0.0);
  //vel_->Update(1.0,velold,0.0);
  //acc_->Update(1.0,accold,0.0);
  dis_->Scale(0.0);
  vel_->Scale(0.0);
  acc_->Scale(0.0);
#endif


  //------ update anything that needs to be updated at the element level
#ifdef STRUGENALPHA_FINTLIKETR
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_update_istep");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    discret_.Evaluate(p,null,null,null,null,null);
  }
#else
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    //p.set("action","calc_struct_update_istep");
    p.set("action","calc_struct_update_genalpha_imrlike");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    discret_.Evaluate(p,null,null,null,null,null);
  }
#endif

  //----------------- update surface stress history variables if present
  if (surf_stress_man_!=null)
  {
    surf_stress_man_->Update();
  }
} // ContactStruGenAlpha::Update()


/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::Output()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double timen         = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);
  int    istep         = params_.get<int>   ("step"                   ,0);
  int    nstep         = params_.get<int>   ("nstep"                  ,5);
  int    numiter       = params_.get<int>   ("num iterations"         ,-1);

  bool   iodisp        = params_.get<bool>  ("io structural disp"     ,true);
  int    updevrydisp   = params_.get<int>   ("io disp every nstep"    ,10);
  string iostress      = params_.get<string>("io structural stress"   ,"none");
  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);
  string iostrain      = params_.get<string>("io structural strain"   ,"none");

  int    writeresevry  = params_.get<int>   ("write restart every"    ,0);

  bool   printscreen   = params_.get<bool>  ("print to screen"        ,true);
  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
  if (!errfile) printerr = false;

  bool isdatawritten = false;
 
  //------------------------------------------------- write restart step
  if (writeresevry and istep%writeresevry==0)
  {
    output_.WriteMesh(istep,timen);
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    isdatawritten = true;
    
    // write restart information for contact
    RCP<Epetra_Vector> zold = contactmanager_->LagrMultOld();
    RCP<Epetra_Vector> activetoggle = contactmanager_->WriteRestart();
    output_.WriteVector("lagrmultold",zold);
    output_.WriteVector("activetoggle",activetoggle);

    if (surf_stress_man_!=null)
    {
      RCP<Epetra_Map> surfrowmap=surf_stress_man_->GetSurfRowmap();
      RCP<Epetra_Vector> A=rcp(new Epetra_Vector(*surfrowmap, true));
      RCP<Epetra_Vector> con=rcp(new Epetra_Vector(*surfrowmap, true));
      surf_stress_man_->GetHistory(A,con);
      output_.WriteVector("Aold", A);
      output_.WriteVector("conquot", con);
    }

    if (constrMan_->HaveConstraint())
    {
      output_.WriteDouble("uzawaparameter",uzawaSolv_->GetUzawaParameter());
    }

    if (discret_.Comm().MyPID()==0 and printscreen)
    {
      cout << "====== Restart written in step " << istep << endl;
      fflush(stdout);
    }
    if (errfile and printerr)
    {
      fprintf(errfile,"====== Restart written in step %d\n",istep);
      fflush(errfile);
    }
  }

  //----------------------------------------------------- output results
  if (iodisp and updevrydisp and istep%updevrydisp==0 and !isdatawritten)
  {
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    output_.WriteElementData();
    isdatawritten = true;
  }

  //------------------------------------- do stress calculation and output
  if (updevrystress and !(istep%updevrystress) and iostress!="none")
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
    Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stress);
    p.set("strain", strain);
    if (iostress == "cauchy")   // output of Cauchy stresses instead of 2PK stresses
    {
      p.set("cauchy", true);
    }
    else
    {
      p.set("cauchy", false);
    }
    p.set("iostrain", iostrain);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
    if (!isdatawritten) output_.NewStep(istep, timen);
    isdatawritten = true;
    if (iostress == "cauchy")
    {
      output_.WriteVector("gauss_cauchy_stresses_xyz",*stress,*discret_.ElementColMap());
    }
    else
    {
      output_.WriteVector("gauss_2PK_stresses_xyz",*stress,*discret_.ElementColMap());
    }
    if (iostrain != "none")
    {
      if (iostrain == "euler_almansi")
      {
        output_.WriteVector("gauss_EA_strains_xyz",*strain,*discret_.ElementColMap());
      }
      else
      {
        output_.WriteVector("gauss_GL_strains_xyz",*strain,*discret_.ElementColMap());
      }
    }
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
} // ContactStruGenAlpha::Output()


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               popp  02/08|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::Integrate()
{
  int    step    = params_.get<int>   ("step" ,0);
  int    nstep   = params_.get<int>   ("nstep",5);
  double maxtime = params_.get<double>("max time",0.0);

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  //in case a constraint is defined, use defined algorithm
  if (constrMan_->HaveConstraint())
  {
    string algo = params_.get<string>("uzawa algorithm","newtonlinuzawa");
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      //Does predicted displacement satisfy constraint?
      double time = params_.get<double>("total time",0.0);
      double dt   = params_.get<double>("delta time",0.01);
      // what algorithm is used?
      // - "newtonlinuzawa":      Potential is linearized wrt displacements and Lagrange multipliers
      //                  Linear problem is solved with Uzawa algorithm
      // - "augmentedlagrange":   Potential is linearized wrt displacements keeping Lagrange multiplier fixed
      //                Until convergence Lagrange multiplier increased by Uzawa_param*(Vol_err)
      if (algo=="newtonlinuzawa")
      {
        FullNewtonLinearUzawa();
      }
      else if (algo=="augmentedlagrange")
      {
         NonLinearUzawaFullNewton(predictor);
      }
      else dserror("Unknown type of algorithm to deal with constraints");
      UpdateandOutput();
      if (time>=maxtime) break;
    }
  }
  
  // Newton as nonlinear iteration scheme
  else if (equil=="full newton")
  {
    //********************************************************************
    // OPTIONS FOR PRIMAL-DUAL ACTIVE SET STRATEGY (PDASS)
    //********************************************************************
    // SEMI-SMOOTH NEWTON
    // The search for the correct active set (=contact nonlinearity) and
    // the large deformstion linearization (=geometrical nonlinearity) are
    // merged into one semi-smooth Newton method and solved within ONE
    // iteration loop
    //********************************************************************
    // FIXED-POINT APPROACH
    // The search for the correct active set (=contact nonlinearity) is
    // represented by a fixed-point approach, whereas the large deformation
    // linearization (=geimetrical nonlinearity) is treated by a standard
    // Newton scheme. This yields TWO nested iteration loops
    //********************************************************************
#ifdef CONTACTSEMISMOOTH
    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      contactmanager_->ActiveSetConverged() = false;
      contactmanager_->ActiveSetSteps() = 1;

      // predictor step
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();

      // LOOP2: nonlinear iteration (Newton)
      FullNewton();
      
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
    
#else
    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      contactmanager_->ActiveSetConverged() = false;
      contactmanager_->ActiveSetSteps() = 1;

      // LOOP2: active set strategy
      while (contactmanager_->ActiveSetConverged()==false)
      {
        // predictor step
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();

        // LOOP3: nonlinear iteration (Newton)
        FullNewton();

        // update of active set
        contactmanager_->UpdateActiveSet(disn_);
      }
      
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
#endif // #ifdef CONTACTSEMISMOOTH
  }
  
  // other types of nonlinear iteration schemes
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void ContactStruGenAlpha::Integrate()


/*----------------------------------------------------------------------*
 |  read restart (public)                                    mwgee 06/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::ReadRestart(int step)
{
  //FIXME
  // currently restart with contact only works for the evaluation of
  // fint IMR-like at the new mid-point. For a TR-like evaluation in
  // the new predictor, one would have to store fint for restart!
#ifdef STRUGENALPHA_FINTLIKETR
  dserror("ERROR: ReadRestart: Not yet implemented for FINTLIKETR!");
#endif // #ifdef STRUGENALPHA_FINTLIKETR
  
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

  // read restart information for contact
  RCP<Epetra_Vector> zold = rcp(new Epetra_Vector(*(contactmanager_->SlaveRowDofs())));
  RCP<Epetra_Vector> activetoggle =rcp(new Epetra_Vector(*(contactmanager_->SlaveRowNodes())));
  reader.ReadVector(zold,"lagrmultold");
  reader.ReadVector(activetoggle,"activetoggle");
  *(contactmanager_->LagrMultOld())=*zold;
  contactmanager_->StoreNodalQuantities("lmold");
  contactmanager_->ReadRestart(activetoggle);
  
  // build restart Mortar matrices D and M
  contactmanager_->Initialize(0);
  contactmanager_->SetState("displacement",dis_);
  contactmanager_->EvaluateMortar();
  contactmanager_->StoreDM("old");
  
  // override current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  if (surf_stress_man_!=null)
  {
    RCP<Epetra_Map> surfmap=surf_stress_man_->GetSurfRowmap();
    RCP<Epetra_Vector> A_old = LINALG::CreateVector(*surfmap,true);
    RCP<Epetra_Vector> con_quot = LINALG::CreateVector(*surfmap,true);
    reader.ReadVector(A_old, "Aold");
    reader.ReadVector(con_quot, "conquot");
    surf_stress_man_->SetHistory(A_old,con_quot);
  }

  if (DRT::Problem::Instance()->ProblemType()=="struct_multi")
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","multi_readrestart");
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
  }

  if (constrMan_->HaveConstraint())
  {
    double uzawatemp = reader.ReadDouble("uzawaparameter");
    uzawaSolv_->SetUzawaParameter(uzawatemp);
  }

  return;
} // void ContactStruGenAlpha::ReadRestart()


#endif  // #ifdef CCADISCRET
