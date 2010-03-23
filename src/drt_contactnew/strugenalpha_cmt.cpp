/*!----------------------------------------------------------------------
\file strugenalpha_cmt.cpp
\brief Gen Alpha time integration for structural problems with meshtying or contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "strugenalpha_cmt.H"
#include "contact_defines.H"
#include "contact_manager.H"
#include "meshtying_defines.H"
#include "meshtying_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"

#include <iostream>

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/07|
 *----------------------------------------------------------------------*/
CONTACT::CmtStruGenAlpha::CmtStruGenAlpha(ParameterList& params,
                                          DRT::Discretization& dis,
                                          LINALG::Solver& solver,
                                          IO::DiscretizationWriter& output,
                                          INPAR::CONTACT::ApplicationType apptype) :
StruGenAlpha(params,dis,solver,output)
{
  //**********************************************************************
  // see whether we have contact boundary conditions
  // and create meshtying OR contact manager if so
  //**********************************************************************
  {
    // check contact conditions
    vector<DRT::Condition*> contactconditions(0);
    discret_.GetCondition("Contact",contactconditions);
    if (!contactconditions.size()) dserror("No contact boundary conditions present");

    // store integration parameter alphaf into cmanager as well
    double alphaf = params_.get<double>("alpha f",0.459);
    
    // decide whether this is meshtying or contact
    if (apptype==INPAR::CONTACT::app_mortarmeshtying)
      cmtmanager_ = rcp(new CONTACT::MtManager(discret_,alphaf));
    else if (apptype==INPAR::CONTACT::app_mortarcontact)
      cmtmanager_ = rcp(new CONTACT::CoManager(discret_,alphaf));
    else
      dserror("ERROR: You should not be here!");
  }

  //**********************************************************************
  // save Dirichlet B.C. status in Contact Manager
  // all nodes on all interfaces then know if D.B.C.s are applied
  //**********************************************************************
  {
    // map containing Dirichlet DOFs
    Teuchos::ParameterList p;
    double time = params_.get<double>("total time",0.0);
    p.set("total time", time);
    RCP<LINALG::MapExtractor> dbcmaps = rcp(new LINALG::MapExtractor());
    discret_.EvaluateDirichlet(p,zeros_,Teuchos::null,Teuchos::null,Teuchos::null,dbcmaps);
    zeros_->PutScalar(0.0); // just in case
    cmtmanager_->GetStrategy().StoreDirichletStatus(dbcmaps);
  }
  
  //********************************************************************
  // print warning messages for multifield problems (e.g FSI)
  //********************************************************************
  const string probtype = DRT::Problem::Instance()->ProblemType(); 
  if (probtype!="structure" && discret_.Comm().MyPID() == 0)
  {
#ifdef CONTACTPSEUDO2D
    cout << RED << "WARNING: The flag CONTACTPSEUDO2D is switched on. If this "
         << "is a real 3D problem, switch it off!" << END_COLOR << endl;
#else
    cout << RED << "WARNING: The flag CONTACTPSEUDO2D is switched off. If this "
         << "is a 2D problem modeled pseudo-3D, switch it on!" << END_COLOR << endl;
#endif // #ifdef CONTACTPSEUDO2D
    cout << RED << "WARNING: Contact and Meshtying are still experimental "
         << "for the chosen problem type \"" << probtype << "\"!\n" << END_COLOR << endl;
  }
    
  //**********************************************************************
  // visualization of initial configuration
  //**********************************************************************
#ifdef CONTACTGMSH3
  cmtmanager_->GetStrategy().VisualizeGmsh(0,0);
#endif // #ifdef CONTACTGMSH2
  
  //**********************************************************************
  // initialization of contact or meshting
  //**********************************************************************
  {
    // FOR MESHTYING (ONLY ONCE), NO FUNCTIONALITY FOR CONTACT CASES
    // (1) Do mortar coupling in reference configuration
    //     (for restart this has already been done in Strategy object)
    // (2) Perform mesh intialization by relocating slave nodes such
    //     that rotational invariance is satisfied.
    if (params_.get<int>("step") == 0)
    {
     cmtmanager_->GetStrategy().MortarCoupling(zeros_);
     cmtmanager_->GetStrategy().MeshInitialization();
    }
    
    // FOR PENALTY CONTACT (ONLY ONCE), NO FUNCTIONALITY FOR OTHER CASES
    // (1) Explicitly store gap-scaling factor kappa
    //     (for restart this need not be done)
    if (params_.get<int>("step") == 0)
     cmtmanager_->GetStrategy().SaveReferenceState(disn_);
  }

  //**********************************************************************
  // output of strategy type to screen
  //**********************************************************************
  {
    // strategy type
    INPAR::CONTACT::SolvingStrategy soltype =
    Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtmanager_->GetStrategy().Params(),"STRATEGY");
    
    // shape function type
    INPAR::MORTAR::ShapeFcn shapefcn =
    Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(cmtmanager_->GetStrategy().Params(),"SHAPEFCN");
    
    // output
    if (discret_.Comm().MyPID() == 0)
    {
      cout << "===== DRT_NEWCONTACT ===========================================" << endl;  
      if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_standard)
        cout << "===== Standard Lagrange multiplier strategy ====================\n" << endl;
      else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        cout << "===== Dual Lagrange multiplier strategy ========================\n" << endl;
      else if (soltype == INPAR::CONTACT::solution_penalty)
        cout << "===== Penalty strategy =========================================\n" << endl;
      else if (soltype == INPAR::CONTACT::solution_auglag)
        cout << "===== Augmented Lagrange strategy ==============================\n" << endl;
    }
  }
  
  return;
} // CmtStruGenAlpha::CmtStruGenAlpha


/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::ConsistentPredictor()
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
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }
  
  // FRICTIONAL CONTACT - evaluate reference state 
  // for frictional contact we need history values (relative velocity) and
  // therefore we store the nodal entries of mortar matrices (reference
  // configuration) before the first time step
  if (params_.get<int>("step") == 0)
    cmtmanager_->GetStrategy().EvaluateReferenceState(disn_);

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
    // in the case of local systems we have to rotate forth and back
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(disn_);
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disn_);
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

  // zerofy velocity and acceleration in case of statics
  if (dynkindstat)
  {
    velm_->PutScalar(0.0);
    accm_->PutScalar(0.0);
    veln_->PutScalar(0.0);
    accn_->PutScalar(0.0);
    vel_->PutScalar(0.0);
    acc_->PutScalar(0.0);
  }

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

    // include potential conditions in fint_ and stiff_
    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute dynamic equilibrium
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     + F_c(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}

  // Please note that due to the contact modifications to the l.h.s. and to
  // the r.h.s below, this dynamic equilibrium duoes NOT play the role of
  // the residual here, as it does in structural dynamics without contact.
  // Of course it is part of the contact residual, but the normal and
  // tangential contact conditions have to be considered as well.

  // build residual
  if (dynkindstat)
  {
    // static residual
    // Res = F_int - F_ext
    fresm_->PutScalar(0.0);
  }
  else
  {
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
  }

  // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
  fresm_->Update(1.0,*fextm_,-1.0);
  fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  // (still without contact, this is just Gen-alpha stuff here)
  if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
  else
  {
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
  }

  // keep a copy of fresm for contact forces / equilibrium check
  RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

  //------------------------- make contact modifications to lhs and rhs
  cmtmanager_->GetStrategy().SetState("displacement",disn_);

  cmtmanager_->GetStrategy().InitEvalInterface();
  cmtmanager_->GetStrategy().InitEvalMortar();
 
  // friction
  // here the relative movement of the contact bodies is evaluated
  // therefore the current configuration and the according mortar
  // matrices are needed
  // it is only evaluated (resetted) for penalty strategy and
  // augmented lagrange strategy
  INPAR::CONTACT::SolvingStrategy strattype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtmanager_->GetStrategy().Params(),"STRATEGY");

  if (strattype == INPAR::CONTACT::solution_penalty
      or strattype == INPAR::CONTACT::solution_auglag)
    cmtmanager_->GetStrategy().EvaluateRelMov();

  cmtmanager_->GetStrategy().Initialize();
  cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

  //---------------------------------------------------- contact forces
  cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);
  
#ifdef CONTACTGMSH2
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  cmtmanager_->GetStrategy().VisualizeGmsh(istep,0);
#endif // #ifdef CONTACTGMSH2

  //------------------------------------ ----complete stiffness matrix
  stiff_->Complete();
      
  // blank residual DOFs that are on Dirichlet BC
  // in the case of local systems we have to rotate forth and back
  {
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
    Epetra_Vector fresmdbc(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
  }

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

  return;
} // CmtStruGenAlpha::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                       popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::ConstantPredictor()
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
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }
  
  // FRICTIONAL CONTACT - evaluate reference state 
  // for frictional contact we need history values (relative velocity) and
  // therefore we store the nodal entries of mortar matrices (reference
  // configuration) before the first time step
  if (params_.get<int>("step") == 0)
    cmtmanager_->GetStrategy().EvaluateReferenceState(disn_);
  
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
    // in the case of local systems we have to rotate forth and back
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(disn_);
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disn_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // zerofy velocity and acceleration in case of statics
  if (dynkindstat)
  {
    vel_->PutScalar(0.0);
    acc_->PutScalar(0.0);
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
  if (dynkindstat)
  {
    fext_->PutScalar(0.0);
  }
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

    // include potential conditions in fint_ and stiff_
    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute dynamic equilibrium
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     + F_c(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}

  // Please note that due to the contact modifications to the l.h.s. and to
  // the r.h.s below, this dynamic equilibrium duoes NOT play the role of
  // the residual here, as it does in structural dynamics without contact.
  // Of course it is part of the contact residual, but the normal and
  // tangential contact conditions have to be considered as well.

  // build residual
  if (dynkindstat)
  {
    // static residual
    // Res = F_int - F_ext
     fresm_->PutScalar(0.0);
  }
  else
  {
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
  }
  // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
  fresm_->Update(1.0,*fextm_,-1.0);
  fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  // (still without contact, this is just Gen-alpha stuff here)
  if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
  else
  {
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
  }

  // keep a copy of fresm for contact forces / equilibrium check
  RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

  //-------------------------- make contact modifications to lhs and rhs
  cmtmanager_->GetStrategy().SetState("displacement",disn_);
  cmtmanager_->GetStrategy().InitEvalInterface();
  cmtmanager_->GetStrategy().InitEvalMortar();
  
  // friction
  // here the relative movement of the contact bodies is evaluated
  // therefore the current configuration and the according mortar
  // matrices are needed
  // it is only evaluated (resetted) for penalty strategy and
  // augmented lagrange strategy
  INPAR::CONTACT::SolvingStrategy strattype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtmanager_->GetStrategy().Params(),"STRATEGY");

  if (strattype == INPAR::CONTACT::solution_penalty
      or strattype == INPAR::CONTACT::solution_auglag)
        cmtmanager_->GetStrategy().EvaluateRelMov();
  
  cmtmanager_->GetStrategy().Initialize();
  cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

  //---------------------------------------------------- contact forces
  cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  cmtmanager_->GetStrategy().VisualizeGmsh(istep,0);
#endif // #ifdef CONTACTGMSH2

  //------------------------------------ ----complete stiffness matrix
  stiff_->Complete();
    
  // blank residual DOFs that are on Dirichlet BC
  // in the case of local systems we have to rotate forth and back
  {
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
    Epetra_Vector fresmdbc(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
  }

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

  return;
} // CmtStruGenAlpha::ConstantPredictor()



/*----------------------------------------------------------------------*
 |  setup equilibrium with additional external forces        u.may 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::ApplyExternalForce(const STR::UTILS::MapExtractor& extractor,
                                                  Teuchos::RCP<Epetra_Vector> iforce)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  double mdamp       = params_.get<double>("damping factor M",0.0);
  int    step        = params_.get<int>   ("step"           ,0);
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
  bool structrobin   = params_.get<bool>  ("structrobin", false);
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // increment time and step
  double timen = time + dt;

  // get new external force vector
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
    // in the case of local systems we have to rotate forth and back
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(disn_);
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disn_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    discret_.SetState("velocity",veln_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // Add iforce to fextn_
  // there might already be (body) forces at the interface nodes
  extractor.AddFSICondVector(iforce,fextn_);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  // add frobin_ directly to fextm_, so we do not change fextn_ which
  // is going to be used later
  if (structrobin)
    fextm_->Update(1.,*frobin_,1.);


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
    discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
    discret_.SetState("displacement",disn_);
    discret_.SetState("velocity",veln_);
#else
    discret_.SetState("displacement",dism_);
    discret_.SetState("velocity",velm_);
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

    if (surf_stress_man_->HaveSurfStress())  dserror("No surface stresses in 'ApplyExternalForce'");

    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later

    // If we have a robin condition we need to modify both the rhs and the
    // matrix diagonal corresponding to the dofs at the robin interface.
    if (structrobin)
    {
      double alphas = params_.get<double>("alpha s",-1.);

      // Add structural part of Robin force
      fsisurface_->AddFSICondVector(alphas/dt,
                                    fsisurface_->ExtractFSICondVector(dism_),
                                    fint_);

      double scale = alphas*(1.-alphaf)/dt;
      const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

  //-------------------------------------------- compute residual forces
  if (dynkindstat)
  {
    // static residual
    // Res = F_int - F_ext
    fresm_->PutScalar(0.0);
  }
  else
  {
    // dynamic residual
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
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
  }

  // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
  fresm_->Update(1.0,*fextm_,-1.0);
  fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  // (still without contact, this is just Gen-alpha stuff here)
  if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
  else
  {
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
  }
  stiff_->Complete();

  // keep a copy of fresm for contact forces / equilibrium check
  RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

  //-------------------------- make contact modifications to lhs and rhs
  cmtmanager_->GetStrategy().SetState("displacement",disn_);
  cmtmanager_->GetStrategy().InitEvalInterface();
  cmtmanager_->GetStrategy().InitEvalMortar();
  cmtmanager_->GetStrategy().Initialize();
  cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

  //---------------------------------------------------- contact forces
  cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
  int istep = step + 1;
  cmtmanager_->GetStrategy().VisualizeGmsh(istep,0);
#endif // #ifdef CONTACTGMSH2

  // blank residual DOFs that are on Dirichlet BC
  // in the case of local systems we have to rotate forth and back
  {
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
    Epetra_Vector fresmdbc(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
  }


  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (step == 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }
}


/*----------------------------------------------------------------------*
 |  compute residual to given state x (public)               u.may 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::computeF(const Epetra_Vector& x, Epetra_Vector& F)
{
  dserror("not yet impemented for cmtstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 |  compute Jacobian at x (public)                           u.may 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::computeJacobian(const Epetra_Vector& x)
{
  dserror("not yet impemented for cmtstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 |  build linear system matrix and rhs (public)              u.kue 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  dserror("not yet impemented for cmtstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 | linear relaxation solve (public)                          u.kue 03/07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::CmtStruGenAlpha::LinearRelaxationSolve(
      Teuchos::RCP<Epetra_Vector> relax)
{
  dserror("not yet impemented for cmtstrugenalpha");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                              popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::FullNewton()
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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // check whether we have a stiffness matrix that is filled
  // and whether mass and damping are present
  // (here you can note the procedural change compared to standard
  // Gen-alpha, where stiff_ must NOT be filled at this point)
  if (!stiff_->Filled()) dserror("stiffness must be filled here");
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
#ifdef CONTACTEIG
  static int globindex = 0;
  ++globindex;
#endif // #ifdef CONTACTEIG

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<maxiter)
  {
    //-----------------------------transform to local coordinate systems
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(SystemMatrix(),fresm_);

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

#ifdef CONTACTEIG
    // print to file in matlab format
    std::ostringstream filename;
    const std::string filebase = "sparsematrix";
    filename << "o/matlab_output/" << filebase << "_" << globindex << "_" << numiter+1 << ".mtl";
    LINALG::PrintMatrixInMatlabFormat(filename.str().c_str(),*(SystemMatrix()->EpetraMatrix()));

    // print sparsity pattern to file
    LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
#endif // #ifdef CONTACTEIG

    //---------------------------------------------- solve linear system
    LinearSolve(numiter,tolres,fresmnorm);
    
    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disi_);

    //--------------------------------------- recover disi and Lag. Mult.
    cmtmanager_->GetStrategy().Recover(disi_);

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

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }

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

      // include potential conditions in fint_ and stiff_
      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
      }

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddFSICondVector(alphas/dt,
                                      fsisurface_->ExtractFSICondVector(dism_),
                                      fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

    //------------------------------------------ compute dynamic equilibrium
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}

    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
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
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    // (again without contact, this is just Gen-alpha stuff here)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
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
    }

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

    //-------------------------make contact modifications to lhs and rhs
    cmtmanager_->GetStrategy().SetState("displacement",disn_);
    cmtmanager_->GetStrategy().InitEvalInterface();
    cmtmanager_->GetStrategy().InitEvalMortar();
    
    // friction
    // here the relative movement of the contact bodies is evaluated
    // therefore the current configuration and the according mortar
    // matrices are needed
    cmtmanager_->GetStrategy().EvaluateRelMov();
 
    cmtmanager_->GetStrategy().Initialize();
    cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

    //--------------------------------------------------- contact forces
    cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
    dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

    //------------------------------------ ----complete stiffness matrix
    stiff_->Complete();
      
    // blank residual DOFs that are on Dirichlet BC
    // in the case of local systems we have to rotate forth and back
    {
      if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;


  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //------------------------------------------------- linear static case
  int nstep = params_.get<int>("nstep",5);
  if (dynkindstat && maxiter==1 && nstep==1)
  {
    if (!myrank_ and printscreen)
    {
      cout << "computed 1 step with 1 iteration: STATIC LINEAR SOLUTION\n";
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }
  //-------------------------------- test whether max iterations was hit
  else if (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter==maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  //--------------------------------------------------- Newton converged
  else
  {
    if (!myrank_ && printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // CmtStruGenAlpha::FullNewton()

/*----------------------------------------------------------------------*
 |  do Newton iteration with line search (public)            mwgee 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::FullNewtonLineSearch()
{
  if (locsysmanager_ != null) dserror("Locsys not yet implemented for LS-Newton!");
  const int myrank = discret_.Comm().MyPID();

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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // check whether we have a stiffness matrix that is filled
  // and whether mass and damping are present
  // (here you can note the procedural change compared to Gen-alpha
  // without contact, where stiff_ must NOT be filled at this point)
  if (!stiff_->Filled()) dserror("stiffness must be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  // storage of old values to be able to repeat a step
  RCP<Epetra_Vector> disno = rcp(new Epetra_Vector(disn_->Map(),false));
  RCP<Epetra_Vector> dismo = rcp(new Epetra_Vector(dism_->Map(),false));
  RCP<Epetra_Vector> velmo = rcp(new Epetra_Vector(velm_->Map(),false));
  RCP<Epetra_Vector> accmo = rcp(new Epetra_Vector(accm_->Map(),false));

  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  double       f0      = fresmnorm;
  const double sigma0  = 0.1;        // lower bound of lambda
  const double sigma1  = 0.5;        // upper bound of lambda
  const double lsalpha = 1.0e-04;    // decrease requirement to be accept a step
  const int    maxarm  = 10;         // max number line search steps

  if (!myrank) printf("Initial    residual      %15.5e\n",fresmnorm);

  //=================================================== equilibrium loop
  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<maxiter)
  {
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //---------------------------------------------- solve linear system
    LinearSolve(numiter,tolres,fresmnorm);
    
    //------------------------------------ -- recover disi and Lag. Mult.
    cmtmanager_->GetStrategy().Recover(disi_);

    //----------------------------------------- store the current values
    disno->Update(1.0,*disn_,0.0);
    dismo->Update(1.0,*dism_,0.0);
    velmo->Update(1.0,*velm_,0.0);
    accmo->Update(1.0,*accm_,0.0);

    //--------------------------------------------- create the step size
    //                                            trying full step first
    double lambda = 1.0;
    double lamm   = 1.0;
    double lamc   = lambda;
    double iarm   = 0;

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(lambda,*disi_,1.0,*disno,0.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
    disn_->Update(lambda,*disi_,1.0,*disno,0.0);
    dism_->Update(lambda*(1.-alphaf),*disi_,1.0,*dismo,0.0);
#endif
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0,*velmo,0.0);
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
    accm_->Update(lambda*(1.-alpham)/(beta*dt*dt),*disi_,1.0,*accmo,0.0);
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

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }

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
      RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
      disim->Scale(lambda);
#else
      // do not touch disi_ here!
      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
      disim->Scale(lambda*(1.-alphaf));
#endif
      discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
#else
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);
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

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddFSICondVector(alphas/dt,
                                      fsisurface_->ExtractFSICondVector(dism_),
                                      fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

    //------------------------------------------ compute dynamic equilibrium
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}

    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
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
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    // (again without contact, this is just Gen-alpha stuff here)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
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
    }

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

    //-------------------------make contact modifications to lhs and rhs
    cmtmanager_->GetStrategy().SetState("displacement",disn_);
    cmtmanager_->GetStrategy().InitEvalInterface();
    cmtmanager_->GetStrategy().InitEvalMortar();
    cmtmanager_->GetStrategy().Initialize();
    cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

    //--------------------------------------------------- contact forces
    cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
    dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

    //------------------------------------ ----complete stiffness matrix
    stiff_->Complete();
      
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);

    double nft = fresmnorm;
    double nf0 = f0;
    double ff0 = nf0*nf0;
    double ffc = nft*nft;
    double ffm = nft*nft;

    //--------------------------------------------------------------------
    // line searching if the step is bad
    //--------------------------------------------------------------------
    bool dolinesearch = true;
    if (dynkindstat && numiter==0) dolinesearch = false;
    while(nft >= (1.0-lsalpha*lambda)*nf0 && dolinesearch)
    {
      //---------------------------- print the rejected step residual norm
      if (!myrank) printf("Bad step                 %15.5e",nft);

      // ---------------------------compute lambda (step length reduction)
      if (!iarm)
        lambda = sigma1*lambda;
      else
      {
        double c2 = lamm*(ffc-ff0)-lamc*(ffm-ff0);
        if (c2>=0.0)
        {
          lambda = sigma1*lamc;
        }
        else
        {
          double c1 = lamc*lamc*(ffm-ff0)-lamm*lamm*(ffc-ff0);
          lambda = -c1*0.5/c2;
          if      (lambda < sigma0*lamc) lambda = sigma0*lamc;
          else if (lambda > sigma1*lamc) lambda = sigma1*lamc;
        }
      }
      //----------------------------- keep track of old lambda values
      lamm = lamc;
      lamc = lambda;

      //---------------------------------- update mid configuration values
      // displacements
      // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
      disn_->Update(lambda,*disi_,1.0,*disno,0.0);
      dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
      disn_->Update(lambda,*disi_,1.0,*disno,0.0);
      dism_->Update(lambda*(1.-alphaf),*disi_,1.0,*dismo,0.0);
#endif
      // velocities
#ifndef STRUGENALPHA_INCRUPDT
      // iterative
      // V_{n+1-alpha_f} := V_{n+1-alpha_f}
      //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
      velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0,*velmo,0.0);
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
      accm_->Update(lambda*(1.-alpham)/(beta*dt*dt),*disi_,1.0,*accmo,0.0);
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

      // zerofy velocity and acceleration in case of statics
      if (dynkindstat)
      {
        velm_->PutScalar(0.0);
        accm_->PutScalar(0.0);
      }

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
        RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
        disim->Scale(lambda);
#else
        // do not touch disi_ here!
        // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
        RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
        disim->Scale(lambda*(1.-alphaf));
#endif
        discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
        discret_.SetState("displacement",disn_);
        discret_.SetState("velocity",veln_);
#else
        discret_.SetState("displacement",dism_);
        discret_.SetState("velocity",velm_);
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

        // do NOT finalize the stiffness matrix to add masses to it later

        // If we have a robin condition we need to modify both the rhs and the
        // matrix diagonal corresponding to the dofs at the robin interface.
        if (structrobin)
        {
          double alphas = params_.get<double>("alpha s",-1.);

          // Add structural part of Robin force
          fsisurface_->AddFSICondVector(alphas/dt,
                                        fsisurface_->ExtractFSICondVector(dism_),
                                        fint_);

          double scale = alphas*(1.-alphaf)/dt;
          const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

      //------------------------------------------ compute dynamic equilibrium
      // Res = M . A_{n+1-alpha_m}
      //     + C . V_{n+1-alpha_f}
      //     + F_int(D_{n+1-alpha_f})
      //     + F_c(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}

      if (dynkindstat)
      {
        // static residual
        // Res = F_int - F_ext
        fresm_->PutScalar(0.0);
      }
      else
      {
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
      }
      // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
      fresm_->Update(1.0,*fextm_,-1.0);
      fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
      fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

      //---------------------------------------------- build effective lhs
      // (using matrix stiff_ as effective matrix)
      // (again without contact, this is just Gen-alpha stuff here)
      if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
      else
      {
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
      }
      
      // keep a copy of fresm for contact forces / equilibrium check
      RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

      //-------------------------make contact modifications to lhs and rhs
      cmtmanager_->GetStrategy().SetState("displacement",disn_);
      cmtmanager_->GetStrategy().InitEvalInterface();
      cmtmanager_->GetStrategy().InitEvalMortar();
      cmtmanager_->GetStrategy().Initialize();
      cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

      //--------------------------------------------------- contact forces
      cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
      dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

      //------------------------------------ ----complete stiffness matrix
      stiff_->Complete();
          
      // blank residual DOFs that are on Dirichlet BC
      {
        Epetra_Vector fresmdbc(*fresm_);
        fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      }

      //---------------------------------------------- build residual norm
      disi_->Norm2(&disinorm);
      fresm_->Norm2(&fresmnorm);

      //---------------------------------------------- keep track of norms
      if (!myrank) printf(" now: %10.5e stepsize %10.5e\n",fresmnorm,lambda); fflush(stdout);
      nft = fresmnorm;
      ffm = ffc;
      ffc = nft*nft;
      iarm++;

      // yes, this can also fail....
      if (iarm>=maxarm) dserror("Line search finally failed");
    } // while(nft >= (1.0-lsalpha*lambda)*nf0)
    //--------------------------------------------------------------------
    //                                                  end of line search
    //--------------------------------------------------------------------
    //------------------------------------- update reference residual norm
    f0 = nft;

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;


  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //------------------------------------------------- linear static case
  int nstep = params_.get<int>("nstep",5);
  if (dynkindstat && maxiter==1 && nstep==1)
  {
    if (!myrank_ and printscreen)
    {
      cout << "computed 1 step with 1 iteration: STATIC LINEAR SOLUTION\n";
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }
  //-------------------------------- test whether max iterations was hit
  else if (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter==maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  //--------------------------------------------------- Newton converged
  else
  {
    if (!myrank_ && printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // CmtStruGenAlpha::FullNewtonLineSearch()

/*----------------------------------------------------------------------*
 |  do semi-smooth Newton iteration (public)                  popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::SemiSmoothNewton()
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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // check whether we have a stiffness matrix that is filled
  // and whether mass and damping are present
  // (here you can note the procedural change compared to standard
  // Gen-alpha, where stiff_ must NOT be filled at this point)
  if (!stiff_->Filled()) dserror("stiffness must be filled here");
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
#ifdef CONTACTEIG
  static int globindex = 0;
  ++globindex;
#endif // #ifdef CONTACTEIG

  // active set search and geometrical nonlinearity are merged into
  // ONE Newton loop, thus we have to check for convergence of the
  // active set here, too!
  while ((!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) ||
          !cmtmanager_->GetStrategy().ActiveSetConverged()) && numiter<maxiter)
  {
#ifdef CONTACTTIME
    const double t_start0 = Teuchos::Time::wallTime();
    const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
    //-----------------------------transform to local coordinate systems
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(SystemMatrix(),fresm_);

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

#ifdef CONTACTEIG
    // print to file in matlab format
    std::ostringstream filename;
    const std::string filebase = "sparsematrix";
    filename << "o/matlab_output/" << filebase << "_" << globindex << "_" << numiter+1 << ".mtl";
    LINALG::PrintMatrixInMatlabFormat(filename.str().c_str(),*(SystemMatrix()->EpetraMatrix()));

    // print sparsity pattern to file
    LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
#endif // #ifdef CONTACTEIG

    //---------------------------------------------- solve linear system
    LinearSolve(numiter,tolres,fresmnorm);
    
    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disi_);
#ifdef CONTACTTIME
    const double t_end = Teuchos::Time::wallTime()-t_start;
    cout << "\n***\nSolve: " << t_end << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

    //--------------------------------------- recover disi and Lag. Mult.
    cmtmanager_->GetStrategy().Recover(disi_);

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

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }
#ifdef CONTACTTIME
    const double t_start2 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME

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

      // include potential conditions in fint_ and stiff_
      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
      }

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddFSICondVector(alphas/dt,
                                      fsisurface_->ExtractFSICondVector(dism_),
                                      fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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
#ifdef CONTACTTIME
    const double t_end2 = Teuchos::Time::wallTime()-t_start2;
    cout << "\n***\nDiscret.Evaluate: " << t_end2 << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

    //------------------------------------------ compute dynamic equilibrium
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}

    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
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
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    // (again without contact, this is just Gen-alpha stuff here)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
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
    }

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));
#ifdef CONTACTTIME
    const double t_start3 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME

    //-------------------------make contact modifications to lhs and rhs
    //-------------------- update active set for semi-smooth Newton case
    {
#ifdef CONTACTTIME
      const double t_start31 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
      cmtmanager_->GetStrategy().SetState("displacement",disn_);
#ifdef CONTACTTIME
      const double t_end31 = Teuchos::Time::wallTime()-t_start31;
      cout << "\n***\nContact.SetState: " << t_end31 << " seconds";
      const double t_start32 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
      cmtmanager_->GetStrategy().InitEvalInterface();
#ifdef CONTACTTIME
      const double t_end321 = Teuchos::Time::wallTime()-t_start32;
#endif // #ifdef CONTACTTIME
      cmtmanager_->GetStrategy().InitEvalMortar();
#ifdef CONTACTTIME
      const double t_end322 = Teuchos::Time::wallTime()-t_start32;
      cout << "\nContact.InitMortar: " << t_end321 << " seconds";
      cout << "\nContact.EvalMortar: " << t_end322 << " seconds\n\n";
#endif // #ifdef CONTACTTIME
      
      // friction
      // here the relative movement of the contact bodies is evaluated
      // therefore the current configuration and the according mortar
      // matrices are needed
      cmtmanager_->GetStrategy().EvaluateRelMov();
      
      // this is the correct place to update the active set!!!
      // (on the one hand we need the new weighted gap vector g, which is
      // computed in InitEvalMortar() above and on the other hand we want to
      // run the Evaluate() routine below with the NEW active set already)
#ifdef CONTACTTIME
      const double t_start33 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
      cmtmanager_->GetStrategy().UpdateActiveSetSemiSmooth();
#ifdef CONTACTTIME
      const double t_end33 = Teuchos::Time::wallTime()-t_start33;
      cout << "\nContact.UpdateActiveSet: " << t_end33 << " seconds";
      const double t_start34 = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
      cmtmanager_->GetStrategy().Initialize();
      cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);
#ifdef CONTACTTIME
      const double t_end34 = Teuchos::Time::wallTime()-t_start34;
      cout << "\nContact.StiffFresm: " << t_end34 << " seconds";
#endif // #ifdef CONTACTTIME
    }

    //--------------------------------------------------- contact forces
    cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);
#ifdef CONTACTTIME
    const double t_end3 = Teuchos::Time::wallTime()-t_start3;
    cout << "\n->Contact.Evaluate: " << t_end3 << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

#ifdef CONTACTGMSH2
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    cmtmanager_->GetStrategy().VisualizeGmsh(istep,numiter+1);
#endif // #ifdef CONTACTGMSH2

    //------------------------------------ ----complete stiffness matrix
    stiff_->Complete();
        
    // blank residual DOFs that are on Dirichlet BC
    // in the case of local systems we have to rotate forth and back
    {
      if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);
#ifdef CONTACTTIME
    const double t_end0 = Teuchos::Time::wallTime()-t_start0;
    cout << "\n***\nIteration Step (overall): " << t_end0 << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;


  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //------------------------------------------------- linear static case
  int nstep = params_.get<int>("nstep",5);
  if (dynkindstat && maxiter==1 && nstep==1)
  {
    dserror("ERROR: Linear Static solution not applicable to semi-smooth Newton case");
  }
  //-------------------------------- test whether max iterations was hit
  else if (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter==maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  //--------------------------------------------------- Newton converged
  else
  {
    if (!myrank_ && printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // CmtStruGenAlpha::SemiSmoothNewton()

/*----------------------------------------------------------------------*
 |  do semismooth Newton iteration with line search (public)  popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::SemiSmoothNewtonLineSearch()
{
  if (locsysmanager_ != null) dserror("Locsys not yet implemented for LS-Newton!");
  const int myrank = discret_.Comm().MyPID();

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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // check whether we have a stiffness matrix that is filled
  // and whether mass and damping are present
  // (here you can note the procedural change compared to Gen-alpha
  // without contact, where stiff_ must NOT be filled at this point)
  if (!stiff_->Filled()) dserror("stiffness must be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  // storage of old values to be able to repeat a step
  RCP<Epetra_Vector> disno = rcp(new Epetra_Vector(disn_->Map(),false));
  RCP<Epetra_Vector> dismo = rcp(new Epetra_Vector(dism_->Map(),false));
  RCP<Epetra_Vector> velmo = rcp(new Epetra_Vector(velm_->Map(),false));
  RCP<Epetra_Vector> accmo = rcp(new Epetra_Vector(accm_->Map(),false));

  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  double       f0      = fresmnorm;
  const double sigma0  = 0.1;        // lower bound of lambda
  const double sigma1  = 0.5;        // upper bound of lambda
  const double lsalpha = 1.0e-04;    // decrease requirement to be accept a step
  const int    maxarm  = 10;         // max number line search steps

  if (!myrank) printf("Initial    residual      %15.5e\n",fresmnorm);

  //=================================================== equilibrium loop
  // active set search and geometrical nonlinearity are merged into
  // ONE Newton loop, thus we have to check for convergence of the
  // active set here, too!
  while ((!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) ||
      !cmtmanager_->GetStrategy().ActiveSetConverged()) && numiter<maxiter)
  {
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //---------------------------------------------- solve linear system
    LinearSolve(numiter,tolres,fresmnorm);
    
    //--------------------------------------- recover disi and Lag. Mult.
    cmtmanager_->GetStrategy().Recover(disi_);

    //----------------------------------------- store the current values
    disno->Update(1.0,*disn_,0.0);
    dismo->Update(1.0,*dism_,0.0);
    velmo->Update(1.0,*velm_,0.0);
    accmo->Update(1.0,*accm_,0.0);

    //--------------------------------------------- create the step size
    //                                            trying full step first
    double lambda = 1.0;
    double lamm   = 1.0;
    double lamc   = lambda;
    double iarm   = 0;

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(lambda,*disi_,1.0,*disno,0.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
    disn_->Update(lambda,*disi_,1.0,*disno,0.0);
    dism_->Update(lambda*(1.-alphaf),*disi_,1.0,*dismo,0.0);
#endif
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0,*velmo,0.0);
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
    accm_->Update(lambda*(1.-alpham)/(beta*dt*dt),*disi_,1.0,*accmo,0.0);
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

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }

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
      RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
      disim->Scale(lambda);
#else
      // do not touch disi_ here!
      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
      disim->Scale(lambda*(1.-alphaf));
#endif
      discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
#else
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);
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

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddFSICondVector(alphas/dt,
                                      fsisurface_->ExtractFSICondVector(dism_),
                                      fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

    //------------------------------------------ compute dynamic equilibrium
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}

    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
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
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    // (again without contact, this is just Gen-alpha stuff here)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
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
    }
    
    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

    //-------------------------make contact modifications to lhs and rhs
    {
      cmtmanager_->GetStrategy().SetState("displacement",disn_);

      cmtmanager_->GetStrategy().InitEvalInterface();
      cmtmanager_->GetStrategy().InitEvalMortar();

      // this is the correct place to update the active set!!!
      // (on the one hand we need the new weighted gap vector g, which is
      // computed in InitEvalMortar() above and on the other hand we want to
      // run the Evaluate() routine below with the NEW active set already)
      cmtmanager_->GetStrategy().UpdateActiveSetSemiSmooth();

      cmtmanager_->GetStrategy().Initialize();
      cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);
    }

    //--------------------------------------------------- contact forces
    cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

    //------------------------------------ ----complete stiffness matrix
    stiff_->Complete();
        
#ifdef CONTACTGMSH2
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    cmtmanager_->GetStrategy().VisualizeGmsh(istep,numiter+1);
#endif // #ifdef CONTACTGMSH2

    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);

    double nft = fresmnorm;
    double nf0 = f0;
    double ff0 = nf0*nf0;
    double ffc = nft*nft;
    double ffm = nft*nft;

    //--------------------------------------------------------------------
    // line searching if the step is bad
    //--------------------------------------------------------------------
    bool dolinesearch = true;
    if (dynkindstat && numiter==0) dolinesearch = false;
    while(nft >= (1.0-lsalpha*lambda)*nf0 && dolinesearch)
    {
      //---------------------------- print the rejected step residual norm
      if (!myrank) printf("Bad step                 %15.5e",nft);

      // ---------------------------compute lambda (step length reduction)
      if (!iarm)
        lambda = sigma1*lambda;
      else
      {
        double c2 = lamm*(ffc-ff0)-lamc*(ffm-ff0);
        if (c2>=0.0)
        {
          lambda = sigma1*lamc;
        }
        else
        {
          double c1 = lamc*lamc*(ffm-ff0)-lamm*lamm*(ffc-ff0);
          lambda = -c1*0.5/c2;
          if      (lambda < sigma0*lamc) lambda = sigma0*lamc;
          else if (lambda > sigma1*lamc) lambda = sigma1*lamc;
        }
      }
      //----------------------------- keep track of old lambda values
      lamm = lamc;
      lamc = lambda;

      //---------------------------------- update mid configuration values
      // displacements
      // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
      disn_->Update(lambda,*disi_,1.0,*disno,0.0);
      dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
      disn_->Update(lambda,*disi_,1.0,*disno,0.0);
      dism_->Update(lambda*(1.-alphaf),*disi_,1.0,*dismo,0.0);
#endif
      // velocities
#ifndef STRUGENALPHA_INCRUPDT
      // iterative
      // V_{n+1-alpha_f} := V_{n+1-alpha_f}
      //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
      velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0,*velmo,0.0);
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
      accm_->Update(lambda*(1.-alpham)/(beta*dt*dt),*disi_,1.0,*accmo,0.0);
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

      // zerofy velocity and acceleration in case of statics
      if (dynkindstat)
      {
        velm_->PutScalar(0.0);
        accm_->PutScalar(0.0);
      }

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
        RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
        disim->Scale(lambda);
#else
        // do not touch disi_ here!
        // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
        RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
        disim->Scale(lambda*(1.-alphaf));
#endif
        discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
        discret_.SetState("displacement",disn_);
        discret_.SetState("velocity",veln_);
#else
        discret_.SetState("displacement",dism_);
        discret_.SetState("velocity",velm_);
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

        // do NOT finalize the stiffness matrix to add masses to it later

        // If we have a robin condition we need to modify both the rhs and the
        // matrix diagonal corresponding to the dofs at the robin interface.
        if (structrobin)
        {
          double alphas = params_.get<double>("alpha s",-1.);

          // Add structural part of Robin force
          fsisurface_->AddFSICondVector(alphas/dt,
                                        fsisurface_->ExtractFSICondVector(dism_),
                                        fint_);

          double scale = alphas*(1.-alphaf)/dt;
          const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

      //------------------------------------------ compute dynamic equilibrium
      // Res = M . A_{n+1-alpha_m}
      //     + C . V_{n+1-alpha_f}
      //     + F_int(D_{n+1-alpha_f})
      //     + F_c(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}

      if (dynkindstat)
      {
        // static residual
        // Res = F_int - F_ext
        fresm_->PutScalar(0.0);
      }
      else
      {
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
      }
      // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
      fresm_->Update(1.0,*fextm_,-1.0);
      fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
      fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

#endif

      //---------------------------------------------- build effective lhs
      // (using matrix stiff_ as effective matrix)
      // (again without contact, this is just Gen-alpha stuff here)
      if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
      else
      {
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
      }

      // keep a copy of fresm for contact forces / equilibrium check
      RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

      //-------------------------make contact modifications to lhs and rhs
      {
        cmtmanager_->GetStrategy().SetState("displacement",disn_);

        cmtmanager_->GetStrategy().InitEvalInterface();
        cmtmanager_->GetStrategy().InitEvalMortar();

        // NO update of the active set here, as this would change the
        // system and thus the residual! During line search the active
        // set has to be kept constant!

        cmtmanager_->GetStrategy().Initialize();
        cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);
      }

      //--------------------------------------------------- contact forces
      cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
      int step  = params_.get<int>("step",0);
      int istep = step + 1;
      cmtmanager_->GetStrategy().VisualizeGmsh(istep,numiter+1);
#endif // #ifdef CONTACTGMSH2

      //------------------------------------ ----complete stiffness matrix
      stiff_->Complete();
          
      // blank residual DOFs that are on Dirichlet BC
      {
        Epetra_Vector fresmdbc(*fresm_);
        fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      }

      //---------------------------------------------- build residual norm
      disi_->Norm2(&disinorm);
      fresm_->Norm2(&fresmnorm);

      //---------------------------------------------- keep track of norms
      if (!myrank) printf(" now: %10.5e stepsize %10.5e\n",fresmnorm,lambda); fflush(stdout);
      nft = fresmnorm;
      ffm = ffc;
      ffc = nft*nft;
      iarm++;

      // yes, this can also fail....
      if (iarm>=maxarm) dserror("Line search finally failed");
    } // while(nft >= (1.0-lsalpha*lambda)*nf0)
    //--------------------------------------------------------------------
    //                                                  end of line search
    //--------------------------------------------------------------------
    //------------------------------------- update reference residual norm
    f0 = nft;

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;


  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //------------------------------------------------- linear static case
  int nstep = params_.get<int>("nstep",5);
  if (dynkindstat && maxiter==1 && nstep==1)
  {
    dserror("ERROR: Linear Static solution not applicable to semi-smooth Newton case");
  }
  //-------------------------------- test whether max iterations was hit
  else if (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter==maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  //--------------------------------------------------- Newton converged
  else
  {
    if (!myrank_ && printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // CmtStruGenAlpha::SemiSmoothNewtonLineSearch()

/*----------------------------------------------------------------------*
 |  do PTC iteration (public)                                mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::PTC()
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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  /// check whether we have a stiffness matrix that is filled
  // and whether mass and damping are present
  // (here you can note the procedural change compared to standard
  // Gen-alpha, where stiff_ must NOT be filled at this point)
  if (!stiff_->Filled()) dserror("stiffness must be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  // hard wired ptc parameters
  double ptcdt = 1.0e-03;
  double nc;
  fresm_->NormInf(&nc);
  double dti = 1/ptcdt;
  double dti0 = dti;
  RCP<Epetra_Vector> x0 = rcp(new Epetra_Vector(*disi_));

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
#if 1 // SER
#else // TTI
    double dtim = dti0;
#endif
    dti0 = dti;
    RCP<Epetra_Vector> xm = rcp(new Epetra_Vector(*x0));
    x0->Update(1.0,*disi_,0.0);

    //------------------------------- do ptc modification to effective LHS
    {
      RCP<Epetra_Vector> tmp = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      tmp->PutScalar(dti);
      RCP<Epetra_Vector> diag = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      SystemMatrix()->ExtractDiagonalCopy(*diag);
      diag->Update(1.0,*tmp,1.0);
      SystemMatrix()->ReplaceDiagonalValues(*diag);
    }

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //---------------------------------------------- solve linear system
    LinearSolve(numiter,tolres,fresmnorm);
    
    //--------------------------------------- recover disi and Lag. Mult.
    cmtmanager_->GetStrategy().Recover(disi_);

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

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }

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

      // include potential conditions in fint_ and stiff_
      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
      }

      // do NOT finalize the stiffness matrix to add masses to it later

      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddFSICondVector(alphas/dt,
                                      fsisurface_->ExtractFSICondVector(dism_),
                                      fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->FSICondMap();
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

    if (surf_stress_man_->HaveSurfStress()) dserror("No surface stresses in case of PTC");

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     + F_c(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}

    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
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
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    // (again without contact, this is just Gen-alpha stuff here)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
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
    }
    
    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

    //-------------------------make contact modifications to lhs and rhs
    cmtmanager_->GetStrategy().SetState("displacement",disn_);
    cmtmanager_->GetStrategy().InitEvalInterface();
    cmtmanager_->GetStrategy().InitEvalMortar();
    cmtmanager_->GetStrategy().Initialize();
    cmtmanager_->GetStrategy().Evaluate(stiff_,fresm_,disn_);

    //--------------------------------------------------- contact forces
    cmtmanager_->GetStrategy().InterfaceForces(fresmcopy);

#ifdef CONTACTGMSH2
    dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

    //------------------------------------ ----complete stiffness matrix
    stiff_->Complete();
        
    // blank residual DOFs that are on Dirichlet BC
    // in the case of local systems we have to rotate forth and back
    {
      if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
    }

    // compute inf norm of residual
    double np;
    fresm_->NormInf(&np);

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,dti);
    }

    //------------------------------------ PTC update of artificial time
#if 1
    // *******************************************************************
    // A modification of the last fresm-infnorm (stored in nc) might be
    // necessary here in the case of an Augmented Lagrangian strategy!
    // This is due to the fresm-infnorm becoming approx. zero for the
    // first Newton in the last augmentation loop shortly before Uzawa
    // convergence. As we do not expect PTC convergence problems in these
    // augmentation steps anyway (as we are very close to the solution
    // already), nc could simply be chosen quite high (e.g. 1000 as during
    // intialization, see above) and PTC stabilization would be reduced!
    // *******************************************************************
    // -> Up to now everything works fine! (popp, 01/2010)
    // -> No need for this modification!
    // *******************************************************************
    
    // SER step size control
    dti *= (np/nc);
    dti = max(dti,0.0);
    nc = np;
#else
    {
      // TTI step size control
      double ttau=0.75;
      RCP<Epetra_Vector> d1 = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      d1->Update(1.0,*disi_,-1.0,*x0,0.0);
      d1->Scale(dti0);
      RCP<Epetra_Vector> d0 = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      d0->Update(1.0,*x0,-1.0,*xm,0.0);
      d0->Scale(dtim);
      double dt0 = 1/dti0;
      double dtm = 1/dtim;
      RCP<Epetra_Vector> xpp = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      xpp->Update(2.0/(dt0+dtm),*d1,-2.0/(dt0+dtm),*d0,0.0);
      RCP<Epetra_Vector> xtt = LINALG::CreateVector(SystemMatrix()->RowMap(),false);
      for (int i=0; i<xtt->MyLength(); ++i) (*xtt)[i] = abs((*xpp)[i])/(1.0+abs((*disi_)[i]));
      double ett;
      xtt->MaxValue(&ett);
      ett = ett / (2.*ttau);
      dti = sqrt(ett);
      nc = np;
    }
#endif

    //--------------------------------- increment equilibrium loop index
    ++numiter;

  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //------------------------------------------------- linear static case
  int nstep = params_.get<int>("nstep",5);
  if (dynkindstat && maxiter==1 && nstep==1)
  {
    if (!myrank_ and printscreen)
    {
      cout << "computed 1 step with 1 iteration: STATIC LINEAR SOLUTION\n";
      PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,dti);
    }
  }
  //-------------------------------- test whether max iterations was hit
  else if (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter==maxiter)
  {
     dserror("PTC unconverged in %d iterations",numiter);
  }
  //--------------------------------------------------- Newton converged
  else
  {
    if (!myrank_ and printscreen)
    {
      PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,dti);
    }
  }

  params_.set<int>("num iterations",numiter);


  return;
} // CmtStruGenAlpha::PTC()

/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::UpdateandOutput()
{
  Update();
  Output();
  UpdateElement();
  return;
} // ContactStruGenAlpha::UpdateandOutput()

/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::Update()
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

  const bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

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

  // zerofy velocity and acceleration in case of statics
  if (dynkindstat)
  {
    vel_->PutScalar(0.0);
    acc_->PutScalar(0.0);
  }

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0,*fextn_,0.0);
  // zerofy external force, such that there is no history from load step to load step
  if (dynkindstat)
  {
    fext_->PutScalar(0.0);
  }
#ifdef STRUGENALPHA_FINTLIKETR
  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0,*fintn_,0.0);
#endif

  // update contact
  cmtmanager_->GetStrategy().Update(istep,dis_);

  /*
  Teuchos::RCP<Epetra_Vector> linmom = LINALG::CreateVector(*(discret_.DofRowMap()), true);
  mass_->Multiply(false, *vel_, *linmom);

  int dim = cmtmanager_->GetStrategy().Dim();
  vector<double> sumlinmom(3);
  vector<double> sumangmom(3);

  for (int k=0; k<(discret_.NodeRowMap())->NumMyElements();++k)
  {
    int gid = (discret_.NodeRowMap())->GID(k);
    DRT::Node* mynode = discret_.gNode(gid);

    vector<double> nodelinmom(3);
    vector<double> position(3);

    for (int d=0;d<dim;++d)
    {
      int dofid = (discret_.Dof(mynode))[d];
      nodelinmom[d] = (*linmom)[dofid];
      sumlinmom[d] += nodelinmom[d];
      position[d] = mynode->X()[d] + (*dis_)[dofid];
    }

    vector<double> nodeangmom(3);
    nodeangmom[0] = position[1]*nodelinmom[2]-position[2]*nodelinmom[1];
    nodeangmom[1] = position[2]*nodelinmom[0]-position[0]*nodelinmom[2];
    nodeangmom[2] = position[0]*nodelinmom[1]-position[1]*nodelinmom[0];

    for (int d=0;d<3;++d)
      sumangmom[d] += nodeangmom[d];

    // vector product position x nodelinmom
  }

  // global calculation of kinetic energy
  double kinergy = 0.0;  // total kinetic energy
  {
    linmom->Dot(*vel_,&kinergy);
    kinergy *= 0.5;
  }

  cout << "\n************************************************************************" << endl;
  cout << "CONTACT CONSERVATION QUANTITIES" << endl;
  cout << "Linear Momentum x-direction:  " << sumlinmom[0] << endl;
  cout << "Linear Momentum y-direction:  " << sumlinmom[1] << endl;
  cout << "Linear Momentum z-direction:  " << sumlinmom[2] << endl << endl;
  cout << "Angular Momentum x-direction: " << sumangmom[0] << endl;
  cout << "Angular Momentum y-direction: " << sumangmom[1] << endl;
  cout << "Angular Momentum z-direction: " << sumangmom[2] << endl << endl;
  cout << "Kinetic Energy:               " << kinergy << endl;
  cout << "************************************************************************" << endl;
  */

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

} // CmtStruGenAlpha::Update()


/*----------------------------------------------------------------------*
 |  update element (public)                                     st 03/10|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::UpdateElement()
{
  double timen         = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);
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
    p.set("action","calc_struct_update_imrlike");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    discret_.Evaluate(p,null,null,null,null,null);
  }
#endif
} // CmtStruGenAlpha::UpdateElement()


/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::Output()
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
  INPAR::STR::StressType iostress = params_.get<INPAR::STR::StressType>("io structural stress",INPAR::STR::stress_none);
  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);
  INPAR::STR::StrainType iostrain = params_.get<INPAR::STR::StrainType>("io structural strain",INPAR::STR::strain_none);

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
    cmtmanager_->WriteRestart(output_);
    
    // evaluate interface tractions for postprocessing
    cmtmanager_->PostprocessTractions(output_);

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
    
    // evaluate interface tractions for postprocessing
    cmtmanager_->PostprocessTractions(output_);
  }

  //------------------------------------- do stress calculation and output
  if (updevrystress and !(istep%updevrystress) and iostress!=INPAR::STR::stress_none)
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
    p.set("iostress", iostress);
    p.set("iostrain", iostrain);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
    if (!isdatawritten) output_.NewStep(istep, timen);
    isdatawritten = true;

    switch (iostress)
    {
    case INPAR::STR::stress_cauchy:
      output_.WriteVector("gauss_cauchy_stresses_xyz",*stress,*discret_.ElementColMap());
      break;
    case INPAR::STR::stress_2pk:
      output_.WriteVector("gauss_2PK_stresses_xyz",*stress,*discret_.ElementColMap());
      break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror ("requested stress type not supported");
    }

    switch (iostrain)
    {
    case INPAR::STR::strain_ea:
      output_.WriteVector("gauss_EA_strains_xyz",*strain,*discret_.ElementColMap());
      break;
    case INPAR::STR::strain_gl:
      output_.WriteVector("gauss_GL_strains_xyz",*strain,*discret_.ElementColMap());
      break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror ("requested strain type not supported");
    }
  }

  // print active set
  cmtmanager_->GetStrategy().PrintActiveSet();

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
} // CmtStruGenAlpha::Output()

/*----------------------------------------------------------------------*
 |  nonlinear solution in one time step                      popp  03/10|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::CmtNonlinearSolve()
{
  //********************************************************************
  // get some parameters
  //********************************************************************
  // application type
  INPAR::CONTACT::ApplicationType apptype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(cmtmanager_->GetStrategy().Params(),"APPLICATION");
  
  // strategy type
  INPAR::CONTACT::SolvingStrategy soltype =
    Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtmanager_->GetStrategy().Params(),"STRATEGY");
  
  // semi-smooth Newton type
  bool semismooth = Teuchos::getIntegralValue<int>(cmtmanager_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");
  
  // iteration type
  string equil = params_.get<string>("equilibrium iteration","full newton");
  if (equil != "full newton" && equil != "line search newton" && equil != "ptc")
    dserror("Unknown type of equilibrium iteration");

  //********************************************************************
  // OPTIONS FOR PRIMAL-DUAL ACTIVE SET STRATEGY (PDASS) FOR CONTACT
  // ONLY ONE OPTION FOR MESHTYING -> FULL NEWTON
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

  //********************************************************************
  // Solving Strategy using Lagrangian Multipliers
  //********************************************************************
  if (soltype == INPAR::CONTACT::solution_lagmult)
  {
    //********************************************************************
    // 1) SEMI-SMOOTH NEWTON FOR CONTACT
    // The search for the correct active set (=contact nonlinearity) and
    // the large deformstion linearization (=geometrical nonlinearity) are
    // merged into one semi-smooth Newton method and solved within ONE
    // iteration loop
    //********************************************************************
    if (apptype == INPAR::CONTACT::app_mortarcontact && semismooth)
    {
#ifdef CONTACTTIME
      const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef CONTACTTIME
      
      // nonlinear iteration
      if (equil=="full newton")              SemiSmoothNewton();
      else if (equil=="line search newton")  SemiSmoothNewtonLineSearch();
      else if (equil=="ptc")                 dserror("ERROR: No semi-smooth PTC available!");
      
#ifdef CONTACTTIME
      const double t_end = Teuchos::Time::wallTime()-t_start;
      cout << "\n***\nTime Step (overall): " << t_end << " seconds\n***\n";
#endif // #ifdef CONTACTTIME
    }
    
    //********************************************************************
    // 2) FIXED-POINT APPROACH FOR CONTACT
    // The search for the correct active set (=contact nonlinearity) is
    // represented by a fixed-point approach, whereas the large deformation
    // linearization (=geimetrical nonlinearity) is treated by a standard
    // Newton scheme. This yields TWO nested iteration loops
    //********************************************************************
    else if (apptype == INPAR::CONTACT::app_mortarcontact && !semismooth)
    {
      // predictor type (needed here because of nested loops)
      string pred  = params_.get<string>("predictor","constant");
      
      // active set strategy
      int activeiter = 0;
      while (cmtmanager_->GetStrategy().ActiveSetConverged()==false)
      {
        // increase active set iteration index
        ++activeiter;
        
        // predictor step (except for first active set step)
        if (activeiter>1 && pred=="constant")        ConstantPredictor();
        else if (activeiter>1 && pred=="consistent") ConsistentPredictor();
        
        // nonlinear iteration
        if (equil=="full newton")             FullNewton();
        else if (equil=="line search newton") FullNewtonLineSearch();
        else if (equil=="ptc")                PTC();

        // update of active set (fixed-point)
        cmtmanager_->GetStrategy().UpdateActiveSet();
      }
    }
    
    //********************************************************************
    // 3) STANDARD NEWTON APPROACH FOR MESHTYING
    // No search for the correct active set has to be resolved for mortar
    // meshtying and mortar coupling is linear in this case. Thus, only
    // the large deformation FE problem remains to be solved as nonlinearity
    // Here, a standard Newton scheme is applied and we have ONLY ONE loop.
    //********************************************************************
    else
    {
      // nonlinear iteration
      if (equil=="full newton")             FullNewton();
      else if (equil=="line search newton") FullNewtonLineSearch();
      else if (equil=="ptc")                PTC();
    }
  }

  //********************************************************************
  // Solving Strategy using Regularization Techniques (Penalty Method)
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_penalty)
  {
    // nonlinear iteration
    if (equil=="full newton")             FullNewton();
    else if (equil=="line search newton") FullNewtonLineSearch();
    else if (equil=="ptc")                PTC();

    // update constraint norm
    cmtmanager_->GetStrategy().UpdateConstraintNorm();
  }

  //********************************************************************
  // Solving Strategy using Augmented Lagrange Techniques (with Uzawa)
  //********************************************************************
  else if (soltype == INPAR::CONTACT::solution_auglag)
  {
    // get tolerance and maximum Uzawa steps
    double eps = cmtmanager_->GetStrategy().Params().get<double>("UZAWACONSTRTOL");
    int maxuzawaiter = cmtmanager_->GetStrategy().Params().get<int>("UZAWAMAXSTEPS");

    // Augmented Lagrangian loop (Uzawa)
    int uzawaiter=0;
    do
    {
      // increase iteration index
      ++uzawaiter;
      if (uzawaiter > maxuzawaiter) dserror("Uzawa unconverged in %d iterations",maxuzawaiter);
      if (discret_.Comm().MyPID() == 0) cout << "Starting Uzawa step No. " << uzawaiter << endl;

      // for second, third,... Uzawa step: out-of-balance force
      if (uzawaiter>1) cmtmanager_->GetStrategy().InitializeUzawa(stiff_,fresm_);

      // nonlinear iteration
      if (equil=="full newton")             FullNewton();
      else if (equil=="line search newton") FullNewtonLineSearch();
      else if (equil=="ptc")                PTC();

      // update constraint norm and penalty parameter
      cmtmanager_->GetStrategy().UpdateConstraintNorm(uzawaiter);

      // store Lagrange multipliers for next Uzawa step
      cmtmanager_->GetStrategy().UpdateAugmentedLagrange();
      cmtmanager_->GetStrategy().StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

    } while (cmtmanager_->GetStrategy().ConstraintNorm() >= eps);
          
    // reset penalty parameter
    cmtmanager_->GetStrategy().ResetPenalty();
  }

  return;
} // void CmtStruGenAlpha::CmtNonlinearSolve()


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               popp  02/08|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::Integrate()
{
  // time parameters
  int    step    = params_.get<int>   ("step" ,0);
  int    nstep   = params_.get<int>   ("nstep",5);
  double maxtime = params_.get<double>("max time",0.0);

  // predictor type
  string pred = params_.get<string>("predictor","constant");
  if (pred != "constant" && pred != "consistent")
    dserror("Unknown type of predictor");

  // time step loop
  for (int i=step; i<nstep; ++i)
  {
    // predictor step
    if (pred=="constant")        ConstantPredictor();
    else if (pred=="consistent") ConsistentPredictor();
   
    // nonlinear solution routine
    CmtNonlinearSolve();
    
    // update and output
    UpdateandOutput();
    
    // check if maxtime reached
    double time = params_.get<double>("total time",0.0);
    if (time>=maxtime) break;
  }
  
  return;
} // void CmtStruGenAlpha::Integrate()


/*----------------------------------------------------------------------*
 |  linear solution in one iteration step                    popp  03/10|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::LinearSolve(int numiter, double wanted, double worst)
{
  // adapt solver tolerance
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);
  if (isadapttol && numiter) solver_.AdaptTolerance(wanted,worst,adaptolbetter);
  
  // shape function and strategy types
  INPAR::MORTAR::ShapeFcn shapefcn        = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(cmtmanager_->GetStrategy().Params(),"SHAPEFCN");  
  INPAR::CONTACT::SolvingStrategy soltype = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtmanager_->GetStrategy().Params(),"STRATEGY");
  
  //**********************************************************************
  // Solving strategy using standard Lagrange multipliers
  //**********************************************************************
  if (shapefcn == INPAR::MORTAR::shape_standard && soltype==INPAR::CONTACT::solution_lagmult)
  {
    // saddle point solver call
    cmtmanager_->GetStrategy().SaddlePointSolve(solver_,stiff_,fresm_,disi_,dirichtoggle_,numiter);
  }
  
  //**********************************************************************
  // All other solving strategies (Dual LM, Penalty, Augmented)
  //**********************************************************************
  else
  {
    // standard solver call
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);  
  }
  
  // reset solver tolerance
  solver_.ResetTolerance();
  
  return;
} // void CmtStruGeanAlpha::LinearSolve()


/*----------------------------------------------------------------------*
 |  read restart (public)                                    mwgee 06/07|
 *----------------------------------------------------------------------*/
void CONTACT::CmtStruGenAlpha::ReadRestart(int step)
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

  //**********************************************************************
  // read restart information for contact or meshtying
  //**********************************************************************
  // NOTE: There is an important difference here between contact and
  // meshtying simulations. In both cases, the current coupling operators
  // have to be re-computed for restart, but in the meshtying case this
  // means evaluating DM in the reference configuration!
  // Thus, dis_ is handed in as current displacement state for contact
  // simulations, but zero_ is handed in for meshtying simulations.
  //**********************************************************************
  INPAR::CONTACT::ApplicationType apptype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(cmtmanager_->GetStrategy().Params(),"APPLICATION");
  
  if (apptype == INPAR::CONTACT::app_mortarcontact)
    cmtmanager_->ReadRestart(reader,dis_);
  else if (apptype == INPAR::CONTACT::app_mortarmeshtying)
    cmtmanager_->ReadRestart(reader,zeros_);

  // override current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  return;
} // void CmtStruGenAlpha::ReadRestart()


#endif  // #ifdef CCADISCRET
