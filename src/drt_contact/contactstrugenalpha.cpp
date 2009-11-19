/*!----------------------------------------------------------------------
\file contactstrugenalpha.cpp
\brief Generalized Alpha time integration for structural problems with contact

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

#include "contactstrugenalpha.H"
#include "contactdefines.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"

#include "drt_contact_abstract_strategy.H"
#include "drt_contact_penalty_strategy.H"

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

    // store integration parameter alphaf into cmanager as well
    double alphaf = params_.get<double>("alpha f",0.459);
    contactmanager_ = rcp(new CONTACT::Manager(discret_,alphaf));
  }

  // map containing Dirichlet DOFs
  Teuchos::ParameterList p;
  double time = params_.get<double>("total time"     ,0.0);
  p.set("total time", time);
  RCP<LINALG::MapExtractor> dbcmaps = rcp(new LINALG::MapExtractor());
  discret_.EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null,
                              Teuchos::null, dbcmaps);
  zeros_->PutScalar(0.0); // just in case of change

  // save Dirichlet B.C. status in Contact Manager
  // all CNodes on all interfaces then know if D.B.C.s are applied on their dofs
  contactmanager_->GetStrategy().StoreDirichletStatus(dbcmaps);

  return;
} // ContactStruGenAlpha::ContactStruGenAlpha


/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                     popp 06/08|
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
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  //-------------friction: calculate quantities of reference configuration
  // for frictional contact we need history values and therefore we store
  // the nodal entries of mortar matrices (reference configuration) before
  // first time step
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(contactmanager_->GetStrategy().Params(),"CONTACT");
  
  INPAR::CONTACT::SolvingStrategy strattype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(contactmanager_->GetStrategy().Params(),"STRATEGY");

  if(params_.get<int>("step") == 0 && ctype != INPAR::CONTACT::contact_normal)
  {
  	// set state and do mortar calculation
    contactmanager_->GetStrategy().SetState("displacement",disn_);
    contactmanager_->GetStrategy().InitEvalInterface();
    contactmanager_->GetStrategy().InitEvalMortar();

    // store contact state to contact nodes (active or inactive)
    contactmanager_->GetStrategy().StoreNodalQuantities(AbstractStrategy::activeold);

  	// store D and M to old ones
    contactmanager_->GetStrategy().StoreDM("old");

   	// store nodal entries from D and M to old ones
    contactmanager_->GetStrategy().StoreDMToNodes(AbstractStrategy::dm);
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
  stiff_->Complete();

  // keep a copy of fresm for contact forces / equilibrium check
  RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

  //------------------------- make contact modifications to lhs and rhs
  contactmanager_->GetStrategy().SetState("displacement",disn_);

  contactmanager_->GetStrategy().InitEvalInterface();
  contactmanager_->GetStrategy().InitEvalMortar();
  
  // friction
  // here the relative movement of the contact bodies is evaluated
  // therefore the current configuration and the according mortar
  // matrices are needed
  // it is only evaluated (resetted) for penalty strategy  
  if (strattype == INPAR::CONTACT::solution_penalty)
  { 
    if(ctype != INPAR::CONTACT::contact_normal)
      contactmanager_->GetStrategy().EvaluateRelMov(disi_);
  }
  contactmanager_->GetStrategy().Initialize();
  contactmanager_->GetStrategy().Evaluate(SystemMatrix(),fresm_);

  //---------------------------------------------------- contact forces
  contactmanager_->GetStrategy().ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  contactmanager_->GetStrategy().VisualizeGmsh(istep,0);
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
} // ContactStruGenAlpha::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                       popp 06/08|
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

  //-------------friction: calculate quantities of reference configuration
  // for frictional contact we need history values and therefore we store
  // the nodal entries of mortar matrices (reference configuration) before
  // first time step
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(contactmanager_->GetStrategy().Params(),"CONTACT");
  
  INPAR::CONTACT::SolvingStrategy strattype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(contactmanager_->GetStrategy().Params(),"STRATEGY");

  if(params_.get<int>("step") == 0 && ctype != INPAR::CONTACT::contact_normal)
  {
  	// set state and do mortar calculation
    contactmanager_->GetStrategy().SetState("displacement",disn_);
    contactmanager_->GetStrategy().InitEvalInterface();
    contactmanager_->GetStrategy().InitEvalMortar();

    // store contact state to contact nodes (active or inactive)
    contactmanager_->GetStrategy().StoreNodalQuantities(AbstractStrategy::activeold);

  	// store D and M to old ones
    contactmanager_->GetStrategy().StoreDM("old");

  	// store nodal entries from D and M to old ones
    contactmanager_->GetStrategy().StoreDMToNodes(AbstractStrategy::dm);
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
  stiff_->Complete();

  // keep a copy of fresm for contact forces / equilibrium check
  RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

  //-------------------------- make contact modifications to lhs and rhs
  contactmanager_->GetStrategy().SetState("displacement",disn_);

  contactmanager_->GetStrategy().InitEvalInterface();
  contactmanager_->GetStrategy().InitEvalMortar();
  
  // friction
  // here the relative movement of the contact bodies is evaluated
  // therefore the current configuration and the according mortar
  // matrices are needed
  // it is only evaluated (resetted) for penalty strategy  
  if (strattype == INPAR::CONTACT::solution_penalty)
  { 
    if(ctype != INPAR::CONTACT::contact_normal)
      contactmanager_->GetStrategy().EvaluateRelMov(disi_);
  }

  contactmanager_->GetStrategy().Initialize();
  contactmanager_->GetStrategy().Evaluate(SystemMatrix(),fresm_);

  //---------------------------------------------------- contact forces
  contactmanager_->GetStrategy().ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  contactmanager_->GetStrategy().VisualizeGmsh(istep,0);
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
} // ContactStruGenAlpha::ConstantPredictor()



/*----------------------------------------------------------------------*
 |  setup equilibrium with additional external forces        u.may 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::ApplyExternalForce(  const STR::UTILS::MapExtractor& extractor,
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
  contactmanager_->GetStrategy().SetState("displacement",disn_);

  contactmanager_->GetStrategy().InitEvalInterface();
  contactmanager_->GetStrategy().InitEvalMortar();
  
  contactmanager_->GetStrategy().Initialize();
  contactmanager_->GetStrategy().Evaluate(SystemMatrix(),fresm_);

  //---------------------------------------------------- contact forces
  contactmanager_->GetStrategy().ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
  int istep = step + 1;
  contactmanager_->GetStrategy().VisualizeGmsh(istep,0);
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
void CONTACT::ContactStruGenAlpha::computeF(const Epetra_Vector& x, Epetra_Vector& F)
{
  dserror("not yet impemented for contactstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 |  compute Jacobian at x (public)                           u.may 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::computeJacobian(const Epetra_Vector& x)
{
  dserror("not yet impemented for contactstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 |  build linear system matrix and rhs (public)              u.kue 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  dserror("not yet impemented for contactstrugenalpha");
  return;
}


/*----------------------------------------------------------------------*
 | linear relaxation solve (public)                          u.kue 03/07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::ContactStruGenAlpha::LinearRelaxationSolve(
      Teuchos::RCP<Epetra_Vector> relax)
{
  dserror("not yet impemented for contactstrugenalpha");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                              popp 06/08|
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
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

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

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disi_);

    //--------------------------------------- recover disi and Lag. Mult.
    contactmanager_->GetStrategy().Recover(disi_);

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
    stiff_->Complete();

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

    //-------------------------make contact modifications to lhs and rhs
    {
      contactmanager_->GetStrategy().SetState("displacement",disn_);

      contactmanager_->GetStrategy().InitEvalInterface();
      contactmanager_->GetStrategy().InitEvalMortar();

      // friction
      // here the relative movement of the contact bodies is evaluated
      // therefore the current configuration and the according mortar
      // matrices are needed
      INPAR::CONTACT::ContactType ctype =
        Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(contactmanager_->GetStrategy().Params(),"CONTACT");
      if(ctype != INPAR::CONTACT::contact_normal)
        contactmanager_->GetStrategy().EvaluateRelMov(disi_);

      contactmanager_->GetStrategy().Initialize();
      contactmanager_->GetStrategy().Evaluate(SystemMatrix(),fresm_);
    }

    //--------------------------------------------------- contact forces
    contactmanager_->GetStrategy().ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
    dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

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
 |  do semi-smooth Newton iteration (public)                  popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::SemiSmoothNewton()
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

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

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
          !contactmanager_->GetStrategy().ActiveSetConverged()) && numiter<maxiter)
  {
#ifdef CONTACTTIME
    const double t_start0 = ds_cputime();
    const double t_start = ds_cputime();
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

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disi_);
#ifdef CONTACTTIME
    const double t_end = ds_cputime()-t_start;
    cout << "\n***\nSolve: " << t_end << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

    //--------------------------------------- recover disi and Lag. Mult.
    contactmanager_->GetStrategy().Recover(disi_);

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
    const double t_start2 = ds_cputime();
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
    const double t_end2 = ds_cputime()-t_start2;
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
    stiff_->Complete();

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));
#ifdef CONTACTTIME
    const double t_start3 = ds_cputime();
#endif // #ifdef CONTACTTIME

    //-------------------------make contact modifications to lhs and rhs
    //-------------------- update active set for semi-smooth Newton case
    {
#ifdef CONTACTTIME
      const double t_start31 = ds_cputime();
#endif // #ifdef CONTACTTIME
      contactmanager_->GetStrategy().SetState("displacement",disn_);
#ifdef CONTACTTIME
      const double t_end31 = ds_cputime()-t_start31;
      cout << "\n***\nContact.SetState: " << t_end31 << " seconds";
      const double t_start32 = ds_cputime();
#endif // #ifdef CONTACTTIME
      contactmanager_->GetStrategy().InitEvalInterface();
#ifdef CONTACTTIME
      const double t_end321 = ds_cputime()-t_start32;
#endif // #ifdef CONTACTTIME
      contactmanager_->GetStrategy().InitEvalMortar();
#ifdef CONTACTTIME
      const double t_end322 = ds_cputime()-t_start32;
      cout << "\nContact.InitMortar: " << t_end321 << " seconds";
      cout << "\nContact.EvalMortar: " << t_end322 << " seconds\n\n";
#endif // #ifdef CONTACTTIME

      // friction
      // here the relative movement of the contact bodies is evaluated
      // therefore the current configuration and the according mortar
      // matrices are needed
      INPAR::CONTACT::ContactType ctype =
        Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(contactmanager_->GetStrategy().Params(),"CONTACT");
      if(ctype != INPAR::CONTACT::contact_normal)
        contactmanager_->GetStrategy().EvaluateRelMov(disi_);

      // this is the correct place to update the active set!!!
      // (on the one hand we need the new weighted gap vector g, which is
      // computed in EvaluateMortar() above and on the other hand we want to
      // run the Evaluate()routine below with the NEW active set already)
#ifdef CONTACTTIME
      const double t_start33 = ds_cputime();
#endif // #ifdef CONTACTTIME
      contactmanager_->GetStrategy().UpdateActiveSetSemiSmooth();
#ifdef CONTACTTIME
      const double t_end33 = ds_cputime()-t_start33;
      cout << "\nContact.UpdateActiveSet: " << t_end33 << " seconds";
      const double t_start34 = ds_cputime();
#endif // #ifdef CONTACTTIME
      contactmanager_->GetStrategy().Initialize();
      contactmanager_->GetStrategy().Evaluate(SystemMatrix(),fresm_);
#ifdef CONTACTTIME
      const double t_end34 = ds_cputime()-t_start34;
      cout << "\nContact.StiffFresm: " << t_end34 << " seconds";
#endif // #ifdef CONTACTTIME
    }

    //--------------------------------------------------- contact forces
    contactmanager_->GetStrategy().ContactForces(fresmcopy);
#ifdef CONTACTTIME
    const double t_end3 = ds_cputime()-t_start3;
    cout << "\n->Contact.Evaluate: " << t_end3 << " seconds\n***\n";
#endif // #ifdef CONTACTTIME

#ifdef CONTACTGMSH2
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    contactmanager_->GetStrategy().VisualizeGmsh(istep,numiter+1);
#endif // #ifdef CONTACTGMSH2

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
    const double t_end0 = ds_cputime()-t_start0;
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
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // ContactStruGenAlpha::SemiSmoothNewton()


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
  contactmanager_->GetStrategy().Update(istep);

  /*
  Teuchos::RCP<Epetra_Vector> linmom = LINALG::CreateVector(*(discret_.DofRowMap()), true);
  mass_->Multiply(false, *vel_, *linmom);

  int dim = contactmanager_->GetStrategy().Dim();
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

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // informations and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information/quantities. These are:
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(contactmanager_->GetStrategy().Params(),"CONTACT");
  if(ctype != INPAR::CONTACT::contact_normal)
  {

  	// store contact state to contact nodes (active or inactive)
    contactmanager_->GetStrategy().StoreNodalQuantities(AbstractStrategy::activeold);

  	// store nodal entries of D and M to old ones
    contactmanager_->GetStrategy().StoreDMToNodes(AbstractStrategy::dm);
    
    // store nodal entries form penalty contact tractions to old ones
    contactmanager_->GetStrategy().StoreDMToNodes(AbstractStrategy::pentrac);

    // store the displacements to contact nodes
    contactmanager_->GetStrategy().SetState("olddisplacement",dis_);
  }

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
    p.set("action","calc_struct_update_imrlike");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    p.set("alpha f",alphaf);
    discret_.Evaluate(p,null,null,null,null,null);
  }
#endif

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
  INPAR::STR::StressType iostress = params_.get<INPAR::STR::StressType>("io structural stress",INPAR::STR::stress_none);
  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);
  INPAR::STR::StrainType iostrain = params_.get<INPAR::STR::StrainType>("io structural strain",INPAR::STR::strain_none);

  int    writeresevry  = params_.get<int>   ("write restart every"    ,0);

  bool   printscreen   = params_.get<bool>  ("print to screen"        ,true);
  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
  if (!errfile) printerr = false;

  bool isdatawritten = false;

  //------------------------------- evaluate contact stresses for output
  contactmanager_->GetStrategy().OutputStresses();
	RCP<Epetra_Map> problem = (contactmanager_->GetStrategy().ProblemRowMap());

	// normal direction
  RCP<Epetra_Vector> normalstresses = contactmanager_->GetStrategy().ContactNorStress();
	RCP<Epetra_Vector> normalstressesexp = rcp(new Epetra_Vector(*problem));
	LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
	RCP<Epetra_Vector> tangentialstresses = contactmanager_->GetStrategy().ContactTanStress();
	RCP<Epetra_Vector> tangentialstressesexp = rcp(new Epetra_Vector(*problem));
	LINALG::Export(*tangentialstresses, *tangentialstressesexp);

  //------------------------------------------------- write restart step
  if (writeresevry and istep%writeresevry==0)
  {
    output_.WriteMesh(istep,timen);
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("norcontactstress",normalstressesexp);
    output_.WriteVector("tancontactstress",tangentialstressesexp);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    isdatawritten = true;

    // write restart information for contact
    contactmanager_->WriteRestart(output_);

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
    output_.WriteVector("norcontactstress",normalstressesexp);
    output_.WriteVector("tancontactstress",tangentialstressesexp);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    output_.WriteElementData();
    isdatawritten = true;
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
  contactmanager_->GetStrategy().PrintActiveSet();

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

#ifdef CONTACTGMSH3
  contactmanager_->GetStrategy().VisualizeGmsh(0,0);
#endif // #ifdef CONTACTGMSH2

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  // unknown types of nonlinear iteration schemes
  if (equil != "full newton" && equil != "line search newton")
    dserror("Unknown type of equilibrium iteration");

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

  INPAR::CONTACT::SolvingStrategy soltype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(contactmanager_->GetStrategy().Params(),"STRATEGY");
  bool semismooth = Teuchos::getIntegralValue<int>(contactmanager_->GetStrategy().Params(),"SEMI_SMOOTH_NEWTON");

  // Solving Strategy using Lagrangian Multipliers
  // here we can do sa SEMI-SMOOTH NEWTON  -or- a FIXED-POINT APPROACH

  if (soltype == INPAR::CONTACT::solution_lagmult)
  {
    if (discret_.Comm().MyPID() == 0)
      cout << "===== Lagrange multiplier strategy =============================" << endl;

    //********************************************************************
    // 1) SEMI-SMOOTH NEWTON
    // The search for the correct active set (=contact nonlinearity) and
    // the large deformstion linearization (=geometrical nonlinearity) are
    // merged into one semi-smooth Newton method and solved within ONE
    // iteration loop
    //********************************************************************
    if (semismooth)
    {
      // LOOP1: time steps
      for (int i=step; i<nstep; ++i)
      {
  #ifdef CONTACTTIME
        const double t_start = ds_cputime();
  #endif // #ifdef CONTACTTIME

        // predictor step
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();

        // LOOP2: nonlinear iteration (Newton)
        if (equil=="full newton") SemiSmoothNewton();
        else if (equil=="line search newton") SemiSmoothNewtonLineSearch();

        UpdateandOutput();

  #ifdef CONTACTTIME
        const double t_end = ds_cputime()-t_start;
        cout << "\n***\nTime Step (overall): " << t_end << " seconds\n***\n";
  #endif // #ifdef CONTACTTIME

        double time = params_.get<double>("total time",0.0);
        if (time>=maxtime) break;
      }
    }

    //********************************************************************
    // FIXED-POINT APPROACH
    // The search for the correct active set (=contact nonlinearity) is
    // represented by a fixed-point approach, whereas the large deformation
    // linearization (=geimetrical nonlinearity) is treated by a standard
    // Newton scheme. This yields TWO nested iteration loops
    //********************************************************************
    else
    {
      // LOOP1: time steps
      for (int i=step; i<nstep; ++i)
      {

        // LOOP2: active set strategy
        while (contactmanager_->GetStrategy().ActiveSetConverged()==false)
        {
          // predictor step
          if      (predictor==1) ConstantPredictor();
          else if (predictor==2) ConsistentPredictor();

          // LOOP3: nonlinear iteration (Newton)
          if (equil=="full newton") FullNewton();
          else if (equil=="line search newton") FullNewtonLineSearch();

          // update of active set (fixed-point)
          contactmanager_->GetStrategy().UpdateActiveSet();
        }

        UpdateandOutput();
        double time = params_.get<double>("total time",0.0);
        if (time>=maxtime) break;
      }
    }
    //********************************************************************
    // END: options for primal-dual active set strategy (PDASS)
    //********************************************************************
  }

  // Solving Strategy using Regularization Techniques (Penalty Method)
  // here we have just one option: the ordinary NEWTON
  else if (soltype == INPAR::CONTACT::solution_penalty)
  {
    if (discret_.Comm().MyPID() == 0)
      cout << "===== Penalty strategy =========================================" << endl;

    // explicitely store gap-scaling
    PenaltyStrategy& strategy = dynamic_cast<PenaltyStrategy&> (contactmanager_->GetStrategy());
    strategy.SaveReferenceState(disn_);

    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      // predictor step
      if (predictor==1)
        ConstantPredictor();
      else if (predictor==2)
        ConsistentPredictor();

      // LOOP2: nonlinear iteration (Newton)
      if (equil=="full newton") FullNewton();
      else if (equil=="line search newton") FullNewtonLineSearch();

      strategy.UpdateConstraintNorm();
      UpdateandOutput();
      double time = params_.get<double>("total time", 0.0);
      if (time>=maxtime) break;
    }

  }

  // Solving Strategy using Augmented Lagrange Techniques (with Uzawa)
  // here we have just one option: the ordinary NEWTON
  else if (soltype == INPAR::CONTACT::solution_auglag)
  {
    if (discret_.Comm().MyPID() == 0)
      cout << "===== Augmented Lagrange strategy ================================" << endl;

    // explicitely store gap-scaling
    PenaltyStrategy& strategy = dynamic_cast<PenaltyStrategy&> (contactmanager_->GetStrategy());
    strategy.SaveReferenceState(disn_);

    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      // reset penalty parameter
      strategy.ResetPenalty();

      // predictor step
       if (predictor==1)
         ConstantPredictor();
       else if (predictor==2)
         ConsistentPredictor();

      // get tolerance and maximum Uzawa steps
      double eps = contactmanager_->GetStrategy().Params().get<double>("UZAWACONSTRTOL");
      int maxuzawaiter = contactmanager_->GetStrategy().Params().get<int>("UZAWAMAXSTEPS");

      // LOOP2: augmented Lagrangian (Uzawa)
      int uzawaiter=0;
      do
      {
        // increase iteration index
        ++uzawaiter;
        if (uzawaiter > maxuzawaiter)
          dserror("Uzawa unconverged in %d iterations",maxuzawaiter);

        if (discret_.Comm().MyPID() == 0)
          cout << "Starting Uzawa step No. " << uzawaiter << endl;

        // LOOP3: nonlinear iteration (Newton)
        if (equil=="full newton") FullNewton();
        else if (equil=="line search newton") FullNewtonLineSearch();

        // update constraint norm and penalty parameter
        strategy.UpdateConstraintNorm(uzawaiter);

        // store Lagrange multipliers for next Uzawa step
        strategy.UpdateAugmentedLagrange();
        strategy.StoreNodalQuantities(AbstractStrategy::lmuzawa);

      } while (strategy.ConstraintNorm() >= eps);

      UpdateandOutput();
      double time = params_.get<double>("total time", 0.0);
      if (time>=maxtime) break;
    }
  }

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
  contactmanager_->ReadRestart(reader,dis_);

  // override current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  return;
} // void ContactStruGenAlpha::ReadRestart()


#endif  // #ifdef CCADISCRET
