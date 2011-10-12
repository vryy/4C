/*!----------------------------------------------------------------------
\file beam3contactstrugenalpha.cpp
\brief Generalized Alpha time integration for structural problems with 3D beam contact

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET


#include "beam3contactstrugenalpha.H"
#include <iostream>
#include "../drt_constraint/constraint_manager.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3ContactStruGenAlpha::Beam3ContactStruGenAlpha(
  ParameterList& params,  DRT::Discretization& dis,
  LINALG::Solver& solver, IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output)
{
  // print welcome message
  if (!myrank_)
  {
    cout << "\n*******************************";
    cout << "\n* Welcome to 3D BEAM CONTACT! *";
    cout << "\n*******************************\n" << endl;
  }
  
  // -------------------------------------------------------------------
  // check again whether we have beam contact and create beam3cmanager
  // -------------------------------------------------------------------
  // Check for beam contact
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->MeshtyingAndContactParams();
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION") == INPAR::CONTACT::app_beamcontact)
    beamcmanager_ = rcp(new CONTACT::Beam3cmanager(dis));
  else
    dserror("ERROR: How did you arrive here...???");
  
  return;
} // Beam3ContactStruGenAlpha::Beam3ContactStruGenAlpha

  
/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::ConstantPredictor()
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
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  bool   loadlin     = params_.get<bool>("LOADLIN",false);

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

  //-------------------------------- zero out effective stiffness matrix
  stiff_->Zero();

  //--------------------------------------------------- predicting state
  // constant predictor : displacement in domain
  disn_->Update(1.0,*dis_,0.0);
  veln_->Update(1.0,*vel_,0.0);

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
    discret_.SetState("velocity",veln_);
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
    if (!loadlin) discret_.EvaluateNeumann(p,fextn_);
    else
    {
      *fextlinold_ = *fextlin_;
      if (!fextlinold_->Filled()) fextlinold_->Complete();
      fextlin_->Zero();
      discret_.EvaluateNeumann(p,fextn_,fextlin_);
      fextlin_->Complete();
    }
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

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,SystemMatrix());
    }

    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
    }

    if (constrMan_->HaveConstraint())
    {
      ParameterList pcon;
      pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
      constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,SystemMatrix(),pcon);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // build residual
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

  //**********************************************************************
  //**********************************************************************
  // evaluate beam contact


  beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);
  


#ifdef GMSHNEWTONSTEPS
  // create gmsh-output to visualize predictor step
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  beamcmanager_->GmshOutput(*disn_,istep,0);
  beamcmanager_->ConsoleOutput();
#endif
  //**********************************************************************
  //**********************************************************************
  
  // blank residual DOFs that are on Dirichlet BC
  // in the case of local systems we have to rotate forth and back
  {
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
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
} // StruGenAlpha::ConstantPredictor()



/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::ConsistentPredictor()
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
  bool   loadlin     = params_.get<bool>("LOADLIN",false);
  if (loadlin) dserror("Pressure load linearization currently only with ConstantPredictor and FullNewton");

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
    discret_.SetState("velocity",veln_);
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

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,SystemMatrix());
    }

    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
    }

    if (constrMan_->HaveConstraint())
    {
      ParameterList pcon;
      pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
      constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,SystemMatrix(),pcon);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // build residual
  if (dynkindstat)
  {
    // static residual
    // Res = F_int - F_ext
    fresm_->PutScalar(0.0);
  }
  else
  {
    // dynamic residual

    // Res = M . A_{n+1-alpha_m} + C . V_{n+1-alpha_f} + F_int(D_{n+1-alpha_f}) - F_{ext;n+1-alpha_f}

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



  //**********************************************************************
  //**********************************************************************
  // evaluate beam contact
  beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);
  

#ifdef GMSHNEWTONSTEPS
  // create gmsh-output to visualize predictor step
  int step  = params_.get<int>("step",0);
  int istep = step + 1;
  beamcmanager_->GmshOutput(*disn_,istep,0);
  beamcmanager_->ConsoleOutput();
#endif
  //**********************************************************************
  //**********************************************************************
      


  // blank residual DOFs that are on Dirichlet BC
  // in the case of local systems we have to rotate forth and back
  {
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
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
} // StruGenAlpha::ConsistentPredictor()

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::FullNewton()
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
  bool  loadlin    = params_.get<bool>("LOADLIN",false);
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  // check whether mass and damping are present
  // note: the stiffness matrix might be filled already
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

  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  InitializeNewtonUzawa();
  
  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    if (dynkindstat) // do nothing, we have the ordinary stiffness matrix ready
    {
      if (loadlin);
    }
    else
    {
#ifdef STRUGENALPHA_BE
      stiff_->Add(*mass_,false,(1.-alpham)/(delta*dt*dt),1.-alphaf);
#else
      stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
#endif
      if(damping)
      {
#ifdef STRUGENALPHA_BE
        stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(delta*dt),1.0);
#else
        stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
#endif
      }

#if 1
      if (loadlin)
      {
//        stiff_->Add(*fextlin_,false,-1.0,1.0); // mid point
        stiff_->Add(*fextlin_,false,-(1.-alphaf),1.0); // TR
        stiff_->Add(*fextlinold_,false,-alphaf,1.0);
      }
#endif
    }

    stiff_->Complete();

    //-----------------------------transform to local coordinate systems
    if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(SystemMatrix(),fresm_);

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
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(disi_);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(1.0,*disi_,1.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
    dism_->Update(1.-alphaf,*disi_,1.0);
    disn_->Update(1.0,*disi_,1.0);
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

    //--------------------------- recompute external forces if nonlinear
    // at state n, the external forces and linearization are interpolated at
    // time 1-alphaf in a TR fashion
    if (loadlin)
    {
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_eleload");
      // other parameters needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      p.set("alpha f",alphaf);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
//      discret_.SetState("displacement",dism_); // mid point
//      discret_.SetState("velocity",velm_);
      fextn_->PutScalar(0.0); // TR
//      fextm_->PutScalar(0.0);
      fextlin_->Zero();
//      discret_.EvaluateNeumann(p,fextm_,fextlin_);
      discret_.EvaluateNeumann(p,fextn_,fextlin_);
      fextlin_->Complete();
      discret_.ClearState();
      fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);
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

      if (surf_stress_man_->HaveSurfStress())
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,SystemMatrix());
      }

      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
      }

      if (constrMan_->HaveConstraint())
      {
        ParameterList pcon;
        pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
        constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,SystemMatrix(),pcon);
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

    //------------------------------------------ compute residual forces
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
    
    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    // Res = Res + (1-alpha_f) * F_c(D_{n+1}) + alpha_f * F_c(D_{n})
    beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);
    
#ifdef GMSHNEWTONSTEPS
    // Create gmsh-output to visualize every step of newton iteration
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    beamcmanager_->GmshOutput(*disn_,istep,numiter+1);
    beamcmanager_->ConsoleOutput();
#endif
    //**********************************************************************
    //**********************************************************************
    
    // blank residual DOFs that are on Dirichlet BC
    // in the case of local systems we have to rotate forth and back
    {
      if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
      Epetra_Vector fresmcopy(*fresm_);

      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

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
  //=============================================== end equilibrium loop
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
} // StruGenAlpha::FullNewton()

/*----------------------------------------------------------------------*
 |  initialize Newton for 2nd, 3rd, ... Uzawa iteration       popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::InitializeNewtonUzawa()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  bool structrobin = params_.get<bool>  ("structrobin"            ,false);
  if (!errfile) printerr = false;
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  bool  loadlin    = params_.get<bool>("LOADLIN",false);

  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  if (beamcmanager_->GetUzawaIter() > 1)
  {
    //--------------------------- recompute external forces if nonlinear
    // at state n, the external forces and linearization are interpolated at
    // time 1-alphaf in a TR fashion
    if (loadlin)
    {
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_eleload");
      // other parameters needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      p.set("alpha f",alphaf);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
//      discret_.SetState("displacement",dism_); // mid point
//      discret_.SetState("velocity",velm_);
      fextn_->PutScalar(0.0); // TR
//      fextm_->PutScalar(0.0);
      fextlin_->Zero();
//      discret_.EvaluateNeumann(p,fextm_,fextlin_);
      discret_.EvaluateNeumann(p,fextn_,fextlin_);
      fextlin_->Complete();
      discret_.ClearState();
      fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);
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

      if (surf_stress_man_->HaveSurfStress())
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,SystemMatrix());
      }

      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,SystemMatrix());
      }

      if (constrMan_->HaveConstraint())
      {
        ParameterList pcon;
        pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
        constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,SystemMatrix(),pcon);
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

    //------------------------------------------ compute residual forces
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
    
    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    // Res = Res + (1-alpha_f) * F_c(D_{n+1}) + alpha_f * F_c(D_{n})
    beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);
    //**********************************************************************
    //**********************************************************************
    
    // blank residual DOFs that are on Dirichlet BC
    // in the case of local systems we have to rotate forth and back
    {
      if (locsysmanager_ != null) locsysmanager_->RotateGlobalToLocal(fresm_);
      Epetra_Vector fresmcopy(*fresm_);

      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

      if (locsysmanager_ != null) locsysmanager_->RotateLocalToGlobal(fresm_);
    }
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  integrate in time                                        popp  04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::Integrate()
{
  // some paramaters
  int    step    = params_.get<int>   ("step" ,0);
  int    nstep   = params_.get<int>   ("nstep",5);
  double maxtime = params_.get<double>("max time",0.0);

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred   = params_.get<string>("predictor","constant");
  int predictor = -1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  // unknown types of nonlinear iteration schemes
  if (equil != "full newton")
    dserror("Unknown type of equilibrium iteration");
  
  // decide which solving strategy
  INPAR::CONTACT::SolvingStrategy soltype =
  DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(beamcmanager_->InputParameters(),"STRATEGY");
  
  // decide wether the original or the modified gapfunction definition is used
  bool newgapfunction =  DRT::INPUT::IntegralValue<int>(beamcmanager_->InputParameters(),"BEAMS_NEWGAP");

  // decide wether the tangent field should be smoothed or not
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->MeshtyingAndContactParams();
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::Smoothing>(scontact,"BEAMS_SMOOTHING") == INPAR::CONTACT::bsm_none)
  {
    //cout << "Test BEAMS_SMOOTHING" << INPAR::CONTACT::bsm_none << endl;
  }

  //**********************************************************************
  // solving strategy using regularization with penalty method
  // (nonlinear solution approach: ordinary NEWTON)
  //**********************************************************************
  if (soltype == INPAR::CONTACT::solution_penalty)
  {
    if (discret_.Comm().MyPID() == 0)
      cout << "===== Penalty strategy =========================================" << endl;
 
    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      // set normal vector of last time "normal_" to old normal vector "normal_old_"
      if (newgapfunction) {beamcmanager_->ShiftAllNormal();}

      if (predictor==1) {

        ConstantPredictor();

      }

      //here

      else if (predictor==2) {

        ConsistentPredictor();

      }
       // LOOP2: nonlinear iteration (Newton)


      FullNewton();
      

      // update constraint norm
      beamcmanager_->UpdateConstrNorm();
              
      // update and output of time step
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  //**********************************************************************
  
  //**********************************************************************
  // solving strategy using regularization with augmented Lagrange method
  // (nonlinear solution approach: nested UZAWA NEWTON)
  //**********************************************************************
  else if (soltype == INPAR::CONTACT::solution_auglag)
  {
    if (discret_.Comm().MyPID() == 0)
      cout << "===== Augmented Lagrange strategy ==============================" << endl;

    // LOOP1: time steps
    for (int i=step; i<nstep; ++i)
    {
      // set normal vector of last time step "normal_" to old normal vector "normal_old_"
      if (newgapfunction) beamcmanager_->ShiftAllNormal();

      // Initialize all lmuzawa to zero at beginning of new time step
      beamcmanager_->ResetAlllmuzawa();

      // predictor step
      if (predictor==1)
        ConstantPredictor();
      else if (predictor==2)
        ConsistentPredictor();
 
      // get tolerance and maximum number of Uzawa steps from input file
      double eps = beamcmanager_->InputParameters().get<double>("UZAWACONSTRTOL");
      int maxuzawaiter = beamcmanager_->InputParameters().get<int>("UZAWAMAXSTEPS");

      // LOOP2: augmented Lagrangian (Uzawa)
      beamcmanager_->ResetUzawaIter();
      do
      {
        // increase iteration index by one
        beamcmanager_->UpdateUzawaIter();
        if (beamcmanager_->GetUzawaIter() > maxuzawaiter)
          dserror("Uzawa unconverged in %d iterations",maxuzawaiter);

        if (discret_.Comm().MyPID() == 0)
          cout << endl << "Starting Uzawa step No. " << beamcmanager_->GetUzawaIter() << endl;
        
        // LOOP3: nonlinear iteration (Newton)
        FullNewton();
        
        // update constraint norm and penalty parameter
        beamcmanager_->UpdateConstrNorm();
        
        // update Uzawa Lagrange multipliers
        beamcmanager_->UpdateAlllmuzawa();
        
      } while (abs(beamcmanager_->GetConstrNorm()) >= eps);
  
      // reset penalty parameter
      beamcmanager_->ResetCurrentpp();
            
      // update and output of time step
      UpdateandOutput();
      double time = params_.get<double>("total time", 0.0);
      if (time>=maxtime) break;
    }
  }
  //**********************************************************************

  return;
} // void beam3contactstrugenalpha::Integrate()

/*----------------------------------------------------------------------*
 |  do update (public)                                       mwgee 03/07|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3ContactStruGenAlpha::Update()
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
  
  //**********************************************************************
  //**********************************************************************
  // update beam contact-specific quantities
  beamcmanager_->Update(*dis_,istep,99);
  //**********************************************************************
  //**********************************************************************
  
  //**********************************************************************
  //**********************************************************************
  // output reaction forces and moments
#ifdef REACTIONFORCES
  beamcmanager_->Reactions(*fint_,*dirichtoggle_,istep);
#endif
  //**********************************************************************
   //**********************************************************************
}

#endif  // #ifdef CCADISCRET
