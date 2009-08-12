/*!----------------------------------------------------------------------
\file strugenalpha_lsnewton.cpp
\brief Newton iteration with parabolic line search

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>

#include "strugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 06/08|
 *----------------------------------------------------------------------*/
void StruGenAlpha::LineSearchNewton()
{
  const int myrank = discret_.Comm().MyPID();

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  //double mdamp     = params_.get<double>("damping factor M"       ,0.0);
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
  if (!errfile) printerr = false;
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  // check whether mass and damping are present
  // note: the stiffness matrix might be filled already
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  // storage of old values to be bla to repeat a step
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
  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    if (dynkindstat); // do nothing
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
    LineSearchUpdateMidConfiguration(lambda,disi_,disno,dismo,velmo,accmo,
                                     dt,alphaf,alpham,beta,gamma);


    // zerofy velocity and acceleration in case of quasistatics
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
      discret_.SetState("residual displacement",disim);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
#else
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);
#endif
#ifdef STRUGENALPHA_FINTLIKETR
      fintn_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
      discret_.ClearState();
      // do NOT finalize the stiffness matrix to add masses to it later
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
      // Res = M . A_{n+1-alpha_m}
      //     + C . V_{n+1-alpha_f}
      //     + F_int(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}
      // add inertia mid-forces
      mass_->Multiply(false,*accm_,*finert_);
      fresm_->Update(1.0,*finert_,0.0);
      // add viscous mid-forces
      if (damping)
      {
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
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
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
    //if (dynkindstat && numiter==0) dolinesearch = false;
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

      //----------------------------- update with new steplength alpha
      LineSearchUpdateMidConfiguration(lambda,disi_,disno,dismo,velmo,accmo,
                                       dt,alphaf,alpham,beta,gamma);

      //----- zerofy velocity and acceleration in case of quasistatics
      if (dynkindstat)
      {
        velm_->PutScalar(0.0);
        accm_->PutScalar(0.0);
      }

      // scale IncD_{n+1} by (1-alphaf) to obtain
      // mid residual displacements IncD_{n+1-alphaf}
      RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*disi_));
      disim->Scale(lambda*(1.-alphaf));

      //-------------------------------------- compute internal forces
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
        discret_.SetState("residual displacement",disim);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
#else
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);
#endif
#ifdef STRUGENALPHA_FINTLIKETR
        fintn_->PutScalar(0.0);  // initialise internal force vector
        discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
        fint_->PutScalar(0.0);  // initialise internal force vector
        discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
        discret_.ClearState();
        // do NOT finalize the stiffness matrix to add masses to it later
      }

      // ----------------------------------------------compute residual
      {
        if (dynkindstat)
        {
          // static residual
          // Res = F_int - F_ext
          fresm_->PutScalar(0.0);
        }
        else
        {
          // Res = M . A_{n+1-alpha_m}
          //     + C . V_{n+1-alpha_f}
          //     + F_int(D_{n+1-alpha_f})
          //     - F_{ext;n+1-alpha_f}
          // add inertia mid-forces
          mass_->Multiply(false,*accm_,*finert_);
          fresm_->Update(1.0,*finert_,0.0);
          // add viscous mid-forces
          if (damping)
          {
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
        // blank residual DOFs that are on Dirichlet BC
        {
          Epetra_Vector fresmcopy(*fresm_);
          fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
        }

        //------------------------------------------------ recompute norms
        disim->Norm2(&disinorm);
        fresm_->Norm2(&fresmnorm);
      }

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
} // StruGenAlpha::LineSearchNewton()



/*----------------------------------------------------------------------*
 |  update mid configuration (prot.)                         mwgee 06/08|
 *----------------------------------------------------------------------*/
void StruGenAlpha::LineSearchUpdateMidConfiguration(double lambda,
                                                    RCP<Epetra_Vector> disi,
                                                    RCP<Epetra_Vector> disno,
                                                    RCP<Epetra_Vector> dismo,
                                                    RCP<Epetra_Vector> velmo,
                                                    RCP<Epetra_Vector> accmo,
                                                    double dt,
                                                    double alphaf,
                                                    double alpham,
                                                    double beta,
                                                    double gamma)
{
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(lambda,*disi_,1.0,*disno,0.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
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
  return;
}





/*----------------------------------------------------------------------*
 |  constant predictor with disp or arclength control        mwgee 08/09|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ControlledConstantPredictor()
{
  //cout << "ControlledConstantPredictor\n";
#ifdef STRUGENALPHA_FINTLIKETR
  dserror("No control method other than load with STRUGENALPHA_FINTLIKETR");
#endif
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  double mdamp       = params_.get<double>("damping factor M",0.0);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  //int controldof     = params_.get<int>("CONTROLDOF",-1);
  //int controlcurve   = params_.get<int>("CONTROLCURVE",-1);
  INPAR::STR::ControlType controltype = params_.get<INPAR::STR::ControlType>("CONTROLTYPE",INPAR::STR::control_load);
  if (!dynkindstat || controltype==INPAR::STR::control_load)
    dserror("this predictor only for static with control other than load control");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
    CalcRefNorms();

  // increment time and step
  double timen = time + dt;
  //int istep = step + 1; 

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
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  //------------------------------ compute interpolated dis, vel and acc
  // constant predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fext_->PutScalar(0.0);
  //fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);
  fextm_->Update(1.0,*fextn_,0.0);

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
    discret_.SetState("displacement",dism_);
    discret_.SetState("velocity",velm_); // do we need this?
    fint_->PutScalar(0.0);  // initialize internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,stiff_);
    }

    if (pot_man_!=null)
    {
      p.set("pot_man", pot_man_);
      pot_man_->EvaluatePotential(p,dism_,fint_,stiff_);
    }

    if (constrMan_->HaveConstraint())
    {
      ParameterList pcon;
      pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
      constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,stiff_,pcon);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // build residual
  // static residual
  // Res = F_int - F_ext
  // add static mid-balance
  //fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
  if (controltype==INPAR::STR::control_disp)
  {
    fresm_->Update(-1.0,*fint_,lambda_+dlambda_,*fextm_,0.0);
  }
  else if (controltype==INPAR::STR::control_arc1)
  {
    dserror("arclength control not implemented");
  }
  else
    dserror("Only disp control implemented");

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
    PrintControlledPredictor(convcheck, fresmnorm, lambda_+dlambda_);
  }

  return;
} // StruGenAlpha::ControlledConstantPredictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration with disp or arclength control       mwgee 08/09|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ControlledFullNewton()
{
  //cout << "ControlledFullNewton\n";
#ifdef STRUGENALPHA_FINTLIKETR
  dserror("No control method other than load with STRUGENALPHA_FINTLIKETR");
#endif
#ifdef STRUGENALPHA_BE
  dserror("No control method other than load with STRUGENALPHA_BE");
#endif
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
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
  int controldof     = params_.get<int>("CONTROLDOF",-1);
  int controlcurve   = params_.get<int>("CONTROLCURVE",-1);
  int controlowner   = params_.get<int>("CONTROLOWNER",-1);
  INPAR::STR::ControlType controltype = params_.get<INPAR::STR::ControlType>("CONTROLTYPE",INPAR::STR::control_load);
  if (!dynkindstat || controltype==INPAR::STR::control_load)
    dserror("ControlledFullNewton only for static with control other than load control");
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;
  
  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    stiff_->Complete();

    //-----------------------------transform to local coordinate systems
    if (locsysmanager_ != null) 
    {
      locsysmanager_->RotateGlobalToLocal(stiff_,fresm_);
      locsysmanager_->RotateGlobalToLocal(vp_);
    }

    //----------------------- apply dirichlet BCs to system of equations
    // in disp control, currently no nonzero dirichlet allowed!
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);
    LINALG::ApplyDirichlettoSystem(vp_,fextn_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,numiter==0);
    solver_.Solve(stiff_->EpetraMatrix(),vp_,fextn_,true,false);
    solver_.ResetTolerance();

    //----------------------- transform back to global coordinate system
    if (locsysmanager_ != null) 
    {
      locsysmanager_->RotateLocalToGlobal(disi_);
      locsysmanager_->RotateLocalToGlobal(vp_);
    }
    
    //---------------------------------------------------------- control 
    double vbar;
    if (controltype==INPAR::STR::control_disp)
    {
      // see Wriggers: Nichtlinear Finite-Elemente Methoden, p. 156ff.
      
      // get curve value at time timen which is prescribed value vbar_i+1
      vbar = DRT::UTILS::TimeCurveManager::Instance().Curve(controlcurve).f(timen);
      
      // get current increment for v_a from disi_
      double va = 0.0;
      double vp = 0.0;
      double vg = 0.0;
      if (myrank_==controlowner) 
      {
        va = (*disn_)[disn_->Map().LID(controldof)];
        vp = (*vp_)[vp_->Map().LID(controldof)];
        vg = (*disi_)[disi_->Map().LID(controldof)];
      }
      disn_->Comm().Broadcast(&va,1,controlowner);
      vp_->Comm().Broadcast(&vp,1,controlowner);
      disi_->Comm().Broadcast(&vg,1,controlowner);
      
      // constraint f = va - vbar
      double fi  = va-vbar;
      //cout << "fi      " << fi << endl;
      // grad_v(f) = 1
      double ft = 1.0; // not so sure
      // grad_lambda(f) = 0.0
      double flambda = 0.0;
      // increment of lambda
      dlambda_ = -(fi+ft*vg)/(flambda+ft*vp);
      
      // do appropriate scaling of increment
      disi_->Update(dlambda_,*vp_,1.0);
      lambda_ += dlambda_;
    }
    else if (controltype==INPAR::STR::control_arc1)
    {
      dserror("arclength control not implemented");
    } 
    else dserror("Currently, only disp/arc1 control implemented");

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.0,*disi_,1.0);
    disn_->Update(1.0,*disi_,1.0);
    double dis;
    if (myrank_==controlowner) 
      dis = (*disn_)[disn_->Map().LID(controldof)];
    disn_->Comm().Broadcast(&dis,1,controlowner);

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

      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      disi_->Scale(1.-alphaf);

      discret_.SetState("residual displacement",disi_);

      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);

      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();

      if (surf_stress_man_->HaveSurfStress())
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fint_,stiff_);
      }

      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,stiff_);
      }

      if (constrMan_->HaveConstraint())
      {
        ParameterList pcon;
        pcon.set("scaleStiffEntries",1.0/(1.0-alphaf));
        constrMan_->StiffnessAndInternalForces(time+dt,dis_,disn_,fint_,stiff_,pcon);
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
    // static residual
    // Res = F_int - F_ext
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,lambda_,*fextm_,0.0);

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
      PrintControlledNewton(lambda_,dis,
                  printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
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
      PrintControlledNewton(lambda_,0.0,
                  printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StruGenAlpha::ControlledFullNewton()



/*----------------------------------------------------------------------*
 |  print to screen and/or error file                          gee 08/09|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PrintControlledNewton(double lambda, double dis,
                               bool printscreen, bool printerr, bool print_unconv,
                               FILE* errfile, Epetra_Time timer, int numiter,
                               int maxiter, double fresmnorm, double disinorm,
                               string convcheck)
{
  bool relres        = (convcheck == "RelRes_And_AbsDis" or convcheck == "RelRes_Or_AbsDis");
  bool relres_reldis = (convcheck == "RelRes_And_RelDis" or convcheck == "RelRes_Or_RelDis");

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
        printf("numiter %2d scaled res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("numiter %2d scaled res-norm %10.5e scaled dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(stdout);
      }
      else
        {
        printf("numiter %2d absolute res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(stdout);
      }
    }
    if (printerr)
    {
      if (relres)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(errfile);
      }
      else if (relres_reldis)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e scaled dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(errfile);
      }
      else
        {
        fprintf(errfile, "numiter %2d absolute res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e d_bar %10.5e\n",numiter+1, fresmnorm, disinorm,lambda,dis);
        fflush(errfile);
      }
    }
  }
  else
  {
    if (constrMan_->HaveMonitor())
    {
      constrMan_->PrintMonitorValues();
    }
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("Newton iteration converged:\nnumiter  %d scaled res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e time %10.5f\n",
             numiter,fresmnorm,disinorm,lambda,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("Newton iteration converged:\nnumiter  %d scaled res-norm %10.5e scaled dis-norm %10.5e lambda %10.5e time %10.5f\n",
             numiter,fresmnorm,disinorm,lambda,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("Newton iteration converged:\nnumiter  %d absolute res-norm %10.5e absolute dis-norm %10.5e lambda %10.5e time %10.5f\n",
             numiter,fresmnorm,disinorm,lambda,timepernlnsolve);
      fflush(stdout);
    }
  }
}

/*----------------------------------------------------------------------*
 |  print to screen                                            gee 08/09|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PrintControlledPredictor(string convcheck, double fresmnorm,
                                            double lambda)
{
  if (convcheck != "AbsRes_And_AbsDis" and convcheck != "AbsRes_Or_AbsDis")
  {
    fresmnorm /= ref_fnorm_;
    cout << "Predictor    scaled res-norm " << fresmnorm;// << endl;
  }
  else
  {
    cout << "Predictor  absolute res-norm " << fresmnorm;// << endl;
  }
  printf(" Lambda %10.5e\n",lambda);
  fflush(stdout);
}


#endif  // #ifdef CCADISCRET
