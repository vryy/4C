/*!----------------------------------------------------------------------
\file contactstrugenalpha_lsnewton.cpp
\brief Line Search Newton for Gen-Alpha time integration with contact

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
 |  do Newton iteration with line search (public)            mwgee 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::FullNewtonLineSearch()
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
  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) && numiter<maxiter)
  {
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

      if (surf_stress_man_!=null)
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
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

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

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

    //-------------------------make contact modifications to lhs and rhs
    {
      contactmanager_->SetState("displacement",disn_);

      contactmanager_->InitializeMortar(numiter+1);
      contactmanager_->EvaluateMortar(numiter+1);

      contactmanager_->Initialize(numiter+1);
      contactmanager_->Evaluate(stiff_,fresm_,numiter+1);
    }

    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    }

    //--------------------------------------------------- contact forces
    contactmanager_->ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
    dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
#endif // #ifdef CONTACTGMSH2

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

        if (surf_stress_man_!=null)
        {
          p.set("surfstr_man", surf_stress_man_);
          surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
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

      // keep a copy of fresm for contact forces / equilibrium check
      RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

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

      //-------------------------make contact modifications to lhs and rhs
      {
        contactmanager_->SetState("displacement",disn_);

        contactmanager_->InitializeMortar(numiter+1);
        contactmanager_->EvaluateMortar(numiter+1);

        contactmanager_->Initialize(numiter+1);
        contactmanager_->Evaluate(stiff_,fresm_,numiter+1);
      }

      // blank residual DOFs that are on Dirichlet BC
      {
        Epetra_Vector fresmdbc(*fresm_);
        fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      }

      //--------------------------------------------------- contact forces
      contactmanager_->ContactForces(fresmcopy);

  #ifdef CONTACTGMSH2
      dserror("Gmsh Output for every iteration only implemented for semi-smooth Newton");
  #endif // #ifdef CONTACTGMSH2

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
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }
  
  params_.set<int>("num iterations",numiter);

  return;
} // ContactStruGenAlpha::FullNewtonLineSearch()


/*----------------------------------------------------------------------*
 |  do semismooth Newton iteration with line search (public)  popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ContactStruGenAlpha::SemiSmoothNewtonLineSearch()
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
  // active set search and geometrical nonlinearity are merged into
  // ONE Newton loop, thus we have to check for convergence of the
  // active set here, too!
  while ((!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) ||
          !contactmanager_->ActiveSetConverged()) && numiter<maxiter)
  {
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

      if (surf_stress_man_!=null)
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
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

    // keep a copy of fresm for contact forces / equilibrium check
    RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

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

    //-------------------------make contact modifications to lhs and rhs
    {
      contactmanager_->SetState("displacement",disn_);

      contactmanager_->InitializeMortar(numiter+1);
      contactmanager_->EvaluateMortar(numiter+1);

      // this is the correct place to update the active set!!!
      // (on the one hand we need the new weighted gap vector g, which is
      // computed in EvaluateMortar() above and on the other hand we want to
      // run the Evaluate()routine below with the NEW active set already)
      contactmanager_->UpdateActiveSetSemiSmooth(disn_);
            
      contactmanager_->Initialize(numiter+1);
      contactmanager_->Evaluate(stiff_,fresm_,numiter+1);
    }

    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmdbc(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
    }

    //--------------------------------------------------- contact forces
    contactmanager_->ContactForces(fresmcopy);

#ifdef CONTACTGMSH2
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    contactmanager_->VisualizeGmsh(istep,numiter+1);
#endif // #ifdef CONTACTGMSH2

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

        if (surf_stress_man_!=null)
        {
          p.set("surfstr_man", surf_stress_man_);
          surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
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

      // keep a copy of fresm for contact forces / equilibrium check
      RCP<Epetra_Vector> fresmcopy= rcp(new Epetra_Vector(*fresm_));

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

      //-------------------------make contact modifications to lhs and rhs
      {
        contactmanager_->SetState("displacement",disn_);

        contactmanager_->InitializeMortar(numiter+1);
        contactmanager_->EvaluateMortar(numiter+1);

        // NO update of the active set here, as this would change the
        // system and thus the residual! During line search the active
        // set has to be kept constant!
        
        contactmanager_->Initialize(numiter+1);
        contactmanager_->Evaluate(stiff_,fresm_,numiter+1);
      }

      // blank residual DOFs that are on Dirichlet BC
      {
        Epetra_Vector fresmdbc(*fresm_);
        fresm_->Multiply(1.0,*invtoggle_,fresmdbc,0.0);
      }

      //--------------------------------------------------- contact forces
      contactmanager_->ContactForces(fresmcopy);

  #ifdef CONTACTGMSH2
      int step  = params_.get<int>("step",0);
      int istep = step + 1;
      contactmanager_->VisualizeGmsh(istep,numiter+1);
  #endif // #ifdef CONTACTGMSH2

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
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // ContactStruGenAlpha::SemiSmoothNewtonLineSearch()


#endif  // #ifdef CCADISCRET
