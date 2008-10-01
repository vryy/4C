/*!----------------------------------------------------------------------
\file statmech.cpp
\brief time integration for structural problems with statistical mechanics

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "statmech_time.H"

#include "../drt_lib/drt_globalproblem.H"

#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
#ifdef D_BEAM2
#include "../drt_beam2/beam2.H"
#endif  // #ifdef D_BEAM2

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 08/08|
 *----------------------------------------------------------------------*/
StatMechTime::StatMechTime(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output)
{
  statmechmanager_ = rcp(new StatMechManager(params,dis));
  return;
} // StatMechTime::StatMechTime


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               cyron 08/08|
 *----------------------------------------------------------------------*/
void StatMechTime::Integrate()
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

  if (equil=="full newton")
  {

  double dt = params_.get<double>("delta time" ,0.01);
  
  int num_dof = (*fext_).GlobalLength();


  Epetra_SerialDenseVector  v0;
  v0.Size(num_dof);
  Epetra_SerialDenseVector d0;
  d0.Size(num_dof);
  Epetra_SerialDenseVector d1_ap;
  d1_ap.Size(num_dof);
  Epetra_SerialDenseVector v1_ap;
  v1_ap.Size(num_dof);
  Epetra_SerialDenseVector fint0;
  fint0.Size(num_dof);
  Epetra_SerialDenseVector relerr_d;
  relerr_d.Size(num_dof);
  Epetra_SerialDenseVector relerr_v;
  relerr_v.Size(num_dof);
  Epetra_SerialDenseVector Delta_d;
  Delta_d.Size(num_dof);
  Epetra_SerialDenseVector Delta_v;
  Delta_v.Size(num_dof);
  int kd = 0;
  double gamma = 0.125663706 / ((num_dof/6) - 1);
  
  {
    //random generator for seeding only (necessary for thermal noise)
    ranlib::Normal<double> seedgenerator(0,1);
    /*by seeding with both current time and processor Id a different random initilization on each processor is made sure;
     * note that each processor will set up its own random generator and they all have to be pairwise independent*/
    int seedvariable = std::time(0)*(discret_.Comm().MyPID() + 1); // +1 in order to make sure random also for 0-th processor!
    //seedvariable = 1;
    seedgenerator.seed((unsigned int)seedvariable);
    
    /*before starting calculations elements have to be initilized in order to make them know viscosity and thermal energy of 
     * their surrounding thermal bath*/
    
    //parameters with respect to statistical mechanics into special variable for later acceess
    Teuchos::ParameterList statisticalparams( DRT::Problem::Instance()->StatisticalMechanicsParams() );
    
    for (int num=0; num<  discret_.NumMyColElements(); ++num)
    {    
      //in case that current element is not a beam3 or beam 2element there is nothing to do and we go back
      //to the head of the loop
      if (discret_.lColElement(num)->Type() != DRT::Element::element_beam3 && discret_.lColElement(num)->Type() != DRT::Element::element_beam2) continue;
      
      #ifdef D_BEAM3
      
      //if we get so far current element is a beam3 or beam2 element and  we get a pointer at it
      if (discret_.lColElement(num)->Type() == DRT::Element::element_beam3)
      {
        DRT::ELEMENTS::Beam3* currele = dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_.lColElement(num));
        if (!currele) dserror("cast to Beam3* failed");
        
        //finally parameters related with thermal bath in beam elements can be set to correct values
        currele->kT_ =  statisticalparams.get<double>("KT",0.0);
        //zeta denotes frictional coefficient per length (approximated by the one for an infinitely long staff)
        currele->zeta_ = 4*PI*currele->lrefe_*statisticalparams.get<double>("ETA",0.0);
        currele->stochasticorder_ = statisticalparams.get<int>("STOCH_ORDER",0); 
      }
      
      #endif  // #ifdef D_BEAM3
      
      #ifdef D_BEAM2
      
      if (discret_.lColElement(num)->Type() == DRT::Element::element_beam2)
      {
        DRT::ELEMENTS::Beam2* currele = dynamic_cast<DRT::ELEMENTS::Beam2*>(discret_.lColElement(num));
        if (!currele) dserror("cast to Beam2* failed");
        
        //finally parameters related with thermal bath in beam elements can be set to correct values
        currele->kT_ =  statisticalparams.get<double>("KT",0.0);
        //zeta denotes frictional coefficient per length (approximated by the one for an infinitely long staff)
        currele->zeta_ = 4*PI*currele->lrefe_*statisticalparams.get<double>("ETA",0.0);
        currele->stochasticorder_ = statisticalparams.get<int>("STOCH_ORDER",0); 
      }
      
      #endif  // #ifdef D_BEAM2
      
    }
  }

  for (int i=step; i<nstep; ++i)
  {
    double time = params_.get<double>("total time",0.0);

    /*
    if      (predictor==1) ConstantPredictor();
    else if (predictor==2) ConsistentPredictor();
    */
    //if(i<10000)
      ConsistentPredictor();
    //else
     // BrownianPredictor3D();
      


    FullNewton();
    UpdateandOutput();
    
    //special update and output for statistical mechanics
    statmechmanager_->StatMechOutput(time,num_dof,i,dt,*dis_);
    statmechmanager_->StatMechUpdate();

    /*
    //Freiheitsgrade längs zur Filamentachse: Da nur geringe axiale Dehnung zu erwarten ist, kann angenommen werden,
    //dass alle Freiheitsgrade in Längsrichtung dieselbe Bewegung Delta_x ausführen, die approximiert werden kann durch:
    // Gamma * Delta_x / dt = fext_axial, wobei fext_axial die Summe der externen Kräfte in Axialrichtung längs des
    //gesamten Filaments ist und Gamma die Gesamtreibung eines Filaments der Länge 10 gegenüber axialer Verschiebung ist
    double fext_axial = 0;
    for(int id = 0; id < num_dof; id = id+6)
    {
      fext_axial += (*fext_)[id];
    }
    //Freiheitsgrade entlang der Filamentachse:
    v1_ap(0) = fext_axial / 0.125663706;
    d1_ap(0) = 0.5*dt*(v0(0) + v1_ap(0)) + d0(0);
    double lrefe = 10.0 / (num_dof/6 - 1);

    for(int jd = 0; jd < (num_dof/6); jd++)
    {
      //Freiheitsgrade entlang der Filamentachse aus Undehnbarkeitsbedingung:
      //v1_ap(jd*6) = fext_axial / 0.125663706;
      //d1_ap(jd*6) = 0.5*dt*(v0(jd*6) + v1_ap(jd*6)) + d0(jd*6);


      //Freiheitsgrade quer zur Filamentachse
      for(int id = 1; id < 3; id++)
      {
        kd = 6*jd + id;
        v1_ap(kd) = ( (*fext_)[kd] - fint0(kd) ) / gamma;
        //an den Randknoten nur jeweils halbes gamma:
        if (jd == 0 || jd == (num_dof/6 -1) )
          v1_ap(kd) = 2*v1_ap(kd);
        d1_ap(kd) = 0.5*dt*(v0(kd) + v1_ap(kd)) + d0(kd);
      }

      if(jd>0)
      {
        double dy = d1_ap(jd*6+1) - d1_ap((jd-1)*6+1);
        double dz = d1_ap(jd*6+2) - d1_ap((jd-1)*6+2);
        d1_ap(jd*6) = pow(lrefe*lrefe - dy*dy - dz*dz  ,0.5) + d1_ap((jd-1)*6) - lrefe;
        v1_ap(jd*6) = 2*(d1_ap(jd*6) - d0(jd*6))/dt - v0(jd*6);
      }

      for(int id = 0; id < 6; id++)
      {
        kd = 6*jd + id;
        //Berechnung des relativen Fehlers im Prädiktorschritt:
        Delta_d(kd) = (*dis_)[kd] - d0(kd);
        Delta_v(kd) = (*velm_)[kd] - v0(kd);
        relerr_d(kd) = ( (d1_ap(kd) - d0(kd) ) - Delta_d(kd) ) / Delta_d(kd);
        relerr_v(kd) = ( (v1_ap(kd) - v0(kd) ) - Delta_v(kd) ) / Delta_v(kd);

        //Zwischenspeichern der Endgrößen im abgeschlossenen Zeitschritt
        d0(kd) = (*dis_)[kd];
        v0(kd) = (*velm_)[kd];
        fint0(kd) = (*fint_)[kd];
      }
    }


    //std::cout<<"\nfext nach update and output"<<*fext_;
    //std::cout<<"\nvelm nach update and output"<<*velm_;
    //std::cout<<"\nfint nach update and output"<<*fint_;
    //std::cout<<"\ndis nach update and output"<< d0;

    //std::cout<<"\n*dis_ nach update and output"<<d0;

    //std::cout<<"\n\nrelerr_d \n" << relerr_d;
    //std::cout<<"\n\nrelerr_v \n" << relerr_v;
    //std::cout<<"\n\nDelta_d \n" << Delta_d;
    //std::cout<<"\n\nDelta_v \n" << Delta_v;

  */

    if (time>=maxtime) break;
  }
#if 0
  for (int i=0; i<discret_.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = discret_.lRowNode(i);
    printf("NODE %d COORD ",actnode->Id()+1);
    for (int j=0; j<discret_.NumDof(actnode); ++j)
    {
      const int gdof = discret_.Dof(actnode,j);
      const int lid  = dis_->Map().LID(gdof);
      printf("%20.15f ",actnode->X()[j]+(*dis_)[lid]);
    }
    printf("\n");
  }
#endif
  }
  else if (equil=="line search newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      LineSearchNewton();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
#if 0
    for (int i=0; i<discret_.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = discret_.lRowNode(i);
      printf("NODE %d COORD ",actnode->Id()+1);
      for (int j=0; j<discret_.NumDof(actnode); ++j)
      {
        const int gdof = discret_.Dof(actnode,j);
        const int lid  = dis_->Map().LID(gdof);
        printf("%20.15f ",actnode->X()[j]+(*dis_)[lid]);
      }
      printf("\n");
    }
#endif
  }
  else if (equil=="modified newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      ModifiedNewton();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else if (equil=="nonlinear cg")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      NonlinearCG();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else if (equil=="ptc")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      PTC();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void StatMechTime::Integrate()

/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void StatMechTime::ConsistentPredictor()
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
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    discret_.SetState("velocity",veln_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
    
    
    
    //declaration of a column and row map Epetra_Vector for evaluation of statistical forces  
    RCP<Epetra_Vector>    fstatcol;
    fstatcol = LINALG::CreateVector(*discret_.DofColMap(),true);
    RCP<Epetra_Vector>    fstatrow;
    fstatrow = LINALG::CreateVector(*discret_.DofRowMap(),true);
    
    //defining parameter list passed down to the elements in order to evalute statistical forces down there
    ParameterList pstat;
    pstat.set("action","calc_stat_forces");
    pstat.set("delta time",dt);
    
    //evaluation of statistical forces on column map vecotor
    discret_.SetState("displacement",dis_); //during evaluation of statistical forces access to current displacement possible
    discret_.Evaluate(pstat,null,null,fstatcol,null,null);
    discret_.ClearState();
    
    /*exporting col map statistical force vector to a row map vector additively, i.e. in such a way that a 
     * vector element with a certain GID in the final row vector is the sum of all elements of the column 
     * vector related to the same GID*/
    Epetra_Export exporter(*discret_.DofColMap(),*discret_.DofRowMap());
    fstatrow->Export(*fstatcol,exporter,Add);
    
    //adding statistical forces to external forces of EvaluateNeumann call in predictor and to residual forces
    fextn_->Update(1.0,*fstatrow,1.0);
    
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

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      p.set("newstep", true);

#ifdef STRUGENALPHA_FINTLIKETR
      p.set("fintliketr", true);
#endif
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

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
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



#endif  // #ifdef CCADISCRET
