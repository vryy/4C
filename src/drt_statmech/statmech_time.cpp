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

#include <Teuchos_Time.hpp>
#include <iostream>
#include <iomanip>
#include "statmech_time.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_io/io_control.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_inpar/inpar_contact.H"

#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
#ifdef D_BEAM3II
#include "../drt_beam3ii/beam3ii.H"
#endif  // #ifdef D_BEAM3II
#ifdef D_BEAM2
#include "../drt_beam2/beam2.H"
#endif  // #ifdef D_BEAM2
#ifdef D_BEAM2R
#include "../drt_beam2r/beam2r.H"
#endif  // #ifdef D_BEAM2R
#ifdef D_TRUSS3
#include "../drt_truss3/truss3.H"
#endif  // #ifdef D_TRUSS3
#ifdef D_TRUSS2
#include "../drt_truss2/truss2.H"
#endif  // #ifdef D_TRUSS2


//#define MEASURETIME
/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 08/08|
 *----------------------------------------------------------------------*/
StatMechTime::StatMechTime(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output),
isconverged_(0),
uzawaiter_(0)
{
	Teuchos::RCP<LINALG::SparseMatrix> stiff = SystemMatrix();
  statmechmanager_ = rcp(new StatMechManager(params,dis));

  //maximal number of random numbers to be generated per time step for any column map element of this processor
  int randomnumbersperlocalelement = 0;

  /*check maximal number of nodes of an element with stochastic forces on this processor*/
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {
    const DRT::ElementType & eot = dis.lColElement(i)->ElementType();
    /*stochastic forces implemented so far only for the following elements:*/
#ifdef D_BEAM3
    if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(i))->HowManyRandomNumbersINeed());

        //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
        if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
          statmechmanager_->PeriodicBoundaryBeam3Init(dis.lColElement(i));
      }
    else
#endif  // #ifdef D_BEAM3
#ifdef D_BEAM3II
      if ( eot == DRT::ELEMENTS::Beam3iiType::Instance() )
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3ii*>(dis.lColElement(i))->HowManyRandomNumbersINeed());

        //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
        if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
          statmechmanager_->PeriodicBoundaryBeam3iiInit(dis.lColElement(i));
      }
    else
#endif  // #ifdef D_BEAM3II
#ifdef D_BEAM2
      if ( eot == DRT::ELEMENTS::Beam2Type::Instance() )
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam2*>(dis.lColElement(i))->HowManyRandomNumbersINeed());
      }
    else
#endif  // #ifdef D_BEAM2
#ifdef D_BEAM2R
      if ( eot == DRT::ELEMENTS::Beam2rType::Instance() )
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam2r*>(dis.lColElement(i))->HowManyRandomNumbersINeed());
      }
    else
#endif  // #ifdef D_BEAM2R
#ifdef D_TRUSS3
      if ( eot == DRT::ELEMENTS::Truss3Type::Instance() )
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Truss3*>(dis.lColElement(i))->HowManyRandomNumbersINeed());

        //in case of periodic boundary conditions truss3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
        if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
          statmechmanager_->PeriodicBoundaryTruss3Init(dis.lColElement(i));
      }
      else
#endif  // #ifdef D_TRUSS3
        continue;
  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  /*so far the maximal number of random numbers required per element has been checked only locally on this processor;
   *now we compare the results of each processor and store the maximal one in maxrandomnumbersperglobalelement_*/
  dis.Comm().MaxAll(&randomnumbersperlocalelement,&maxrandomnumbersperglobalelement_ ,1);

  //suppress all output printed to screen in case of single filament studies in order not to generate too much output on the cluster
  if( Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechmanager_->statmechparams_, "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_endtoendlog ||
      Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechmanager_->statmechparams_, "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_endtoendconst ||
      Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechmanager_->statmechparams_, "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_orientationcorrelation ||
      Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechmanager_->statmechparams_, "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_anisotropic)
  {
    params_.set("print to screen",false);
    std::cout<<"\n\nPay Attention: from no on regular output to screen suppressed !!!\n\n";
  }


  //in case that beam contact is activated by respective input parameter, a Beam3cmanager object is created
  if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
  {
    //check wheter appropriate parameters are set in the parameter list "CONTACT & MESHTYING"
    const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->MeshtyingAndContactParams();
    if (Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION") == INPAR::CONTACT::app_beamcontact)
      beamcmanager_ = rcp(new CONTACT::Beam3cmanager(dis));
    else
      dserror("beam contact switched on in parameter list STATISTICAL MECHANICS, but not in in parameter list MESHTYING AND CONTACT!!!");
  }


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


  double dt = params_.get<double>("delta time" ,0.01);

  //getting number of dimensions for diffusion coefficient calculation
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  int ndim=psize.get<int>("DIM");

  // debug cout
	Epetra_Vector discol(*discret_.DofColMap(), true);

  //defining solution strategz for beam contact
	INPAR::CONTACT::SolvingStrategy soltype;
  if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
    soltype = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(beamcmanager_->InputParameters(),"STRATEGY");

  for (int i=step; i<nstep; ++i)
  {
    /*in the very first step and in case that special output for statistical mechanics is requested we have
     * to initialize the related output method*/
    if(i == 0)
    {
      statmechmanager_->InitOutput(ndim,dt);
      // handling gmsh output seperately
      if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"GMSHOUTPUT"))
    	{
    		std::ostringstream filename;
    			filename << "./GmshOutput/networkInit.pos";
    		//calling method for writing Gmsh output
    		statmechmanager_->GmshOutput(*dis_,filename,step);
    	}
    }

    //time_ is time at the end of this time step
    double time = params_.get<double>("total time",0.0);
    if(time + statmechmanager_->statmechparams_.get<double>("DELTA_T_NEW",dt) > statmechmanager_->statmechparams_.get<double>("STARTTIMEACT", 0.0))
    {
      if(statmechmanager_->statmechparams_.get<double>("DELTA_T_NEW",dt)>0.0)
      {
        dt = statmechmanager_->statmechparams_.get<double>("DELTA_T_NEW",dt);
        params_.set("delta time", dt);
      }
    }
    statmechmanager_->time_ = time + dt;

    //save relevant class variables at the beginning of this time step
    statmechmanager_->WriteConv();

    //seed random generators of statmechmanager_ to generate the same random numbers even if the simulation was interrupted by a restart
    statmechmanager_->SeedRandomGenerators(i);


    //processor 0 indicates beginning of new time step on console
    if(!discret_.Comm().MyPID() && params_.get<bool>  ("print to screen",true))
      std::cout<<"\nbegin time step "<<i+1<<":";

    //redo time step in case of bad random configuration
    do
    {
      //set and delete crosslinkers compared to converged configuration of last time step
      const double t_admin = Teuchos::Time::wallTime();
      statmechmanager_->Update(i, dt, *dis_, stiff_,ndim);

      //processor 0 write total number of elements at the beginning of time step i to console as well as how often a time step had to be restarted due to bad random numbers
      if(!discret_.Comm().MyPID() && params_.get<bool>  ("print to screen",true))
      {
        std::cout<<"\ntime for update of crosslinkers: " << Teuchos::Time::wallTime() - t_admin<< " seconds";
        std::cout<<"\nTotal number of elements after crosslinker update: "<<discret_.NumGlobalElements();
        std::cout<<"\nNumber of unconverged steps since simulation start: "<<statmechmanager_->unconvergedsteps_<<"\n";
      }

      //assuming that iterations will converge
      isconverged_ = 1;

      /*multivector for stochastic forces evaluated by each element; the numbers of vectors in the multivector equals the maximal
       *number of random numbers required by any element in the discretization per time step; therefore this multivector is suitable
       *for synchrinisation of these random numbers in parallel computing*/
      RCP<Epetra_MultiVector> randomnumbers = rcp( new Epetra_MultiVector(*(discret_.ElementColMap()),maxrandomnumbersperglobalelement_) );

      //pay attention: for a constant predictor an incremental velocity update is necessary, which has
      //been deleted out of the code in oder to simplify it

      //generate gaussian random numbers for parallel use with mean value 0 and standard deviation (2KT / dt)^0.5
      statmechmanager_->GenerateGaussianRandomNumbers(randomnumbers,0,pow(2.0 * (statmechmanager_->statmechparams_).get<double>("KT",0.0) / dt,0.5));



      //in case that beam contact is activated special solution strategies are required
      if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
      {
        switch (soltype)
        {
          //solving strategy using regularization with penalty method (nonlinear solution approach: ordinary NEWTON)
          case INPAR::CONTACT::solution_penalty:
          {
            ConsistentPredictor(randomnumbers);

            if(ndim ==3)
              PTC(randomnumbers);
            else
              FullNewton(randomnumbers);

            // update constraint norm
            beamcmanager_->UpdateConstrNorm();

            break;
          }
          //solving strategy using regularization with augmented Lagrange method (nonlinear solution approach: nested UZAWA NEWTON)
          case INPAR::CONTACT::solution_auglag:
          {
            // Initialize all lmuzawa to zero at beginning of new time step
            beamcmanager_->ResetAlllmuzawa();

            ConsistentPredictor(randomnumbers);

            // get tolerance and maximum number of Uzawa steps from input file
            double eps = beamcmanager_->InputParameters().get<double>("UZAWACONSTRTOL");
            int maxuzawaiter = beamcmanager_->InputParameters().get<int>("UZAWAMAXSTEPS");

            // LOOP2: augmented Lagrangian (Uzawa)
            uzawaiter_=0;
            do
            {
              // increase iteration index
              ++uzawaiter_;
              if (uzawaiter_ > maxuzawaiter)
                dserror("Uzawa unconverged in %d iterations",maxuzawaiter);

              if (discret_.Comm().MyPID() == 0)
                cout << endl << "Starting Uzawa step No. " << uzawaiter_ << endl;

              // LOOP3: nonlinear iteration (Newton)
              if(ndim ==3)
                PTC(randomnumbers);
              else
                FullNewton(randomnumbers);

              // update constraint norm and penalty parameter
              beamcmanager_->UpdateConstrNorm(uzawaiter_);

              // update Uzawa Lagrange multipliers
              beamcmanager_->UpdateAlllmuzawa(uzawaiter_);

            } while (abs(beamcmanager_->GetConstrNorm()) >= eps);

            // reset penalty parameter
            beamcmanager_->ResetCurrentpp();

            break;
          }
          default:
            dserror("Only penaltz and augmented Lagrange implemented in statmech_time.cpp for beam contact");
          }
      }
      else
      {

        ConsistentPredictor(randomnumbers);

        if(ndim ==3)
          PTC(randomnumbers);
        else
          FullNewton(randomnumbers);
      }


      /*if iterations have not converged a new trial requires setting all intern element variables, statmechmanager class variables
       *and the state of the discretization to status at the beginning of this time step*/
      if(isconverged_ == 0)
      {
        ParameterList p;
        p.set("action","calc_struct_reset_istep");
        discret_.Evaluate(p,null,null,null,null,null);
        statmechmanager_->RestoreConv(stiff_);
      }

    }
    while(isconverged_ == 0);

    //periodic shift of configuration at the end of the time step in order to avoid improper output
    statmechmanager_->PeriodicBoundaryShift(*dism_, ndim, dt);

    UpdateandOutput();

    //**********************************************************************
    //**********************************************************************
    // update beam contact-specific quantities
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
    {
      beamcmanager_->Update(*dis_,params_.get<int>("step",0),99,99);
      // output reaction forces and moments
      #ifdef REACTIONFORCES
      beamcmanager_->Reactions(*fint_,*dirichtoggle_,params_.get<int>("step",0));
      #endif
    }
    //**********************************************************************
     //**********************************************************************

    //special output for statistical mechanics
    statmechmanager_->Output(params_,ndim,time,i,dt,*dis_,*fint_);

    if (time>=maxtime) break;
  }


  return;
} // void StatMechTime::Integrate()


/*----------------------------------------------------------------------*
 |do consistent predictor step for Brownian dynamics (public)cyron 10/09|
 *----------------------------------------------------------------------*/
void StatMechTime::ConsistentPredictor(RCP<Epetra_MultiVector> randomnumbers)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (!firststep_ and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
    CalcRefNorms();

  // increment time and step
  double timen = time + dt;  // t_{n+1}
  //int istep = step + 1;  // n+1

  /*special part for STATMECH: initialize disn_ and veln_ with zero; this is necessary only for the following case:
   * assume that an iteration step did not converge and is repeated with new random numbers; if the failure of conver
   * gence lead to disn_ = NaN and veln_ = NaN this would affect also the next trial as e.g. disn_->Update(1.0,*dis_,0.0);
   * would set disn_ to NaN as even 0*NaN = NaN!; this would defeat the purpose of the repeated iterations with new
   * random numbers and has thus to be avoided; therefore we initialized disn_ and veln_ with zero which has no effect
   * in any other case*/
  disn_->PutScalar(0.0);
  veln_->PutScalar(0.0);
  dism_->PutScalar(0.0);
  velm_->PutScalar(0.0);
  fresm_->PutScalar(0.0);

  //consistent predictor for backward Euler time integration scheme
  disn_->Update(1.0,*dis_,0.0);
  veln_->Update(1.0/dt,*disn_,-1.0/dt,*dis_,0.0);

  //evaluate deterministic external forces
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
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements

    // determine evaluation mode
    if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) <= 0.0 &&
    	 Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"PERIODICDBC"))
    	dserror("Set PeriodLength > 0.0 if periodic DBCs are to be applied");
    if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0 &&
    	 Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"CONVENTIONALDBC"))
    	dserror("Set PeriodLength to Zero if conventional DBCs are to be applied");
    // in case of activated periodic boundary conditions
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"PERIODICDBC"))
    	statmechmanager_->EvaluateDirichletPeriodic(p, disn_, dirichtoggle_, invtoggle_);
		// "common" case without periodic boundary conditions
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"CONVENTIONALDBC") )
    	discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);


    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    discret_.SetState("velocity",veln_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }


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


  //------------------------------ compute interpolated dis, vel and acc
  // consistent predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
  velm_->Update(1.-alphaf,*veln_,alphaf,*vel_,0.0);


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


    //passing statistical mechanics parameters to elements
    p.set("ETA",(statmechmanager_->statmechparams_).get<double>("ETA",0.0));
    p.set("THERMALBATH",Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmechmanager_->statmechparams_,"THERMALBATH"));
    p.set("FRICTION_MODEL",Teuchos::getIntegralValue<INPAR::STATMECH::FrictionModel>(statmechmanager_->statmechparams_,"FRICTION_MODEL"));
    p.set("RandomNumbers",randomnumbers);
    p.set("SHEARAMPLITUDE",(statmechmanager_->statmechparams_).get<double>("SHEARAMPLITUDE",0.0));
    p.set("CURVENUMBER",(statmechmanager_->statmechparams_).get<int>("CURVENUMBER",-1));
    p.set("OSCILLDIR",(statmechmanager_->statmechparams_).get<int>("OSCILLDIR",-1));
    p.set("STARTTIMEACT",(statmechmanager_->statmechparams_).get<double>("STARTTIMEACT",0.0));
    p.set("DELTA_T_NEW",(statmechmanager_->statmechparams_).get<double>("DELTA_T_NEW",0.0));
    p.set("PeriodLength",(statmechmanager_->statmechparams_).get<double>("PeriodLength",0.0));


    // set vector values needed by elements
    discret_.ClearState();
    disi_->PutScalar(0.0);
    discret_.SetState("residual displacement",disi_);

    discret_.SetState("displacement",dism_);
    discret_.SetState("velocity",velm_);

    //discret_.SetState("velocity",velm_); // not used at the moment

    fint_->PutScalar(0.0);  // initialise internal force vector

    p.set("action","calc_struct_nlnstiff");
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();
  }

  //-------------------------------------------- compute residual forces
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);

  //**********************************************************************
  //**********************************************************************
  if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
  {
    // evaluate beam contact
    beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);

#ifdef GMSHNEWTONSTEPS
    // create gmsh-output to visualize predictor step
    int step  = params_.get<int>("step",0);
    int istep = step + 1;
    beamcmanager_->GmshOutput(*disn_,istep,0,0);
    beamcmanager_->ConsoleOutput();
#endif
  }
  //**********************************************************************
  //**********************************************************************

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
} //StatMechTime::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StatMechTime::FullNewton(RCP<Epetra_MultiVector> randomnumbers)
{
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
  if (!errfile) printerr = false;

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);


#ifndef STRUGENALPHA_BE
  //double delta = beta;
#endif

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

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {


    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    //stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(delta*dt),1.0);
    //stiff_->Complete();

    //backward Euler
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
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();


    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}

    dism_->Update(1.-alphaf,*disi_,1.0);
    disn_->Update(1.0,*disi_,1.0);

    // velocities

    // incremental (required for constant predictor)

    //backward Euler
    velm_->Update(1.0/dt,*dism_,-1.0/dt,*dis_,0.0);

    //velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    //velm_->Update((delta-(1.0-alphaf)*gamma)/delta,*vel_,gamma/(delta*dt));


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

      //passing statistical mechanics parameters to elements
      p.set("ETA",(statmechmanager_->statmechparams_).get<double>("ETA",0.0));
      p.set("THERMALBATH",Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmechmanager_->statmechparams_,"THERMALBATH"));
      p.set("FRICTION_MODEL",Teuchos::getIntegralValue<INPAR::STATMECH::FrictionModel>(statmechmanager_->statmechparams_,"FRICTION_MODEL"));
      p.set("RandomNumbers",randomnumbers);
      p.set("SHEARAMPLITUDE",(statmechmanager_->statmechparams_).get<double>("SHEARAMPLITUDE",0.0));
      p.set("CURVENUMBER",(statmechmanager_->statmechparams_).get<int>("CURVENUMBER",-1));
      p.set("STARTTIMEACT",(statmechmanager_->statmechparams_).get<double>("STARTTIMEACT",0.0));
      p.set("DELTA_T_NEW",(statmechmanager_->statmechparams_).get<double>("DELTA_T_NEW",0.0));
      p.set("OSCILLDIR",(statmechmanager_->statmechparams_).get<int>("OSCILLDIR",-1));
      p.set("PeriodLength",(statmechmanager_->statmechparams_).get<double>("PeriodLength",0.0));


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

      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces

    // dynamic residual
    // Res =  C . V_{n+1-alpha_f}
    //        + F_int(D_{n+1-alpha_f})
    //        - F_{ext;n+1-alpha_f}
    // add mid-inertial force


    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);

    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
    {
      beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);

#ifdef GMSHNEWTONSTEPS
      // Create gmsh-output to visualize every step of newton iteration
      int step  = params_.get<int>("step",0);
      int istep = step + 1;
      beamcmanager_->GmshOutput(*disn_,istep,uzawaiter_,numiter+1);
      beamcmanager_->ConsoleOutput();
#endif
    }
    //**********************************************************************
    //**********************************************************************


    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);
    fresm_->Norm2(&fresmnorm);



    //if code is compiled with DEBUG flag each iteration is written into file for Gmsh visualization
#ifdef DEBUG
    // first index = time step index
    std::ostringstream filename;

    //creating complete file name dependent on step number with 5 digits and leading zeros
    if (numiter<100000)
      filename << "./GmshOutput/konvergenz"<< std::setw(5) << setfill('0') << numiter <<".pos";
    else
      dserror("Gmsh output implemented for a maximum of 99999 steps");

    //statmechmanager_->GmshOutput(*dism_,filename,numiter);
#endif  // #ifdef DEBUG


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
  //if on convergence arises within maxiter iterations the time step is restarted with new random numbers
  if (numiter>=maxiter)
  {
    isconverged_ = 0;
    statmechmanager_->unconvergedsteps_++;
    if(discret_.Comm().MyPID() == 0)
      std::cout<<"\n\niteration unconverged - new trial with new random numbers!\n\n";
     //dserror("PTC unconverged in %d iterations",numiter);
  }
  else if(!myrank_ and printscreen)
  {
    PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                fresmnorm,disinorm,convcheck);
  }


  params_.set<int>("num iterations",numiter);

  return;
} // StatMechTime::FullNewton()

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             cyron 12/10|
 *----------------------------------------------------------------------*/
void StatMechTime::PTC(RCP<Epetra_MultiVector> randomnumbers)
{
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

  #ifndef STRUGENALPHA_BE
    //double delta = beta;
  #endif


  double sumsolver     = 0;
  double sumevaluation = 0;
  double sumptc = 0;
  const double tbegin = Teuchos::Time::wallTime();

  if (!errfile) printerr = false;
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  const bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  if (dynkindstat) dserror("Static case not implemented");


  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  InitializeNewtonUzawa();

  // hard wired ptc parameters
  double ctransptc = (statmechmanager_->statmechparams_).get<double>("CTRANSPTC0",0.0);
  double crotptc   = (statmechmanager_->statmechparams_).get<double>("CROTPTC0",0.145);
  double alphaptc  = (statmechmanager_->statmechparams_).get<double>("ALPHAPTC",6.0);

  double nc;
  fresm_->NormInf(&nc);
  RCP<Epetra_Vector> x0 = rcp(new Epetra_Vector(*disi_));

  double resinit = nc;


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
    // debug cout
  	Epetra_Vector discol(*discret_.DofColMap(), true);

    RCP<Epetra_Vector> xm = rcp(new Epetra_Vector(*x0));
    x0->Update(1.0,*disi_,0.0);

    //backward Euler
    stiff_->Complete();

    //the following part was especially introduced for Brownian dynamics
    {
      const double t_ptc = Teuchos::Time::wallTime();

      // create the parameters for the discretization
      ParameterList p;

      p.set("action","calc_struct_ptcstiff");
      p.set("delta time",dt);
      p.set("crotptc",crotptc);
      p.set("ctransptc",ctransptc);

      //add statistical vector to parameter list for statistical forces and damping matrix computation
      p.set("ETA",(statmechmanager_->statmechparams_).get<double>("ETA",0.0));
      p.set("THERMALBATH",Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmechmanager_->statmechparams_,"THERMALBATH"));
      p.set("FRICTION_MODEL",Teuchos::getIntegralValue<INPAR::STATMECH::FrictionModel>(statmechmanager_->statmechparams_,"FRICTION_MODEL"));
      p.set("SHEARAMPLITUDE",(statmechmanager_->statmechparams_).get<double>("SHEARAMPLITUDE",0.0));
      p.set("CURVENUMBER",(statmechmanager_->statmechparams_).get<int>("CURVENUMBER",-1));
      p.set("OSCILLDIR",(statmechmanager_->statmechparams_).get<int>("OSCILLDIR",-1));
      p.set("PeriodLength",(statmechmanager_->statmechparams_).get<double>("PeriodLength",0.0));

      //evaluate ptc stiffness contribution in all the elements
      discret_.Evaluate(p,stiff_,null,null,null,null);
      sumptc += Teuchos::Time::wallTime() - t_ptc;
    }


    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    const double t_solver = Teuchos::Time::wallTime();
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraOperator(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    sumsolver += Teuchos::Time::wallTime() - t_solver;



    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}

    dism_->Update(1.-alphaf,*disi_,1.0);
    disn_->Update(1.0,*disi_,1.0);

    // velocities

    //backward Euler
    // incremental (required for constant predictor)
    velm_->Update(1.0/dt,*dism_,-1.0/dt,*dis_,0.0);
    //velm_->Update((delta-(1.0-alphaf)*gamma)/delta,*vel_,gamma/(delta*dt));



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

      //passing statistical mechanics parameters to elements
      p.set("ETA",(statmechmanager_->statmechparams_).get<double>("ETA",0.0));
      p.set("THERMALBATH",Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmechmanager_->statmechparams_,"THERMALBATH"));
      p.set("FRICTION_MODEL",Teuchos::getIntegralValue<INPAR::STATMECH::FrictionModel>(statmechmanager_->statmechparams_,"FRICTION_MODEL"));
      p.set("RandomNumbers",randomnumbers);
      p.set("SHEARAMPLITUDE",(statmechmanager_->statmechparams_).get<double>("SHEARAMPLITUDE",0.0));
      p.set("CURVENUMBER",(statmechmanager_->statmechparams_).get<int>("CURVENUMBER",-1));
      p.set("STARTTIMEACT",(statmechmanager_->statmechparams_).get<double>("STARTTIMEACT",0.0));
      p.set("DELTA_T_NEW",(statmechmanager_->statmechparams_).get<double>("DELTA_T_NEW",0.0));
      p.set("OSCILLDIR",(statmechmanager_->statmechparams_).get<int>("OSCILLDIR",-1));
      p.set("PeriodLength",(statmechmanager_->statmechparams_).get<double>("PeriodLength",0.0));

      // set vector values needed by elements
      discret_.ClearState();

      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      disi_->Scale(1.-alphaf);

      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);

      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector

      const double t_evaluate = Teuchos::Time::wallTime();

      discret_.Evaluate(p,stiff_,null,fint_,null,null);

      sumevaluation += Teuchos::Time::wallTime() - t_evaluate;

      discret_.ClearState();

    }

    //------------------------------------------ compute residual forces

    // dynamic residual
    // Res =  C . V_{n+1-alpha_f}
    //        + F_int(D_{n+1-alpha_f})
    //        - F_{ext;n+1-alpha_f}
    // add mid-inertial force


    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);

    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
    {
      beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);

#ifdef GMSHNEWTONSTEPS
      // Create gmsh-output to visualize every step of newton iteration
      int step  = params_.get<int>("step",0);
      int istep = step + 1;
      beamcmanager_->GmshOutput(*disn_,istep,uzawaiter_,numiter+1);
      beamcmanager_->ConsoleOutput();
#endif
    }
    //**********************************************************************
    //**********************************************************************


    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // compute inf norm of residual
    double np;
    fresm_->NormInf(&np);

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    if (!myrank_ and (printscreen or printerr))
    {
      PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,crotptc);
    }

    //------------------------------------ PTC update of artificial time
    // SER step size control
    crotptc *= pow((np/nc),alphaptc);
    ctransptc *= pow((np/nc),alphaptc);
    nc = np;


    //Modifikation: sobald Residuum klein, PTC ausgeschaltet
    if(np < 0.001*resinit || numiter > 5) //
    {
      ctransptc = 0.0;
      crotptc = 0.0;
    }


    //--------------------------------- increment equilibrium loop index
    ++numiter;
  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  //if on convergence arises within maxiter iterations the time step is restarted with new random numbers
  if (numiter>=maxiter)
  {
    isconverged_ = 0;
    statmechmanager_->unconvergedsteps_++;
    if(discret_.Comm().MyPID()==0 and printscreen)
      std::cout<<"\n\niteration unconverged - new trial with new random numbers!\n\n";
     //dserror("FullNewton unconverged in %d iterations",numiter);
  }
  else
  {
     if (!myrank_ and printscreen)
     {
       PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                   fresmnorm,disinorm,convcheck,crotptc);
     }
  }

  params_.set<int>("num iterations",numiter);

  if(!discret_.Comm().MyPID() and printscreen)
  std::cout << "\n***\nevaluation time: " << sumevaluation<< " seconds\nptc time: "<< sumptc <<" seconds\nsolver time: "<< sumsolver <<" seconds\ntotal solution time: "<<Teuchos::Time::wallTime() - tbegin<<" seconds\n***\n";

  return;
} // StatMechTime::PTC()

/*----------------------------------------------------------------------*
 |  initialize Newton for 2nd, 3rd, ... Uzawa iteration      cyron 12/10|
 *----------------------------------------------------------------------*/
void StatMechTime::InitializeNewtonUzawa()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  bool  loadlin    = params_.get<bool>("LOADLIN",false);

  // create out-of-balance force for 2nd, 3rd, ... Uzawa iteration
  if (uzawaiter_>1)
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

      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      disi_->Scale(1.-alphaf);

      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);

      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);

      discret_.ClearState();
    }

    //------------------------------------------ compute residual forces
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);


    //**********************************************************************
    //**********************************************************************
    // evaluate beam contact
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"BEAMCONTACT"))
      beamcmanager_->Evaluate(*SystemMatrix(),*fresm_,*disn_,alphaf);
    //**********************************************************************
    //**********************************************************************

    // blank residual DOFs that are on Dirichlet BC
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

  }

  return;
}


/*----------------------------------------------------------------------*
 |  do output including statistical mechanics data(public)    cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechTime::Output()
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
  INPAR::STR::StrainType iostrain      = params_.get<INPAR::STR::StrainType>("io structural strain",INPAR::STR::strain_none);
  bool   iosurfactant  = params_.get<bool>  ("io surfactant"          ,false);

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

#ifdef INVERSEDESIGNCREATE // indicate that this restart is from INVERSEDESIGCREATE phase
    output_.WriteInt("InverseDesignRestartFlag",0);
#endif
#ifdef INVERSEDESIGNUSE // indicate that this restart is from INVERSEDESIGNUSE phase
    output_.WriteInt("InverseDesignRestartFlag",1);
#endif

    isdatawritten = true;

//____________________________________________________________________________________________________________
//note:the following block is the only difference to Output() in strugenalpha.cpp-----------------------------
/* write restart information for statistical mechanics problems; all the information is saved as class variables
 * of StatMechManager*/
    statmechmanager_->WriteRestart(output_);
//------------------------------------------------------------------------------------------------------------
//____________________________________________________________________________________________________________

    if (surf_stress_man_->HaveSurfStress())
      surf_stress_man_->WriteRestart(istep, timen);

    if (constrMan_->HaveConstraint())
    {
      output_.WriteDouble("uzawaparameter",constrSolv_->GetUzawaParameter());
      output_.WriteVector("lagrmultiplier",constrMan_->GetLagrMultVector());
      output_.WriteVector("refconval",constrMan_->GetRefBaseValues());
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

    if (surf_stress_man_->HaveSurfStress() and iosurfactant)
      surf_stress_man_->WriteResults(istep,timen);

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
    p.set("iostress", iostress);
    p.set("strain", strain);
    p.set("iostrain", iostrain);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    discret_.SetState("velocity",vel_);
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
    if (!isdatawritten) output_.NewStep(istep, timen);
    isdatawritten = true;

    switch (iostress)
    {
    case INPAR::STR::stress_cauchy:
      output_.WriteVector("gauss_cauchy_stresses_xyz",*stress,*discret_.ElementRowMap());
      break;
    case INPAR::STR::stress_2pk:
      output_.WriteVector("gauss_2PK_stresses_xyz",*stress,*discret_.ElementRowMap());
      break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror ("requested stress type not supported");
    }

    switch (iostrain)
    {
    case INPAR::STR::strain_ea:
      output_.WriteVector("gauss_EA_strains_xyz",*strain,*discret_.ElementRowMap());
      break;
    case INPAR::STR::strain_gl:
      output_.WriteVector("gauss_GL_strains_xyz",*strain,*discret_.ElementRowMap());
      break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not supported");
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
}//StatMechTime::Output()

/*----------------------------------------------------------------------*
 |  read restart (public)                                    cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechTime::ReadRestart(int step)
{
  RCP<DRT::Discretization> rcpdiscret = rcp(&discret_,false);
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
  statmechmanager_->ReadRestart(reader);

#ifdef INVERSEDESIGNUSE
  int idrestart = -1;
  idrestart = reader.ReadInt("InverseDesignRestartFlag");
  if (idrestart==-1) dserror("expected inverse design restart flag not on file");
  // if idrestart==0 then the file is from a INVERSEDESIGCREATE phase
  // and we have to zero out the inverse design displacements.
  // The stored reference configuration is on record at the element level
  if (!idrestart)
  {
    dis_->PutScalar(0.0);
    vel_->PutScalar(0.0);
    acc_->PutScalar(0.0);
  }
#endif

  // override current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  if (surf_stress_man_->HaveSurfStress())
    surf_stress_man_->ReadRestart(rstep, DRT::Problem::Instance()->InputControlFile()->FileName());

  if (constrMan_->HaveConstraint())
  {
    double uzawatemp = reader.ReadDouble("uzawaparameter");
    constrSolv_->SetUzawaParameter(uzawatemp);
    RCP<Epetra_Map> constrmap=constrMan_->GetConstraintMap();
    RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*constrmap,true);
    reader.ReadVector(tempvec, "lagrmultiplier");
    constrMan_->SetLagrMultVector(tempvec);
    reader.ReadVector(tempvec, "refconval");
    constrMan_->SetRefBaseValues(tempvec,time);
  }

  return;
}//StatMechTime::ReadRestart()

#endif  // #ifdef CCADISCRET
