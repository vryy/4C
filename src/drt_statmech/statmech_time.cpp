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

#include "statmech_time.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_io/io_control.H"

#include <random/normal.h>

#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
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


#define MEASURETIME
/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 08/08|
 *----------------------------------------------------------------------*/
StatMechTime::StatMechTime(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output),
isconverged_(0),
unconvergedsteps_(0),
isinit_(false),
timecurve_(DRT::Problem::Instance()->Curve(0)),
dbcswitch_(LINALG::CreateVector(*(discret_.DofRowMap()),true)),
drefnew_(LINALG::CreateVector(*(discret_.DofRowMap()),true))
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff = SystemMatrix();
  statmechmanager_ = rcp(new StatMechManager(params,dis));

  //maximal number of random numbers to be generated per time step for any column map element of this processor
  int randomnumbersperlocalelement = 0;

  /*check maximal number of nodes of an element with stochastic forces on this processor*/
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {
    /*stochastic forces implemented so far only for the following elements:*/
    switch(dis.lColElement(i)->Type())
    {
#ifdef D_BEAM3
      case DRT::Element::element_beam3:
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(i))->HowManyRandomNumbersINeed());

        //in case of periodic boundary conditions beam3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
        if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
          statmechmanager_->PeriodicBoundaryBeam3Init(dis.lColElement(i));
        break;
      }
#endif  // #ifdef D_BEAM3
#ifdef D_BEAM2
      case DRT::Element::element_beam2:
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam2*>(dis.lColElement(i))->HowManyRandomNumbersINeed());
        break;
      }
#endif  // #ifdef D_BEAM2
#ifdef D_BEAM2R
      case DRT::Element::element_beam2r:
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Beam2r*>(dis.lColElement(i))->HowManyRandomNumbersINeed());
        break;
      }
#endif  // #ifdef D_BEAM2R
#ifdef D_TRUSS3
      case DRT::Element::element_truss3:
      {
        //see whether current element needs more random numbers per time step than any other before
        randomnumbersperlocalelement = max(randomnumbersperlocalelement,dynamic_cast<DRT::ELEMENTS::Truss3*>(dis.lColElement(i))->HowManyRandomNumbersINeed());

        //in case of periodic boundary conditions truss3 elements require a special initialization if they are broken by the periodic boundaries in the initial configuration
        if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
          statmechmanager_->PeriodicBoundaryTruss3Init(dis.lColElement(i));
        break;
      }
#endif  // #ifdef D_BEAM2
      default:
        continue;
    }
  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  /*so far the maximal number of random numbers required per element has been checked only locally on this processor;
   *now we compare the results of each processor and store the maximal one in maxrandomnumbersperglobalelement_*/
  dis.Comm().MaxAll(&randomnumbersperlocalelement,&maxrandomnumbersperglobalelement_ ,1);

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


  double dt = params_.get<double>("delta time" ,0.01);

  //getting number of dimensions for diffusion coefficient calculation
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  int ndim=psize.get<int>("DIM");


  for (int i=step; i<nstep; ++i)
  {

    //if input flag is set random number seed should be the same for all realizations
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"FIXEDSEED"))
    {
      //random generator for seeding only (necessary for thermal noise)
      ranlib::Normal<double> seedgenerator(0,1);
      //seeding random generator
      int seedvariable = i;
      seedgenerator.seed((unsigned int)seedvariable);
    }


    /*in the very first step and in case that special output for statistical mechanics is requested we have
     * to initialized the related output method*/
    if(i == 0)
      statmechmanager_->StatMechInitOutput(ndim,dt);

    //processor 0 write total number of elements at the beginning of time step i to console as well as how often a time step had to be restarted due to bad random numbers
    if(!discret_.Comm().MyPID())
    {
      std::cout<<"\nNumber of elements at the beginning of time step "<<i<<" : "<<discret_.NumGlobalElements()<<"\n";
      std::cout<<"\nNumber of unconverged steps "<<unconvergedsteps_<<"\n";
    }


      double time = params_.get<double>("total time",0.0);
      statmechmanager_->time_ = time;


      do
      {
        //assuming that iterations will converge
        isconverged_ = 1;

        /*multivector for stochastic forces evaluated by each element; the numbers of vectors in the multivector equals the maximal
         *number of random numbers required by any element in the discretization per time step; therefore this multivector is suitable
         *for synchrinisation of these random numbers in parallel computing*/
        RCP<Epetra_MultiVector> randomnumbers = rcp( new Epetra_MultiVector(*(discret_.ElementColMap()),maxrandomnumbersperglobalelement_) );

        //pay attention: for a constant predictor an incremental velocity update is necessary, which has
        //been deleted out of the code in oder to simplify it

        //generate gaussian random numbers for parallel use with mean value 0 and standard deviation (2KT / dt)0.5
        statmechmanager_->GenerateGaussianRandomNumbers(randomnumbers,0,pow(2.0 * (statmechmanager_->statmechparams_).get<double>("KT",0.0) / dt,0.5));

        ConsistentPredictor(randomnumbers);


        if(ndim ==3)
          PTC(randomnumbers);
        else
          FullNewton(randomnumbers);

        /*if iterations have not converged a new trial requires setting all intern element variables to
         * status at the beginning of this time step*/
        if(isconverged_ == 0)
        {
          ParameterList p;
          p.set("action","calc_struct_reset_istep");
          discret_.Evaluate(p,null,null,null,null,null);
        }

      }
      while(isconverged_ == 0);

        const double t_admin = Teuchos::Time::wallTime();

    UpdateandOutput();

    /*special update for statistical mechanics; this output has to be handled seperately from the time integration scheme output
     * as it may take place independently on writing geometric output data in a specific time step or not*/
    statmechmanager_->StatMechUpdate(dt,*dis_,stiff_,ndim);

    statmechmanager_->StatMechOutput(params_,ndim,time,i,dt,*dis_,*fint_);

    if(!discret_.Comm().MyPID())
    cout << "\n***\ntotal administration time: " << Teuchos::Time::wallTime() - t_admin<< " seconds\n***\n";

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
  bool   dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");

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
    	dserror("Proper PeriodLength is required if periodic DBCs are to be applied");

    // in case of activated periodic boundary conditions
    if(statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0) > 0.0 &&
			 Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"PERIODICDBC"))
    {
    	/* Reinitialize disn_ and dirichtoggle_ once.
    	 * Now, why is this done? For t==0, disn_ and dirichtoggle_ are initialized in strugenalpha.cpp.
    	 * Especially dirichtoggle_ and invtoggle_ contain information that is incorrect if DBC DOFs are
    	 * selected anew for each timestep and periodic boundary conditions are to be applied.
    	 * Also, initialize drefnew_ with large abs. values
    	 * The "incorrect" initialization occurs due to DBCs defined in the input file.
    	 * One may interject, that they've been defined for a reason, so why cancel them here?
    	 * answer: reasons of flexibility.
    	 * If those conventional DBCs are really needed they can be readded later.
    	 */
    	if(!isinit_)
    	{
    		for(int i=0;i<disn_->MyLength();i++)
    		  (*disn_)[i] = 0.0;
    		for(int i=0;i<dirichtoggle_->MyLength();i++)
    		  (*dirichtoggle_)[i] = 0.0;
    		for(int i=0;i<invtoggle_->MyLength();i++)
    		    		  (*invtoggle_)[i] = 1.0;
    		for(int i=0;i<drefnew_->MyLength();i++)
    			(*drefnew_)[i] = 9e99;
    	}
    	EvaluateDirichletPeriodic(p);
    	//cout<<"Periodic: "<<*disn_<<endl;
    }
		// "common" case without periodic boundary conditions
    if(Teuchos::getIntegralValue<int>(statmechmanager_->statmechparams_,"CONVENTIONALDBC"))
    {
    	discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
			//cout<<"Conventional: "<<*disn_<<endl;
    }

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


  // zerofy velocity and acceleration in case of statics
  if (dynkindstat)
  {
    velm_->PutScalar(0.0);
    veln_->PutScalar(0.0);
    vel_->PutScalar(0.0);
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


    //passing statistical mechanics parameters to elements
    p.set("ETA",(statmechmanager_->statmechparams_).get<double>("ETA",0.0));
    p.set("THERMALBATH",Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmechmanager_->statmechparams_,"THERMALBATH"));
    p.set("FRICTION_MODEL",Teuchos::getIntegralValue<INPAR::STATMECH::FrictionModel>(statmechmanager_->statmechparams_,"FRICTION_MODEL"));
    p.set("RandomNumbers",randomnumbers);
    p.set("SHEARAMPLITUDE",(statmechmanager_->statmechparams_).get<double>("SHEARAMPLITUDE",0.0));
    p.set("CURVENUMBER",(statmechmanager_->statmechparams_).get<int>("CURVENUMBER",-1));
    p.set("OSCILLDIR",(statmechmanager_->statmechparams_).get<int>("OSCILLDIR",-1));
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

  }
  // add static mid-balance

  fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);



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
			//test output to determine source of divergence
      // Test 1: Check residuals (show only those that surpass a certain treshold)
      /*if(fresmnorm>10)
      {
      	vector<int> dbcnodes;
      	for(int i=0; i<dirichtoggle_->MyLength(); i++)
      	{
      		if((*fresm_)[i]>0.5)
      		{
      			cout<<"Node "<<(int)floor((double)i/6.0)<<" DOFabs "<<i<<" DOF "<<i%6<<" : fresm@"<<i<<": "<<(*fresm_)[i]<<" DBC-DOF: ";
      			if((*dirichtoggle_)[i]==1)
      				cout<<" Yes"<<endl;
      			else
      				cout<<" No"<<endl;
      		}
    			if((*dirichtoggle_)[i]==1)
    			{
    				if(!dbcnodes.empty() && dbcnodes.back()==(int)floor((double)i/6.0))
    					continue;
    				bool dirichlet;
    				int tmpid = (int)floor((double)i/6.0);
    				for(int j=0; j<6; j++)
    				{
    					if((*fresm_)[6*tmpid+j]>0.5)
    					{
    						dirichlet=true;
    						break;
    					}
    					else
    						dirichlet=false;
    				}

    				if(dirichlet==true)
    					dbcnodes.push_back(tmpid);
    			}
      	}
      	cout<<"DBC-Nodes with high non-DBC DOF residuals: "<<endl;
      	for(int i=0; i<(int)dbcnodes.size(); i++)
      	{
      		bool done=false;
      		cout<<dbcnodes.at(i)<<" ";
      		DRT::Node* tmpnode = discret_.gNode(dbcnodes.at(i));
      		for(int j=0; j<tmpnode->NumElement(); j++)
      		{
      			DRT::Element* tmpelement = tmpnode->Elements()[j];
      			LINALG::SerialDenseMatrix coord(3,(int)discret_.lColElement(0)->NumNode(), true);
      		  LINALG::SerialDenseMatrix cut(3,(int)discret_.lColElement(0)->NumNode()-1,true);
      			bool broken;
      			vector<int> lids;
      			statmechmanager_->GetElementNodeCoords(tmpelement,disn_,coord, &lids);
      			statmechmanager_->CheckForBrokenElement(coord, cut, &broken);
      			if(cut(2,0)==0)
      				continue;
      			else if(cut(2,0)!=0 && !done)
      			{
      				for(int k=0; k<coord.N(); k++)
      					if((int)floor((double)lids.at(3*k)/6.0)==tmpnode->LID())
      					{
      						if(coord(2,k)>(statmechmanager_->statmechparams_.get<double>("PeriodLength", 0.0))/2.0)
      							cout<<"oscillating node  ";
      						else
      							cout<<"fixed node  ";
      					}
      				done=true;
      			}
      		}

      		if(i+1<(int)dbcnodes.size())
      			if(dbcnodes.at(i+1)-dbcnodes.at(i)>1)
      				cout<<"\n"<<endl;
      	}
      	cout<<"\n"<<endl;
      }*/
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
    unconvergedsteps_++;
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
 |  do Newton iteration (public)                             mwgee 03/07|
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


  // hard wired ptc parameters
  double ptcdt = 1.3e1; //1.3e1
  double nc;
  fresm_->NormInf(&nc);
  double dti = 1/ptcdt;
  double dti0 = dti;
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


#if 1 // SER
#else // TTE
    double dtim = dti0;
#endif

    dti0 = dti;
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
      p.set("dti",dti);

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

      // do NOT finalize the stiffness matrix to add damping to it later
    }

    //------------------------------------------ compute residual forces

    // dynamic residual
    // Res =  C . V_{n+1-alpha_f}
    //        + F_int(D_{n+1-alpha_f})
    //        - F_{ext;n+1-alpha_f}
    // add mid-inertial force


    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,0.0);


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
                  fresmnorm,disinorm,convcheck,dti);
    }

    //------------------------------------ PTC update of artificial time
#if 1
    // SER step size control
    dti *= (np/nc);
    dti = max(dti,0.0);
    nc = np;


    //Modifikation: sobald Residuum klein, PTC ausgeschaltet
    if(np < 0.01*resinit)
      dti = 0.0;


#else
    {
      // TTE step size control
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

      //Modifikation: sobald Residuum klein, PTC ausgeschaltet
      if(np < 0.01*resinit)
        dti = 0.0;

    }
#endif


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
    unconvergedsteps_++;
    std::cout<<"\n\niteration unconverged - new trial with new random numbers!\n\n";
     //dserror("FullNewton unconverged in %d iterations",numiter);
  }
  else
  {
     if (!myrank_ and printscreen)
     {
       PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                   fresmnorm,disinorm,convcheck,dti);
     }
  }

  params_.set<int>("num iterations",numiter);

  if(!discret_.Comm().MyPID())
  std::cout << "\n***\nevaluation time: " << sumevaluation<< " seconds\nptc time: "<< sumptc <<" seconds\nsolver time: "<< sumsolver <<" seconds\ntotal solution time: "<<Teuchos::Time::wallTime() - tbegin<<" seconds\n***\n";

  return;
} // StatMechTime::PTC()


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
    statmechmanager_->StatMechWriteRestart(output_);
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
  statmechmanager_->StatMechReadRestart(reader);

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

/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (public)         mueller  2/10|
 *----------------------------------------------------------------------*/
void StatMechTime::EvaluateDirichletPeriodic(ParameterList& params)
/*The idea behind EvaluateDirichletPeriodic() is simple:
 * Give Dirichlet values to nodes of an element that is broken in
 * z-direction due to the application of periodic boundary conditions.
 * The motion of the node close to z=0.0 in the cubic volume of edge
 * length l (==PeriodLength in this case) is inhibited in direction of the
 * oscillatory motion. The oscillation is imposed on the node close to z=l.
 * This method is triggered in case of PeriodLength>0.0 (i.e. Periodic BCs
 * exist). Since the DBC setup happens dynamically by checking element
 * positions with each new time step, the static definition of DBCs in the input
 * file is only used to get the direction of the oscillatory motion as well as
 * the time curve. Therefore, only one DBC needs to be specified.
 *
 * How it works:
 * Each time this method is called, the system vector and the toggle vector are
 * modified to fit the current geometric situation.
 * DOFs holdings Dirichlet values are marked by setting the corresponding toggle
 * vector component to 1.0. In case of an element which was broken the step before
 * and is now whole again,just the toggle vector components in question are
 * reset to 0.0.
 * A position vector drefnew_ is needed in order to calculate the correct Dirichlet
 * values to be imposed on nodes of an element which has drifted over the boundaries
 * and thus has been broken.
 * These positions are used to calculate the zero position of the oscillation which then can
 * be added to the time curve value in DoDirichletConditionPeriodic().
 */
{
#ifdef MEASURETIME
  const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

	if (!(discret_.Filled())) dserror("FillComplete() was not called");
	if (!(discret_.HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");

  //----------------------------------- some variables
	// indicates broken element
	bool 											broken;
	// time trigger
	bool 											usetime = true;
	// to avoid redundant or wrong actions when filling vectors or deleting the last element of the free nodes vector
	bool											alreadydone = false;
	// store node Id of previously handled node
  int 											tmpid = -1;
  // store LIDs of element nodes
  vector<int>								lids;
  // store LIDs of DOFs subjected to oscillation in order to identify force sensor locations
  vector<int>								fsensorlids;
  // vectors to manipulate DBC properties
  vector<int> 							oscillnodes;
  vector<int> 							fixednodes;
  vector<int> 							freenodes;

  // get the current time
	const double time = params.get("total time",-1.0);
	if (time<0.0) usetime = false;

  // init
  if(!isinit_)
  {
  	// get the apmlitude of the oscillation
  	amp_= statmechmanager_->statmechparams_.get<double>("SHEARAMPLITUDE",0.0)*
					statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0);
		// retrieve direction of oscillatory motion
		oscdir_ = statmechmanager_->statmechparams_.get<int>("OSCILLDIR",0)-1;
		// retrieve number of time curve that is to be applied
		curvenumber_ = statmechmanager_->statmechparams_.get<int>("CURVENUMBER",0)-1;

		isinit_=true;
	} // init

//---------------------------------------------------------- loop through elements
	for(int i=0;i<discret_.NumMyRowElements(); i++)
	{
		// An element used to browse through Row Elements
	  DRT::Element* 						element = discret_.lRowElement(i);
	  // positions of nodes of an element with n nodes
	  LINALG::SerialDenseMatrix coord(3,(int)discret_.lRowElement(i)->NumNode(), true);
	  // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
	  LINALG::SerialDenseMatrix cut(3,(int)discret_.lRowElement(i)->NumNode()-1,true);
//-------------------------------- obtain nodal coordinates of the current element
	  // get nodal coordinates and LIDs of the nodal DOFs
	  statmechmanager_->GetElementNodeCoords(element, disn_, coord, &lids);
//-----------------------detect broken/fixed/free elements and fill position vector
	  // determine existence and location of broken element
	  statmechmanager_->CheckForBrokenElement(coord, cut, &broken);

	  for(int n=0; n<cut.N(); n++)
	  {
			// case: broken element (in z-dir); node_n+1 oscillates, node_n is fixed in dir. of oscillation
			if(broken && cut(2,n)==1.0)
			{
				// indicates beginning of a new filament (in the very special case that this is needed)
				bool newfilament = false;
				// check for case: last element of filament I as well as first element of filament I+1 broken
				if(tmpid!=element->Nodes()[n]->Id() && alreadydone)
				{
					// in this case, reset alreadydone...
					alreadydone = false;
					// ...and set newfilament to true. Otherwise the last free nodes vector element will be deleted
					newfilament = true;
				}

				// add GID of fixed node to fixed-nodes-vector (to be added to condition later)
				if(!alreadydone)
					fixednodes.push_back(element->Nodes()[n]->Id());
				// add GID of oscillating node to osc.-nodes-vector
				oscillnodes.push_back(element->Nodes()[n+1]->Id());
				// add displacement of the node to drefnew_
				// (at this stage, this is not the displacement of the new reference position for oscillating nodes yet:
				// drefnew_ is set to current displacement (without Dirichlet) if its lid-th entry is set (or reset) to its initial value (9e99),
				// i.e. the element in question was not broken in the preceding timestep.)
				for(int j = n; j < n+2; j++)
					for(int k=0; k<cut.M(); k++)
						if((*drefnew_)[lids.at(3*j+k)] > statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0))
							(*drefnew_)[lids.at(3*j+k)] = (*disn_)[lids.at(3*j+k)];

				// add DOF LID where a force sensor is to be set
				if(Teuchos::getIntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"DYN_CROSSLINKERS"))
				{
					int lid = lids.at(3*(n+1)+oscdir_);
					fsensorlids.push_back(lid);
				}

				// delete last Id of freenodes if it was previously and falsely added
				if(element->Nodes()[n]->Id()==tmpid && !alreadydone && !newfilament)
					freenodes.pop_back();
				// store gid of the "n+1"-node to avoid overwriting during the following iteration,
				// e.g. oscillating node becomes free if the following CheckForBrokenElement() call yields "!broken".
				tmpid = element->Nodes()[n+1]->Id();
				// Set to true to initiate certain actions if the following element is also broken.
				// If the following element isn't broken, alreadydone will be reset to false (see case: !broken)
				alreadydone=true;
			}
			// case: broken element (in z-dir); node_n oscillates, node_n+1 is fixed in dir. of oscillation
			if(broken && cut(2,n)==2.0)
			{
				bool newfilament = false;

				if(tmpid!=element->Nodes()[n]->Id() && alreadydone)
				{
					alreadydone = false;
					newfilament = true;
				}

				if(!alreadydone)
					oscillnodes.push_back(element->Nodes()[n]->Id());
				fixednodes.push_back(element->Nodes()[n+1]->Id());

				for(int j=n; j<n+2; j++)
					for(int k=0; k<cut.M(); k++)
						if((*drefnew_)[lids.at(3*j+k)] > statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0))
							(*drefnew_)[lids.at(3*j+k)] = (*disn_)[lids.at(3*j+k)];

				if(Teuchos::getIntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"DYN_CROSSLINKERS"))
				{
					int lid = lids.at(3*n+oscdir_);
					fsensorlids.push_back(lid);
				}

				if(element->Nodes()[n]->Id()==tmpid && !alreadydone && !newfilament)
					freenodes.pop_back();

				tmpid = element->Nodes()[n+1]->Id();
				alreadydone = true;
			}
			// case: unbroken element or broken in another than z-direction
			if(cut(2,n)!=2)
			{
				if(element->Nodes()[n]->Id()!=tmpid)
				{
					freenodes.push_back(element->Nodes()[n]->Id());
					freenodes.push_back(element->Nodes()[n+1]->Id());
				}
				else
					freenodes.push_back(element->Nodes()[n+1]->Id());
				tmpid=element->Nodes()[n+1]->Id();
				// set to false to handle annoying special cases
				alreadydone = false;
			}
	  }
	}

//---------check/set force sensors anew for each time step
  if(Teuchos::getIntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"DYN_CROSSLINKERS"))
  	statmechmanager_->UpdateForceSensors(fsensorlids);

//------------------------------------set Dirichlet values
	// preliminary
	DRT::Node* node = discret_.gNode(discret_.NodeRowMap()->GID(0));
	int numdof = (int)discret_.Dof(node).size();

	vector<int> 	 addcurve(numdof, -1);
	vector<int> 	 addfunct(numdof, 0);
	vector<int>  	 addonoff(numdof, 0);
	vector<double> addval(numdof, 0.0);

  // set condition for oscillating nodes
	addcurve.at(oscdir_) = curvenumber_;
	addval.at(oscdir_) = amp_;
	// inhibit planar DOFs, leave free only z-direction (for now, inhibiting only oscdir_ seems to cause instabilities)
	for(int i=0; i<3; i++)
		if(i==oscdir_ || i!=2)
			addonoff.at(i) = 1;

  // do not do anything if vector is empty
  if(!oscillnodes.empty())
  	DoDirichletConditionPeriodic(usetime, time, &oscillnodes, &addcurve, &addfunct, &addonoff, &addval);

  // set condition for fixed nodes
	addcurve.at(oscdir_) = -1;
	addval.at(oscdir_) = 0.0;

	if(!fixednodes.empty())
  	DoDirichletConditionPeriodic(usetime, time, &fixednodes, &addcurve, &addfunct, &addonoff, &addval);

  // set condition for free or recently set free nodes
  for(int i=0; i<3; i++)
  	if(i==oscdir_ || i!=2)
  		addonoff.at(i) = 0;

	if(!freenodes.empty())
  	DoDirichletConditionPeriodic(usetime, time, &freenodes, &addcurve, &addfunct, &addonoff, &addval);

#ifdef MEASURETIME
  const double t_end = Teuchos::Time::wallTime();
  cout<<"DBC Evaluation time: "<<t_end-t_start<<endl;
#endif // #ifdef MEASURETIME
	return;
}

/*----------------------------------------------------------------------*
 |  fill system vector and toggle vector (public)          mueller  3/10|
 *----------------------------------------------------------------------*/
void StatMechTime::DoDirichletConditionPeriodic(const bool usetime,
                                                const double time,
																								vector<int>* nodeids,
																								vector<int>* curve,
																								vector<int>* funct,
																								vector<int>* onoff,
																								vector<double>* val)
/*
 * This basically does the same thing as DoDirichletCondition() (to be found in drt_discret_evaluate.cpp),
 * but with the slight difference of taking current displacements into account.
 * Time curve values aren't added to the reference position(s) of the discretization as usual,
 * but to the latest known 0-position(s). These positions are calculated using the drefnew_
 * vector holding the latest node positions.
 */
{
	/*/ test output
	cout<<"NODES: ";
	for(int i=0; i<(int)nodeids->size(); i++)
		cout<<nodeids->at(i)<<" ";
	cout<<endl;
	cout<<"curve: "<<endl;
	for(int i=0; i<6; i++)
		cout<<curve->at(i)<<" ";
	cout<<endl;
	cout<<"funct: "<<endl;
	for(int i=0; i<6; i++)
		cout<<funct->at(i)<<" ";
	cout<<endl;
	cout<<"onoff: "<<endl;
	for(int i=0; i<6; i++)
		cout<<onoff->at(i)<<" ";
	cout<<endl;
	cout<<"val: "<<endl;
	for(int i=0; i<6; i++)
		cout<<val->at(i)<<" ";
	cout<<endl;*/
	// highest degree of requested time derivative (may be needed in the future(?))
	unsigned deg = 0;
	// some checks for errors
	if (!nodeids) dserror("No Node IDs were handed over!");
	if(disn_==Teuchos::null) dserror("Displacement vector must be unequal to null");
	if(dbcswitch_==Teuchos::null || drefnew_==Teuchos::null || dirichtoggle_==Teuchos::null)
		dserror("dbcwitch_, drefnew_ and dirichtoggle_ must be non-empty");

	// get the condition properties
	const int nnode = nodeids->size();


	// loop over all nodes in condition
	for (int i=0; i<nnode; ++i)
	{
		// do only nodes in my row map
		if (!discret_.NodeRowMap()->MyGID(nodeids->at(i))) continue;
		DRT::Node* actnode = discret_.gNode(nodeids->at(i));
		if (!actnode) dserror("Cannot find global node %d",nodeids->at(i));
		// call explicitly the main dofset, i.e. the first column
		vector<int> dofs = discret_.Dof(0,actnode);
		const unsigned numdf = dofs.size();

		// loop over DOFs
		for (unsigned j=0; j<numdf; ++j)
		{
			// get the LID and GID of the currently handled DOF
			const int lid = (*disn_).Map().LID(dofs[j]);
      const int gid = dofs[j];

			// if DOF in question is not subject to DBCs (anymore)
      if (onoff->at(j)==0)
      {
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (dirichtoggle_!=Teuchos::null)
        {
        	// turn off application of Dirichlet value
          (*dirichtoggle_)[lid] = 0.0;
          // in addition, modify the inverse vector (needed for manipulation of the residual vector)
          (*invtoggle_)[lid] = 1.0;
        }
        // if formerly subject to Dirichlet values, turn dbcswitch_ off for this LID
        if((*dbcswitch_)[lid]!=0.0)
        	(*dbcswitch_)[lid] = 0.0;
        // reset reference position if lid bleongs to node formerly subjected to a DBC
        if((*drefnew_)[lid] <= statmechmanager_->statmechparams_.get<double>("PeriodLength",0.0))
        	(*drefnew_)[lid] = 9e99;
        continue;
      }

//---------------------------------------------Dirichlet Value Assignment

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct) funct_num = funct->at(j);
      {
         if (funct_num>0)
           functfac =
             DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(j,
                                                                  actnode->X(),
                                                                  time,
                                                                  &discret_);
      }

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = curve->at(j);
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      // Calculate displacement of new reference position.
      // Especially important when node is added to condition afterwards
			if((*dbcswitch_)[lid]==0.0 && val->at(j)!=0.0)
			{
				double dt = params_.get<double>("delta time" ,0.01);
				// calculate the new reference value for this DOF ("-dt" because the value preceding the current one is needed)
				double drefnm = DRT::Problem::Instance()->Curve(curvenum).f(time-dt);
				drefnm *= val->at(j)*functfac;
				// calculate new displacement d_ref,new,j = d_node,j - d_dbc,j of j-th DOF
				(*drefnew_)[lid] -= drefnm;
				// turn dbcswitch_ "on" to make sure this offset value is not changed until the the node,
				// to which the DOF belongs, is set free and then by chance is subjected to a DBC again.
				(*dbcswitch_)[lid] = 1.0;
			}

			vector<double> value(deg+1,val->at(j));

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
        // add offset to Dirichlet value
        value[i] += (*drefnew_)[lid];
      }

	    // assign value
	    if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
	    if (disn_ != Teuchos::null)
	      (*disn_)[lid] = value[0];
	    // set toggle vector and the inverse vector
	    if (dirichtoggle_ != Teuchos::null)
	    {
	      (*dirichtoggle_)[lid] = 1.0;
	      (*invtoggle_)[lid] = 0.0;
	    }
	  }  // loop over nodal DOFs
	}  // loop over nodes
	return;
}

#endif  // #ifdef CCADISCRET
