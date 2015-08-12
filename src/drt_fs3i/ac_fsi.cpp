/*!----------------------------------------------------------------------
\file ac_fsi.cpp

\brief cpp-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches. This file thereby holds
       all functions related with the small time scales simulation and all the basic control structures.

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

\level 3
*----------------------------------------------------------------------*/


#include "ac_fsi.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_ac_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../linalg/linalg_utils.H"

//FOR WSS CALCULATIONS
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"

/*!
 * What does the problem type Atherosclerosis_Fluid_Structure_Interaction?
 * Short answer: cool stuff!
 * And here is the long answer:
 * Its doing a multiscale (in time) approach with an full FS3I simulation at the small time scale (seconds) and
 * a scalar transport simulation at the larger time scale (days). The solving strategy is as follows:
 *   We start with the full small time scale FS3i simulation (including fluid Windkessel and WSS permeability).
 * After each FSI cycle we check if the FSI problem is periodic, by looking if the Windkessel does produce
 * peridodic results. Afterwards we continue the small time scale but do not solve the fsi subproblem anymore. Instead
 * we peridically repeat it by calling suitible restarts. When the fluid scatra subproblem gets periodic
 * at the FS3I interface we stop the small time scale and switch to the large time scale
 *   Now we higer dt_ and only solve the structure scatra problem. We thereby use the WSS and interface
 * conentrtions of the small time scale. Each time when there as been 'created' enough grothw inducing mass we
 * do a growth update. If we have finally grew to much, we go back to the small time scale. And so on, and so on,...
 */

/*----------------------------------------------------------------------*
 | constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
FS3I::ACFSI::ACFSI(const Epetra_Comm& comm)
  :PartFS3I(comm),
   structureincrement_(LINALG::CreateVector(*fsi_->StructureField()->DofRowMap(0),true)),
   fluidincrement_(LINALG::CreateVector(*fsi_->FluidField()->DofRowMap(0),true)),
   aleincrement_(LINALG::CreateVector(*fsi_->AleField()->DofRowMap(),true)),
   fluidphinp_lp_(LINALG::CreateVector(*scatravec_[0]->ScaTraField()->DofRowMap()/**scatrafieldexvec_[0]->Map(1)*/,true)),
   structurephinp_bltsl_(LINALG::CreateVector(*scatravec_[1]->ScaTraField()->DofRowMap(),true)),
   structurephinp_blgu_(LINALG::CreateVector(*scatravec_[1]->ScaTraField()->DofRowMap(),true)),
   WallShearStress_lp_(LINALG::CreateVector(*fsi_->FluidField()->DofRowMap(0),true)),
   fsiperiod_(DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("PERIODICITY") ),
   dt_large_(DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("LARGE_TIMESCALE_TIMESTEP")),
   fsiisperiodic_(false),
   scatraisperiodic_(false),
   fsineedsupdate_(false),
   extractjthscalar_(BuildMapExtractor())
{
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = scatravec_[1]->ScaTraField()->Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that two-way coupling is possible
  if (fsi_->StructureField()->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");


  //Some AC FSI specific testings:

  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  const Teuchos::ParameterList& fs3idynpart = fs3idyn.sublist("PARTITIONED");
  const Teuchos::ParameterList& fs3idynac = fs3idyn.sublist("AC");

  if ( fs3idyn.get<int>("UPRES") != 1 )
    dserror("If you want the fsi problem to be periodically repeated from some point, you have to have UPRES set to 1!");
  if ( fs3idyn.get<int>("RESTARTEVRY") != 1 )
    dserror("If you want the fsi problem to be periodically repeated from some point, you have to have RESTARTEVRY set to 1!");

  if ( fsiperiod_ <= 1e-14 )
    dserror("You need to specify the PERIODICITY. If you don't want that change your PROBLEMTYPE to gas_fluid_structure_interaction!");

  if ( (dt_large_+1e-14) <= dt_ )
    dserror("You need to specify a LARGE_TIMESCALE_TIMESTEP and it must be larger than the 'normal' fs3i TIMESTEP!");

  if ( (bool)DRT::INPUT::IntegralValue<int>( DRT::Problem::Instance()->IOParams(),"FLUID_WALL_SHEAR_STRESS") )
    dserror("You need FLUID_WALL_SHEAR_STRESS in section IO to be NO,"
            " otherwise GetWallShearStresses() will be called twice, which is bad if you use WSS_TYPE Mean."
            " Nevertheless we do write the wss output!");

  std::vector<DRT::Condition*> ImpCond;
  DRT::Problem::Instance()->GetDis("fluid")->GetCondition("ImpedanceCond",ImpCond);
  for (unsigned int i=0; i<ImpCond.size(); i++)
  {
    const double thisperiod = (ImpCond[i])->GetDouble("timeperiod");

    if (thisperiod != fsiperiod_)
    {
      dserror("your impedance TIMEPERIOD and your fs3i PERIODICITY do not match!");
    }
  }

  if ( not ModuloIsRealtiveZero(fsiperiod_,dt_,fsiperiod_) )
    dserror("Choose a time step such that TIMESTEP = PERIODIC / n with n being an integer!");

  if ( not ModuloIsRealtiveZero(dt_large_,fsiperiod_,dt_large_) )
      dserror("Choose LARGE_TIMESCALE_TIMESTEP as a multiple of PERIODICITY!");

  if (infperm_)
    dserror("AC-FS3I does have a finite interface permeability. So set INF_PERM to NO!");

  const int fsiperssisteps = fs3idynac.get<int>("FSI_STEPS_PER_SCATRA_STEP");
  const INPAR::FS3I::SolutionSchemeOverFields couplingalgo = DRT::INPUT::IntegralValue<INPAR::FS3I::SolutionSchemeOverFields>(fs3idynpart,"COUPALGO");
  const double wk_rel_tol = fs3idynac.get<double>("WINDKESSEL_REL_TOL");

  if (couplingalgo == INPAR::FS3I::fs3i_IterStagg and fsiperssisteps != 1)
    dserror("In an iteratively staggered FS3I scheme subcycling is not supported!");

  if (couplingalgo == INPAR::FS3I::fs3i_IterStagg and not IsRealtiveEqualTo(wk_rel_tol,-1.0,1.0) )
    dserror("In an iteratively staggered FS3I scheme periodical repetition of the fsi problem is not supported!");

  if ( not IsRealtiveEqualTo(dt_,fsiperssisteps*fsi_->FluidField()->Dt(),1.0) )
    dserror("Your fluid time step does not match!");
  if ( not IsRealtiveEqualTo(dt_,fsiperssisteps*fsi_->StructureField()->Dt(),1.0) )
    dserror("Your structure time step does not match!");
  if ( not IsRealtiveEqualTo(dt_,fsiperssisteps*fsi_->AleField()->Dt(),1.0) )
    dserror("Your ale time step does not match!");
  if ( not IsRealtiveEqualTo(dt_,scatravec_[0]->ScaTraField()->Dt(),1.0) )
    dserror("Your fluid scatra time step does not match!");
  if ( not IsRealtiveEqualTo(dt_,scatravec_[1]->ScaTraField()->Dt(),1.0) )
    dserror("Your structure scatra time step does not match!");
}

/*----------------------------------------------------------------------*
 | Read restart                                              Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    const Teuchos::ParameterList& fs3idynac = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC");

    const bool restartfromfsi = DRT::INPUT::IntegralValue<int>(fs3idynac,"RESTART_FROM_PART_FSI");

    if (not restartfromfsi) //standard restart
    {
      fsi_->ReadRestart(restart);

      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->ReadRestart(restart);
      }

      //AC-FSI specific input
      IO::DiscretizationReader reader = IO::DiscretizationReader(fsi_->FluidField()->Discretization(),restart);
      reader.ReadVector(WallShearStress_lp_, "wss");
    }
    else //we do not want to read the scatras values and the lagrange multiplyer, since we start from a partitioned FSI
    {
      //restart FSI problem
      fsi_->StructureField()->ReadRestart(restart);
      fsi_->FluidField()->ReadRestart(restart);

      fsi_->AleField()->ReadRestart(restart);

      fsi_->SetTimeStep(fsi_->FluidField()->Time(),fsi_->FluidField()->Step());

      //AC-FSI specific input
      IO::DiscretizationReader reader = IO::DiscretizationReader(fsi_->FluidField()->Discretization(),restart);
      reader.ReadVector(WallShearStress_lp_, "wss");

      // We may want to calculate the initial time deriv, since the scatra field is new.
      // This has to be done before setting the step in the scatra fields, because step_ it needs to be zero.
      const Teuchos::ParameterList& scatradynparams = DRT::Problem::Instance()->ScalarTransportDynamicParams();
      const bool skiptimederiv = DRT::INPUT::IntegralValue<int>(scatradynparams,"SKIPINITDER");

      if (not skiptimederiv)
      {
        SetMeshDisp();
        SetVelocityFields();
        //SetFSISolution(); //Before we have to set the velocity fields

        for (unsigned i=0; i<scatravec_.size(); ++i)
          {
            Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
            currscatra->ScaTraField()->PrepareFirstTimeStep();
          }
      }

      //we need to set time and step in scatra to have it matching with the fsi ones
      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->SetTimeStep(fsi_->FluidField()->Time(),fsi_->FluidField()->Step());
      }
    }

    time_ = fsi_->FluidField()->Time();
    step_ = fsi_->FluidField()->Step();

    //AC-FSI specific input
    IO::DiscretizationReader reader = IO::DiscretizationReader(fsi_->FluidField()->Discretization(),restart);

    reader.ReadVector(WallShearStress_lp_, "wss");

    //recover IsPeriodic flags by recalculating them
    IsFsiPeriodic();
    //Note: this does not work for the scatra test  since we don't save phinp_lp for the sake of less storage
    //IsScatraSteady();
  }
}

/*----------------------------------------------------------------------*
 | Timeloop                                                  Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::Timeloop()
{
  // output of initial state
//  fsi_->PrepareOutput();
//  FsiOutput();
  ScatraOutput();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetVelocityFields(); //doing this we can use the flag SKIPINITDER

  // do time loop
  while (NotFinished())
  {
    SmallTimeScaleLoop();
    LargeTimeScaleLoop();
  }
}

/*----------------------------------------------------------------------*
 | timeloop for small time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleLoop()
{
  //print info
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                         SMALL TIME SCALE LOOP"
               "\n************************************************************************\n"<<std::endl;
  }

  while (SmallTimeScaleLoopNotFinished())
  {
  //Do a time step
  SmallTimeScalePrepareTimeStep();
  SmallTimeScaleOuterLoop();
  SmallTimeScaleUpdateAndOutput();

  //Update the isperiodic_ flags
  IsSmallTimeScalePeriodic();
  }
}

/*--------------------------------------------------------------------------*
 | flag whether small time scale time loop should be finished    Thon 07/15 |
 *--------------------------------------------------------------------------*/
bool FS3I::ACFSI::SmallTimeScaleLoopNotFinished()
{
  return ( NotFinished() and not (fsiisperiodic_ and scatraisperiodic_) ) ;
}

/*----------------------------------------------------------------------*
 | Prepare time step                                         Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScalePrepareTimeStep()
{
  IncrementTimeAndStep();

  //Print to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "\n\n"<< "TIME:  "    << std::scientific <<std::setprecision(8)<< time_ << "/" << std::scientific << timemax_
             << "     DT = " << std::scientific << dt_
             << "     STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << numstep_
             << "\n\n";
  }

  fsi_->PrepareTimeStep();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->PrepareTimeStep();
  }
}

/*----------------------------------------------------------------------*
 | OuterLoop                                                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleOuterLoop()
{
  const Teuchos::ParameterList& fs3idynpart = DRT::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED");
  // get coupling algorithm from input file
  const INPAR::FS3I::SolutionSchemeOverFields couplingalgo = DRT::INPUT::IntegralValue<INPAR::FS3I::SolutionSchemeOverFields>(fs3idynpart,"COUPALGO");

  switch (couplingalgo)
  {
  case INPAR::FS3I::fs3i_SequStagg:
    SmallTimeScaleOuterLoopSequStagg();
    break;
  case INPAR::FS3I::fs3i_IterStagg:
    SmallTimeScaleOuterLoopIterStagg();
    break;
  default:
    dserror("partitioned FS3I coupling scheme not implemented!");
    break;
  }
}

/*----------------------------------------------------------------------*
 | OuterLoop for sequentially staggered FS3I scheme          Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleOuterLoopSequStagg()
{
  SetStructScatraSolution();

  DoFSIStep();

  SetFSISolution();

  SmallTimeScaleDoScatraStep();
}

/*----------------------------------------------------------------------*
 | OuterLoop for iterative staggered FS3I scheme             Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleOuterLoopIterStagg()
{
  int itnum = 0;

  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************\n"
                 "                         OUTER ITERATION START"
             <<"\n************************************************************************"<<std::endl;
  }

  while ( stopnonliniter==false )
  {
    itnum++;

    structureincrement_->Update(1.0,*fsi_->StructureField()->Dispnp(),0.0);
    fluidincrement_->Update(1.0,*fsi_->FluidField()->Velnp(),0.0);
    aleincrement_->Update(1.0,*fsi_->AleField()->Dispnp(),0.0);

    SetStructScatraSolution();

    DoFSIStepStandard(itnum);
    //subcycling is not allowed, since we use this function for the growth update. Nevertheless it should work..
    //periodical repetition is not allowed, since we want to converge the problems

    SetFSISolution();

    SmallTimeScaleDoScatraStep();

    stopnonliniter = PartFs3iConvergenceCkeck(itnum);
  }
}

/*----------------------------------------------------------------------*
 | Do a single fsi step                                                 |
 | (including subcycling and periodic repetition )           Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStep()
{
  if ( not fsiisperiodic_ )
  {
    const int fsiperssisteps = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<int>("FSI_STEPS_PER_SCATRA_STEP");
    if (fsiperssisteps == 1) //no subcycling
    {
      DoFSIStepStandard();
    }
    else //subcycling
    {
      DoFSIStepSubcycled(fsiperssisteps);
    }
  }
  else //fsi problem is periodic
  {
    DoFSIStepPeriodic();
  }
}

void FS3I::ACFSI::IsSmallTimeScalePeriodic()
{
  if ( step_>1 and ModuloIsRealtiveZero(time_,fsiperiod_,time_) )
  {
    // Check fsi periodicity
    IsFsiPeriodic();

    // Check if scatra is periodic
    IsScatraPeriodic();
  }
}

/*----------------------------------------------------------------------*
 | Decide if fsi problem is already periodic                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::IsFsiPeriodic()
{
  fsiisperiodic_ = true;

  Teuchos::RCP<ADAPTER::FluidACFSI> fluid = Teuchos::rcp_dynamic_cast<ADAPTER::FluidACFSI>( fsi_->FluidField() );
  std::vector<double> wk_rel_errors = fluid->GetWindkesselErrors();

  const double wk_rel_tol = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("WINDKESSEL_REL_TOL");

  for (unsigned int i=0; i<wk_rel_errors.size(); i++)
  {
    if (Comm().MyPID()==0)
      std::cout<<std::scientific<<std::setprecision(2)<<"The "<<i+1<<"-th Windkessel has an relative error of "<<abs(wk_rel_errors[i])<<std::endl;

    if ( abs(wk_rel_errors[i]) > wk_rel_tol)
      fsiisperiodic_ = false;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Do a standard fsi step                                    Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepStandard(const int callnum)
{
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                               FSI SOLVER "
               "\n************************************************************************"<<std::endl;
  }
  fsi_->TimeStep(fsi_,callnum);
}

/*----------------------------------------------------------------------*
 | Do a fsi step with subcycling                             Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepSubcycled( const int subcyclingsteps)
{
  int subcyclingiter = 0;
  while (subcyclingiter < subcyclingsteps)
  {
    subcyclingiter++;

    if (subcyclingiter != 1) //for the first subcycling step we...
    {
      fsi_->PrepareOutput();    //... will do this in UpdateAndOutput()
      fsi_->Update();           //... will do this in UpdateAndOutput()
      fsi_->PrepareTimeStep();  //... have already done this in PrepareTimeStep()
      //now fix the step_ counter. When subcycling the fsi subproblem we do not want to proceed the step_ AND the time_, but just the time_.
      SetTimeStepInFSI(fsi_->FluidField()->Time(),step_);
    }

    if (Comm().MyPID()==0)
    {
      std::cout<<"\n************************************************************************"
                 "\n                     FSI SUBCYCLING SOLVER " << subcyclingiter <<"/"<<subcyclingsteps
               <<"\n************************************************************************"<<std::endl;
    }

    fsi_->TimeStep(fsi_); //all necessary changes for the fsi problem (i.e. adapting dt) has already been done in PartFS3I::ManipulateDt()
  }

  // do time and step in fsi (including subfields) and fs3i match after we will have called SmallTimeScaleUpdateAndOutput()
  if ( not IsRealtiveEqualTo(fsi_->FluidField()->Time(),time_,time_) )
    dserror("After the subcycling the fsi time and fs3i time do not match anymore!");

  if ( not IsRealtiveEqualTo( fsi_->FluidField()->Step(),step_,step_) )
    dserror("After the subcycling the fsi step and fs3i step do not match anymore!");
}

/*----------------------------------------------------------------------*
 | Get fsi solution from one period before                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepPeriodic()
{
  //instead of solving the FSI problem for the present time step, we take the
  //solution from the FSI period before. We do this by replacing all values
  //in FSI via ReadRestart(). Afterwards we just have to repair everything we
  //destroyed by calling ReadRestart().

  //this is the related step of the last period
  int previousperiodstep = GetStepOfPreviousPeriodAndPrepareReading(step_,time_);

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                         PERIODICAL FSI STEP"
               "\n************************************************************************\n"
               "\n Using results from timestep "<<previousperiodstep<<" as solution of the current FSI step"<<std::endl;
  }

  //do the reading
  fsi_->ReadRestart(previousperiodstep ); //ReadRestartfromTime()

  //we first fix the grid velocity of the fluid. This calculation is normally done
  //in ADAPTER::FluidFSI::ApplyMeshDisplacement(), but since we never call this function
  //we have to call it ourself
  fsi_->FluidField()->UpdateGridv(); //calculate grid velocity via FD approximation

  //update time and step in FSI and all subproblems
  SetTimeStepInFSI(time_,step_);
}

/*--------------------------------------------------------------------------------------*
 | Get step number of on cycle ago                                           Thon 07/15 |
 *--------------------------------------------------------------------------------------*/
double FS3I::ACFSI::GetStepOfPreviousPeriodAndPrepareReading(const int actstep,const double acttime)
{
  if ( not ModuloIsRealtiveZero(fsiperiod_,dt_,fsiperiod_) )
    dserror("PERIODICITY should be an multiple of TIMESTEP!");

  const int previousperiodstep = actstep- round(fsiperiod_/dt_ ); //here we assume a constant timestep over the last period

  //Is this the right step? Let's check:
  {
    //we have to clean the mapstack, otherwise it would fill with each iteration
    DRT::Problem::Instance()->GetDis("structure")->Writer()->ClearMapCache();

    //get filename in which the equivalent step of the last period is written
    std::string filename = GetFileName(actstep);
    //we always have to recreate the InputControl() since our Inputfile (=Outputfile) has changed in since the last reading (new timestep written)
    Teuchos::RCP<IO::InputControl> inputreader = Teuchos::rcp( new IO::InputControl(filename, Comm()) );
    //overwrite existing InputControl()
    DRT::Problem::Instance()->SetInputControlFile(inputreader);

    //AC-FSI specific input
    IO::DiscretizationReader reader = IO::DiscretizationReader(fsi_->FluidField()->Discretization(),previousperiodstep);

    double previousperiodtime = reader.ReadDouble("time");

    //Now check if the candidate is right
    if ( not IsRealtiveEqualTo(previousperiodtime + fsiperiod_,acttime,acttime) )
      dserror("You can't change your TIMESTEP when you are within the FSI PERIODIC cycle!");
  }

  return previousperiodstep;
}


/*--------------------------------------------------------------------------------------*
 | Get filename in which the equivalent step of the last period is written   Thon 12/14 |
 *--------------------------------------------------------------------------------------*/
std::string FS3I::ACFSI::GetFileName(const int actstep)
{
  //we have to be careful since in the fist steps after a restart the last fsi period is written
  //in the input and not in the output file

  std::string filename;

  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    const int crit_step = (int)(restart+fsiperiod_/dt_+1e-14);

    if (actstep <= crit_step) //the last period is written in the file we have restarted from
    {
      filename = DRT::Problem::Instance()->InputControlFile()->FileName();
    }
    else //the last period is written in the newly written output file
    {
      filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
    }
  }
  else
  {
    filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  }
  return filename;
}

/*----------------------------------------------------------------------*
 | Set time and step in FSI and all subfields                Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SetTimeStepInFSI(const double time, const int step)
{
  // Set time and step in FSI. This looks a bit strange but in case of subcyling
  // we want a 'proper' screen output. And since the fsi time and step are only used
  // for this purpose, it's not really important what we do here
  fsi_->SetTimeStep(time,step);
  //Note: The last function did not touch the subfields, so we have to do it yourself

  // Set time and step in structure field
  fsi_->StructureField()->SetTime(time-fsi_->Dt());
  fsi_->StructureField()->SetTimen(time);
  fsi_->StructureField()->SetStep(step-1);
  fsi_->StructureField()->SetStepn(step);

  // Set time and step in fluid field
  fsi_->FluidField()->SetTimeStep(time,step);

  // Set time and step in ale field
  fsi_->AleField()->SetTimeStep(time,step);
}

/*----------------------------------------------------------------------*
 | Do a single scatra step                                   Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleDoScatraStep()
{
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                        AC COUPLED SCATRA SOLVER"
               "\n************************************************************************\n"<<std::endl;

    std::cout<<"+- step/max -+- tol ---- [norm] -+-- scal-res --+-- scal-inc --+"<<std::endl;
  }

  bool stopnonliniter=false;
  int itnum = 0;

  while (stopnonliniter==false)
  {
    ScatraEvaluateSolveIterUpdate();
    itnum++;
    if (ScatraConvergenceCheck(itnum))
      break;
  }

}

/*----------------------------------------------------------------------*
 | Decide if fluid scatra problem is steady state            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::IsScatraPeriodic( )
{
  const double rel_scatra_error = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("FLUID_SCATRA_REL_TOL");

  //Note: we compare the difference of phinp and phin only on the FS3I interface, since this is the interesting part for the multiscale problem
  const Teuchos::RCP<const Epetra_Vector> phinp_boundary = scatrafieldexvec_[0]->ExtractVector( scatravec_[0]->ScaTraField()->Phinp() ,1); //fluidscatra
  const Teuchos::RCP<const Epetra_Vector> phin_boundary =  scatrafieldexvec_[0]->ExtractVector( fluidphinp_lp_  ,1); //fluidscatra
  const Teuchos::RCP<Epetra_Vector> phidiff_boundary = scatrafieldexvec_[0]->ExtractVector( LINALG::CreateVector(*scatravec_[0]->ScaTraField()->DofRowMap(0),true),1);

  phidiff_boundary->Update(1.0,*phinp_boundary,-1.0,*phin_boundary,0.0);

  double phin_boundary_norm(0.0);
  phin_boundary->Norm2(&phin_boundary_norm);
  if (phin_boundary_norm < 1e-05) //e.g. at the first time step this will happen...
  {
    phin_boundary_norm = 1.0;   //... so we simply check the absolut error.
  }

  double phidiff_boundary_norm(0.0);
  phidiff_boundary->Norm2(&phidiff_boundary_norm);

  //compute the actual rel difference per timestep
  const int numscal = scatravec_[0]->ScaTraField()->NumScal(); //numscal of fluid scatra field

  double fluid_scatra_rel_error = phidiff_boundary_norm/phin_boundary_norm/numscal;

  if (Comm().MyPID()==0)
  {
    std::cout<<std::scientific<<std::setprecision(2)<<"The relative norm of the fluid scatra time period increment is "<<fluid_scatra_rel_error<<std::endl;
  }

  *fluidphinp_lp_=*scatravec_[0]->ScaTraField()->Phinp();
  scatraisperiodic_= (fluid_scatra_rel_error < rel_scatra_error);

  return;
}

/*----------------------------------------------------------------------*
 | Update and output the small time scale                    Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleUpdateAndOutput()
{
  //It is important to first update the fluid field (including the impedance BC)
  //before checking the fsi periodicity. So be careful with the order here!!

  //FSI update and output
  fsi_->PrepareOutput();
  fsi_->Update();
  FsiOutput();

  //Scatra update and output
  UpdateScatraFields();
  ScatraOutput();
}

/*----------------------------------------------------------------------*
 | Write FSI output                                          Thon 03/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::FsiOutput()
{
  int coupling = Teuchos::getIntegralValue<int>(  DRT::Problem::Instance()->FSIDynamicParams() ,"COUPALGO");

  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  const int uprestart = fs3idyn.get<int>("RESTARTEVRY");
  const int upres = fs3idyn.get<int>("UPRES");

  //    /* Note: The order is important here! In here control file entries are
  //     * written. And these entries define the order in which the filters handle
  //     * the Discretizations, which in turn defines the dof number ordering of the
  //     * Discretizations.

  //structure output
  fsi_->StructureField()->Output();
  if (coupling == fsi_iter_monolithicstructuresplit)
    fsi_->OutputLambda();

  //fluid output
  fsi_->FluidField()->    Output();
  if (coupling == fsi_iter_monolithicfluidsplit)
    fsi_->OutputLambda();
  //AC specific fluid output
  if ((uprestart != 0 && step_ % uprestart == 0) || step_ % upres == 0)
  {
    fsi_->FluidField()->DiscWriter()->WriteVector("wss", WallShearStress_lp_); //TODO: is it enough to this vector only at the beginning of a new fsi period?
  }

  //ale output
  fsi_->AleField()->Output();

  fsi_->FluidField()->LiftDrag();
}

/*----------------------------------------------------------------------*
 | check convergence of scatra fields                        Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::ScatraConvergenceCheck(const int itnum)
{

  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int scatraitemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double scatraittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double scatraabstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");


  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatraincrement_->Norm2(&incconnorm);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> con = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp(); //fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp(); //structurescatra
  SetupCoupledScatraVector(con,scatra1,scatra2);

  double connorm(0.0);
  con->Norm2(&connorm);

  // care for the case that nothing really happens in the concentration field
  if (connorm < 1e-5)
    connorm = 1.0;

  // print the screen info
  if (Comm().MyPID()==0)
  {
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |\n",
           itnum,scatraitemax,scatraittol,conresnorm,incconnorm/connorm);
  }

  // this is the convergence check
  // We always require at least one solve. We test the L_2-norm of the
  // current residual. Norm of residual is just printed for information
  if (conresnorm <= scatraittol and incconnorm/connorm <= scatraittol)
  {
    if (Comm().MyPID()==0)
    {
      // print 'finish line'
      printf("+------------+-------------------+--------------+--------------+\n\n");
    }
    return true;
  }

  // abort iteration, when there's nothing more to do! -> more robustness
  else if (conresnorm < scatraabstolres)
  {
    // print 'finish line'
    if (Comm().MyPID()==0)
    {
      printf("+------------+-------------------+--------------+--------------+\n\n");
    }
    return true;
  }

  // if itemax is reached without convergence stop the simulation
  else if (itnum == scatraitemax)
  {
    if (Comm().MyPID()==0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|    scalar-scalar field did not converge in itemax steps!     |\n");
      printf("+---------------------------------------------------------------+\n");
    }
    // yes, we stop!
    dserror("Scatra not converged in itemax steps!");
    return true;
  }
  else
    return false;

}

/*----------------------------------------------------------------------*
 | Convergence check for iterative staggered FS3I scheme     Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::PartFs3iConvergenceCkeck( const int itnum)
{
  const Teuchos::ParameterList& fs3idynpart = DRT::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED");
  // get control parameters from input file
  double ittol = fs3idynpart.get<double>("CONVTOL");
  int itmax = fs3idynpart.get<int>("ITEMAX");

  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //calculate fsi increments. scatra increment is already done in the scatra fields convergence check
  structureincrement_->Update(1.0,*fsi_->StructureField()->Dispnp(),-1.0);
  fluidincrement_->Update(1.0,*fsi_->FluidField()->Velnp(),-1.0);
  aleincrement_->Update(1.0,*fsi_->AleField()->Dispnp(),-1.0);

  //L2-norms of increment vectors
  double scatraincconnorm_L2(0.0);
  scatraincrement_->Norm2(&scatraincconnorm_L2);
  double structureincconnorm_L2(0.0);
  structureincrement_->Norm2(&structureincconnorm_L2);
  double fluidincconnorm_L2(0.0);
  fluidincrement_->Norm2(&fluidincconnorm_L2);
  double aleincconnorm_L2(0.0);
  aleincrement_->Norm2(&aleincconnorm_L2);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> scatra = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp(); //fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp(); //structurescatra
  SetupCoupledScatraVector(scatra,scatra1,scatra2);

  //norms of solution vectors
  double scatranorm_L2(0.0);
  scatra->Norm2(&scatranorm_L2);
  if (scatranorm_L2 < 1e-05) scatranorm_L2=1.0;
  double structurenorm_L2(0.0);
  fsi_->StructureField()->Dispnp()->Norm2(&structurenorm_L2);
  if (structurenorm_L2 < 1e-05) structurenorm_L2=1.0;
  double fluidnorm_L2(0.0);
  fsi_->FluidField()->Velnp()->Norm2(&fluidnorm_L2);
  if (fluidnorm_L2 < 1e-05) fluidnorm_L2=1.0;
  double alenorm_L2(0.0);
  fsi_->AleField()->Dispnp()->Norm2(&alenorm_L2);
  if (alenorm_L2 < 1e-05) alenorm_L2=1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n************************************************************************\n"
                 "                       OUTER ITERATION STEP "<<itnum<<"/"<<itmax<<std::endl;

    printf("+--------------+---------------------+----------------+---------------+---------------+--------------+\n");
    printf("|   step/max   |   tol      [norm]   |   scalar-inc   |   disp-inc    |   vel-inc     |   ale-inc    |\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E   |  %10.3E   | %10.3E   |\n",
         itnum,itmax,ittol,scatraincconnorm_L2/scatranorm_L2,structureincconnorm_L2/structurenorm_L2,fluidincconnorm_L2/fluidnorm_L2,aleincconnorm_L2/alenorm_L2);
    printf("+--------------+---------------------+----------------+---------------+---------------+--------------+\n");
    std::cout<<"************************************************************************\n";
  }

  if ( scatraincconnorm_L2/scatranorm_L2 <= ittol and structureincconnorm_L2/structurenorm_L2 <= ittol
     and fluidincconnorm_L2/fluidnorm_L2 <= ittol and aleincconnorm_L2/alenorm_L2 <= ittol ) //converged!
   {
     stopnonliniter = true;
     if ((Comm().MyPID()==0) )
     {
       std::cout<<"\n************************************************************************\n"
                    "                     OUTER ITERATION STEP CONVERGED"
                  "\n************************************************************************\n"<<std::endl;
     }
   }
  else if (itnum == itmax)
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      std::cout<<"\n************************************************************************\n"
                   "           OUTER ITERATION STEP NOT CONVERGED IN ITEMAX STEPS"
                 "\n************************************************************************\n"<<std::endl;
    }
    dserror("The partitioned FS3I solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 |  Extract wall shear stresses                              Thon 03/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss)
{
  Teuchos::RCP<Epetra_Vector> WallShearStress;

  //############ Fluid Field ###############
  if ( not fsiisperiodic_ )
  {
    Teuchos::RCP<ADAPTER::FluidACFSI> fluid = Teuchos::rcp_dynamic_cast<ADAPTER::FluidACFSI>( fsi_->FluidField() );
    if (fluid == Teuchos::null)
      dserror("Dynamic cast to ADAPTER::FluidFSI failed!");

    if ( step_>1 and ModuloIsRealtiveZero(time_-dt_,fsiperiod_,time_) )
      {
        fluid->ResetHistoryVectors( ); //Reset StressManager
        if (Comm().MyPID()==0)
        {
        std::cout<<"\nReseting stress manager"<<std::endl;
        }
      }

    WallShearStress = fluid->CalculateWallShearStresses(); //Instantaneous wss vector

    switch ( DRT::INPUT::IntegralValue<INPAR::FLUID::WSSType>( DRT::Problem::Instance()->FluidDynamicParams() ,"WSS_TYPE") )
    {
    case INPAR::FLUID::wss_standard: //we simply use the instantaneous WSS
      WallShearStress_lp_->Update(1.0,*WallShearStress,0.0); //Update
      break;
    case INPAR::FLUID::wss_mean: //we use the mean WSS from the last period..
      if ( step_>1 and ModuloIsRealtiveZero(time_,fsiperiod_,time_) )
      {
        WallShearStress_lp_->Update(1.0,*WallShearStress,0.0); //Update
        if (Comm().MyPID()==0)
        {
          std::cout<<"\nUpdating mean wss vector in ac-fs3i"<<std::endl;
        }
        //Note: we can not reset the stress manager here, since we first need to write its data, to be able to restart.
        //Hence the correct order is: 1. Calculate WSS; 2. Write fluid output; 3. Reset stress Manager(in next time step)
      }
      break;
    default:
      dserror("WSS_TYPE not supported for AC-FS3I!");
      break;
    }
  }
  //else //nothing to do in this case since the WallShearStress_ vector has been read from on period before

  wss.push_back(WallShearStress_lp_);

  //############ Structure Field ###############

  // extract FSI-Interface from fluid field
  WallShearStress = fsi_->FluidField()->Interface()->ExtractFSICondVector(WallShearStress_lp_);

  // replace global fluid interface dofs through structure interface dofs
  WallShearStress = fsi_->FluidToStruct(WallShearStress);

  // insert structure interface entries into vector with full structure length
  Teuchos::RCP<Epetra_Vector> structure = LINALG::CreateVector(*(fsi_->StructureField()->Interface()->FullMap()),true);

  //Parameter int block of function InsertVector: (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
  fsi_->StructureField()->Interface()->InsertVector(WallShearStress,1,structure);
  wss.push_back(structure);
}
