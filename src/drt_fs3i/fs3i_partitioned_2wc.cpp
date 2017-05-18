/*!----------------------------------------------------------------------
\file fs3i_partitioned_2wc.cpp
\brief Algorithmic routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I) specifically
       related to two-way-coupled problem configurations

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


*----------------------------------------------------------------------*/


#include "fs3i_partitioned_2wc.H"

#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/fluid_timint_loma.H"

#include "../drt_scatra/scatra_timint_loma.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I_2WC::PartFS3I_2WC(const Epetra_Comm& comm)
  : PartFS3I(comm),
    itmax_(DRT::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED").get<int>("ITEMAX")),
    ittol_(DRT::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED").get<double>("CONVTOL")),
    consthermpress_(DRT::Problem::Instance()->FS3IDynamicParams().get<std::string>("CONSTHERMPRESS"))
{
  // constructor is supposed to stay empty
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Init()
{
  // call Init() in base class
  FS3I::PartFS3I::Init();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Setup()
{
  // call Setup() in base class
  FS3I::PartFS3I::Setup();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Timeloop()
{
  CheckIsInit();
  CheckIsSetup();

  InitialCalculations();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetFSISolution();

  //calculate inital time derivative, when restart was done from a part. FSI simulation
  if ( DRT::Problem::Instance()->Restart() and DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->FS3IDynamicParams(),"RESTART_FROM_PART_FSI") )
  {
    scatravec_[0]->ScaTraField()->PrepareFirstTimeStep();
    scatravec_[1]->ScaTraField()->PrepareFirstTimeStep();
  }

  // output of initial state
  fsi_->PrepareOutput();
  fsi_->Output();
  ScatraOutput();

  while (NotFinished())
  {
    IncrementTimeAndStep();

    PrepareTimeStep();

    OuterLoop();

    TimeUpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::InitialCalculations()
{
  // set initial fluid velocity field for evaluation of initial scalar
  // time derivative in fluid-based scalar transport
  scatravec_[0]->ScaTraField()->SetVelocityField(fsi_->FluidField()->Velnp(),
                                               Teuchos::null,
                                               Teuchos::null,
                                               Teuchos::null,
                                               1);

  // set initial value of thermodynamic pressure in fluid-based scalar
  // transport
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->SetInitialThermPressure();

  // energy conservation: compute initial time derivative of therm. pressure
  // mass conservation: compute initial mass (initial time deriv. assumed zero)
  if (consthermpress_=="No_energy")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ComputeInitialThermPressureDeriv();
  else if (consthermpress_=="No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ComputeInitialMass();

  // set initial scalar field and thermodynamic pressure for evaluation of
  // Neumann boundary conditions in fluid at beginning of first time step
  fsi_->FluidField()->SetScalarFields(scatravec_[0]->ScaTraField()->Phinp(),
                                       Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressNp(),
                                       Teuchos::null,
                                       scatravec_[0]->ScaTraField()->Discretization());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::PrepareTimeStep()
{
  CheckIsInit();
  CheckIsSetup();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->PrepareTimeStep();
  }

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->PredictThermPressure();

  // prepare time step for fluid, structure and ALE fields
  fsi_->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::OuterLoop()
{
  // set iteration number and stop criterion
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           scatravec_[0]->ScaTraField()->Time(),timemax_,dt_,scatravec_[0]->ScaTraField()->MethodTitle().c_str(),scatravec_[0]->ScaTraField()->Step(),numstep_);
  }

  // set mesh displacement and velocity fields
  // TO DO: temporally consistent transfer of velocity and other fields,
  // for the time being, zero velocity field from structure
  SetFSISolution();

  // initially solve coupled scalar transport equation system
  // (values for intermediate time steps were calculated at the end of PrepareTimeStep)
  if (Comm().MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
  ScatraEvaluateSolveIterUpdate();

  while (stopnonliniter==false)
  {
    itnum++;

    // in case of non-constant thermodynamic pressure: compute
    // (either based on energy conservation or based on mass conservation)
    if (consthermpress_=="No_energy")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ComputeThermPressure();
    else if (consthermpress_=="No_mass")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ComputeThermPressureFromMassCons();

    // set fluid- and structure-based scalar transport values required in FSI
    SetScaTraValuesInFSI();

    // solve FSI system
    if (Comm().MyPID()==0) std::cout<<"\n****************************************\n               FSI SOLVER\n****************************************\n";
    fsi_->TimeStep(fsi_);

    // set mesh displacement and velocity fields
    // TO DO: temporally consistent transfer of velocity and other fields,
    // for the time being, zero velocity field from structure
    SetFSISolution();

    // solve scalar transport equation
    if (Comm().MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
    ScatraEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::SetScaTraValuesInFSI()
{
  // set scalar and thermodynamic pressure values as well as time
  // derivatives and discretization from fluid scalar in fluid
  // based on time-integration scheme
  // TO DO: check dynamic cast with "Iter" routine
  switch(fsi_->FluidField()->TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  {
    fsi_->FluidField()->SetLomaIterScalarFields(FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiaf()),
                                                FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiam()),
                                                FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtam()),
                                                Teuchos::null,
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressAf(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressAm(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressDtAf(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressDtAm(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->Discretization());
  }
  break;
  case INPAR::FLUID::timeint_one_step_theta:
  {
    fsi_->FluidField()->SetLomaIterScalarFields(FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phinp()),
                                                FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phin()),
                                                FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtnp()),
                                                Teuchos::null,
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressNp(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressN(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressDtNp(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressDtNp(),
                                                Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->Discretization());
  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }

  // set structure-scalar field in structure
  // (Note potential inconsistencies related to this call in case of generalized-alpha time integration!)
  fsi_->StructureField()->Discretization()->SetState( 1,"temperature",StructureScalarToStructure(scatravec_[1]->ScaTraField()->Phinp()) );

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_2WC::ConvergenceCheck(int itnum)
{
  // define flags for fluid and scatra convergence check
  bool fluidstopnonliniter  = false;
  bool scatrastopnonliniter = false;

  // dump on screen
  if (Comm().MyPID()==0) std::cout<<"\n****************************************\n  CONVERGENCE CHECK FOR ITERATION STEP\n****************************************\n";

  // fsi convergence check
  if (fsi_->NoxStatus() == NOX::StatusTest::Converged) fluidstopnonliniter = true;

  // scatra convergence check
  scatrastopnonliniter = ScatraConvergenceCheck(itnum);

  // warn if itemax is reached without convergence of FSI solver,
  // but proceed to next timestep
  if ((itnum == itmax_) and (fluidstopnonliniter == false))
  {
    fluidstopnonliniter=true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf(">>>>>> FSI solver not converged in itemax steps!\n");
      printf("\n");
    }
  }

  if (fluidstopnonliniter == true and scatrastopnonliniter == true) return true;
  else                                                              return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_2WC::ScatraConvergenceCheck(int itnum)
{
  // define flags for convergence check for scatra fields
  bool scatra1stopnonliniter = false;
  bool scatra2stopnonliniter = false;

  // convergence check of scatra fields
  if (Comm().MyPID() == 0)
  {
    std::cout<<"\n****************************************\n         SCALAR TRANSPORT CHECK\n****************************************\n";
    std::cout<<"\n****************************************\n   FLUID-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
  }
  scatra1stopnonliniter = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ConvergenceCheck(itnum,itmax_,ittol_);

  if (Comm().MyPID() == 0) std::cout<<"\n****************************************\n STRUCTURE-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
  //dserror("ConvergenceCheck in scatra currently only for loma scatra!Fix this!");
  //scatra2stopnonliniter = scatravec_[1]->ScaTraField()->ConvergenceCheck(itnum,itmax_,ittol_);
  scatra2stopnonliniter = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[1]->ScaTraField())->ConvergenceCheck(itnum,itmax_,ittol_);

  if (scatra1stopnonliniter == true and scatra2stopnonliniter == true) return true;
  else                                                                 return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::TimeUpdateAndOutput()
{
  // prepare output for FSI
  fsi_->PrepareOutput();

  // update fluid- and structure-based scalar transport
  UpdateScatraFields();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->UpdateThermPressure();

  // update structure, fluid and ALE
  fsi_->Update();

  // set scalar and thermodynamic pressure at n+1 and SCATRA trueresidual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  fsi_->FluidField()->SetScalarFields(FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phinp()),
                                      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())->ThermPressNp(),
                                      FluidScalarToFluid(scatravec_[0]->ScaTraField()->TrueResidual()),
                                      scatravec_[0]->ScaTraField()->Discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  fsi_->Output();

  // output of fluid- and structure-based scalar transport
  ScatraOutput();

  return;
}


