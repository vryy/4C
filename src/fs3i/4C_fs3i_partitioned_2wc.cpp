/*----------------------------------------------------------------------*/
/*! \file
\brief Algorithmic routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I) specifically
       related to two-way-coupled problem configurations

\level 2



*----------------------------------------------------------------------*/


#include "4C_fs3i_partitioned_2wc.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fluid_timint_loma.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_loma.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I2Wc::PartFS3I2Wc(const Epetra_Comm& comm)
    : PartFS3I(comm),
      itmax_(GLOBAL::Problem::Instance()
                 ->FS3IDynamicParams()
                 .sublist("PARTITIONED")
                 .get<int>("ITEMAX")),
      ittol_(GLOBAL::Problem::Instance()
                 ->FS3IDynamicParams()
                 .sublist("PARTITIONED")
                 .get<double>("CONVTOL")),
      consthermpress_(
          GLOBAL::Problem::Instance()->FS3IDynamicParams().get<std::string>("CONSTHERMPRESS"))
{
  // constructor is supposed to stay empty
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::Init()
{
  // call Init() in base class
  FS3I::PartFS3I::Init();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::Setup()
{
  // call Setup() in base class
  FS3I::PartFS3I::Setup();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::Timeloop()
{
  check_is_init();
  check_is_setup();

  initial_calculations();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetFSISolution();

  // calculate inital time derivative, when restart was done from a part. FSI simulation
  if (GLOBAL::Problem::Instance()->Restart() and
      CORE::UTILS::IntegralValue<int>(
          GLOBAL::Problem::Instance()->FS3IDynamicParams(), "RESTART_FROM_PART_FSI"))
  {
    scatravec_[0]->ScaTraField()->prepare_first_time_step();
    scatravec_[1]->ScaTraField()->prepare_first_time_step();
  }

  // output of initial state
  constexpr bool force_prepare = true;
  fsi_->prepare_output(force_prepare);
  fsi_->output();
  ScatraOutput();

  while (NotFinished())
  {
    increment_time_and_step();

    prepare_time_step();

    outer_loop();

    TimeUpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::initial_calculations()
{
  // set initial fluid velocity field for evaluation of initial scalar
  // time derivative in fluid-based scalar transport
  scatravec_[0]->ScaTraField()->set_velocity_field(
      fsi_->fluid_field()->Velnp(), Teuchos::null, Teuchos::null, Teuchos::null);

  // set initial value of thermodynamic pressure in fluid-based scalar
  // transport
  // for constant thermodynamic pressure in low-Mach-number flow and
  // for temperature-dependent water, thermodynamic pressure is set
  // to a constant here and never touched again
  if ((fsi_->fluid_field()->PhysicalType() == INPAR::FLUID::tempdepwater) &&
      (consthermpress_ != "Yes"))
    FOUR_C_THROW(
        "Constant thermodynamic pressure required if TFSI algorithm is used with "
        "temperature-dependent water!");
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
      ->set_initial_therm_pressure();

  // mass conservation: compute initial mass (initial time deriv. assumed zero)
  if (consthermpress_ == "No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
        ->ComputeInitialMass();

  // set initial scalar field and thermodynamic pressure for evaluation of
  // Neumann boundary conditions in fluid at beginning of first time step
  fsi_->fluid_field()->SetScalarFields(scatravec_[0]->ScaTraField()->Phinp(),
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
          ->ThermPressNp(),
      Teuchos::null, scatravec_[0]->ScaTraField()->discretization());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::prepare_time_step()
{
  check_is_init();
  check_is_setup();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->prepare_time_step();
  }

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_ == "No_energy")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
        ->predict_therm_pressure();

  // prepare time step for fluid, structure and ALE fields
  fsi_->prepare_time_step();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::outer_loop()
{
  // set iteration number and stop criterion
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
        scatravec_[0]->ScaTraField()->Time(), timemax_, dt_,
        scatravec_[0]->ScaTraField()->MethodTitle().c_str(), scatravec_[0]->ScaTraField()->Step(),
        numstep_);
  }

  // set mesh displacement and velocity fields
  // TO DO: temporally consistent transfer of velocity and other fields,
  // for the time being, zero velocity field from structure
  SetFSISolution();

  // initially solve coupled scalar transport equation system
  // (values for intermediate time steps were calculated at the end of prepare_time_step)
  if (Comm().MyPID() == 0)
    std::cout << "\n****************************************\n        SCALAR TRANSPORT "
                 "SOLVER\n****************************************\n";
  scatra_evaluate_solve_iter_update();

  while (stopnonliniter == false)
  {
    itnum++;

    // in case of non-constant thermodynamic pressure: compute
    // (either based on energy conservation or based on mass conservation)
    if (consthermpress_ == "No_energy")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
          ->compute_therm_pressure();
    else if (consthermpress_ == "No_mass")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
          ->compute_therm_pressure_from_mass_cons();

    // set fluid- and structure-based scalar transport values required in FSI
    set_sca_tra_values_in_fsi();

    // solve FSI system
    if (Comm().MyPID() == 0)
      std::cout << "\n****************************************\n               FSI "
                   "SOLVER\n****************************************\n";
    fsi_->TimeStep(fsi_);

    // set mesh displacement and velocity fields
    // TO DO: temporally consistent transfer of velocity and other fields,
    // for the time being, zero velocity field from structure
    SetFSISolution();

    // solve scalar transport equation
    if (Comm().MyPID() == 0)
      std::cout << "\n****************************************\n        SCALAR TRANSPORT "
                   "SOLVER\n****************************************\n";
    scatra_evaluate_solve_iter_update();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::set_sca_tra_values_in_fsi()
{
  // set scalar and thermodynamic pressure values as well as time
  // derivatives and discretization from fluid scalar in fluid
  // based on time-integration scheme
  // TO DO: check dynamic cast with "Iter" routine
  switch (fsi_->fluid_field()->TimIntScheme())
  {
    case INPAR::FLUID::timeint_afgenalpha:
    {
      if (fsi_->fluid_field()->PhysicalType() == INPAR::FLUID::tempdepwater)
        fsi_->fluid_field()->SetIterScalarFields(
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiaf()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiam()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtam()),
            scatravec_[0]->ScaTraField()->discretization());
      else
        fsi_->fluid_field()->set_loma_iter_scalar_fields(
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiaf()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phiam()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtam()), Teuchos::null,
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressAf(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressAm(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressDtAf(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressDtAm(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->discretization());
    }
    break;
    case INPAR::FLUID::timeint_one_step_theta:
    {
      if (fsi_->fluid_field()->PhysicalType() == INPAR::FLUID::tempdepwater)
        fsi_->fluid_field()->SetIterScalarFields(
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phinp()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phin()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtnp()),
            scatravec_[0]->ScaTraField()->discretization());
      else
        fsi_->fluid_field()->set_loma_iter_scalar_fields(
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phinp()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phin()),
            FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phidtnp()), Teuchos::null,
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressNp(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressN(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressDtNp(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->ThermPressDtNp(),
            Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
                ->discretization());
    }
    break;
    default:
      FOUR_C_THROW("Time integration scheme not supported");
      break;
  }

  // set structure-scalar field in structure
  // (Note potential inconsistencies related to this call in case of generalized-alpha time
  // integration!)
  fsi_->structure_field()->discretization()->set_state(
      1, "scalarfield", structure_scalar_to_structure(scatravec_[1]->ScaTraField()->Phinp()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I2Wc::convergence_check(int itnum)
{
  // define flags for fluid and scatra convergence check
  bool fluidstopnonliniter = false;
  bool scatrastopnonliniter = false;

  // dump on screen
  if (Comm().MyPID() == 0)
    std::cout << "\n****************************************\n  CONVERGENCE CHECK FOR ITERATION "
                 "STEP\n****************************************\n";

  // fsi convergence check
  if (fsi_->NoxStatus() == ::NOX::StatusTest::Converged) fluidstopnonliniter = true;

  // scatra convergence check
  scatrastopnonliniter = scatra_convergence_check(itnum);

  // warn if itemax is reached without convergence of FSI solver,
  // but proceed to next timestep
  if ((itnum == itmax_) and (fluidstopnonliniter == false))
  {
    fluidstopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf(">>>>>> FSI solver not converged in itemax steps!\n");
      printf("\n");
    }
  }

  if (fluidstopnonliniter == true and scatrastopnonliniter == true)
    return true;
  else
    return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I2Wc::scatra_convergence_check(int itnum)
{
  // define flags for convergence check for scatra fields
  bool scatra1stopnonliniter = false;
  bool scatra2stopnonliniter = false;

  // convergence check of scatra fields
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n         SCALAR TRANSPORT "
                 "CHECK\n****************************************\n";
    std::cout << "\n****************************************\n   FLUID-BASED SCALAR TRANSPORT "
                 "CHECK\n****************************************\n";
  }
  scatra1stopnonliniter =
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
          ->convergence_check(itnum, itmax_, ittol_);

  if (Comm().MyPID() == 0)
    std::cout << "\n****************************************\n STRUCTURE-BASED SCALAR TRANSPORT "
                 "CHECK\n****************************************\n";
  // FOUR_C_THROW("convergence_check in scatra currently only for loma scatra!Fix this!");
  // scatra2stopnonliniter = scatravec_[1]->ScaTraField()->convergence_check(itnum,itmax_,ittol_);
  scatra2stopnonliniter =
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[1]->ScaTraField())
          ->convergence_check(itnum, itmax_, ittol_);

  if (scatra1stopnonliniter == true and scatra2stopnonliniter == true)
    return true;
  else
    return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I2Wc::TimeUpdateAndOutput()
{
  // prepare output for FSI
  constexpr bool force_prepare = false;
  fsi_->prepare_output(force_prepare);

  // update fluid- and structure-based scalar transport
  UpdateScatraFields();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_ == "No_energy" or consthermpress_ == "No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
        ->UpdateThermPressure();

  // update structure, fluid and ALE
  fsi_->update();

  // set scalar and thermodynamic pressure at n+1 and SCATRA trueresidual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  fsi_->fluid_field()->SetScalarFields(FluidScalarToFluid(scatravec_[0]->ScaTraField()->Phinp()),
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(scatravec_[0]->ScaTraField())
          ->ThermPressNp(),
      FluidScalarToFluid(scatravec_[0]->ScaTraField()->TrueResidual()),
      scatravec_[0]->ScaTraField()->discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  fsi_->output();

  // output of fluid- and structure-based scalar transport
  ScatraOutput();

  return;
}

FOUR_C_NAMESPACE_CLOSE
