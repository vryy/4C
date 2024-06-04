/*----------------------------------------------------------------------*/
/*! \file

\brief Transport of passive scalars in Navier-Stokes velocity field

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_algorithm.hpp"

#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fluid_turbulence_statistic_manager.hpp"
#include "4C_global_data.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraAlgorithm::ScaTraAlgorithm(const Epetra_Comm& comm,  ///< communicator
    const Teuchos::ParameterList& scatradyn,                       ///< scatra parameter list
    const Teuchos::ParameterList& fdyn,                            ///< fluid parameter list
    const std::string scatra_disname,                              ///< scatra discretization name
    const Teuchos::ParameterList& solverparams                     ///< solver parameter list
    )
    : ScaTraFluidCouplingAlgorithm(comm, scatradyn, false, scatra_disname, solverparams),
      natconv_(CORE::UTILS::IntegralValue<int>(scatradyn, "NATURAL_CONVECTION")),
      natconvitmax_(scatradyn.sublist("NONLINEAR").get<int>("ITEMAX_OUTER")),
      natconvittol_(scatradyn.sublist("NONLINEAR").get<double>("CONVTOL_OUTER")),
      velincnp_(Teuchos::null),
      phiincnp_(Teuchos::null),
      samstart_(fdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_START")),
      samstop_(fdyn.sublist("TURBULENCE MODEL").get<int>("SAMPLING_STOP"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Setup()
{
  // call setup in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  // create vectors
  velincnp_ = Teuchos::rcp(
      new Epetra_Vector(*(fluid_field()->ExtractVelocityPart(fluid_field()->Velnp()))));
  phiincnp_ = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())));

  if (velincnp_ == Teuchos::null) FOUR_C_THROW("velincnp_ == Teuchos::null");
  if (phiincnp_ == Teuchos::null) FOUR_C_THROW("phiincnp_ == Teuchos::null");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::Init()
{
  // call init in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::TimeLoop()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // provide information about initial field
  prepare_time_loop();

  if (!natconv_)
    time_loop_one_way();  // one-way coupling
  else
    time_loop_two_way();  // two-way coupling (natural convection)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::time_loop_one_way()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // write velocities since they may be needed in prepare_first_time_step() of the scatra field
  set_velocity_field();
  // write initial output to file
  if (Step() == 0) output();

  // time loop (no-subcycling at the moment)
  while (NotFinished())
  {
    // prepare next time step
    prepare_time_step();

    // solve nonlinear Navier-Stokes system
    do_fluid_step();

    // solve transport equations
    do_transport_step();

    // update all field solvers
    update();

    // compute error for problems with analytical solution
    ScaTraField()->evaluate_error_compared_to_analytical_sol();

    // write output to screen and files
    output();
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::time_loop_two_way()
{
  // safety checks
  check_is_init();
  check_is_setup();

  prepare_time_loop_two_way();

  // time loop
  while (NotFinished())
  {
    increment_time_and_step();

    // prepare next time step
    prepare_time_step_convection();

    // Outer iteration loop
    outer_iteration_convection();

    // update all single field solvers and density
    update_convection();

    // write output to screen and files
    output();
  }
}

/*----------------------------------------------------------------------*
 | Initial calculations for two-way coupling time loop       fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::prepare_time_loop_two_way()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // safety check
  switch ((fluid_field()->TimIntScheme()))
  {
    case INPAR::FLUID::timeint_afgenalpha:
    case INPAR::FLUID::timeint_bdf2:
    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_stationary:
      break;
    default:
    {
      FOUR_C_THROW("Selected time integration scheme is not available!");
      break;
    }
  }

  // compute initial mean concentrations and load densification coefficients
  ScaTraField()->SetupNatConv();

  // compute initial density
  ScaTraField()->ComputeDensity();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::prepare_time_step()
{
  // safety checks
  check_is_init();
  check_is_setup();

  increment_time_and_step();

  fluid_field()->prepare_time_step();

  // prepare time step
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField()->prepare_time_step();

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n******************\n   TIME STEP     \n******************\n";
    std::cout << "\nStep:   " << Step() << " / " << n_step() << "\n";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::prepare_time_step_convection()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // OST/BDF2 time integration schemes are implemented
  // pass density to fluid discretization at time steps n+1 and n
  // density time derivative is not used for OST and BDF2 (pass zero vector)
  // thermodynamic pressure values are set to 1.0 and its derivative to 0.0

  switch ((fluid_field()->TimIntScheme()))
  {
    case INPAR::FLUID::timeint_afgenalpha:
    case INPAR::FLUID::timeint_bdf2:
    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_stationary:
    {
      fluid_field()->SetIterScalarFields(ScaTraField()->Densafnp(),
          ScaTraField()->Densafnp(),  // not needed, provided as dummy vector
          Teuchos::null, ScaTraField()->discretization());
      break;
    }
    default:
    {
      FOUR_C_THROW("Selected time integration scheme is not available!");
      break;
    }
  }

  fluid_field()->prepare_time_step();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step() == 1)
    ScaTraField()->set_velocity_field(
        fluid_field()->Velnp(), fluid_field()->Hist(), Teuchos::null, Teuchos::null);

  // prepare time step (+ initialize one-step-theta scheme correctly with
  // velocity given above)
  ScaTraField()->prepare_time_step();
}

/*----------------------------------------------------------------------*
 | Print scatra solver type to screen                        fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::print_sca_tra_solver()
{
  if (Comm().MyPID() == 0)
    std::cout
        << "\n****************************\n      TRANSPORT SOLVER\n****************************\n";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::do_fluid_step()
{
  // solve nonlinear Navier-Stokes system
  if (Comm().MyPID() == 0)
    std::cout << "\n************************\n      FLUID SOLVER\n************************\n";

  // currently only required for forced homogeneous isotropic turbulence with
  // passive scalar transport; does nothing otherwise
  fluid_field()->calc_intermediate_solution();

  fluid_field()->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::do_transport_step()
{
  print_sca_tra_solver();

  // transfer velocities to scalar transport field solver
  set_velocity_field();

  // solve the transport equation(s)
  ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::set_velocity_field()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // transfer velocities to scalar transport field solver
  // NOTE: so far, the convective velocity is chosen to equal the fluid velocity
  //       since it is not yet clear how the grid velocity should be interpolated
  //       properly -> hence, ScaTraAlgorithm does not support moving
  //       meshes yet

  // this is ugly, but FsVel() may give a Null pointer which we canNOT give to the volmortart
  // framework
  // TODO (thon): make this somehow prettier..
  Teuchos::RCP<const Epetra_Vector> fsvel = fluid_field()->FsVel();
  if (fsvel != Teuchos::null) fsvel = fluid_to_scatra(fsvel);

  switch (fluid_field()->TimIntScheme())
  {
    case INPAR::FLUID::timeint_npgenalpha:
    case INPAR::FLUID::timeint_afgenalpha:
    {
      ScaTraField()->set_velocity_field(fluid_to_scatra(fluid_field()->Velaf()),
          fluid_to_scatra(fluid_field()->Accam()), fluid_to_scatra(fluid_field()->Velaf()), fsvel);
      break;
    }
    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_bdf2:
    case INPAR::FLUID::timeint_stationary:
    {
      ScaTraField()->set_velocity_field(fluid_to_scatra(fluid_field()->Velnp()),
          fluid_to_scatra(fluid_field()->Hist()), fluid_to_scatra(fluid_field()->Velnp()), fsvel);
      break;
    }
    default:
    {
      FOUR_C_THROW("Time integration scheme not supported");
      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::outer_iteration_convection()
{
  // safety checks
  check_is_init();
  check_is_setup();

  int natconvitnum = 0;
  bool stopnonliniter = false;

  // Outer Iteration loop starts
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "**************************************************************\n";
    std::cout << "      OUTER ITERATION LOOP (" << ScaTraField()->MethodTitle() << ")\n";
    printf("      Time Step %3d/%3d \n", Step(), ScaTraField()->NStep());
    std::cout << "**************************************************************\n";
  }

#ifdef OUTPUT
  // Output after each Outer Iteration step
  const int numdim = 3;
  // create output file name
  std::stringstream temp;
  temp << GLOBAL::Problem::Instance()->OutputControlFile()->FileName() << "_nonliniter_step"
       << Step();
  std::string outname = temp.str();
  std::string probtype = GLOBAL::Problem::Instance()->ProblemName();

  Teuchos::RCP<CORE::IO::OutputControl> myoutputcontrol =
      Teuchos::rcp(new CORE::IO::OutputControl(ScaTraField().discretization()->Comm(), probtype,
          CORE::FE::ShapeFunctionType::polynomial, "myinput", outname, numdim, 0, 1000));
  // create discretization writer with my own control settings
  Teuchos::RCP<CORE::IO::DiscretizationWriter> myoutput = ScaTraField().discretization()->Writer();
  myoutput->SetOutput(myoutputcontrol);
  // write mesh at step 0
  myoutput->WriteMesh(0, 0.0);
#endif

  while (!stopnonliniter)
  {
    natconvitnum++;

    phiincnp_->Update(1.0, *ScaTraField()->Phinp(), 0.0);
    velincnp_->Update(1.0, *fluid_field()->ExtractVelocityPart(fluid_field()->Velnp()), 0.0);

    // solve nonlinear Navier-Stokes system with body forces
    do_fluid_step();

    // solve scalar transport equation
    do_transport_step();

    // compute new density field and pass it to the fluid discretization
    ScaTraField()->ComputeDensity();
    fluid_field()->SetScalarFields(
        ScaTraField()->Densafnp(), 0.0, Teuchos::null, ScaTraField()->discretization());

    // convergence check based on incremental values
    stopnonliniter = convergence_check(natconvitnum, natconvitmax_, natconvittol_);

    // Test: Output of the flux across the boundary into the output text file
    // after each outer Iteration step
    // print mean concentration
#ifdef OUTPUT
    if (stopnonliniter == false)
    {
      printf("\n");
      printf("Flux: Outer Iterations step: %3d \n", natconvitnum);
      ScaTraField().output_problem_specific();
      ScaTraField().output_total_and_mean_scalars();
      printf("\n");
    }

    myoutput->NewStep(natconvitnum, natconvitnum);
    myoutput->WriteVector("phinp", ScaTraField().Phinp());
    myoutput->WriteVector("convec_velocity", ScaTraField().ConVel());
    // myoutput->WriteVector("velnp", fluid_field()->Velnp());
#endif
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::update()
{
  fluid_field()->Update();
  ScaTraField()->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::update_convection()
{
  // update scatra and fluid fields
  ScaTraField()->Update();
  fluid_field()->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.

  if ((Step() >= samstart_) and (Step() <= samstop_))
  {
    // if statistics for one-way coupled problems is performed, provide
    // the field for the first scalar!
    fluid_field()->SetScalarFields(
        ScaTraField()->Phinp(), 0.0, Teuchos::null, ScaTraField()->discretization(),
        0  // do statistics for FIRST dof at every node!!
    );
  }

  fluid_field()->StatisticsAndOutput();
  ScaTraField()->check_and_write_output_and_restart();

  // we have to call the output of averaged fields for scatra separately
  if (fluid_field()->turbulence_statistic_manager() != Teuchos::null)
    fluid_field()->turbulence_statistic_manager()->DoOutputForScaTra(
        *ScaTraField()->DiscWriter(), ScaTraField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScaTraAlgorithm::convergence_check(
    int natconvitnum, const int natconvitmax, const double natconvittol)
{
  // convergence check based on phi and velocity increment

  //   | phi increment |_2
  //  -------------------------- < Tolerance
  //     | phi_n+1 |_2

  bool stopnonliniter = false;
  Teuchos::RCP<CORE::LINALG::MapExtractor> phisplitter = ScaTraField()->Splitter();
  // Variables to save different L2 - Norms

  double velincnorm_L2(0.0);
  double velnorm_L2(0.0);
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);

  // Calculate velocity increment and velocity L2 - Norm
  // velincnp_ = 1.0 * convelnp_ - 1.0 * conveln_

  velincnp_->Update(1.0, *fluid_field()->ExtractVelocityPart(fluid_field()->Velnp()), -1.0);
  velincnp_->Norm2(&velincnorm_L2);  // Estimation of the L2 - norm save values to both variables
                                     // (velincnorm_L2 and velnorm_L2)
  fluid_field()->ExtractVelocityPart(fluid_field()->Velnp())->Norm2(&velnorm_L2);

  // Calculate phi increment and phi L2 - Norm
  // tempincnp_ includes the concentration and the potential increment
  // tempincnp_ = 1.0 * phinp_ - 1.0 * phin_

  phiincnp_->Update(1.0, *ScaTraField()->Phinp(), -1.0);
  phiincnp_->Norm2(&phiincnorm_L2);
  ScaTraField()->Phinp()->Norm2(&phinorm_L2);

  // care for the case that there is (almost) zero temperature or velocity
  // (usually not required for temperature)
  if (velnorm_L2 < 1e-6) velnorm_L2 = 1.0;
  if (phinorm_L2 < 1e-6) phinorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (natconvitnum != 1)
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << "\n";
      std::cout
          << "*****************************************************************************\n";
      std::cout << "                          OUTER ITERATION STEP\n";
      std::cout
          << "*****************************************************************************\n";
      printf("+------------+-------------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E  | %10.3E  |", natconvitnum, natconvitmax,
          natconvittol, phiincnorm_L2 / phinorm_L2, velincnorm_L2 / velnorm_L2);
      printf("\n");
      printf("+------------+-------------------+-------------+-------------+\n");
    }

    // Converged or Not
    if ((phiincnorm_L2 / phinorm_L2 <= natconvittol) &&
        (velincnorm_L2 / velnorm_L2 <= natconvittol))
    // if ((incconnorm_L2/connorm_L2 <= natconvittol))
    {
      stopnonliniter = true;
      if (Comm().MyPID() == 0)
      {
        printf("| Outer Iteration loop converged after iteration %3d/%3d                    |\n",
            natconvitnum, natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
      }
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        printf("| Outer Iteration loop is not converged after iteration %3d/%3d             |\n",
            natconvitnum, natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
      }
    }
  }
  else
  {
    // first outer iteration loop: fluid solver has not got the new density yet
    // => minimum two outer iteration loops
    stopnonliniter = false;
    if (Comm().MyPID() == 0)
    {
      std::cout << "\n";
      std::cout
          << "*****************************************************************************\n";
      std::cout << "                          OUTER ITERATION STEP\n";
      std::cout
          << "*****************************************************************************\n";
      printf("+------------+-------------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |       -     |      -      |", natconvitnum,
          natconvitmax, natconvittol);
      printf("\n");
      printf("+------------+-------------------+-------------+-------------+\n");
    }
  }

  // warn if natconvitemax is reached without convergence, but proceed to next timestep
  // natconvitemax = 1 is also possible for segregated coupling approaches (not fully implicit)
  if (natconvitnum == natconvitmax)
  {
    if (((phiincnorm_L2 / phinorm_L2 > natconvittol) ||
            (velincnorm_L2 / velnorm_L2 > natconvittol)) or
        (natconvitmax == 1))
    {
      stopnonliniter = true;
      if ((Comm().MyPID() == 0))
      {
        printf("|     >>>>>> not converged in itemax steps!     |\n");
        printf("+-----------------------------------------------+\n");
        printf("\n");
        printf("\n");
      }
    }
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::ReadInflowRestart(int restart)
{
  // in case a inflow generation in the inflow section has been performed,
  // there are not any scatra results available and the initial field is used
  fluid_field()->read_restart(restart);
  // as read_restart is only called for the fluid_field
  // time and step have not been set in the superior class and the ScaTraField
  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
  ScaTraField()->SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraAlgorithm::TestResults()
{
  GLOBAL::Problem::Instance()->AddFieldTest(fluid_field()->CreateFieldTest());
  GLOBAL::Problem::Instance()->AddFieldTest(create_sca_tra_field_test());
  GLOBAL::Problem::Instance()->TestAll(Comm());
}
FOUR_C_NAMESPACE_CLOSE
