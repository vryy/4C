/*----------------------------------------------------------------------*/
/*! \file

\brief Main control routine for electromagnetic simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_elemag_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_comm_utils.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_elemag_timeint.hpp"
#include "4C_elemag_utils_clonestrategy.hpp"
#include "4C_fem_discretization_hdg.hpp"
#include "4C_fem_dofset_independent.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_timint_stat.hpp"
#include "4C_scatra_timint_stat_hdg.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

void electromagnetics_drt()
{
  // declare abbreviation
  Global::Problem* problem = Global::Problem::instance();

  // The function num_dof_per_element_auxiliary() of the electromagnetic elements return nsd_*2.
  // This does not assure that the code will work in any case (more spatial dimensions might give
  // problems)
  if (problem->n_dim() != 3)
  {
    FOUR_C_THROW(
        "The implementation of electromagnetic propagation only supports 3D problems.\n"
        "It is necessary to change the spatial dimension of your problem.");
  }

  // declare problem-specific parameter list for electromagnetics
  const Teuchos::ParameterList& elemagparams = problem->electromagnetic_params();

  // declare discretization and check their existence
  Teuchos::RCP<Core::FE::DiscretizationHDG> elemagdishdg =
      Teuchos::rcp_dynamic_cast<Core::FE::DiscretizationHDG>(problem->get_dis("elemag"));
  if (elemagdishdg == Teuchos::null)
    FOUR_C_THROW(
        "Failed to cast Core::FE::Discretization to "
        "Core::FE::DiscretizationHDG.");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  elemagdishdg->print_faces(std::cout);
#endif

  // declare communicator and print module information to screen
  const Epetra_Comm& comm = elemagdishdg->get_comm();
  if (comm.MyPID() == 0)
  {
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "---------- You are now about to enter the module for electromagnetics! ----------"
              << std::endl;
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
  }

  // call fill complete on discretization
  if (not elemagdishdg->filled() || not elemagdishdg->have_dofs()) elemagdishdg->fill_complete();
  // Asking the discretization how many internal DOF the elements have and creating the additional
  // DofSet
  int eledofs = dynamic_cast<Discret::ELEMENTS::Elemag*>(elemagdishdg->l_col_element(0))
                    ->num_dof_per_element_auxiliary();
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(0, eledofs, 0, false));
  elemagdishdg->add_dof_set(dofsetaux);

  // call fill complete on discretization
  elemagdishdg->fill_complete();

  // create solver
  const int linsolvernumber_elemag = elemagparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber_elemag == (-1))
    FOUR_C_THROW(
        "There is not any linear solver defined for electromagnetic problem. Please set "
        "LINEAR_SOLVER in ELECTROMAGNETIC DYNAMIC to a valid number!");

  Teuchos::RCP<Core::LinAlg::Solver> solver =
      Teuchos::rcp(new Core::LinAlg::Solver(problem->solver_params(linsolvernumber_elemag), comm,
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY")));

  // declare output writer
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = elemagdishdg->writer();

  // declare electromagnetic parameter list
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList(elemagparams));

  // set restart step if required
  int restart = problem->restart();
  params->set<int>("restart", restart);

  // create algorithm depending on time-integration scheme
  Inpar::EleMag::DynamicType elemagdyna =
      Core::UTILS::IntegralValue<Inpar::EleMag::DynamicType>(elemagparams, "TIMEINT");
  Teuchos::RCP<EleMag::ElemagTimeInt> elemagalgo;
  switch (elemagdyna)
  {
    case Inpar::EleMag::elemag_ost:
    {
      FOUR_C_THROW("One step theta not yet implemented.");
      // elemagalgo = Teuchos::rcp(new EleMag::TimIntOST(elemagdishdg,solver,params,output));
      break;
    }
    case Inpar::EleMag::elemag_bdf1:
    case Inpar::EleMag::elemag_bdf2:
    case Inpar::EleMag::elemag_bdf4:
    {
      elemagalgo = Teuchos::rcp(new EleMag::ElemagTimeInt(elemagdishdg, solver, params, output));
      break;
    }
    case Inpar::EleMag::elemag_genAlpha:
    {
      FOUR_C_THROW("Generalized-alpha method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new EleMag::ElemagGenAlpha(elemagdishdg, solver, params,
      // output));
      break;
    }
    case Inpar::EleMag::elemag_explicit_euler:
    {
      FOUR_C_THROW("Explicit euler method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new EleMag::TimeIntExplEuler(elemagdishdg,solver,params,output));
      break;
    }
    case Inpar::EleMag::elemag_rk:
    {
      FOUR_C_THROW("Runge-Kutta methods not yet implemented.");
      // elemagalgo = Teuchos::rcp(new EleMag::TimeIntRK(elemagdishdg,solver,params,output));
      break;
    }
    case Inpar::EleMag::elemag_cn:
    {
      FOUR_C_THROW("Crank-Nicolson method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new EleMag::TimeIntCN(elemagdishdg,solver,params,output));
      break;
    }
    default:
      FOUR_C_THROW("Unknown time-integration scheme for problem type electromagnetics");
      break;
  }

  // Initialize the evolution algorithm
  elemagalgo->init();

  // set initial field
  if (restart)
    elemagalgo->read_restart(restart);
  else
  {
    Inpar::EleMag::InitialField init =
        Core::UTILS::IntegralValue<Inpar::EleMag::InitialField>(elemagparams, "INITIALFIELD");

    bool ishdg = false;

    switch (init)
    {
      case Inpar::EleMag::initfield_scatra_hdg:
        ishdg = true;
        [[fallthrough]];
      case Inpar::EleMag::initfield_scatra:
      {
        Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(elemagdishdg->get_comm().Clone());

        Teuchos::RCP<Core::FE::Discretization> scatradis;

        if (ishdg)
        {
          scatradis = Teuchos::rcp(
              new Core::FE::DiscretizationHDG((std::string) "scatra", newcomm, problem->n_dim()));

          scatradis->fill_complete();

          Core::FE::CloneDiscretization<
              EleMag::UTILS::ScatraCloneStrategy<Core::FE::ShapeFunctionType::hdg>>(
              elemagdishdg, scatradis, Global::Problem::instance()->cloning_material_map());
        }
        else
        {
          scatradis = Teuchos::rcp(
              new Core::FE::Discretization((std::string) "scatra", newcomm, problem->n_dim()));
          scatradis->fill_complete();

          Core::FE::CloneDiscretization<
              EleMag::UTILS::ScatraCloneStrategy<Core::FE::ShapeFunctionType::polynomial>>(
              elemagdishdg, scatradis, Global::Problem::instance()->cloning_material_map());
        }

        // call fill complete on discretization
        scatradis->fill_complete();

        Teuchos::RCP<Core::IO::DiscretizationWriter> output_scatra = scatradis->writer();

        // This is necessary to have the dirichlet conditions done also in the scatra problmem. It
        // might be necessary to rethink how things are handled inside the
        // Core::FE::UTILS::DbcHDG::do_dirichlet_condition.
        problem->set_problem_type(Core::ProblemType::scatra);

        // access the problem-specific parameter list
        const Teuchos::ParameterList& scatradyn =
            Global::Problem::instance()->scalar_transport_dynamic_params();

        // do the scatra
        const Inpar::ScaTra::VelocityField veltype =
            Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(scatradyn, "VELOCITYFIELD");
        switch (veltype)
        {
          case Inpar::ScaTra::velocity_zero:  // zero  (see case 1)
          {
            // we directly use the elements from the scalar transport elements section
            if (scatradis->num_global_nodes() == 0)
              FOUR_C_THROW("No elements in the ---TRANSPORT ELEMENTS section");

            // add proxy of velocity related degrees of freedom to scatra discretization
            Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
                Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
                    Global::Problem::instance()->n_dim() + 1, 0, 0, true));
            if (scatradis->add_dof_set(dofsetaux) != 1)
              FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");

            // finalize discretization
            scatradis->fill_complete(true, true, true);

            // create scatra output
            // access the problem-specific parameter list
            Teuchos::RCP<Teuchos::ParameterList> scatraparams =
                Teuchos::rcp(new Teuchos::ParameterList(
                    Global::Problem::instance()->scalar_transport_dynamic_params()));

            // TODO (berardocco) Might want to add the scatra section in the input file to avoid
            // adding params to the elemag or using existing ones for scatra purposes
            scatraparams->set("TIMEINTEGR", "Stationary");
            scatraparams->set("NUMSTEP", 1);
            // This way we avoid writing results and restart
            scatraparams->set("RESULTSEVRY", 1000);
            scatraparams->set("RESTARTEVRY", 1000);
            // This has to be changed accordingly to the intial time
            // As of now elemag simulation can only start at 0.

            // The solver type still has to be fixed as the problem is linear but the steady state
            // does not always behave correctly with linear solvers.
            scatraparams->set("SOLVERTYPE", "nonlinear");

            // create necessary extra parameter list for scatra
            Teuchos::RCP<Teuchos::ParameterList> scatraextraparams;
            scatraextraparams = Teuchos::rcp(new Teuchos::ParameterList());
            scatraextraparams->set<bool>("isale", false);
            const Teuchos::ParameterList& fdyn =
                Global::Problem::instance()->fluid_dynamic_params();
            scatraextraparams->sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");
            scatraextraparams->sublist("SUBGRID VISCOSITY") = fdyn.sublist("SUBGRID VISCOSITY");
            scatraextraparams->sublist("MULTIFRACTAL SUBGRID SCALES") =
                fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
            scatraextraparams->sublist("TURBULENT INFLOW") = fdyn.sublist("TURBULENT INFLOW");

            scatraextraparams->set("ELECTROMAGNETICDIFFUSION", true);
            scatraextraparams->set("EMDSOURCE", elemagparams.get<int>("SOURCEFUNCNO"));

            // In case the scatra solver is not defined just use the elemag one
            if (scatraparams->get<int>("LINEAR_SOLVER") == -1)
              scatraparams->set<int>("LINEAR_SOLVER", elemagparams.get<int>("LINEAR_SOLVER"));

            // create solver
            Teuchos::RCP<Core::LinAlg::Solver> scatrasolver = Teuchos::rcp(new Core::LinAlg::Solver(
                Global::Problem::instance()->solver_params(scatraparams->get<int>("LINEAR_SOLVER")),
                scatradis->get_comm(), Global::Problem::instance()->solver_params_callback(),
                Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                    Global::Problem::instance()->io_params(), "VERBOSITY")));

            // create instance of scalar transport basis algorithm (empty fluid discretization)
            Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatraalgo;
            if (ishdg)
            {
              // Add parameters for HDG
              scatraparams->sublist("STABILIZATION").set("STABTYPE", "centered");
              scatraparams->sublist("STABILIZATION").set("DEFINITION_TAU", "Numerical_Value");
              // If the input file does not specify a tau parameter then use the one given to the
              // elemag discretization
              if (scatraparams->sublist("STABILIZATION").get<double>("TAU_VALUE") == 0.0)
                scatraparams->sublist("STABILIZATION")
                    .set("TAU_VALUE", elemagparams.get<double>("TAU"));

              scatraalgo = Teuchos::rcp(new ScaTra::TimIntStationaryHDG(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }
            else
            {
              // Add parameters for CG
              // There is no need for stabilization as the problem is a pure diffusion problem
              scatraparams->sublist("STABILIZATION").set("STABTYPE", "no_stabilization");
              scatraalgo = Teuchos::rcp(new ScaTra::TimIntStationary(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }

            // scatraparams->print(std::cout);

            scatraalgo->init();
            scatraalgo->set_number_of_dof_set_velocity(1);
            scatraalgo->setup();
            scatraalgo->set_velocity_field();
            scatraalgo->time_loop();

            // scatraalgo->compute_interior_values();

            // Create a vector that is going to be filled differently depending on the
            // discretization. If HDG we already have the information about the gradient, otherwise
            // the gradient has to be computed.
            Teuchos::RCP<Epetra_Vector> phi;

            // If it is an HDG discretization return the interior variables else return the nodal
            // values
            if (ishdg)
              phi = Teuchos::rcp_dynamic_cast<ScaTra::TimIntStationaryHDG>(scatraalgo)
                        ->return_int_phinp();
            else
              phi = scatraalgo->phinp();

            // This is a shortcut for output reason
            // TODO (berardocco) Fix the output
            output->create_new_result_and_mesh_file();

            // Given the results of the scatra solver obtain the initial value of the electric field
            elemagalgo->set_initial_electric_field(phi, scatradis);

            // Once work is done change back to problem elemag
            problem->set_problem_type(Core::ProblemType::elemag);

            break;
          }
          default:
            FOUR_C_THROW(
                "Does not make sense to have a velocity field to initialize the electric potential "
                "field.\nCheck your input file.");
            break;
        }
        break;
      }
      default:
      {
        int startfuncno = elemagparams.get<int>("STARTFUNCNO");
        elemagalgo->set_initial_field(init, startfuncno);
        break;
      }
    }
  }

  // print information to screen
  elemagalgo->print_information_to_screen();

  elemagalgo->integrate();

  // Computing the error at the las time step (the conditional stateme nt is inside for now)
  if (Core::UTILS::IntegralValue<bool>(elemagparams, "CALCERR"))
  {
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> errors = elemagalgo->compute_error();
    elemagalgo->print_errors(errors);
  }

  // print computing time
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // do result test if required
  problem->add_field_test(elemagalgo->create_field_test());
  problem->test_all(comm);

  return;
}

FOUR_C_NAMESPACE_CLOSE
