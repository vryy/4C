/*----------------------------------------------------------------------*/
/*! \file

\brief Main control routine for electromagnetic simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_elemag_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_comm_utils.hpp"
#include "4C_discretization_dofset_independent.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_elemag_timeint.hpp"
#include "4C_elemag_utils_clonestrategy.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret_hdg.hpp"
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
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // The function num_dof_per_element_auxiliary() of the electromagnetic elements return nsd_*2.
  // This does not assure that the code will work in any case (more spatial dimensions might give
  // problems)
  if (problem->NDim() != 3)
  {
    FOUR_C_THROW(
        "The implementation of electromagnetic propagation only supports 3D problems.\n"
        "It is necessary to change the spatial dimension of your problem.");
  }

  // declare problem-specific parameter list for electromagnetics
  const Teuchos::ParameterList& elemagparams = problem->electromagnetic_params();

  // declare discretization and check their existence
  Teuchos::RCP<DRT::DiscretizationHDG> elemagdishdg =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationHDG>(problem->GetDis("elemag"));
  if (elemagdishdg == Teuchos::null)
    FOUR_C_THROW("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  elemagdishdg->PrintFaces(std::cout);
#endif

  // declare communicator and print module information to screen
  const Epetra_Comm& comm = elemagdishdg->Comm();
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
  if (not elemagdishdg->Filled() || not elemagdishdg->HaveDofs()) elemagdishdg->fill_complete();
  // Asking the discretization how many internal DOF the elements have and creating the additional
  // DofSet
  int eledofs = dynamic_cast<DRT::ELEMENTS::Elemag*>(elemagdishdg->lColElement(0))
                    ->num_dof_per_element_auxiliary();
  Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(0, eledofs, 0, false));
  elemagdishdg->AddDofSet(dofsetaux);

  // call fill complete on discretization
  elemagdishdg->fill_complete();

  // create solver
  const int linsolvernumber_elemag = elemagparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber_elemag == (-1))
    FOUR_C_THROW(
        "There is not any linear solver defined for electromagnetic problem. Please set "
        "LINEAR_SOLVER in ELECTROMAGNETIC DYNAMIC to a valid number!");

  Teuchos::RCP<CORE::LINALG::Solver> solver =
      Teuchos::rcp(new CORE::LINALG::Solver(problem->SolverParams(linsolvernumber_elemag), comm));

  // declare output writer
  Teuchos::RCP<IO::DiscretizationWriter> output = elemagdishdg->Writer();

  // declare electromagnetic parameter list
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList(elemagparams));

  // set restart step if required
  int restart = problem->Restart();
  params->set<int>("restart", restart);

  // create algorithm depending on time-integration scheme
  INPAR::ELEMAG::DynamicType elemagdyna =
      CORE::UTILS::IntegralValue<INPAR::ELEMAG::DynamicType>(elemagparams, "TIMEINT");
  Teuchos::RCP<ELEMAG::ElemagTimeInt> elemagalgo;
  switch (elemagdyna)
  {
    case INPAR::ELEMAG::elemag_ost:
    {
      FOUR_C_THROW("One step theta not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimIntOST(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_bdf1:
    case INPAR::ELEMAG::elemag_bdf2:
    case INPAR::ELEMAG::elemag_bdf4:
    {
      elemagalgo = Teuchos::rcp(new ELEMAG::ElemagTimeInt(elemagdishdg, solver, params, output));
      break;
    }
    case INPAR::ELEMAG::elemag_genAlpha:
    {
      FOUR_C_THROW("Generalized-alpha method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::ElemagGenAlpha(elemagdishdg, solver, params,
      // output));
      break;
    }
    case INPAR::ELEMAG::elemag_explicit_euler:
    {
      FOUR_C_THROW("Explicit euler method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntExplEuler(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_rk:
    {
      FOUR_C_THROW("Runge-Kutta methods not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntRK(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_cn:
    {
      FOUR_C_THROW("Crank-Nicolson method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntCN(elemagdishdg,solver,params,output));
      break;
    }
    default:
      FOUR_C_THROW("Unknown time-integration scheme for problem type electromagnetics");
      break;
  }

  // Initialize the evolution algorithm
  elemagalgo->Init();

  // set initial field
  if (restart)
    elemagalgo->read_restart(restart);
  else
  {
    INPAR::ELEMAG::InitialField init =
        CORE::UTILS::IntegralValue<INPAR::ELEMAG::InitialField>(elemagparams, "INITIALFIELD");

    bool ishdg = false;

    switch (init)
    {
      case INPAR::ELEMAG::initfield_scatra_hdg:
        ishdg = true;
        [[fallthrough]];
      case INPAR::ELEMAG::initfield_scatra:
      {
        Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(elemagdishdg->Comm().Clone());

        Teuchos::RCP<DRT::Discretization> scatradis;

        if (ishdg)
        {
          scatradis = Teuchos::rcp(new DRT::DiscretizationHDG((std::string) "scatra", newcomm));

          scatradis->fill_complete();

          CORE::FE::CloneDiscretization<
              ELEMAG::UTILS::ScatraCloneStrategy<CORE::FE::ShapeFunctionType::hdg>>(
              elemagdishdg, scatradis, GLOBAL::Problem::Instance()->CloningMaterialMap());
        }
        else
        {
          scatradis = Teuchos::rcp(new DRT::Discretization((std::string) "scatra", newcomm));
          scatradis->fill_complete();

          CORE::FE::CloneDiscretization<
              ELEMAG::UTILS::ScatraCloneStrategy<CORE::FE::ShapeFunctionType::polynomial>>(
              elemagdishdg, scatradis, GLOBAL::Problem::Instance()->CloningMaterialMap());
        }

        // call fill complete on discretization
        scatradis->fill_complete();

        Teuchos::RCP<IO::DiscretizationWriter> output_scatra = scatradis->Writer();

        // This is necessary to have the dirichlet conditions done also in the scatra problmem. It
        // might be necessary to rethink how things are handled inside the
        // DRT::UTILS::DbcHDG::do_dirichlet_condition.
        problem->SetProblemType(GLOBAL::ProblemType::scatra);

        // access the problem-specific parameter list
        const Teuchos::ParameterList& scatradyn =
            GLOBAL::Problem::Instance()->scalar_transport_dynamic_params();

        // do the scatra
        const INPAR::SCATRA::VelocityField veltype =
            CORE::UTILS::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");
        switch (veltype)
        {
          case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
          {
            // we directly use the elements from the scalar transport elements section
            if (scatradis->NumGlobalNodes() == 0)
              FOUR_C_THROW("No elements in the ---TRANSPORT ELEMENTS section");

            // add proxy of velocity related degrees of freedom to scatra discretization
            Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux =
                Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(
                    GLOBAL::Problem::Instance()->NDim() + 1, 0, 0, true));
            if (scatradis->AddDofSet(dofsetaux) != 1)
              FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");

            // finalize discretization
            scatradis->fill_complete(true, true, true);

            // create scatra output
            // access the problem-specific parameter list
            Teuchos::RCP<Teuchos::ParameterList> scatraparams =
                Teuchos::rcp(new Teuchos::ParameterList(
                    GLOBAL::Problem::Instance()->scalar_transport_dynamic_params()));

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
            const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();
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
            Teuchos::RCP<CORE::LINALG::Solver> scatrasolver = Teuchos::rcp(new CORE::LINALG::Solver(
                GLOBAL::Problem::Instance()->SolverParams(scatraparams->get<int>("LINEAR_SOLVER")),
                scatradis->Comm()));

            // create instance of scalar transport basis algorithm (empty fluid discretization)
            Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatraalgo;
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

              scatraalgo = Teuchos::rcp(new SCATRA::TimIntStationaryHDG(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }
            else
            {
              // Add parameters for CG
              // There is no need for stabilization as the problem is a pure diffusion problem
              scatraparams->sublist("STABILIZATION").set("STABTYPE", "no_stabilization");
              scatraalgo = Teuchos::rcp(new SCATRA::TimIntStationary(
                  scatradis, scatrasolver, scatraparams, scatraextraparams, output));
            }

            // scatraparams->print(std::cout);

            scatraalgo->Init();
            scatraalgo->set_number_of_dof_set_velocity(1);
            scatraalgo->Setup();
            scatraalgo->set_velocity_field();
            scatraalgo->TimeLoop();

            // scatraalgo->compute_interior_values();

            // Create a vector that is going to be filled differently depending on the
            // discretization. If HDG we already have the information about the gradient, otherwise
            // the gradient has to be computed.
            Teuchos::RCP<Epetra_Vector> phi;

            // If it is an HDG discretization return the interior variables else return the nodal
            // values
            if (ishdg)
              phi = Teuchos::rcp_dynamic_cast<SCATRA::TimIntStationaryHDG>(scatraalgo)
                        ->ReturnIntPhinp();
            else
              phi = scatraalgo->Phinp();

            // This is a shortcut for output reason
            // TODO (berardocco) Fix the output
            output->create_new_result_and_mesh_file();

            // Given the results of the scatra solver obtain the initial value of the electric field
            elemagalgo->set_initial_electric_field(phi, scatradis);

            // Once work is done change back to problem elemag
            problem->SetProblemType(GLOBAL::ProblemType::elemag);

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
        elemagalgo->SetInitialField(init, startfuncno);
        break;
      }
    }
  }

  // print information to screen
  elemagalgo->print_information_to_screen();

  elemagalgo->Integrate();

  // Computing the error at the las time step (the conditional stateme nt is inside for now)
  if (CORE::UTILS::IntegralValue<bool>(elemagparams, "CALCERR"))
  {
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> errors = elemagalgo->compute_error();
    elemagalgo->PrintErrors(errors);
  }

  // print computing time
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = CORE::COMM::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // do result test if required
  problem->AddFieldTest(elemagalgo->CreateFieldTest());
  problem->TestAll(comm);

  return;
}

FOUR_C_NAMESPACE_CLOSE
