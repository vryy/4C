/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for Electrochemistry module.

\level 2


*----------------------------------------------------------------------*/
#include "4C_elch_dyn.hpp"

#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_elch_algorithm.hpp"
#include "4C_elch_moving_boundary_algorithm.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void elch_dyn(int restart)
{
  // pointer to problem
  auto* problem = GLOBAL::Problem::Instance();

  // access the communicator
  const auto& comm = problem->GetDis("fluid")->Comm();

  // print ELCH-Logo to screen
  if (comm.MyPID() == 0) printlogo();

  // access the fluid discretization
  auto fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  auto scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra/elch dof
  fluiddis->fill_complete();
  scatradis->fill_complete();

  // access the problem-specific parameter list
  const auto& elchcontrol = problem->ELCHControlParams();

  // print default parameters to screen
  if (comm.MyPID() == 0) INPUT::PrintDefaultParameters(IO::cout, elchcontrol);

  // access the scalar transport parameter list
  const auto& scatradyn = problem->scalar_transport_dynamic_params();
  const auto veltype =
      CORE::UTILS::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:      // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // spatial function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes() == 0)
        FOUR_C_THROW("No elements in the ---TRANSPORT ELEMENTS section");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
      {
        FOUR_C_THROW(
            "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
            "TRANSPORT DYNAMIC to a valid number!");
      }

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      auto scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(
          scatradyn, scatradyn, GLOBAL::Problem::Instance()->SolverParams(linsolvernumber)));

      // add proxy of velocity related degrees of freedom to scatra discretization
      auto dofsetaux = Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(
          GLOBAL::Problem::Instance()->NDim() + 1, 0, 0, true));
      if (scatradis->AddDofSet(dofsetaux) != 1)
        FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      scatraonly->ScaTraField()->set_number_of_dof_set_velocity(1);

      // now me may redistribute or ghost the scatra discretization
      // finalize discretization
      scatradis->fill_complete(true, true, true);

      // now we can call Init() on the base algorithm
      // scatra time integrator is constructed and initialized inside
      scatraonly->Init();

      // now me may redistribute or ghost the scatra discretization
      // finalize discretization
      scatradis->fill_complete(true, true, true);

      // only now we must call Setup() on the base algorithm.
      // all objects relying on the parallel distribution are
      // created and pointers are set.
      // calls Setup() on time integrator inside.
      scatraonly->Setup();

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->read_restart(restart);

      // set velocity field
      // note: The order read_restart() before SetVelocityField() is important here!!
      // for time-dependent velocity fields, SetVelocityField() is additionally called in each
      // prepare_time_step()-call
      scatraonly->ScaTraField()->SetVelocityField();

      // enter time loop to solve problem with given convective velocity
      scatraonly->ScaTraField()->TimeLoop();

      // perform the result test if required
      scatraonly->ScaTraField()->TestResults();

      break;
    }
    case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
    {
      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes() == 0)
      {
        // fill scatra discretization by cloning fluid discretization
        CORE::FE::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(
            fluiddis, scatradis, GLOBAL::Problem::Instance()->CloningMaterialMap());
        scatradis->fill_complete();
        // determine implementation type of cloned scatra elements
        INPAR::SCATRA::ImplType impltype = INPAR::SCATRA::impltype_undefined;
        if (CORE::UTILS::IntegralValue<int>(elchcontrol, "DIFFCOND_FORMULATION"))
          impltype = INPAR::SCATRA::impltype_elch_diffcond;
        else
          impltype = INPAR::SCATRA::impltype_elch_NP;

        // set implementation type
        for (int i = 0; i < scatradis->NumMyColElements(); ++i)
        {
          auto* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if (element == nullptr)
            FOUR_C_THROW("Invalid element type!");
          else
            element->SetImplType(impltype);
        }
      }

      else
        FOUR_C_THROW("Fluid AND ScaTra discretization present. This is not supported.");

      // support for turbulent flow statistics
      const auto& fdyn = (problem->FluidDynamicParams());

      Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
      if (!aledis->Filled()) aledis->fill_complete(false, false, false);
      // is ALE needed or not?
      const auto withale = CORE::UTILS::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(
          elchcontrol, "MOVINGBOUNDARY");

      if (withale != INPAR::ELCH::elch_mov_bndry_no)
      {
        // create ale elements only if the ale discretization is empty
        if (aledis->NumGlobalNodes() == 0)
        {
          // clone ALE discretization from fluid discretization
          CORE::FE::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
              fluiddis, aledis, GLOBAL::Problem::Instance()->CloningMaterialMap());

          aledis->fill_complete(true, true, false);
          // setup material in every ALE element
          Teuchos::ParameterList params;
          params.set<std::string>("action", "setup_material");
          aledis->Evaluate(params);
        }
        else
          FOUR_C_THROW("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

        // get linear solver id from SCALAR TRANSPORT DYNAMIC
        const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
        if (linsolvernumber == -1)
        {
          FOUR_C_THROW(
              "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
              "TRANSPORT DYNAMIC to a valid number!");
        }

        // create an ELCH::MovingBoundaryAlgorithm instance
        // NOTE: elch reads time parameters from scatra dynamic section!
        auto elch = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(
            comm, elchcontrol, scatradyn, problem->SolverParams(linsolvernumber)));

        // add proxy of fluid degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
        elch->ScaTraField()->set_number_of_dof_set_velocity(1);

        // add proxy of ALE degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(aledis->GetDofSetProxy()) != 2)
          FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
        elch->ScaTraField()->set_number_of_dof_set_displacement(2);

        // now we must call Init()
        elch->Init();

        // NOTE : At this point we may redistribute and/or
        //        ghost our discretizations at will.
        scatradis->fill_complete();
        fluiddis->fill_complete();
        aledis->fill_complete();

        // now we can call Setup() on the scatra time integrator
        elch->Setup();

        // read the restart information, set vectors and variables
        if (restart) elch->read_restart(restart);

        // solve the whole electrochemistry problem
        elch->TimeLoop();

        // summarize the performance measurements
        Teuchos::TimeMonitor::summarize();

        // perform the result test
        elch->TestResults();
      }
      else
      {
        // get linear solver id from SCALAR TRANSPORT DYNAMIC
        const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
        if (linsolvernumber == -1)
        {
          FOUR_C_THROW(
              "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
              "TRANSPORT DYNAMIC to a valid number!");
        }

        // create an ELCH::Algorithm instance
        // NOTE: elch reads time parameters from scatra dynamic section!
        auto elch = Teuchos::rcp(new ELCH::Algorithm(
            comm, elchcontrol, scatradyn, fdyn, problem->SolverParams(linsolvernumber)));

        // add proxy of fluid degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
        elch->ScaTraField()->set_number_of_dof_set_velocity(1);

        // now we must call Init()
        elch->Init();

        // NOTE : At this point we may redistribute and/or
        //        ghost our discretizations at will.
        scatradis->fill_complete();
        fluiddis->fill_complete();
        aledis->fill_complete();

        // discretizations are done, now we can call Setup() on the algorithm
        elch->Setup();

        // read the restart information, set vectors and variables
        if (restart) elch->read_restart(restart);

        // solve the whole electrochemistry problem
        elch->TimeLoop();

        // summarize the performance measurements
        Teuchos::TimeMonitor::summarize();

        // perform the result test
        elch->TestResults();
      }

      break;
    }
    default:
      FOUR_C_THROW("Unknown velocity field type for transport of passive scalar: %d", veltype);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void printlogo()
{
  // more at http://www.ascii-art.de under entry "moose" (or "elk")
  std::cout << "     ___            ___    " << '\n';
  std::cout << "    /   \\          /   \\ " << '\n';
  std::cout << "    \\_   \\        /  __/ " << '\n';
  std::cout << "     _\\   \\      /  /__  "
            << "     _____ _     _____  _   _   " << '\n';
  std::cout << "     \\___  \\____/   __/  "
            << "    |  ___| |   /  __ \\| | | |  " << '\n';
  std::cout << "         \\_       _/      "
            << "   | |__ | |   | /  \\/| |_| |  " << '\n';
  std::cout << "           | @ @  \\_      "
            << "   |  __|| |   | |    |  _  |   " << '\n';
  std::cout << "           |               "
            << "  | |___| |___| \\__/\\| | | | " << '\n';
  std::cout << "         _/     /\\        "
            << "   \\____/\\_____/\\____/\\_| |_/ " << '\n';
  std::cout << "        /o)  (o/\\ \\_     " << '\n';
  std::cout << "        \\_____/ /         " << '\n';
  std::cout << "          \\____/          " << '\n';
  std::cout << "                           " << '\n';
}

FOUR_C_NAMESPACE_CLOSE
