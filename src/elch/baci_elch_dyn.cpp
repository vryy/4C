/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for Electrochemistry module.

\level 2


*----------------------------------------------------------------------*/
#include "baci_elch_dyn.H"

#include "baci_ale_utils_clonestrategy.H"
#include "baci_elch_algorithm.H"
#include "baci_elch_moving_boundary_algorithm.H"
#include "baci_inpar_elch.H"
#include "baci_inpar_validparameters.H"
#include "baci_lib_dofset_predefineddofnumber.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils_createdis.H"
#include "baci_scatra_ele.H"
#include "baci_scatra_resulttest_elch.H"
#include "baci_scatra_timint_elch.H"
#include "baci_scatra_utils_clonestrategy.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn(int restart)
{
  // pointer to problem
  auto* problem = DRT::Problem::Instance();

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
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // access the problem-specific parameter list
  const auto& elchcontrol = problem->ELCHControlParams();

  // print default parameters to screen
  if (comm.MyPID() == 0) DRT::INPUT::PrintDefaultParameters(IO::cout, elchcontrol);

  // access the scalar transport parameter list
  const auto& scatradyn = problem->ScalarTransportDynamicParams();
  const auto veltype =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:      // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // spatial function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes() == 0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // add proxy of velocity related degrees of freedom to scatra discretization
      auto dofsetaux = Teuchos::rcp(
          new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim() + 1, 0, 0, true));
      if (scatradis->AddDofSet(dofsetaux) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");


      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
      {
        dserror(
            "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
            "TRANSPORT DYNAMIC to a valid number!");
      }

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      auto scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

      // now we can call Init() on the base algorithm
      // scatra time integrator is constructed and initialized inside
      scatraonly->Init(
          scatradyn, scatradyn, DRT::Problem::Instance()->SolverParams(linsolvernumber));
      scatraonly->ScaTraField()->SetNumberOfDofSetVelocity(1);

      // now me may redistribute or ghost the scatra discretization
      // finalize discretization
      scatradis->FillComplete(true, true, true);

      // only now we must call Setup() on the base algorithm.
      // all objects relying on the parallel distribution are
      // created and pointers are set.
      // calls Setup() on time integrator inside.
      scatraonly->Setup();

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->ReadRestart(restart);

      // set velocity field
      // note: The order ReadRestart() before SetVelocityField() is important here!!
      // for time-dependent velocity fields, SetVelocityField() is additionally called in each
      // PrepareTimeStep()-call
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
      if (fluiddis->NumGlobalNodes() == 0) dserror("Fluid discretization is empty!");

      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes() == 0)
      {
        // fill scatra discretization by cloning fluid discretization
        DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis, scatradis);
        scatradis->FillComplete();
        // determine implementation type of cloned scatra elements
        INPAR::SCATRA::ImplType impltype = INPAR::SCATRA::impltype_undefined;
        if (DRT::INPUT::IntegralValue<int>(elchcontrol, "DIFFCOND_FORMULATION"))
          impltype = INPAR::SCATRA::impltype_elch_diffcond;
        else
          impltype = INPAR::SCATRA::impltype_elch_NP;

        // set implementation type
        for (int i = 0; i < scatradis->NumMyColElements(); ++i)
        {
          auto* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if (element == nullptr)
            dserror("Invalid element type!");
          else
            element->SetImplType(impltype);
        }
      }

      else
        dserror("Fluid AND ScaTra discretization present. This is not supported.");

      // support for turbulent flow statistics
      const auto& fdyn = (problem->FluidDynamicParams());

      Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
      if (!aledis->Filled()) aledis->FillComplete(false, false, false);
      // is ALE needed or not?
      const auto withale =
          DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol, "MOVINGBOUNDARY");

      if (withale != INPAR::ELCH::elch_mov_bndry_no)
      {
        // create ale elements only if the ale discretization is empty
        if (aledis->NumGlobalNodes() == 0)
        {
          // clone ALE discretization from fluid discretization
          DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);

          aledis->FillComplete(true, true, false);
          // setup material in every ALE element
          Teuchos::ParameterList params;
          params.set<std::string>("action", "setup_material");
          aledis->Evaluate(params);
        }
        else
          dserror("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

        // add proxy of fluid degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          dserror("Scatra discretization has illegal number of dofsets!");

        // add proxy of ALE degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(aledis->GetDofSetProxy()) != 2)
          dserror("Scatra discretization has illegal number of dofsets!");

        // get linear solver id from SCALAR TRANSPORT DYNAMIC
        const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
        if (linsolvernumber == -1)
        {
          dserror(
              "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
              "TRANSPORT DYNAMIC to a valid number!");
        }

        // create an ELCH::MovingBoundaryAlgorithm instance
        auto elch = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(
            comm, elchcontrol, scatradyn, problem->SolverParams(linsolvernumber)));

        // now we must call Init()
        // NOTE : elch reads time parameters from scatra dynamic section !
        elch->Init(scatradyn, scatradyn, problem->SolverParams(linsolvernumber), "scatra", true);
        elch->ScaTraField()->SetNumberOfDofSetVelocity(1);

        // NOTE : At this point we may redistribute and/or
        //        ghost our discretizations at will.
        scatradis->FillComplete();
        fluiddis->FillComplete();
        aledis->FillComplete();

        // now we can call Setup() on the scatra time integrator
        elch->Setup();

        // read the restart information, set vectors and variables
        if (restart) elch->ReadRestart(restart);

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
          dserror(
              "no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR "
              "TRANSPORT DYNAMIC to a valid number!");
        }

        // create an ELCH::Algorithm instance
        auto elch = Teuchos::rcp(new ELCH::Algorithm(
            comm, elchcontrol, scatradyn, fdyn, problem->SolverParams(linsolvernumber)));

        // add proxy of fluid degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          dserror("Scatra discretization has illegal number of dofsets!");

        // now we must call Init()
        // NOTE : elch reads time parameters from scatra dynamic section !
        elch->Init(scatradyn, scatradyn, problem->SolverParams(linsolvernumber));
        elch->ScaTraField()->SetNumberOfDofSetVelocity(1);

        // NOTE : At this point we may redistribute and/or
        //        ghost our discretizations at will.
        scatradis->FillComplete();
        fluiddis->FillComplete();
        aledis->FillComplete();

        // discretizations are done, now we can call Setup() on the algorithm
        elch->Setup();

        // read the restart information, set vectors and variables
        if (restart) elch->ReadRestart(restart);

        // solve the whole electrochemistry problem
        elch->TimeLoop();

        // summarize the performance measurements
        Teuchos::TimeMonitor::summarize();

        // perform the result test
        elch->TestResults();
      }

      break;
    }  // case 2
    default:
      dserror("Unknown velocity field type for transport of passive scalar: %d", veltype);
      break;
  }
}


/*----------------------------------------------------------------------*/
// print ELCH-Module logo
/*----------------------------------------------------------------------*/
void printlogo()
{
  // more at http://www.ascii-art.de under entry "moose" (or "elk")
  std::cout << "     ___            ___    " << std::endl;
  std::cout << "    /   \\          /   \\ " << std::endl;
  std::cout << "    \\_   \\        /  __/ " << std::endl;
  std::cout << "     _\\   \\      /  /__  "
            << "     _____ _     _____  _   _   " << std::endl;
  std::cout << "     \\___  \\____/   __/  "
            << "    |  ___| |   /  __ \\| | | |  " << std::endl;
  std::cout << "         \\_       _/      "
            << "   | |__ | |   | /  \\/| |_| |  " << std::endl;
  std::cout << "           | @ @  \\_      "
            << "   |  __|| |   | |    |  _  |   " << std::endl;
  std::cout << "           |               "
            << "  | |___| |___| \\__/\\| | | | " << std::endl;
  std::cout << "         _/     /\\        "
            << "   \\____/\\_____/\\____/\\_| |_/ " << std::endl;
  std::cout << "        /o)  (o/\\ \\_     " << std::endl;
  std::cout << "        \\_____/ /         " << std::endl;
  std::cout << "          \\____/          " << std::endl;
  std::cout << "                           " << std::endl;
}
