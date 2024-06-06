/*----------------------------------------------------------------------*/
/*! \file

\brief entry point for scalar transport problems

\level 1


*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_resulttest.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, incl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditions)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *        o scalar transport in velocity field given by Navier-Stokes with natural convection
 *          (two-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = Global::Problem::Instance()->GetDis("fluid")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################" << '\n';
    std::cout << "# YOUR PROBLEM TYPE: " << Global::Problem::Instance()->ProblemName() << '\n';
    std::cout << "###################################################" << '\n';
  }

  // access the problem-specific parameter list
  const auto& scatradyn = Global::Problem::Instance()->scalar_transport_dynamic_params();

  // access the fluid discretization
  auto fluiddis = Global::Problem::Instance()->GetDis("fluid");
  // access the scatra discretization
  auto scatradis = Global::Problem::Instance()->GetDis("scatra");

  // ensure that all dofs are assigned in the right order;
  // this creates dof numbers with fluid dof < scatra dof
  fluiddis->fill_complete(true, true, true);
  scatradis->fill_complete(true, true, true);

  // determine coupling type
  const auto fieldcoupling = Core::UTILS::IntegralValue<Inpar::ScaTra::FieldCoupling>(
      Global::Problem::Instance()->scalar_transport_dynamic_params(), "FIELDCOUPLING");

  // determine velocity type
  const auto veltype =
      Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(scatradyn, "VELOCITYFIELD");

  if (scatradis->NumGlobalNodes() == 0)
  {
    if (fieldcoupling != Inpar::ScaTra::coupling_match and
        veltype != Inpar::ScaTra::velocity_Navier_Stokes)
    {
      FOUR_C_THROW(
          "If you want matching fluid and scatra meshes, do clone you fluid mesh and use "
          "FIELDCOUPLING match!");
    }
  }
  else
  {
    if (fieldcoupling != Inpar::ScaTra::coupling_volmortar and
        veltype == Inpar::ScaTra::velocity_Navier_Stokes)
    {
      FOUR_C_THROW(
          "If you want non-matching fluid and scatra meshes, "
          "you need to use FIELDCOUPLING volmortar!");
    }
  }

  switch (veltype)
  {
    case Inpar::ScaTra::velocity_zero:      // zero  (see case 1)
    case Inpar::ScaTra::velocity_function:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes() == 0)
        FOUR_C_THROW("No elements in the ---TRANSPORT ELEMENTS section");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
      {
        FOUR_C_THROW(
            "no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in "
            "SCALAR TRANSPORT DYNAMIC to a valid number!");
      }

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      auto scatraonly = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
          scatradyn, scatradyn, Global::Problem::Instance()->SolverParams(linsolvernumber)));

      // add proxy of velocity related degrees of freedom to scatra discretization
      auto dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
          Global::Problem::Instance()->NDim() + 1, 0, 0, true));
      if (scatradis->AddDofSet(dofsetaux) != 1)
        FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      scatraonly->ScaTraField()->set_number_of_dof_set_velocity(1);

      // allow TRANSPORT conditions, too
      // NOTE: we can not use the conditions given by 'conditions_to_copy =
      // clonestrategy.conditions_to_copy()' since we than may have some scatra condition twice. So
      // we only copy the Dirichlet and Neumann conditions:
      const std::map<std::string, std::string> conditions_to_copy = {
          {"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
          {"TransportLineNeumann", "LineNeumann"}, {"TransportSurfaceNeumann", "SurfaceNeumann"},
          {"TransportVolumeNeumann", "VolumeNeumann"}};

      Core::FE::DiscretizationCreatorBase creator;
      creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

      // finalize discretization
      scatradis->fill_complete(true, false, true);

      // now we can call Init() on the base algo.
      // time integrator is initialized inside
      scatraonly->Init();

      // redistribution between Init(...) and Setup()
      // redistribute scatra elements in case of heterogeneous reactions
      if (scatradis->GetCondition("ScatraHeteroReactionSlave") != nullptr)
      {
        // create vector of discr.
        std::vector<Teuchos::RCP<Discret::Discretization>> dis;
        dis.push_back(scatradis);

        Core::Rebalance::RebalanceDiscretizationsByBinning(dis, false);
      }

      // assign degrees of freedom and rebuild geometries
      scatradis->fill_complete(true, false, true);

      // now we must call Setup()
      scatraonly->Setup();

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->read_restart(restart);

      // set initial velocity field
      // note: The order read_restart() before set_velocity_field() is important here!!
      // for time-dependent velocity fields, set_velocity_field() is additionally called in each
      // prepare_time_step()-call
      scatraonly->ScaTraField()->set_velocity_field();

      // set external force
      if (scatraonly->ScaTraField()->HasExternalForce())
        scatraonly->ScaTraField()->SetExternalForce();

      // enter time loop to solve problem with given convective velocity
      scatraonly->ScaTraField()->TimeLoop();

      // perform the result test if required
      scatraonly->ScaTraField()->TestResults();
      break;
    }
    case Inpar::ScaTra::velocity_Navier_Stokes:  // Navier_Stokes
    {
      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      // create scatra elements by cloning from fluid dis in matching case
      if (fieldcoupling == Inpar::ScaTra::coupling_match)
      {
        // fill scatra discretization by cloning fluid discretization
        Core::FE::CloneDiscretization<ScaTra::ScatraFluidCloneStrategy>(
            fluiddis, scatradis, Global::Problem::Instance()->CloningMaterialMap());

        // set implementation type of cloned scatra elements
        for (int i = 0; i < scatradis->NumMyColElements(); ++i)
        {
          auto* element = dynamic_cast<Discret::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if (element == nullptr)
            FOUR_C_THROW("Invalid element type!");
          else
            element->SetImplType(Inpar::ScaTra::impltype_std);
        }
      }

      // support for turbulent flow statistics
      const Teuchos::ParameterList& fdyn = (Global::Problem::Instance()->FluidDynamicParams());

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
      {
        FOUR_C_THROW(
            "no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in "
            "SCALAR TRANSPORT DYNAMIC to a valid number!");
      }

      // create a scalar transport algorithm instance
      auto algo = Teuchos::rcp(new ScaTra::ScaTraAlgorithm(comm, scatradyn, fdyn, "scatra",
          Global::Problem::Instance()->SolverParams(linsolvernumber)));

      // create scatra elements by cloning from fluid dis in matching case
      if (fieldcoupling == Inpar::ScaTra::coupling_match)
      {
        // add proxy of fluid transport degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
        algo->ScaTraField()->set_number_of_dof_set_velocity(1);
      }

      // we create  the aux dofsets before Init(...)
      // volmortar adapter Init(...) relies on this
      if (fieldcoupling == Inpar::ScaTra::coupling_volmortar)
      {
        // allow TRANSPORT conditions, too
        ScaTra::ScatraFluidCloneStrategy clonestrategy;
        const auto conditions_to_copy = clonestrategy.conditions_to_copy();
        Core::FE::DiscretizationCreatorBase creator;
        creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

        // build the element and node maps
        scatradis->fill_complete(false, false, false);
        fluiddis->fill_complete(false, false, false);

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
        const int ndofperelement_scatra = 0;
        const int ndofpernode_fluid = fluiddis->NumDof(0, fluiddis->lRowNode(0));
        const int ndofperelement_fluid = 0;

        // add proxy of velocity related degrees of freedom to scatra discretization
        Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (fluiddis->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in fluid field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        algo->ScaTraField()->set_number_of_dof_set_velocity(1);

        // call assign_degrees_of_freedom also for auxiliary dofsets
        // note: the order of fill_complete() calls determines the gid numbering!
        // 1. fluid dofs
        // 2. scatra dofs
        // 3. fluid auxiliary dofs
        // 4. scatra auxiliary dofs
        fluiddis->fill_complete(true, false, false);
        scatradis->fill_complete(true, false, false);
      }

      // init algo (init fluid time integrator and scatra time integrator inside)
      algo->Init();

      // redistribution between Init(...) and Setup()
      // redistribute scatra elements if the scatra discretization is not empty
      if (fieldcoupling == Inpar::ScaTra::coupling_volmortar)
      {
        // create vector of discr.
        std::vector<Teuchos::RCP<Discret::Discretization>> dis;
        dis.push_back(fluiddis);
        dis.push_back(scatradis);

        Core::Rebalance::RebalanceDiscretizationsByBinning(dis, false);
      }

      // ensure that all dofs are assigned in the right order;
      // this creates dof numbers with fluid dof < scatra dof
      fluiddis->fill_complete(true, false, true);
      scatradis->fill_complete(true, false, true);

      // setup algo
      //(setup fluid time integrator and scatra time integrator inside)
      algo->Setup();

      // read restart information
      // in case a inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        if (Core::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") and
            restart == fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP"))
          algo->ReadInflowRestart(restart);
        else
          algo->read_restart(restart);
      }
      else if (Core::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW"))
      {
        FOUR_C_THROW(
            "Turbulent inflow generation for passive scalar transport should be performed as fluid "
            "problem!");
      }

      // solve the whole scalar transport problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      algo->TestResults();

      break;
    }
    default:
    {
      FOUR_C_THROW("unknown velocity field type for transport of passive scalar");
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
