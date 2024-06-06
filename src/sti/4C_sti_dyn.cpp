/*----------------------------------------------------------------------*/
/*! \file

\brief main control routine for monolithic scalar-thermo interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_sti_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_scatra_resulttest_elch.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_sti_clonestrategy.hpp"
#include "4C_sti_monolithic.hpp"
#include "4C_sti_partitioned.hpp"
#include "4C_sti_resulttest.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------------*
 | entry point for simulations of scalar-thermo interaction problems   fang 04/15 |
 *--------------------------------------------------------------------------------*/
void sti_dyn(const int& restartstep  //! time step for restart
)
{
  // access global problem
  Global::Problem* problem = Global::Problem::Instance();

  // access communicator
  const Epetra_Comm& comm = problem->GetDis("scatra")->Comm();

  // access scatra discretization
  Teuchos::RCP<Discret::Discretization> scatradis = problem->GetDis("scatra");

  // add dofset for velocity-related quantities to scatra discretization
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(problem->NDim() + 1, 0, 0, true));
  if (scatradis->AddDofSet(dofsetaux) != 1)
    FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");

  // finalize scatra discretization
  scatradis->fill_complete();

  // safety check
  if (scatradis->NumGlobalNodes() == 0)
    FOUR_C_THROW(
        "The scatra discretization must not be empty, since the thermo discretization needs to be "
        "cloned from it!");

  // access thermo discretization
  Teuchos::RCP<Discret::Discretization> thermodis = problem->GetDis("thermo");

  // add dofset for velocity-related quantities to thermo discretization
  dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(problem->NDim() + 1, 0, 0, true));
  if (thermodis->AddDofSet(dofsetaux) != 1)
    FOUR_C_THROW("Thermo discretization has illegal number of dofsets!");

  // equip thermo discretization with noderowmap for subsequent safety check
  // final fill_complete() is called at the end of discretization cloning
  thermodis->fill_complete(false, false, false);

  // safety check
  if (thermodis->NumGlobalNodes() != 0)
    FOUR_C_THROW(
        "The thermo discretization must be empty, since it is cloned from the scatra "
        "discretization!");

  // clone thermo discretization from scatra discretization, using clone strategy for scatra-thermo
  // interaction
  Core::FE::CloneDiscretization<STI::ScatraThermoCloneStrategy>(
      scatradis, thermodis, Global::Problem::Instance()->CloningMaterialMap());
  thermodis->fill_complete(false, true, true);

  // add proxy of scalar transport degrees of freedom to thermo discretization and vice versa
  if (thermodis->AddDofSet(scatradis->GetDofSetProxy()) != 2)
    FOUR_C_THROW("Thermo discretization has illegal number of dofsets!");
  if (scatradis->AddDofSet(thermodis->GetDofSetProxy()) != 2)
    FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
  thermodis->fill_complete(true, false, false);
  scatradis->fill_complete(true, false, false);

  // add material of scatra elements to thermo elements and vice versa
  for (int i = 0; i < scatradis->NumMyColElements(); ++i)
  {
    Core::Elements::Element* scatraele = scatradis->lColElement(i);
    Core::Elements::Element* thermoele = thermodis->gElement(scatraele->Id());

    thermoele->AddMaterial(scatraele->Material());
    scatraele->AddMaterial(thermoele->Material());
  }

  // access parameter lists for scatra-thermo interaction and scalar transport field
  const Teuchos::ParameterList& stidyn = problem->STIDynamicParams();
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  // extract and check ID of linear solver for scatra field
  const int solver_id_scatra = scatradyn.get<int>("LINEAR_SOLVER");
  if (solver_id_scatra == -1)
    FOUR_C_THROW(
        "No linear solver for scalar transport field was specified in input file section 'SCALAR "
        "TRANSPORT DYNAMIC'!");

  // extract and check ID of linear solver for thermo field
  const int solver_id_thermo = stidyn.get<int>("THERMO_LINEAR_SOLVER");
  if (solver_id_thermo == -1)
    FOUR_C_THROW(
        "No linear solver for temperature field was specified in input file section 'STI "
        "DYNAMIC'!");

  // instantiate coupling algorithm for scatra-thermo interaction
  Teuchos::RCP<STI::Algorithm> sti_algorithm(Teuchos::null);

  switch (Teuchos::getIntegralValue<Inpar::STI::CouplingType>(stidyn, "COUPLINGTYPE"))
  {
    // monolithic algorithm
    case Inpar::STI::CouplingType::monolithic:
    {
      // extract and check ID of monolithic linear solver
      const int solver_id = stidyn.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER");
      if (solver_id == -1)
        FOUR_C_THROW(
            "No global linear solver was specified in input file section 'STI "
            "DYNAMIC/MONOLITHIC'!");

      sti_algorithm = Teuchos::rcp(new STI::Monolithic(comm, stidyn, scatradyn,
          Global::Problem::Instance()->SolverParams(solver_id),
          Global::Problem::Instance()->SolverParams(solver_id_scatra),
          Global::Problem::Instance()->SolverParams(solver_id_thermo)));

      break;
    }

    // partitioned algorithm
    case Inpar::STI::CouplingType::oneway_scatratothermo:
    case Inpar::STI::CouplingType::oneway_thermotoscatra:
    case Inpar::STI::CouplingType::twoway_scatratothermo:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit:
    case Inpar::STI::CouplingType::twoway_thermotoscatra:
    case Inpar::STI::CouplingType::twoway_thermotoscatra_aitken:
    {
      sti_algorithm = Teuchos::rcp(new STI::Partitioned(comm, stidyn, scatradyn,
          Global::Problem::Instance()->SolverParams(solver_id_scatra),
          Global::Problem::Instance()->SolverParams(solver_id_thermo)));

      break;
    }

    // unknown algorithm
    default:
    {
      FOUR_C_THROW("Unknown coupling algorithm for scatra-thermo interaction!");
    }
  }

  // read restart data if necessary
  if (restartstep) sti_algorithm->read_restart(restartstep);

  // provide scatra and thermo fields with velocities
  sti_algorithm->ScaTraField()->set_velocity_field();
  sti_algorithm->ThermoField()->set_velocity_field();

  // enter time loop and solve scatra-thermo interaction problem
  sti_algorithm->TimeLoop();

  // summarize performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform result tests
  problem->AddFieldTest(
      Teuchos::rcp<Core::UTILS::ResultTest>(new STI::STIResultTest(sti_algorithm)));
  if (Teuchos::getIntegralValue<Inpar::STI::ScaTraTimIntType>(
          problem->STIDynamicParams(), "SCATRATIMINTTYPE") == Inpar::STI::ScaTraTimIntType::elch)
    problem->AddFieldTest(Teuchos::rcp<Core::UTILS::ResultTest>(new ScaTra::ElchResultTest(
        Teuchos::rcp_dynamic_cast<ScaTra::ScaTraTimIntElch>(sti_algorithm->ScaTraField()))));
  else
    FOUR_C_THROW(
        "Scatra-thermo interaction is currently only available for thermodynamic electrochemistry, "
        "but not for other kinds of thermodynamic scalar transport!");
  problem->AddFieldTest(Teuchos::rcp<Core::UTILS::ResultTest>(
      new ScaTra::ScaTraResultTest(sti_algorithm->ThermoField())));
  problem->TestAll(comm);

  return;
}  // sti_dyn()

FOUR_C_NAMESPACE_CLOSE
