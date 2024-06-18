/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure interaction

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "4C_ssi_dyn.hpp"

#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_rebalance_print.hpp"
#include "4C_ssi_monolithic.hpp"
#include "4C_ssi_partitioned_1wc.hpp"
#include "4C_ssi_partitioned_2wc.hpp"
#include "4C_ssi_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssi_drt()
{
  // 1.- Initialization
  Teuchos::RCP<SSI::SSIBase> ssi = Teuchos::null;
  Global::Problem* problem = Global::Problem::Instance();
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  {
    TEUCHOS_FUNC_TIME_MONITOR("SSI: setup");

    // 2.- Parameter reading
    auto& ssiparams = const_cast<Teuchos::ParameterList&>(problem->SSIControlParams());
    // access scatra params list
    auto& scatradyn =
        const_cast<Teuchos::ParameterList&>(problem->scalar_transport_dynamic_params());
    // access structural dynamic params list which will be possibly modified while creating the time
    // integrator
    auto& sdyn = const_cast<Teuchos::ParameterList&>(
        Global::Problem::Instance()->structural_dynamic_params());

    // introduce additional scatra field on manifold?
    const bool is_scatra_manifold =
        Core::UTILS::IntegralValue<bool>(ssiparams.sublist("MANIFOLD"), "ADD_MANIFOLD");

    // Modification of time parameter list
    SSI::UTILS::ChangeTimeParameter(comm, ssiparams, scatradyn, sdyn);

    const auto coupling =
        Teuchos::getIntegralValue<Inpar::SSI::SolutionSchemeOverFields>(ssiparams, "COUPALGO");


    // 3.- Creation of Structure + Scalar_Transport problem.
    // by default we employ an ale formulation for our scatra
    bool isale = true;

    // 3.1 choose algorithm depending on solution type
    switch (coupling)
    {
      case Inpar::SSI::SolutionSchemeOverFields::ssi_OneWay_ScatraToSolid:
      {
        ssi = Teuchos::rcp(new SSI::SSIPart1WCScatraToSolid(comm, ssiparams));
        isale = false;
      }
      break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_OneWay_SolidToScatra:
        ssi = Teuchos::rcp(new SSI::SSIPart1WCSolidToScatra(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_IterStagg:
        ssi = Teuchos::rcp(new SSI::SSIPart2WC(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_IterStaggFixedRel_ScatraToSolid:
        ssi = Teuchos::rcp(new SSI::SSIPart2WCScatraToSolidRelax(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_IterStaggFixedRel_SolidToScatra:
        ssi = Teuchos::rcp(new SSI::SSIPart2WCSolidToScatraRelax(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_IterStaggAitken_ScatraToSolid:
        ssi = Teuchos::rcp(new SSI::SSIPart2WCScatraToSolidRelaxAitken(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_IterStaggAitken_SolidToScatra:
        ssi = Teuchos::rcp(new SSI::SSIPart2WCSolidToScatraRelaxAitken(comm, ssiparams));
        break;
      case Inpar::SSI::SolutionSchemeOverFields::ssi_Monolithic:
        ssi = Teuchos::rcp(new SSI::SsiMono(comm, ssiparams));
        break;
      default:
        FOUR_C_THROW("unknown coupling algorithm for SSI!");
        break;
    }

    // 3.1.1 initial fill_complete
    problem->GetDis("structure")->fill_complete(true, true, true);
    problem->GetDis("scatra")->fill_complete(true, true, true);
    if (is_scatra_manifold) problem->GetDis("scatra_manifold")->fill_complete(true, true, true);

    // 3.1.2 init the chosen ssi algorithm
    // Construct time integrators of subproblems inside.
    ssi->init(comm, ssiparams, scatradyn, sdyn, "structure", "scatra", isale);

    // now we can finally fill our discretizations
    // reinitialization of the structural elements is
    // vital for parallelization here!
    problem->GetDis("structure")->fill_complete(true, true, true);
    problem->GetDis("scatra")->fill_complete(true, false, true);
    if (is_scatra_manifold) problem->GetDis("scatra_manifold")->fill_complete(true, false, true);

    Core::Rebalance::UTILS::print_parallel_distribution(*problem->GetDis("structure"));
    Core::Rebalance::UTILS::print_parallel_distribution(*problem->GetDis("scatra"));
    if (is_scatra_manifold)
      Core::Rebalance::UTILS::print_parallel_distribution(*problem->GetDis("scatra_manifold"));

    // 3.1.4 Setup the coupled problem
    // now as we redistributed our discretizations we can construct all
    // objects relying on the parallel distribution
    ssi->setup();

    // 3.2- Read restart if needed. (discretization called inside)
    if (ssi->IsRestart()) ssi->read_restart(problem->restart());

    // 3.3 AFTER restart: reset input filename of the problem so that results from other runs can be
    // read
    bool flag_readscatra = Core::UTILS::IntegralValue<bool>(ssiparams, "SCATRA_FROM_RESTART_FILE");
    if (coupling == Inpar::SSI::SolutionSchemeOverFields::ssi_OneWay_ScatraToSolid and
        flag_readscatra)
    {
      std::string filename = Teuchos::getNumericStringParameter(ssiparams, "SCATRA_FILENAME");
      auto inputscatra = Teuchos::rcp(new Core::IO::InputControl(filename, comm));
      problem->SetInputControlFile(inputscatra);
    }

    // 4.- Run of the actual problem.
    // 4.1.- Some setup needed for the subproblems.
    ssi->SetupSystem();
  }

  // 4.2.- Solve the whole problem
  ssi->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  ssi->TestResults(comm);
}

FOUR_C_NAMESPACE_CLOSE
