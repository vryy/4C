/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure interaction

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssi_dyn.H"
#include "ssi_monolithic.H"
#include "ssi_partitioned_1wc.H"
#include "ssi_partitioned_2wc.H"
#include "ssi_utils.H"
#include "io_control.H"
#include "globalproblem.H"
#include "utils_parallel.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssi_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // 1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // 2.- Parameter reading
  auto& ssiparams = const_cast<Teuchos::ParameterList&>(problem->SSIControlParams());
  // access scatra params list
  auto& scatradyn = const_cast<Teuchos::ParameterList&>(problem->ScalarTransportDynamicParams());
  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  auto& sdyn =
      const_cast<Teuchos::ParameterList&>(DRT::Problem::Instance()->StructuralDynamicParams());

  // introduce additional scatra field on manifold?
  const bool is_scatra_manifold =
      DRT::INPUT::IntegralValue<bool>(ssiparams.sublist("MANIFOLD"), "ADD_MANIFOLD");

  // Modification of time parameter list
  SSI::UTILS::ChangeTimeParameter(comm, ssiparams, scatradyn, sdyn);

  const auto coupling =
      Teuchos::getIntegralValue<INPAR::SSI::SolutionSchemeOverFields>(ssiparams, "COUPALGO");


  // 3.- Creation of Structure + Scalar_Transport problem.
  Teuchos::RCP<SSI::SSIBase> ssi = Teuchos::null;

  // by default we employ an ale formulation for our scatra
  bool isale = true;

  // 3.1 choose algorithm depending on solution type
  switch (coupling)
  {
    case INPAR::SSI::SolutionSchemeOverFields::ssi_OneWay_ScatraToSolid:
    {
      ssi = Teuchos::rcp(new SSI::SSIPart1WCScatraToSolid(comm, ssiparams));
      isale = false;
    }
    break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_OneWay_SolidToScatra:
      ssi = Teuchos::rcp(new SSI::SSIPart1WCSolidToScatra(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_IterStagg:
      ssi = Teuchos::rcp(new SSI::SSIPart2WC(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_IterStaggFixedRel_ScatraToSolid:
      ssi = Teuchos::rcp(new SSI::SSIPart2WCScatraToSolidRelax(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_IterStaggFixedRel_SolidToScatra:
      ssi = Teuchos::rcp(new SSI::SSIPart2WCSolidToScatraRelax(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_IterStaggAitken_ScatraToSolid:
      ssi = Teuchos::rcp(new SSI::SSIPart2WCScatraToSolidRelaxAitken(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_IterStaggAitken_SolidToScatra:
      ssi = Teuchos::rcp(new SSI::SSIPart2WCSolidToScatraRelaxAitken(comm, ssiparams));
      break;
    case INPAR::SSI::SolutionSchemeOverFields::ssi_Monolithic:
      ssi = Teuchos::rcp(new SSI::SSIMono(comm, ssiparams));
      break;
    default:
      dserror("unknown coupling algorithm for SSI!");
      break;
  }

  // 3.1.1 initial FillComplete
  problem->GetDis("structure")->FillComplete(true, true, true);
  problem->GetDis("scatra")->FillComplete(true, true, true);
  if (is_scatra_manifold) problem->GetDis("scatra_manifold")->FillComplete(true, true, true);

  // 3.1.2 init the chosen ssi algorithm
  // Construct time integrators of subproblems inside.
  ssi->Init(comm, ssiparams, scatradyn, sdyn, "structure", "scatra", isale);

  // now we can finally fill our discretizations
  // reinitialization of the structural elements is
  // vital for parallelization here!
  problem->GetDis("structure")->FillComplete(true, true, true);
  problem->GetDis("scatra")->FillComplete(true, false, true);
  if (is_scatra_manifold) problem->GetDis("scatra_manifold")->FillComplete(true, false, true);

  DRT::UTILS::PrintParallelDistribution(*problem->GetDis("structure"));
  DRT::UTILS::PrintParallelDistribution(*problem->GetDis("scatra"));
  if (is_scatra_manifold)
    DRT::UTILS::PrintParallelDistribution(*problem->GetDis("scatra_manifold"));

  // 3.1.4 Setup the coupled problem
  // now as we redistributed our discretizations we can construct all
  // objects relying on the parallel distribution
  ssi->Setup();

  // 3.2- Read restart if needed. (Discretization called inside)
  if (ssi->IsRestart())
  {
    const int restart = problem->Restart();
    const double restarttime = problem->RestartTime();
    if (restarttime > 0.0)
      ssi->ReadRestartfromTime(restarttime);
    else if (restart)
      ssi->ReadRestart(restart);
  }

  // 3.3 AFTER restart: reset input filename of the problem so that results from other runs can be
  // read
  bool flag_readscatra = DRT::INPUT::IntegralValue<bool>(ssiparams, "SCATRA_FROM_RESTART_FILE");
  if (coupling == INPAR::SSI::SolutionSchemeOverFields::ssi_OneWay_ScatraToSolid and
      flag_readscatra)
  {
    std::string filename = Teuchos::getNumericStringParameter(ssiparams, "SCATRA_FILENAME");
    Teuchos::RCP<IO::InputControl> inputscatra = Teuchos::rcp(new IO::InputControl(filename, comm));
    problem->SetInputControlFile(inputscatra);
  }

  // 4.- Run of the actual problem.
  // 4.1.- Some setup needed for the subproblems.
  ssi->SetupSystem();

  // 4.2.- Solve the whole problem
  ssi->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  ssi->TestResults(comm);
}
