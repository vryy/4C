/*---------------------------------------------------------------------------*/
/*!
\brief control routine for particle structure interaction problems

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_dyn.H"
#include "pasi_partitioned_onewaycoup.H"
#include "pasi_partitioned_twowaycoup.H"
#include "pasi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"

#include "../drt_inpar/inpar_pasi.H"

/*---------------------------------------------------------------------------*
 | control routine for particle structure interaction         sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
void pasi_dyn()
{
  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // create a communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // print pasi logo to screen
  if (comm.MyPID() == 0) PASI::UTILS::Logo();

  // get parameter lists
  const Teuchos::ParameterList& pasi_params = problem->PASIDynamicParams();
  const Teuchos::ParameterList& particle_params = problem->ParticleParamsOld();
  const Teuchos::ParameterList& struct_params = problem->StructuralDynamicParams();

  // set time parameters from pasi parameter list for subproblems
  PASI::UTILS::ChangeTimeParameter(comm, pasi_params,
      const_cast<Teuchos::ParameterList&>(particle_params),
      const_cast<Teuchos::ParameterList&>(struct_params));

  // create particle structure interaction algorithm
  Teuchos::RCP<PASI::PartitionedAlgo> algo = Teuchos::null;

  // choose algorithm
  int coupling = DRT::INPUT::IntegralValue<int>(pasi_params, "COUPALGO");

  // query algorithm
  switch (coupling)
  {
    case INPAR::PASI::partitioned_onewaycoup:
    {
      algo = Teuchos::rcp(new PASI::PASI_PartOneWayCoup(comm, pasi_params));
      break;
    }
    case INPAR::PASI::partitioned_twowaycoup:
    {
      algo = Teuchos::rcp(new PASI::PASI_PartTwoWayCoup(comm, pasi_params));
      break;
    }
    case INPAR::PASI::partitioned_twowaycoup_forcerelax:
    {
      algo = Teuchos::rcp(new PASI::PASI_PartTwoWayCoup_ForceRelax(comm, pasi_params));
      break;
    }
    case INPAR::PASI::partitioned_twowaycoup_forcerelaxaitken:
    {
      algo = Teuchos::rcp(new PASI::PASI_PartTwoWayCoup_ForceRelaxAitken(comm, pasi_params));
      break;
    }
    default:
    {
      dserror("no valid coupling type for particle structure interaction specified!");
      break;
    }
  }

  // init pasi algorithm
  algo->Init();

  // setup pasi algorithm
  algo->Setup();

  // read restart information
  const int restart = problem->Restart();
  if (restart) algo->ReadRestart(restart);

  // solve partitioned particle structure interaction
  algo->Timeloop();

  // perform result tests
  algo->TestResults(comm);

  // print summary statistics for all timers
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}
