/*----------------------------------------------------------------------*/
/*! \file

\brief Main control routine for particle simulations

\level 1

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/


#include "particle_dyn.H"
#include "particle_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
// entry point for particle simulations
/*----------------------------------------------------------------------*/
void particle_old_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("particle")->Comm();

  const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParamsOld();
  /// algorithm is created
  Teuchos::RCP<PARTICLE::Algorithm> particlesimulation =
      Teuchos::rcp(new PARTICLE::Algorithm(comm, params));

  /// init particle simulation
  particlesimulation->Init(false);

  /// read the restart information, set vectors and variables ---
  const int restart = problem->Restart();
  if (restart)
  {
    particlesimulation->ReadRestart(restart);
  }

  /// setup particle simulation
  particlesimulation->SetupSystem();

  /// solve the whole problem
  particlesimulation->Timeloop();

  /// perform the result test
  particlesimulation->TestResults(comm);

  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);

  return;
}
