/*!----------------------------------------------------------------------
\file particle_dyn.cpp
\brief Main control routine for particle simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*----------------------------------------------------------------------*/


#include "particle_dyn.H"
#include "particle_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
// entry point for particle simulations
/*----------------------------------------------------------------------*/
void particle_drt()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("particle")->Comm();

  /// algorithm is created
  Teuchos::RCP<PARTICLE::Algorithm> particlesimulation = Teuchos::rcp(new PARTICLE::Algorithm(comm));

  /// read the restart information, set vectors and variables ---
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    particlesimulation->ReadRestart(restart);
  }

  /// ???
  particlesimulation->SetupSystem();

  /// solve the whole problem
  particlesimulation->Timeloop();

  /// perform the result test
  particlesimulation->TestResults(comm);

  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);

  return;

}
