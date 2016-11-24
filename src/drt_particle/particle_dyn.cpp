/*!----------------------------------------------------------------------
\file particle_dyn.cpp
\brief Main control routine for particle simulations

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237

\level 1
</pre>
*----------------------------------------------------------------------*/


#include "particle_dyn.H"
#include "particle_algorithm.H"
#include "cavitation_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
// entry point for particle simulations
/*----------------------------------------------------------------------*/
void particle_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("particle")->Comm();

  // what's the current problem type?
  PROBLEM_TYP probtype = problem->ProblemType();

  switch (probtype)
  {
    case prb_meshfree:
    {
      const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParams();

      problem->GetDis("rendering")->FillComplete();
      /// algorithm is created
      Teuchos::RCP<PARTICLE::Algorithm> particlesimulation = Teuchos::rcp(new PARTICLE::Algorithm(comm,params));

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
    }
    break;
    case prb_particle:
    {
      const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParams();
      /// algorithm is created
      Teuchos::RCP<PARTICLE::Algorithm> particlesimulation = Teuchos::rcp(new PARTICLE::Algorithm(comm,params));

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
    }
    break;
    case prb_cavitation:
    {
      const Teuchos::ParameterList& params = DRT::Problem::Instance()->CavitationParams();

      problem->GetDis("fluid")->FillComplete();

      /// algorithm is created
      Teuchos::RCP<CAVITATION::Algorithm> cavitation = Teuchos::rcp(new CAVITATION::Algorithm(comm,params));

      /// init cavitation simulation
      cavitation->InitCavitation();

      /// read the restart information, set vectors and variables ---
      const int restart = problem->Restart();
      if (restart)
      {
        cavitation->ReadRestart(restart);
      }

      /// setup cavitation simulation
      cavitation->SetupSystem();

      /// solve the whole problem
      cavitation->Timeloop();

      /// perform the result test
      cavitation->TestResults(comm);
    }
    break;
    default:
      dserror("solution of unknown problemtyp %d requested", probtype);
    break;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);

  return;

}
