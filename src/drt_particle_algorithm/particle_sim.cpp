/*---------------------------------------------------------------------------*/
/*!

\brief main control routine for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_sim.H"

#include "../drt_particle_algorithm/particle_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_resulttest.H"

#include "../drt_comm/comm_utils.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

/*---------------------------------------------------------------------------*
 | main control routine for particle simulations              sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void particle_drt()
{
  // get instance of global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get local communicator
  const Epetra_Comm& comm = *problem->GetNPGroup()->LocalComm().get();

  // get parameter lists
  const Teuchos::ParameterList& params = problem->ParticleParams();

  // reference to vector of initial particles
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles = problem->Particles();

  // create and init particle algorithm
  auto particlealgorithm = std::unique_ptr<PARTICLEALGORITHM::ParticleAlgorithm>(
      new PARTICLEALGORITHM::ParticleAlgorithm(comm, params));
  particlealgorithm->Init(initialparticles);

  // read restart information
  const int restart = problem->Restart();
  if (restart) particlealgorithm->ReadRestart(restart);

  // setup particle algorithm
  particlealgorithm->Setup();

  // solve particle problem
  particlealgorithm->Timeloop();

  // perform result tests
  {
    // create particle field specific result test objects
    std::vector<std::shared_ptr<DRT::ResultTest>> allresulttests =
        particlealgorithm->CreateResultTests();

    // add particle field specific result test objects
    for (auto& resulttest : allresulttests)
      if (resulttest) problem->AddFieldTest(Teuchos::rcp(resulttest));

    // perform all tests
    problem->TestAll(comm);
  }

  // print summary statistics for all timers
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}
