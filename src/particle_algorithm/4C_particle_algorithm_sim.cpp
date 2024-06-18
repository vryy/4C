/*---------------------------------------------------------------------------*/
/*! \file
\brief main control routine for particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_algorithm_sim.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_particle_algorithm.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
void particle_drt()
{
  // get instance of global problem
  Global::Problem* problem = Global::Problem::Instance();

  // get local communicator
  const Epetra_Comm& comm = *problem->GetCommunicators()->LocalComm().get();

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
  if (restart) particlealgorithm->read_restart(restart);

  // setup particle algorithm
  particlealgorithm->setup();

  // solve particle problem
  particlealgorithm->Timeloop();

  // perform result tests
  {
    // create particle field specific result test objects
    std::vector<std::shared_ptr<Core::UTILS::ResultTest>> allresulttests =
        particlealgorithm->CreateResultTests();

    // add particle field specific result test objects
    for (auto& resulttest : allresulttests)
      if (resulttest) problem->AddFieldTest(Teuchos::rcp(resulttest));

    // perform all tests
    problem->TestAll(comm);
  }

  // print summary statistics for all timers
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
