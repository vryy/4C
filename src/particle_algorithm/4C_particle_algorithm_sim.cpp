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
  Global::Problem* problem = Global::Problem::instance();

  // get local communicator
  const Epetra_Comm& comm = *problem->get_communicators()->local_comm().get();

  // get parameter lists
  const Teuchos::ParameterList& params = problem->particle_params();

  // reference to vector of initial particles
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles = problem->particles();

  // create and init particle algorithm
  auto particlealgorithm = std::unique_ptr<PARTICLEALGORITHM::ParticleAlgorithm>(
      new PARTICLEALGORITHM::ParticleAlgorithm(comm, params));
  particlealgorithm->init(initialparticles);

  // read restart information
  const int restart = problem->restart();
  if (restart) particlealgorithm->read_restart(restart);

  // setup particle algorithm
  particlealgorithm->setup();

  // solve particle problem
  particlealgorithm->timeloop();

  // perform result tests
  {
    // create particle field specific result test objects
    std::vector<std::shared_ptr<Core::UTILS::ResultTest>> allresulttests =
        particlealgorithm->create_result_tests();

    // add particle field specific result test objects
    for (auto& resulttest : allresulttests)
      if (resulttest) problem->add_field_test(Teuchos::rcp(resulttest));

    // perform all tests
    problem->test_all(comm);
  }

  // print summary statistics for all timers
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
