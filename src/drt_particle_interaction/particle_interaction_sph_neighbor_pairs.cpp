/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_neighbor_pairs.cpp

\brief neighbor pair handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_sph_kernel.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHNeighborPairs::SPHNeighborPairs()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init neighbor pair handler                                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHNeighborPairs::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup neighbor pair handler                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHNeighborPairs::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set kernel handler
  kernel_ = kernel;

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold neighbor pairs
  neighborpairdata_.resize(typevectorsize);
}

/*---------------------------------------------------------------------------*
 | write restart of neighbor pair handler                     sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHNeighborPairs::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of neighbor pair handler                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHNeighborPairs::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | evaluate neighbor pairs                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHNeighborPairs::EvaluateNeighborPairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHNeighborPairs::EvaluateNeighborPairs");

  // get reference to particle neighbors
  const PARTICLEENGINE::ParticleNeighbors& particleneighbors =
      particleengineinterface_->GetParticleNeighbors();

  // iterate over particle types
  for (auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    int particlestored = container_i->ParticlesStored();

    // allocate memory for neighbor pairs of owned particles of current particle type
    neighborpairdata_[type_i].assign(particlestored, std::vector<SPHNeighborPair>(0));

    // check for owned particles of current type
    if (particleneighbors[type_i].empty()) continue;

    // iterate over particles of current type
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get reference to vector of neighbors of current particle
      const std::vector<PARTICLEENGINE::LocalIndexTuple>& currentNeighbors =
          (particleneighbors[type_i])[particle_i];

      // check for neighbors of owned particles of current type
      if (currentNeighbors.empty()) continue;

      // get reference to vector of neighbor pairs of current particle
      std::vector<SPHNeighborPair>& currentNeighborPairs = (neighborpairdata_[type_i])[particle_i];

      // allocate memory for neighbor pairs of current particle
      currentNeighborPairs.reserve(currentNeighbors.size());

      // declare pointer variables for particle i
      const double *pos_i, *rad_i;

      // get pointer to particle states
      pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

      // iterate over neighboring particles
      for (auto& neighborParticleIt : currentNeighbors)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborParticleIt;

        // no evaluation for neighboring boundary and rigid particles
        if (type_i == PARTICLEENGINE::BoundaryPhase and
            (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase))
          continue;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // get pointer to particle position
        const double* pos_j =
            container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);

        // vector from particle i to j
        double r_ji[3];

        // distance between particles considering periodic boundaries
        particleengineinterface_->DistanceBetweenParticles(pos_i, pos_j, r_ji);

        // absolute distance between particles
        const double absdist = std::sqrt(r_ji[0] * r_ji[0] + r_ji[1] * r_ji[1] + r_ji[2] * r_ji[2]);

        // neighboring particle out of support radius
        if (absdist > rad_i[0]) continue;

        // initialize particle pair
        currentNeighborPairs.push_back(std::make_pair(neighborParticleIt, SPHParticlePair()));

        // get reference to current particle pair
        SPHParticlePair& particlepair = (currentNeighborPairs.back()).second;

        // set absolute distance between particles
        particlepair.absdist_ = absdist;

        // versor from particle j to i
        for (int i = 0; i < 3; ++i) particlepair.e_ij_[i] = -r_ji[i] / absdist;

        // evaluate kernel
        particlepair.Wij_ = kernel_->W(absdist, rad_i[0]);

        // evaluate first derivative of kernel
        particlepair.dWdrij_ = kernel_->dWdrij(absdist, rad_i[0]);
      }

      // free superfluous allocated memory
      currentNeighborPairs.resize(currentNeighborPairs.size());
    }
  }
}
