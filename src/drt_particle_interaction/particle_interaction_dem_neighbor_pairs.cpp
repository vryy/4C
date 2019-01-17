/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_neighbor_pairs.cpp

\brief neighbor pair handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_neighbor_pairs.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMNeighborPairs::DEMNeighborPairs()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init neighbor pair handler                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMNeighborPairs::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup neighbor pair handler                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMNeighborPairs::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold neighbor pairs
  neighborpairdata_.resize(typevectorsize);
}

/*---------------------------------------------------------------------------*
 | write restart of neighbor pair handler                     sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMNeighborPairs::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of neighbor pair handler                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMNeighborPairs::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | evaluate neighbor pairs                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateNeighborPairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMNeighborPairs::EvaluateNeighborPairs");

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
    neighborpairdata_[type_i].assign(particlestored, std::vector<DEMNeighborPair>(0));

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
      std::vector<DEMNeighborPair>& currentNeighborPairs = (neighborpairdata_[type_i])[particle_i];

      // allocate memory for neighbor pairs of current particle
      currentNeighborPairs.reserve(currentNeighbors.size());

      // declare pointer variables for particle i
      const double *pos_i, *rad_i, *mass_i;

      // get pointer to particle states
      pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

      // iterate over neighboring particles
      for (auto& neighborParticleIt : currentNeighbors)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborParticleIt;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // declare pointer variables for particle j
        const double *pos_j, *rad_j, *mass_j;

        // get pointer to particle states
        pos_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);
        rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
        mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

        // vector from particle i to j
        double r_ji[3];

        // distance between particles considering periodic boundaries
        particleengineinterface_->DistanceBetweenParticles(pos_i, pos_j, r_ji);

        // absolute distance between particles
        const double absdist = std::sqrt(r_ji[0] * r_ji[0] + r_ji[1] * r_ji[1] + r_ji[2] * r_ji[2]);

        // gap between particles
        const double gap = absdist - rad_i[0] - rad_j[0];

        // neighboring particle out of interaction distance
        if (gap > 0.0) continue;

        // initialize particle pair
        currentNeighborPairs.push_back(std::make_pair(neighborParticleIt, DEMParticlePair()));

        // get reference to current particle pair
        DEMParticlePair& particlepair = (currentNeighborPairs.back()).second;

        // set gap between particles
        particlepair.gap_ = gap;

        // versor from particle i to j
        for (int i = 0; i < 3; ++i) particlepair.e_ji_[i] = r_ji[i] / absdist;

        // set effective mass of particles i and j
        particlepair.m_eff_ = mass_i[0] * mass_j[0] / (mass_i[0] + mass_j[0]);
      }

      // free superfluous allocated memory
      currentNeighborPairs.resize(currentNeighborPairs.size());
    }
  }
}
