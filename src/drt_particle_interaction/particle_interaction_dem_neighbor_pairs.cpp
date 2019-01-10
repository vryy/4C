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
  }

  // iterate over potential particle neighbors
  for (auto& neighborpair : particleengineinterface_->GetPotentialParticleNeighbors())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = neighborpair.first;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = neighborpair.second;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainerShrdPtr container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // declare pointer variables for particle i and j
    const double *pos_i, *pos_j, *rad_i, *rad_j, *mass_i, *mass_j;

    // get pointer to particle states
    pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

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
    (neighborpairdata_[type_i])[particle_i].push_back(
        std::make_pair(neighborpair.second, DEMParticlePair()));

    // get reference to current particle pair
    DEMParticlePair& particlepair_i = ((neighborpairdata_[type_i])[particle_i].back()).second;

    // set gap between particles
    particlepair_i.gap_ = gap;

    // versor from particle i to j
    for (int i = 0; i < 3; ++i) particlepair_i.e_ji_[i] = r_ji[i] / absdist;

    // set effective mass of particles i and j
    particlepair_i.m_eff_ = mass_i[0] * mass_j[0] / (mass_i[0] + mass_j[0]);

    if (status_j == PARTICLEENGINE::Owned)
    {
      // initialize particle pair
      (neighborpairdata_[type_j])[particle_j].push_back(
          std::make_pair(neighborpair.first, DEMParticlePair()));

      // get reference to current particle pair
      DEMParticlePair& particlepair_j = ((neighborpairdata_[type_j])[particle_j].back()).second;

      // set gap between particles
      particlepair_j.gap_ = gap;

      // versor from particle j to i
      for (int i = 0; i < 3; ++i) particlepair_j.e_ji_[i] = -particlepair_i.e_ji_[i];

      // set effective mass of particles i and j
      particlepair_j.m_eff_ = particlepair_i.m_eff_;
    }
  }
}
