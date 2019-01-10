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
    const double *pos_i, *pos_j, *rad_i, *rad_j;

    // get pointer to particle states
    pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    pos_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);

    // vector from particle i to j
    double r_ji[3];

    // distance between particles considering periodic boundaries
    particleengineinterface_->DistanceBetweenParticles(pos_i, pos_j, r_ji);

    // absolute distance between particles
    const double absdist = std::sqrt(r_ji[0] * r_ji[0] + r_ji[1] * r_ji[1] + r_ji[2] * r_ji[2]);

    // particle j within support radius of particle i
    if (absdist < rad_i[0])
    {
      // initialize particle pair
      (neighborpairdata_[type_i])[particle_i].push_back(
          std::make_pair(neighborpair.second, SPHParticlePair()));

      // get reference to current particle pair
      SPHParticlePair& particlepair_i = ((neighborpairdata_[type_i])[particle_i].back()).second;

      // set absolute distance between particles
      particlepair_i.absdist_ = absdist;

      // versor from particle j to i
      for (int i = 0; i < 3; ++i) particlepair_i.e_ij_[i] = -r_ji[i] / absdist;

      // evaluate kernel
      particlepair_i.Wij_ = kernel_->W(absdist, rad_i[0]);

      // evaluate first derivative of kernel
      particlepair_i.dWdrij_ = kernel_->dWdrij(absdist, rad_i[0]);
    }

    // particle i within support radius of owned particle j
    if (absdist < rad_j[0] and (status_j == PARTICLEENGINE::Owned))
    {
      // initialize particle pair
      (neighborpairdata_[type_j])[particle_j].push_back(
          std::make_pair(neighborpair.first, SPHParticlePair()));

      // get reference to current particle pair
      SPHParticlePair& particlepair_j = ((neighborpairdata_[type_j])[particle_j].back()).second;

      // set absolute distance between particles
      particlepair_j.absdist_ = absdist;

      // equal support radius for particle i and j
      if (rad_i[0] == rad_j[0])
      {
        // get reference to current particle pair
        SPHParticlePair& particlepair_i = ((neighborpairdata_[type_i])[particle_i].back()).second;

        // versor from particle i to j
        for (int i = 0; i < 3; ++i) particlepair_j.e_ij_[i] = -particlepair_i.e_ij_[i];

        // evaluate kernel
        particlepair_j.Wij_ = particlepair_i.Wij_;

        // evaluate first derivative of kernel
        particlepair_j.dWdrij_ = particlepair_i.dWdrij_;
      }
      else
      {
        // versor from particle i to j
        for (int i = 0; i < 3; ++i) particlepair_j.e_ij_[i] = r_ji[i] / absdist;

        // evaluate kernel
        particlepair_j.Wij_ = kernel_->W(absdist, rad_j[0]);

        // evaluate first derivative of kernel
        particlepair_j.dWdrij_ = kernel_->dWdrij(absdist, rad_j[0]);
      }
    }
  }
}
