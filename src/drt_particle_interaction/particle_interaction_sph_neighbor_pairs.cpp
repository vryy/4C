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

#include "particle_interaction_utils.H"

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

  // allocate memory to hold index of neighbor pairs for each type
  indexofneighborpairs_.resize(typevectorsize);
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

  // clear neighbor pair data
  neighborpairdata_.clear();

  // clear index of neighbor pairs for each type
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
    indexofneighborpairs_[type_i].clear();

  // index of neighbor pairs
  int neighborpairindex = 0;

  // iterate over potential particle neighbors
  for (auto& potentialneighbors : particleengineinterface_->GetPotentialParticleNeighbors())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = potentialneighbors.first;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = potentialneighbors.second;

    if (type_i == PARTICLEENGINE::BoundaryPhase and type_j == PARTICLEENGINE::BoundaryPhase)
      continue;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // declare pointer variables for particle i and j
    const double *pos_i, *rad_i;
    const double *pos_j, *rad_j;

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
    const double absdist = UTILS::vec_norm2(r_ji);

    if (absdist < rad_i[0] or (absdist < rad_j[0] and status_j == PARTICLEENGINE::Owned))
    {
      // initialize neighbor pair
      neighborpairdata_.push_back(SPHNeighborPair());

      // get reference to current neighbor pair
      SPHNeighborPair& neighborpair = neighborpairdata_[neighborpairindex];

      // store index of neighbor pairs for each type (owned and ghosted status)
      indexofneighborpairs_[type_i].push_back(neighborpairindex);
      if (type_i != type_j) indexofneighborpairs_[type_j].push_back(neighborpairindex);

      // increase index
      ++neighborpairindex;

      // set local index tuple of particles i and j
      neighborpair.tuple_i_ = potentialneighbors.first;
      neighborpair.tuple_j_ = potentialneighbors.second;

      // set absolute distance between particles
      neighborpair.absdist_ = absdist;

      // versor from particle j to i
      UTILS::vec_setscale(neighborpair.e_ij_, -1.0 / absdist, r_ji);

      // particle j within support radius of particle i
      if (absdist < rad_i[0])
      {
        // evaluate kernel
        neighborpair.Wij_ = kernel_->W(absdist, rad_i[0]);

        // evaluate first derivative of kernel
        neighborpair.dWdrij_ = kernel_->dWdrij(absdist, rad_i[0]);
      }

      // particle i within support radius of owned particle j
      if (absdist < rad_j[0] and status_j == PARTICLEENGINE::Owned)
      {
        // equal support radius for particle i and j
        if (rad_i[0] == rad_j[0])
        {
          // evaluate kernel
          neighborpair.Wji_ = neighborpair.Wij_;

          // evaluate first derivative of kernel
          neighborpair.dWdrji_ = neighborpair.dWdrij_;
        }
        else
        {
          // evaluate kernel
          neighborpair.Wji_ = kernel_->W(absdist, rad_j[0]);

          // evaluate first derivative of kernel
          neighborpair.dWdrji_ = kernel_->dWdrij(absdist, rad_j[0]);
        }
      }
    }
  }
}
