/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_boundary_particle.cpp

\brief boundary particle handler for smoothed particle hydrodynamics (SPH) interactions

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
#include "particle_interaction_sph_boundary_particle.H"

#include "particle_interaction_sph_neighbor_pairs.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHBoundaryParticleBase::SPHBoundaryParticleBase(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init boundary particle handler                             sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup boundary particle handler                            sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // check if boundary and/or rigid particles are present
  if ((particlecontainerbundle_->GetParticleTypes()).count(PARTICLEENGINE::BoundaryPhase))
    typestoconsider_.insert(PARTICLEENGINE::BoundaryPhase);
  if ((particlecontainerbundle_->GetParticleTypes()).count(PARTICLEENGINE::RigidPhase))
    typestoconsider_.insert(PARTICLEENGINE::RigidPhase);

  // safety check
  if (typestoconsider_.size() == 0)
    dserror("no boundary or rigid particles defined but a boundary particle formulation is set!");
}

/*---------------------------------------------------------------------------*
 | write restart of boundary particle handler                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of boundary particle handler                  sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHBoundaryParticleAdami::SPHBoundaryParticleAdami(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHBoundaryParticleBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | setup boundary particle handler                            sfuchs 01/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleAdami::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHBoundaryParticleBase::Setup(particleengineinterface, neighborpairs);

  // setup modified states of ghosted boundary particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{
        PARTICLEENGINE::BoundaryPressure, PARTICLEENGINE::BoundaryVelocity};

    for (auto& typeEnum : typestoconsider_)
      boundarystatestorefresh_.push_back(std::make_pair(typeEnum, states));
  }
}

/*---------------------------------------------------------------------------*
 | initialize modified boundary particle states               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleAdami::InitBoundaryParticles(
    std::vector<double>& gravity) const
{
  // get reference to neighbor pair data
  const SPHNeighborPairData& neighborpairdata = neighborpairs_->GetRefToNeighborPairData();

  // iterate over particle types
  for (auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // evaluation only for boundary and rigid particles
    if (type_i != PARTICLEENGINE::BoundaryPhase and type_i != PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear modified boundary particle states
    container_i->ClearState(PARTICLEENGINE::BoundaryPressure);
    container_i->ClearState(PARTICLEENGINE::BoundaryVelocity);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get reference to vector of neighbor pairs of current particle
      const std::vector<SPHNeighborPair>& currentNeighborPairs =
          (neighborpairdata[type_i])[particle_i];

      // check for neighbor pairs of current particle
      if (currentNeighborPairs.empty()) continue;

      // initialize sum of evaluated kernel values for particle i due to neighbor particles j
      double sumj_Wij(0.0);
      double sumj_press_j_Wij(0.0);
      double sumj_dens_j_r_ij_Wij[3];
      double sumj_vel_j_Wij[3];
      for (int i = 0; i < 3; ++i)
      {
        sumj_dens_j_r_ij_Wij[i] = 0.0;
        sumj_vel_j_Wij[i] = 0.0;
      }

      // iterate over neighbor pairs
      for (auto& neighborIt : currentNeighborPairs)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborIt.first;

        // get reference to current particle pair
        const SPHParticlePair& particlepair = neighborIt.second;

        // no evaluation for neighboring boundary or rigid particles
        if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
          continue;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // declare pointer variables for neighbor particle j
        const double *vel_j, *dens_j, *press_j;

        // get pointer to particle states
        vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
        dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
        press_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_j);

        // sum up contributions
        sumj_Wij += particlepair.Wij_;
        sumj_press_j_Wij += press_j[0] * particlepair.Wij_;

        for (int i = 0; i < 3; ++i)
        {
          sumj_dens_j_r_ij_Wij[i] +=
              dens_j[i] * particlepair.absdist_ * particlepair.e_ij_[i] * particlepair.Wij_;
          sumj_vel_j_Wij[i] += vel_j[i] * particlepair.Wij_;
        }
      }

      // set modified boundary particle states
      if (sumj_Wij > 0.0)
      {
        // declare pointer variables for boundary particle i
        const double *vel_i, *acc_i;
        double *boundarypress_i, *boundaryvel_i;

        // get pointer to particle states
        vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
        acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);
        boundarypress_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::BoundaryPressure, particle_i);
        boundaryvel_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::BoundaryVelocity, particle_i);

        // get relative acceleration of boundary particle
        double relacc[3];
        for (int i = 0; i < 3; ++i) relacc[i] = gravity[i] - acc_i[i];

        // set modified boundary pressure
        boundarypress_i[0] =
            (sumj_press_j_Wij + relacc[0] * sumj_dens_j_r_ij_Wij[0] +
                relacc[1] * sumj_dens_j_r_ij_Wij[1] + relacc[2] * sumj_dens_j_r_ij_Wij[2]) /
            sumj_Wij;

        // set modified boundary velocity
        for (int i = 0; i < 3; ++i)
          boundaryvel_i[i] = 2.0 * vel_i[i] - (sumj_vel_j_Wij[i] / sumj_Wij);
      }
    }
  }

  // refresh modified states of ghosted boundary particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(boundarystatestorefresh_);
}
