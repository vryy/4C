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
 | initialize modified boundary particle states               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleAdami::InitBoundaryParticles(
    std::vector<double>& gravity) const
{
  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // evaluation only for boundary and rigid particles
    if (type_i != PARTICLEENGINE::BoundaryPhase and type_i != PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear modified boundary particle states
    container_i->ClearState(PARTICLEENGINE::BoundaryPressure);
    container_i->ClearState(PARTICLEENGINE::BoundaryVelocity);

    // boundary particles with neighbors
    const auto& currboundaryparticles = typeIt.second;

    // iterate over boundary particles
    for (auto& boundaryParticleIt : currboundaryparticles)
    {
      // get local index of boundary particle i
      const int boundaryparticle_i = boundaryParticleIt.first;

      // initialize sum of evaluated kernel values for boundary particle i due to neighbor particles
      // j
      double sumj_Wij(0.0);
      double sumj_press_j_Wij(0.0);
      double sumj_dens_j_r_ij_Wij[3];
      double sumj_vel_j_Wij[3];
      for (int i = 0; i < 3; ++i)
      {
        sumj_dens_j_r_ij_Wij[i] = 0.0;
        sumj_vel_j_Wij[i] = 0.0;
      }

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : boundaryParticleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for neighboring boundary or rigid particles
        if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
          continue;

        // iterate over particle status of neighboring particles
        for (auto& neighborStatusIt : neighborTypeIt.second)
        {
          // get status of neighboring particles of current type
          PARTICLEENGINE::StatusEnum status_j = neighborStatusIt.first;

          // get container of neighboring particles of current particle type and state
          PARTICLEENGINE::ParticleContainerShrdPtr container_j =
              particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

          // iterate over neighboring particles of current type and status
          for (auto& neighborParticleIt : neighborStatusIt.second)
          {
            // get local index of neighbor particle j
            const int particle_j = neighborParticleIt.first;

            // get reference to particle pair
            const ParticlePairSPH& particlepair = neighborParticleIt.second;

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
        }
      }

      // set modified boundary particle states
      if (sumj_Wij > 0.0)
      {
        // declare pointer variables for boundary particle i
        const double *vel_i, *acc_i;
        double *boundarypress_i, *boundaryvel_i;

        // get pointer to particle states
        vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, boundaryparticle_i);
        acc_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, boundaryparticle_i);
        boundarypress_i = container_i->GetPtrToParticleState(
            PARTICLEENGINE::BoundaryPressure, boundaryparticle_i);
        boundaryvel_i = container_i->GetPtrToParticleState(
            PARTICLEENGINE::BoundaryVelocity, boundaryparticle_i);

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
  RefreshModifiedStatesOfBoundaryParticles();
}

/*---------------------------------------------------------------------------*
 | refresh modified states of ghosted boundary particles      sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleAdami::RefreshModifiedStatesOfBoundaryParticles() const
{
  // init map
  std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;

  // set state enums to map
  for (auto& type : typestoconsider_)
    particlestatestotypes[type] = {
        PARTICLEENGINE::BoundaryPressure, PARTICLEENGINE::BoundaryVelocity};

  // refresh specific states of particles of specific types
  particleengineinterface_->RefreshSpecificStatesOfParticlesOfSpecificTypes(particlestatestotypes);
}
