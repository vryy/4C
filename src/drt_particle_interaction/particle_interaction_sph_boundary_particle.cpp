/*---------------------------------------------------------------------------*/
/*!
\brief boundary particle handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_boundary_particle.H"

#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

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
    boundarytypes_.insert(PARTICLEENGINE::BoundaryPhase);
  if ((particlecontainerbundle_->GetParticleTypes()).count(PARTICLEENGINE::RigidPhase))
    boundarytypes_.insert(PARTICLEENGINE::RigidPhase);

  // safety check
  if (boundarytypes_.size() == 0)
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

    for (const auto& typeEnum : boundarytypes_)
      boundarystatestorefresh_.push_back(std::make_pair(typeEnum, states));
  }

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--boundarytypes_.end()) + 1;

  // allocate memory to hold contributions of neighboring particles
  sumj_Wij_.resize(typevectorsize);
  sumj_press_j_Wij_.resize(typevectorsize);
  sumj_dens_j_r_ij_Wij_.resize(typevectorsize);
  sumj_vel_j_Wij_.resize(typevectorsize);
}

/*---------------------------------------------------------------------------*
 | init boundary particle states                              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHBoundaryParticleAdami::InitBoundaryParticleStates(
    std::vector<double>& gravity)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHBoundaryParticleAdami::InitBoundaryParticleStates");

  // iterate over boundary particle types
  for (const auto& type_i : boundarytypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // allocate memory
    sumj_Wij_[type_i].assign(particlestored, 0.0);
    sumj_press_j_Wij_[type_i].assign(particlestored, 0.0);
    sumj_dens_j_r_ij_Wij_[type_i].assign(particlestored, std::vector<double>(3, 0.0));
    sumj_vel_j_Wij_[type_i].assign(particlestored, std::vector<double>(3, 0.0));
  }

  // get relevant particle pair indices for particle types
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndices(boundarytypes_, relindices);

  // iterate over relevant particle pairs
  for (const int particlepairindex : relindices)
  {
    const SPHParticlePair& particlepair =
        neighborpairs_->GetRefToParticlePairData()[particlepairindex];

    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlepair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = particlepair.tuple_j_;

    // check for boundary or rigid particles
    bool isboundaryrigid_i = boundarytypes_.count(type_i);
    bool isboundaryrigid_j = boundarytypes_.count(type_j);

    // no evaluation for both boundary or rigid particles
    if (isboundaryrigid_i and isboundaryrigid_j) continue;

    // evaluate contribution of neighboring particle j
    if (isboundaryrigid_i)
    {
      // get container of owned particles
      PARTICLEENGINE::ParticleContainer* container_j =
          particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

      // declare pointer variables for particle j
      const double *vel_j, *dens_j, *press_j;

      // get pointer to particle states
      vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
      dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
      press_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_j);

      // sum contribution of neighboring particle j
      sumj_Wij_[type_i][particle_i] += particlepair.Wij_;
      sumj_press_j_Wij_[type_i][particle_i] += press_j[0] * particlepair.Wij_;

      const double fac = dens_j[0] * particlepair.absdist_ * particlepair.Wij_;
      UTILS::vec_addscale(&sumj_dens_j_r_ij_Wij_[type_i][particle_i][0], fac, particlepair.e_ij_);

      UTILS::vec_addscale(&sumj_vel_j_Wij_[type_i][particle_i][0], particlepair.Wij_, vel_j);
    }

    // evaluate contribution of neighboring particle i
    if (isboundaryrigid_j and status_j == PARTICLEENGINE::Owned)
    {
      // get container of owned particles
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

      // declare pointer variables for particle i
      const double *vel_i, *dens_i, *press_i;

      // get pointer to particle states
      vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      press_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);

      // sum contribution of neighboring particle i
      sumj_Wij_[type_j][particle_j] += particlepair.Wji_;
      sumj_press_j_Wij_[type_j][particle_j] += press_i[0] * particlepair.Wji_;

      const double fac = -dens_i[0] * particlepair.absdist_ * particlepair.Wji_;
      UTILS::vec_addscale(&sumj_dens_j_r_ij_Wij_[type_j][particle_j][0], fac, particlepair.e_ij_);

      UTILS::vec_addscale(&sumj_vel_j_Wij_[type_j][particle_j][0], particlepair.Wji_, vel_i);
    }
  }

  // iterate over boundary particle types
  for (const auto& type_i : boundarytypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear modified boundary particle states
    container_i->ClearState(PARTICLEENGINE::BoundaryPressure);
    container_i->ClearState(PARTICLEENGINE::BoundaryVelocity);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // set modified boundary particle states
      if (sumj_Wij_[type_i][particle_i] > 0.0)
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
        UTILS::vec_set(relacc, &gravity[0]);
        UTILS::vec_sub(relacc, acc_i);

        const double inv_sumj_Wij = 1.0 / sumj_Wij_[type_i][particle_i];

        // set modified boundary pressure
        boundarypress_i[0] =
            (sumj_press_j_Wij_[type_i][particle_i] +
                UTILS::vec_dot(relacc, &sumj_dens_j_r_ij_Wij_[type_i][particle_i][0])) *
            inv_sumj_Wij;

        // set modified boundary velocity
        UTILS::vec_setscale(boundaryvel_i, 2.0, vel_i);
        UTILS::vec_addscale(boundaryvel_i, -inv_sumj_Wij, &sumj_vel_j_Wij_[type_i][particle_i][0]);
      }
    }
  }

  // refresh modified states of ghosted boundary particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(boundarystatestorefresh_);
}
