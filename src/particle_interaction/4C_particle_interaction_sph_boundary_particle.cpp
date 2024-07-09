/*---------------------------------------------------------------------------*/
/*! \file
\brief boundary particle handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph_boundary_particle.hpp"

#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_interaction_sph_neighbor_pairs.hpp"
#include "4C_particle_interaction_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHBoundaryParticleBase::SPHBoundaryParticleBase(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

void ParticleInteraction::SPHBoundaryParticleBase::init()
{
  // init with potential fluid particle types
  fluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::DirichletPhase,
      PARTICLEENGINE::NeumannPhase};

  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};
}

void ParticleInteraction::SPHBoundaryParticleBase::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<ParticleInteraction::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->get_particle_container_bundle();

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // update with actual fluid particle types
  const auto fluidtypes = fluidtypes_;
  for (const auto& type_i : fluidtypes)
    if (not particlecontainerbundle_->get_particle_types().count(type_i)) fluidtypes_.erase(type_i);

  // update with actual boundary particle types
  const auto boundarytypes = boundarytypes_;
  for (const auto& type_i : boundarytypes)
    if (not particlecontainerbundle_->get_particle_types().count(type_i))
      boundarytypes_.erase(type_i);

  // safety check
  if (boundarytypes_.empty())
    FOUR_C_THROW(
        "no boundary or rigid particles defined but a boundary particle formulation is set!");
}

ParticleInteraction::SPHBoundaryParticleAdami::SPHBoundaryParticleAdami(
    const Teuchos::ParameterList& params)
    : ParticleInteraction::SPHBoundaryParticleBase(params)
{
  // empty constructor
}

void ParticleInteraction::SPHBoundaryParticleAdami::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<ParticleInteraction::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHBoundaryParticleBase::setup(particleengineinterface, neighborpairs);

  // setup modified states of ghosted boundary particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{
        PARTICLEENGINE::BoundaryPressure, PARTICLEENGINE::BoundaryVelocity};

    for (const auto& type_i : boundarytypes_)
      boundarystatestorefresh_.push_back(std::make_pair(type_i, states));
  }

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--boundarytypes_.end()) + 1;

  // allocate memory to hold contributions of neighboring particles
  sumj_wij_.resize(typevectorsize);
  sumj_press_j_wij_.resize(typevectorsize);
  sumj_dens_j_r_ij_wij_.resize(typevectorsize);
  sumj_vel_j_wij_.resize(typevectorsize);
}

void ParticleInteraction::SPHBoundaryParticleAdami::init_boundary_particle_states(
    std::vector<double>& gravity)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "ParticleInteraction::SPHBoundaryParticleAdami::init_boundary_particle_states");

  // iterate over boundary particle types
  for (const auto& type_i : boundarytypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->particles_stored();

    // allocate memory
    sumj_wij_[type_i].assign(particlestored, 0.0);
    sumj_press_j_wij_[type_i].assign(particlestored, 0.0);
    sumj_dens_j_r_ij_wij_[type_i].assign(particlestored, std::vector<double>(3, 0.0));
    sumj_vel_j_wij_[type_i].assign(particlestored, std::vector<double>(3, 0.0));
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->get_relevant_particle_pair_indices_for_disjoint_combination(
      boundarytypes_, fluidtypes_, relindices);

  // iterate over relevant particle pairs
  for (const int particlepairindex : relindices)
  {
    const SPHParticlePair& particlepair =
        neighborpairs_->get_ref_to_particle_pair_data()[particlepairindex];

    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlepair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = particlepair.tuple_j_;

    // evaluate contribution of neighboring fluid particle j
    if (boundarytypes_.count(type_i))
    {
      // get container of owned particles
      PARTICLEENGINE::ParticleContainer* container_j =
          particlecontainerbundle_->get_specific_container(type_j, status_j);

      // get pointer to particle states
      const double* vel_j = container_j->get_ptr_to_state(PARTICLEENGINE::Velocity, particle_j);
      const double* dens_j = container_j->get_ptr_to_state(PARTICLEENGINE::Density, particle_j);
      const double* press_j = container_j->get_ptr_to_state(PARTICLEENGINE::Pressure, particle_j);

      // sum contribution of neighboring particle j
      sumj_wij_[type_i][particle_i] += particlepair.Wij_;
      sumj_press_j_wij_[type_i][particle_i] += press_j[0] * particlepair.Wij_;

      const double fac = dens_j[0] * particlepair.absdist_ * particlepair.Wij_;
      UTILS::VecAddScale(sumj_dens_j_r_ij_wij_[type_i][particle_i].data(), fac, particlepair.e_ij_);

      UTILS::VecAddScale(sumj_vel_j_wij_[type_i][particle_i].data(), particlepair.Wij_, vel_j);
    }

    // evaluate contribution of neighboring fluid particle i
    if (boundarytypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // get container of owned particles
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->get_specific_container(type_i, status_i);

      // get pointer to particle states
      const double* vel_i = container_i->get_ptr_to_state(PARTICLEENGINE::Velocity, particle_i);
      const double* dens_i = container_i->get_ptr_to_state(PARTICLEENGINE::Density, particle_i);
      const double* press_i = container_i->get_ptr_to_state(PARTICLEENGINE::Pressure, particle_i);

      // sum contribution of neighboring particle i
      sumj_wij_[type_j][particle_j] += particlepair.Wji_;
      sumj_press_j_wij_[type_j][particle_j] += press_i[0] * particlepair.Wji_;

      const double fac = -dens_i[0] * particlepair.absdist_ * particlepair.Wji_;
      UTILS::VecAddScale(sumj_dens_j_r_ij_wij_[type_j][particle_j].data(), fac, particlepair.e_ij_);

      UTILS::VecAddScale(sumj_vel_j_wij_[type_j][particle_j].data(), particlepair.Wji_, vel_i);
    }
  }

  // iterate over boundary particle types
  for (const auto& type_i : boundarytypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // clear modified boundary particle states
    container_i->clear_state(PARTICLEENGINE::BoundaryPressure);
    container_i->clear_state(PARTICLEENGINE::BoundaryVelocity);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->particles_stored(); ++particle_i)
    {
      // set modified boundary particle states
      if (sumj_wij_[type_i][particle_i] > 0.0)
      {
        // get pointer to particle states
        const double* vel_i = container_i->get_ptr_to_state(PARTICLEENGINE::Velocity, particle_i);
        const double* acc_i =
            container_i->get_ptr_to_state(PARTICLEENGINE::Acceleration, particle_i);
        double* boundarypress_i =
            container_i->get_ptr_to_state(PARTICLEENGINE::BoundaryPressure, particle_i);
        double* boundaryvel_i =
            container_i->get_ptr_to_state(PARTICLEENGINE::BoundaryVelocity, particle_i);

        // get relative acceleration of boundary particle
        double relacc[3];
        UTILS::VecSet(relacc, gravity.data());
        UTILS::VecSub(relacc, acc_i);

        const double inv_sumj_Wij = 1.0 / sumj_wij_[type_i][particle_i];

        // set modified boundary pressure
        boundarypress_i[0] =
            (sumj_press_j_wij_[type_i][particle_i] +
                UTILS::VecDot(relacc, sumj_dens_j_r_ij_wij_[type_i][particle_i].data())) *
            inv_sumj_Wij;

        // set modified boundary velocity
        UTILS::VecSetScale(boundaryvel_i, 2.0, vel_i);
        UTILS::VecAddScale(
            boundaryvel_i, -inv_sumj_Wij, sumj_vel_j_wij_[type_i][particle_i].data());
      }
    }
  }

  // refresh modified states of ghosted boundary particles
  particleengineinterface_->refresh_particles_of_specific_states_and_types(
      boundarystatestorefresh_);
}

FOUR_C_NAMESPACE_CLOSE
