/*---------------------------------------------------------------------------*/
/*! \file
\brief barrier force handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_surface_tension_barrier_force.H"

#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHBarrierForce::SPHBarrierForce(const Teuchos::ParameterList& params)
    : params_sph_(params),
      liquidtype_(PARTICLEENGINE::Phase1),
      gastype_(PARTICLEENGINE::Phase2),
      dist_(params_sph_.get<double>("BARRIER_FORCE_DISTANCE")),
      stiff_(params_sph_.get<double>("BARRIER_FORCE_STIFF")),
      damp_(params_sph_.get<double>("BARRIER_FORCE_DAMP"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHBarrierForce::Init()
{
  // init fluid particle types
  fluidtypes_ = {liquidtype_, gastype_};

  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};
}

void PARTICLEINTERACTION::SPHBarrierForce::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // safety check
  for (const auto& type_i : fluidtypes_)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());

  // update with actual boundary particle types
  const auto boundarytypes = boundarytypes_;
  for (const auto& type_i : boundarytypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      boundarytypes_.erase(type_i);

  // safety check
  if (not(dist_ > 0.0)) dserror("barrier force distance not positive!");
  if (not(stiff_ > 0.0)) dserror("barrier force stiffness not positive!");
  if (damp_ < 0.0) dserror("barrier force damping parameter not positive or zero!");
}

void PARTICLEINTERACTION::SPHBarrierForce::ComputeBarrierForceContribution() const
{
  // compute barrier force contribution (particle contribution)
  ComputeBarrierForceParticleContribution();

  // compute barrier force contribution (particle-boundary contribution)
  ComputeBarrierForceParticleBoundaryContribution();
}

void PARTICLEINTERACTION::SPHBarrierForce::ComputeBarrierForceParticleContribution() const
{
  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(fluidtypes_, relindices);

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

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    double* acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

    const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    const double* vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
    double* acc_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_j);

    if (particlepair.absdist_ < dist_)
    {
      const double gap = particlepair.absdist_ - dist_;
      const double gapdot =
          UTILS::vec_dot(vel_i, particlepair.e_ij_) - UTILS::vec_dot(vel_j, particlepair.e_ij_);

      // magnitude of barrier force
      const double fac = (stiff_ * gap + damp_ * std::abs(gap) * gapdot);

      // sum contribution of neighboring particle j
      UTILS::vec_addscale(acc_i, -fac / mass_i[0], particlepair.e_ij_);

      // sum contribution of neighboring particle i
      if (status_j == PARTICLEENGINE::Owned)
        UTILS::vec_addscale(acc_j, fac / mass_j[0], particlepair.e_ij_);
    }
  }
}

void PARTICLEINTERACTION::SPHBarrierForce::ComputeBarrierForceParticleBoundaryContribution() const
{
  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForDisjointCombination(
      fluidtypes_, boundarytypes_, relindices);

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

    // swap fluid particle and boundary particle
    const bool swapparticles = boundarytypes_.count(type_i);
    if (swapparticles)
    {
      std::tie(type_i, status_i, particle_i) = particlepair.tuple_j_;
      std::tie(type_j, status_j, particle_j) = particlepair.tuple_i_;
    }

    // absolute distance between particles
    const double absdist = particlepair.absdist_;

    // versor from particle j to i
    double e_ij[3];
    UTILS::vec_set(e_ij, particlepair.e_ij_);
    if (swapparticles) UTILS::vec_scale(e_ij, -1.0);

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

    double* acc_i = nullptr;
    if (status_i == PARTICLEENGINE::Owned)
      acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

    // get pointer to boundary particle states
    const double* vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

    if (absdist < dist_)
    {
      const double gap = absdist - dist_;
      const double gapdot = UTILS::vec_dot(vel_i, e_ij) - UTILS::vec_dot(vel_j, e_ij);

      // magnitude of barrier force
      const double fac = (stiff_ * gap + damp_ * std::abs(gap) * gapdot);

      // sum contribution of neighboring particle j
      if (acc_i) UTILS::vec_addscale(acc_i, -fac / mass_i[0], e_ij);
    }
  }
}
