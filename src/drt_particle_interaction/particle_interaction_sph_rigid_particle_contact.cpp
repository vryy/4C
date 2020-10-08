/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid particle contact handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_rigid_particle_contact.H"

#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHRigidParticleContactBase::SPHRigidParticleContactBase(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHRigidParticleContactBase::Init()
{
  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};
}

void PARTICLEINTERACTION::SPHRigidParticleContactBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // update with actual boundary particle types
  for (const auto& type_i : boundarytypes_)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      boundarytypes_.erase(type_i);

  // safety check
  if (not boundarytypes_.count(PARTICLEENGINE::RigidPhase))
    dserror("no rigid particles defined but a rigid particle contact formulation is set!");
}

PARTICLEINTERACTION::SPHRigidParticleContactElastic::SPHRigidParticleContactElastic(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHRigidParticleContactBase(params),
      stiff_(params_sph_.get<double>("RIGIDPARTICLECONTACTSTIFF")),
      damp_(params_sph_.get<double>("RIGIDPARTICLECONTACTDAMP"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHRigidParticleContactElastic::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHRigidParticleContactBase::Setup(particleengineinterface, neighborpairs);

  // safety check
  if (not(stiff_ > 0.0)) dserror("rigid particle contact stiffness not positive!");
  if (damp_ < 0.0) dserror("rigid particle contact damping parameter not positive or zero!");
}

void PARTICLEINTERACTION::SPHRigidParticleContactElastic::AddForceContribution()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHRigidParticleContactElastic::AddForceContribution");

  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(boundarytypes_, relindices);

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
    const double* vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

    double* force_i = nullptr;
    if (container_i->HaveStoredState(PARTICLEENGINE::Force))
      force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    // get pointer to particle states
    const double* vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

    double* force_j = nullptr;
    if (container_j->HaveStoredState(PARTICLEENGINE::Force) and status_j == PARTICLEENGINE::Owned)
      force_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Force, particle_j);

    if (particlepair.absdist_ < initialparticlespacing)
    {
      const double gap = initialparticlespacing - particlepair.absdist_;
      const double gapdot =
          UTILS::vec_dot(vel_j, particlepair.e_ij_) - UTILS::vec_dot(vel_i, particlepair.e_ij_);

      // magnitude of rigid particle contact force
      const double fac = stiff_ * gap + damp_ * gapdot;

      // add contributions
      if (force_i) UTILS::vec_addscale(force_i, fac, particlepair.e_ij_);
      if (force_j) UTILS::vec_addscale(force_j, -fac, particlepair.e_ij_);
    }
  }
}
