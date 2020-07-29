/*---------------------------------------------------------------------------*/
/*! \file
\brief density handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_density.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_virtual_wall_particle.H"
#include "particle_interaction_sph_density_correction.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityBase::SPHDensityBase(const Teuchos::ParameterList& params)
    : params_sph_(params), dt_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHDensityBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs,
    const std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // set kernel handler
  kernel_ = kernel;

  // set particle material handler
  particlematerial_ = particlematerial;

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // set virtual wall particle handler
  virtualwallparticle_ = virtualwallparticle;

  // setup density of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Density};

    for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

      if (container_i->HaveStoredState(PARTICLEENGINE::DensitySum) or
          container_i->HaveStoredState(PARTICLEENGINE::DensityDot))
        densitytorefresh_.push_back(std::make_pair(type_i, states));
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

void PARTICLEINTERACTION::SPHDensityBase::SumWeightedMass() const
{
  // clear density sum state
  ClearDensitySumState();

  // sum weighted mass (self contribution)
  SumWeightedMassSelfContribution();

  // sum weighted mass (particle contribution)
  SumWeightedMassParticleContribution();

  // sum weighted mass (particle-wall contribution)
  if (virtualwallparticle_) SumWeightedMassParticleWallContribution();
}

void PARTICLEINTERACTION::SPHDensityBase::ClearDensitySumState() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear density sum state
    if (container_i->HaveStoredState(PARTICLEENGINE::DensitySum))
      container_i->ClearState(PARTICLEENGINE::DensitySum);
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassSelfContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassSelfContribution");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    if (not container_i->HaveStoredState(PARTICLEENGINE::DensitySum)) continue;

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // declare pointer variables for particle i
      const double *rad_i, *mass_i;
      double* denssum_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      denssum_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensitySum, particle_i);

      // evaluate kernel
      const double Wii = kernel_->W0(rad_i[0]);

      // add self contribution
      denssum_i[0] += Wii * mass_i[0];
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassParticleContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassParticleContribution");

  // iterate over particle pairs
  for (auto& particlepair : neighborpairs_->GetRefToParticlePairData())
  {
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

    // declare pointer variables for particle i and j
    const double* mass_i;
    double* denssum_i = nullptr;

    const double* mass_j;
    double* denssum_j = nullptr;

    // get pointer to particle states
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    if (container_i->HaveStoredState(PARTICLEENGINE::DensitySum))
      denssum_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensitySum, particle_i);

    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    if (container_j->HaveStoredState(PARTICLEENGINE::DensitySum))
      denssum_j = container_j->GetPtrToParticleState(PARTICLEENGINE::DensitySum, particle_j);

    // sum contribution of neighboring particle j
    if (denssum_i) denssum_i[0] += particlepair.Wij_ * mass_i[0];

    // sum contribution of neighboring particle i
    if (denssum_j and status_j == PARTICLEENGINE::Owned)
      denssum_j[0] += particlepair.Wji_ * mass_j[0];
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassParticleWallContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassParticleWallContribution");

  // iterate over particle-wall pairs
  for (auto& particlewallpair : neighborpairs_->GetRefToParticleWallPairData())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // no evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // no evaluation for open boundary particles
    if (type_i == PARTICLEENGINE::DirichletPhase or type_i == PARTICLEENGINE::NeumannPhase)
      continue;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // declare pointer variables for particle i
    const double *rad_i, *mass_i;
    double* denssum_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    denssum_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensitySum, particle_i);

    // compute vector from wall contact point j to particle i
    double r_ij[3];
    UTILS::vec_setscale(r_ij, particlewallpair.absdist_, particlewallpair.e_ij_);

    // unit surface tangent vectors in wall contact point j
    double t_j_1[3];
    double t_j_2[3];
    UTILS::unitsurfacetangents(particlewallpair.e_ij_, t_j_1, t_j_2);

    // iterate over virtual particles
    for (const std::vector<double>& virtualparticle :
        virtualwallparticle_->GetRelativePositionsOfVirtualParticles())
    {
      // vector from virtual particle k to particle i
      double r_ik[3];
      UTILS::vec_set(r_ik, r_ij);
      UTILS::vec_addscale(r_ik, virtualparticle[0], particlewallpair.e_ij_);
      UTILS::vec_addscale(r_ik, virtualparticle[1], t_j_1);
      UTILS::vec_addscale(r_ik, virtualparticle[2], t_j_2);

      // absolute distance between virtual particle k and particle i
      const double absdist = UTILS::vec_norm2(r_ik);

      // virtual particle within interaction distance
      if (absdist < rad_i[0])
      {
        // evaluate kernel
        const double Wik = kernel_->W(absdist, rad_i[0]);

        // sum contribution of virtual particle k
        denssum_i[0] += Wik * mass_i[0];
      }
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumColorfield() const
{
  // clear colorfield state
  ClearColorfieldState();

  // sum colorfield (self contribution)
  SumColorfieldSelfContribution();

  // sum colorfield (particle contribution)
  SumColorfieldParticleContribution();

  // sum colorfield (particle-wall contribution)
  if (virtualwallparticle_) SumColorfieldParticleWallContribution();
}

void PARTICLEINTERACTION::SPHDensityBase::ClearColorfieldState() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear colorfield state
    if (container_i->HaveStoredState(PARTICLEENGINE::Colorfield))
      container_i->ClearState(PARTICLEENGINE::Colorfield);
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumColorfieldSelfContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityBase::SumColorfieldSelfContribution");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    if (not container_i->HaveStoredState(PARTICLEENGINE::Colorfield)) continue;

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *dens_i;
      double* colorfield_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      colorfield_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Colorfield, particle_i);

      // evaluate kernel
      const double Wii = kernel_->W0(rad_i[0]);

      // add self contribution
      colorfield_i[0] += (Wii / dens_i[0]) * mass_i[0];
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumColorfieldParticleContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::SumColorfieldParticleContribution");

  // iterate over particle pairs
  for (auto& particlepair : neighborpairs_->GetRefToParticlePairData())
  {
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

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    const MAT::PAR::ParticleMaterialBase* material_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    // declare pointer variables for particle i and j
    const double *mass_i, *dens_i;
    double* colorfield_i = nullptr;

    const double *mass_j, *dens_j;
    double* colorfield_j = nullptr;

    // get pointer to particle states
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    dens_i = (container_i->HaveStoredState(PARTICLEENGINE::Density))
                 ? container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i)
                 : &(material_j->initDensity_);

    if (container_i->HaveStoredState(PARTICLEENGINE::Colorfield))
      colorfield_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Colorfield, particle_i);

    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    dens_j = (container_j->HaveStoredState(PARTICLEENGINE::Density))
                 ? container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j)
                 : &(material_i->initDensity_);

    if (container_j->HaveStoredState(PARTICLEENGINE::Colorfield))
      colorfield_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Colorfield, particle_j);

    // sum contribution of neighboring particle j
    if (colorfield_i) colorfield_i[0] += (particlepair.Wij_ / dens_j[0]) * mass_j[0];

    // sum contribution of neighboring particle i
    if (colorfield_j and status_j == PARTICLEENGINE::Owned)
      colorfield_j[0] += (particlepair.Wji_ / dens_i[0]) * mass_i[0];
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SumColorfieldParticleWallContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::SumColorfieldParticleWallContribution");

  // iterate over particle-wall pairs
  for (auto& particlewallpair : neighborpairs_->GetRefToParticleWallPairData())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // no evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // no evaluation for open boundary particles
    if (type_i == PARTICLEENGINE::DirichletPhase or type_i == PARTICLEENGINE::NeumannPhase)
      continue;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // declare pointer variables for particle i
    const double* rad_i;
    double* colorfield_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    colorfield_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Colorfield, particle_i);

    // declare pointer variables for virtual particle k
    const double *mass_k, *dens_k;

    // get pointer to virtual particle states
    mass_k = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    dens_k = &(material_i->initDensity_);

    // (current) volume of virtual particle k
    const double V_k = mass_k[0] / dens_k[0];

    // compute vector from wall contact point j to particle i
    double r_ij[3];
    UTILS::vec_setscale(r_ij, particlewallpair.absdist_, particlewallpair.e_ij_);

    // unit surface tangent vectors in wall contact point j
    double t_j_1[3];
    double t_j_2[3];
    UTILS::unitsurfacetangents(particlewallpair.e_ij_, t_j_1, t_j_2);

    // iterate over virtual particles
    for (const std::vector<double>& virtualparticle :
        virtualwallparticle_->GetRelativePositionsOfVirtualParticles())
    {
      // vector from virtual particle k to particle i
      double r_ik[3];
      UTILS::vec_set(r_ik, r_ij);
      UTILS::vec_addscale(r_ik, virtualparticle[0], particlewallpair.e_ij_);
      UTILS::vec_addscale(r_ik, virtualparticle[1], t_j_1);
      UTILS::vec_addscale(r_ik, virtualparticle[2], t_j_2);

      // absolute distance between virtual particle k and particle i
      const double absdist = UTILS::vec_norm2(r_ik);

      // virtual particle within interaction distance
      if (absdist < rad_i[0])
      {
        // evaluate kernel
        const double Wik = kernel_->W(absdist, rad_i[0]);

        // sum contribution of virtual particle k
        colorfield_i[0] += V_k * Wik;
      }
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::ContinuityEquation() const
{
  // clear density dot state
  ClearDensityDotState();

  // continuity equation (particle contribution)
  ContinuityEquationParticleContribution();

  // continuity equation (particle-wall contribution)
  if (virtualwallparticle_) ContinuityEquationParticleWallContribution();
}

void PARTICLEINTERACTION::SPHDensityBase::ClearDensityDotState() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear density dot state
    if (container_i->HaveStoredState(PARTICLEENGINE::DensityDot))
      container_i->ClearState(PARTICLEENGINE::DensityDot);
  }
}

void PARTICLEINTERACTION::SPHDensityBase::ContinuityEquationParticleContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::ContinuityEquationParticleContribution");

  // iterate over particle pairs
  for (auto& particlepair : neighborpairs_->GetRefToParticlePairData())
  {
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

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    const MAT::PAR::ParticleMaterialBase* material_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    // declare pointer variables for particle i and j
    const double *vel_i, *mass_i, *dens_i;
    double* densdot_i = nullptr;

    const double *vel_j, *mass_j, *dens_j;
    double* densdot_j = nullptr;

    // get pointer to particle states
    vel_i = (container_i->HaveStoredState(PARTICLEENGINE::ModifiedVelocity))
                ? container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_i)
                : container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    dens_i = (container_i->HaveStoredState(PARTICLEENGINE::Density))
                 ? container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i)
                 : &(material_j->initDensity_);

    if (container_i->HaveStoredState(PARTICLEENGINE::DensityDot))
      densdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensityDot, particle_i);

    vel_j = (container_j->HaveStoredState(PARTICLEENGINE::ModifiedVelocity))
                ? container_j->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_j)
                : container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    dens_j = (container_j->HaveStoredState(PARTICLEENGINE::Density))
                 ? container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j)
                 : &(material_i->initDensity_);

    if (container_j->HaveStoredState(PARTICLEENGINE::DensityDot))
      densdot_j = container_j->GetPtrToParticleState(PARTICLEENGINE::DensityDot, particle_j);

    // relative velocity (use modified velocities in case of transport velocity formulation)
    double vel_ij[3];
    UTILS::vec_set(vel_ij, vel_i);
    UTILS::vec_sub(vel_ij, vel_j);

    const double e_ij_vel_ij = UTILS::vec_dot(particlepair.e_ij_, vel_ij);

    // sum contribution of neighboring particle j
    if (densdot_i)
      densdot_i[0] += dens_i[0] * (mass_j[0] / dens_j[0]) * particlepair.dWdrij_ * e_ij_vel_ij;

    // sum contribution of neighboring particle i
    if (densdot_j and status_j == PARTICLEENGINE::Owned)
      densdot_j[0] += dens_j[0] * (mass_i[0] / dens_i[0]) * particlepair.dWdrji_ * e_ij_vel_ij;
  }
}

void PARTICLEINTERACTION::SPHDensityBase::ContinuityEquationParticleWallContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHDensityBase::ContinuityEquationParticleWallContribution");

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // iterate over particle-wall pairs
  for (auto& particlewallpair : neighborpairs_->GetRefToParticleWallPairData())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // no evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // no evaluation for open boundary particles
    if (type_i == PARTICLEENGINE::DirichletPhase or type_i == PARTICLEENGINE::NeumannPhase)
      continue;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // declare pointer variables for particle i
    const double *vel_i, *rad_i, *dens_i;
    double* densdot_i;

    // get pointer to particle states
    vel_i = (container_i->HaveStoredState(PARTICLEENGINE::ModifiedVelocity))
                ? container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_i)
                : container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
    densdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensityDot, particle_i);

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // number of nodes of wall element
    const int numnodes = ele->NumNode();

    // shape functions and location vector of wall element
    Epetra_SerialDenseVector funct(numnodes);
    std::vector<int> lmele;

    if (walldatastate->GetVelCol() != Teuchos::null)
    {
      // evaluate shape functions of element at wall contact point
      DRT::UTILS::shape_function_2D(
          funct, particlewallpair.elecoords_[0], particlewallpair.elecoords_[1], ele->Shape());

      // get location vector of wall element
      lmele.reserve(numnodes * 3);
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      ele->LocationVector(
          *particlewallinterface_->GetWallDiscretization(), lmele, lmowner, lmstride);
    }

    // velocity of wall contact point j
    double vel_j[3] = {0.0};

    if (walldatastate->GetVelCol() != Teuchos::null)
    {
      // get nodal velocities
      std::vector<double> nodal_vel(numnodes * 3);
      DRT::UTILS::ExtractMyValues(*walldatastate->GetVelCol(), nodal_vel, lmele);

      // determine velocity of wall contact point j
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim) vel_j[dim] += funct[node] * nodal_vel[node * 3 + dim];
    }

    // declare pointer variables for virtual particle k
    const double *vel_k, *mass_k, *dens_k;

    // get pointer to virtual particle states
    mass_k = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    dens_k = &(material_i->initDensity_);
    vel_k = &vel_j[0];

    // (current) volume of virtual particle k
    const double V_k = mass_k[0] / dens_k[0];

    // compute vector from wall contact point j to particle i
    double r_ij[3];
    UTILS::vec_setscale(r_ij, particlewallpair.absdist_, particlewallpair.e_ij_);

    // relative velocity (use modified velocities in case of transport velocity formulation)
    double vel_ik[3];
    UTILS::vec_set(vel_ik, vel_i);
    UTILS::vec_sub(vel_ik, vel_k);

    // unit surface tangent vectors in wall contact point j
    double t_j_1[3];
    double t_j_2[3];
    UTILS::unitsurfacetangents(particlewallpair.e_ij_, t_j_1, t_j_2);

    // iterate over virtual particles
    for (const std::vector<double>& virtualparticle :
        virtualwallparticle_->GetRelativePositionsOfVirtualParticles())
    {
      // vector from virtual particle k to particle i
      double r_ik[3];
      UTILS::vec_set(r_ik, r_ij);
      UTILS::vec_addscale(r_ik, virtualparticle[0], particlewallpair.e_ij_);
      UTILS::vec_addscale(r_ik, virtualparticle[1], t_j_1);
      UTILS::vec_addscale(r_ik, virtualparticle[2], t_j_2);

      // absolute distance between virtual particle k and particle i
      const double absdist = UTILS::vec_norm2(r_ik);

      // virtual particle within interaction distance
      if (absdist < rad_i[0])
      {
        const double e_ik_vel_ik = UTILS::vec_dot(r_ik, vel_ik) / absdist;

        // evaluate first derivative of kernel
        const double dWdrik = kernel_->dWdrij(absdist, rad_i[0]);

        // sum contribution of virtual particle k
        densdot_i[0] += dens_i[0] * V_k * dWdrik * e_ik_vel_ik;
      }
    }
  }
}

void PARTICLEINTERACTION::SPHDensityBase::SetDensitySum() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // update density of all particles
    if (container_i->HaveStoredState(PARTICLEENGINE::DensitySum))
      container_i->UpdateState(0.0, PARTICLEENGINE::Density, 1.0, PARTICLEENGINE::DensitySum);
  }
}

void PARTICLEINTERACTION::SPHDensityBase::AddTimeStepScaledDensityDot() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // update density of all particles
    if (container_i->HaveStoredState(PARTICLEENGINE::DensityDot))
      container_i->UpdateState(1.0, PARTICLEENGINE::Density, dt_, PARTICLEENGINE::DensityDot);
  }
}

PARTICLEINTERACTION::SPHDensitySummation::SPHDensitySummation(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHDensitySummation::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // no states for open boundary particles
    if (type == PARTICLEENGINE::DirichletPhase or type == PARTICLEENGINE::NeumannPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(PARTICLEENGINE::DensitySum);
  }
}

void PARTICLEINTERACTION::SPHDensitySummation::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensitySummation::ComputeDensity");

  // evaluate sum of weighted mass
  SumWeightedMass();

  // set density sum to density field
  SetDensitySum();

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

PARTICLEINTERACTION::SPHDensityIntegration::SPHDensityIntegration(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHDensityIntegration::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // no states for open boundary particles
    if (type == PARTICLEENGINE::DirichletPhase or type == PARTICLEENGINE::NeumannPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(PARTICLEENGINE::DensityDot);
  }
}

void PARTICLEINTERACTION::SPHDensityIntegration::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityIntegration::ComputeDensity");

  // evaluate continuity equation
  ContinuityEquation();

  // add time step scaled density dot to density field
  AddTimeStepScaledDensityDot();

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

PARTICLEINTERACTION::SPHDensityPredictCorrect::SPHDensityPredictCorrect(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

PARTICLEINTERACTION::SPHDensityPredictCorrect::~SPHDensityPredictCorrect() = default;

void PARTICLEINTERACTION::SPHDensityPredictCorrect::Init()
{
  // call base class init
  SPHDensityBase::Init();

  // init density correction handler
  InitDensityCorrectionHandler();
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs,
    const std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle)
{
  // call base class setup
  SPHDensityBase::Setup(particleengineinterface, particlewallinterface, kernel, particlematerial,
      equationofstatebundle, neighborpairs, virtualwallparticle);

  // setup density correction handler
  densitycorrection_->Setup();
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::WriteRestart(
    const int step, const double time) const
{
  // call base class function
  SPHDensityBase::WriteRestart(step, time);

  // write restart of density correction handler
  densitycorrection_->WriteRestart(step, time);
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  SPHDensityBase::ReadRestart(reader);

  // read restart of density correction handler
  densitycorrection_->ReadRestart(reader);
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // no states for open boundary particles
    if (type == PARTICLEENGINE::DirichletPhase or type == PARTICLEENGINE::NeumannPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(
        {PARTICLEENGINE::DensityDot, PARTICLEENGINE::DensitySum, PARTICLEENGINE::Colorfield});
  }
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityPredictCorrect::ComputeDensity");

  // evaluate continuity equation
  ContinuityEquation();

  // add time step scaled density dot to density field
  AddTimeStepScaledDensityDot();

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);

  // evaluate sum of weighted mass
  SumWeightedMass();

  // evaluate sum of colorfield
  SumColorfield();

  // correct density of interior/surface particles
  CorrectDensity();

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::InitDensityCorrectionHandler()
{
  // get type of density correction scheme
  INPAR::PARTICLE::DensityCorrectionScheme densitycorrectionscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityCorrectionScheme>(
          params_sph_, "DENSITYCORRECTION");

  // create density correction handler
  switch (densitycorrectionscheme)
  {
    case INPAR::PARTICLE::InteriorCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionInterior>(
          new PARTICLEINTERACTION::SPHDensityCorrectionInterior());
      break;
    }
    case INPAR::PARTICLE::NormalizedCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionNormalized>(
          new PARTICLEINTERACTION::SPHDensityCorrectionNormalized());
      break;
    }
    case INPAR::PARTICLE::RandlesCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionRandles>(
          new PARTICLEINTERACTION::SPHDensityCorrectionRandles());
      break;
    }
    default:
    {
      dserror("no density correction scheme set via parameter 'DENSITYCORRECTION'!");
      break;
    }
  }

  // init density correction handler
  densitycorrection_->Init();
}

void PARTICLEINTERACTION::SPHDensityPredictCorrect::CorrectDensity() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    if (not container_i->HaveStoredState(PARTICLEENGINE::DensitySum) or
        not container_i->HaveStoredState(PARTICLEENGINE::Colorfield))
      continue;

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // declare pointer variables for particle
    const double *denssum, *colorfield;
    double* dens;

    // get pointer to particle state
    dens = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, 0);
    denssum = container_i->GetPtrToParticleState(PARTICLEENGINE::DensitySum, 0);
    colorfield = container_i->GetPtrToParticleState(PARTICLEENGINE::Colorfield, 0);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // get equation of state for current particle type
    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate =
        equationofstatebundle_->GetPtrToSpecificEquationOfState(type_i);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      if (colorfield[i] >= 1.0)
      {
        // set corrected density of interior particles
        densitycorrection_->CorrectedDensityInterior(&denssum[i], &dens[i]);
      }
      else
      {
        double dens_bc = 0.0;
        if (densitycorrection_->ComputeDensityBC())
        {
          double press_bc = 0.0;
          dens_bc = equationofstate->PressureToDensity(press_bc, material->initDensity_);
        }

        // set corrected density of free surface particles
        densitycorrection_->CorrectedDensityFreeSurface(
            &denssum[i], &colorfield[i], &dens_bc, &dens[i]);
      }
    }
  }
}
