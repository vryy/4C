/*---------------------------------------------------------------------------*/
/*! \file
\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_sph_momentum.H"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_io_visualization_manager.H"
#include "baci_lib_element.H"
#include "baci_lib_utils.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_interface.H"
#include "baci_particle_interaction_material_handler.H"
#include "baci_particle_interaction_runtime_writer.H"
#include "baci_particle_interaction_sph_artificialviscosity.H"
#include "baci_particle_interaction_sph_equationofstate.H"
#include "baci_particle_interaction_sph_equationofstate_bundle.H"
#include "baci_particle_interaction_sph_kernel.H"
#include "baci_particle_interaction_sph_momentum_formulation.H"
#include "baci_particle_interaction_sph_neighbor_pairs.H"
#include "baci_particle_interaction_sph_virtual_wall_particle.H"
#include "baci_particle_interaction_utils.H"
#include "baci_particle_wall_datastate.H"
#include "baci_particle_wall_interface.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentum::SPHMomentum(const Teuchos::ParameterList& params)
    : params_sph_(params),
      boundaryparticleinteraction_(
          DRT::INPUT::IntegralValue<INPAR::PARTICLE::BoundaryParticleInteraction>(
              params_sph_, "BOUNDARYPARTICLEINTERACTION")),
      transportvelocityformulation_(
          DRT::INPUT::IntegralValue<INPAR::PARTICLE::TransportVelocityFormulation>(
              params_sph_, "TRANSPORTVELOCITYFORMULATION")),
      writeparticlewallinteraction_(
          DRT::INPUT::IntegralValue<int>(params_sph_, "WRITE_PARTICLE_WALL_INTERACTION"))
{
  // empty constructor
}

PARTICLEINTERACTION::SPHMomentum::~SPHMomentum() = default;

void PARTICLEINTERACTION::SPHMomentum::Init()
{
  // init momentum formulation handler
  InitMomentumFormulationHandler();

  // init artificial viscosity handler
  InitArtificialViscosityHandler();

  // init with potential fluid particle types
  allfluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::DirichletPhase,
      PARTICLEENGINE::NeumannPhase};
  intfluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::NeumannPhase};
  purefluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2};

  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};
}

void PARTICLEINTERACTION::SPHMomentum::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
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

  // set particle interaction writer
  particleinteractionwriter_ = particleinteractionwriter;

  // setup particle interaction writer
  SetupParticleInteractionWriter();

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // set virtual wall particle handler
  virtualwallparticle_ = virtualwallparticle;

  // setup momentum formulation handler
  momentumformulation_->Setup();

  // setup artificial viscosity handler
  artificialviscosity_->Setup();

  // update with actual fluid particle types
  const auto allfluidtypes = allfluidtypes_;
  for (const auto& type_i : allfluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      allfluidtypes_.erase(type_i);

  const auto intfluidtypes = intfluidtypes_;
  for (const auto& type_i : intfluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      intfluidtypes_.erase(type_i);

  const auto purefluidtypes = purefluidtypes_;
  for (const auto& type_i : purefluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      purefluidtypes_.erase(type_i);

  // update with actual boundary particle types
  const auto boundarytypes = boundarytypes_;
  for (const auto& type_i : boundarytypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      boundarytypes_.erase(type_i);

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  fluidmaterial_.resize(typevectorsize);

  // iterate over all fluid particle types
  for (const auto& type_i : allfluidtypes_)
  {
    fluidmaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));
  }
}

void PARTICLEINTERACTION::SPHMomentum::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // current particle type is not a pure fluid particle type
    if (not purefluidtypes_.count(type_i)) continue;

    // additional states for transport velocity formulation
    if (transportvelocityformulation_ !=
        INPAR::PARTICLE::TransportVelocityFormulation::NoTransportVelocity)
      particlestates.insert(
          {PARTICLEENGINE::ModifiedVelocity, PARTICLEENGINE::ModifiedAcceleration});
  }
}

void PARTICLEINTERACTION::SPHMomentum::AddAccelerationContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHMomentum::AddAccelerationContribution");

  // momentum equation (particle contribution)
  MomentumEquationParticleContribution();

  // momentum equation (particle-boundary contribution)
  MomentumEquationParticleBoundaryContribution();

  // momentum equation (particle-wall contribution)
  if (virtualwallparticle_) MomentumEquationParticleWallContribution();
}

void PARTICLEINTERACTION::SPHMomentum::InitMomentumFormulationHandler()
{
  // get type of smoothed particle hydrodynamics momentum formulation
  INPAR::PARTICLE::MomentumFormulationType momentumformulationtype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::MomentumFormulationType>(
          params_sph_, "MOMENTUMFORMULATION");

  // create momentum formulation handler
  switch (momentumformulationtype)
  {
    case INPAR::PARTICLE::AdamiMomentumFormulation:
    {
      momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationAdami>(
          new PARTICLEINTERACTION::SPHMomentumFormulationAdami());
      break;
    }
    case INPAR::PARTICLE::MonaghanMomentumFormulation:
    {
      momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationMonaghan>(
          new PARTICLEINTERACTION::SPHMomentumFormulationMonaghan());
      break;
    }
    default:
    {
      dserror("unknown acceleration formulation type!");
      break;
    }
  }

  // init momentum formulation handler
  momentumformulation_->Init();
}

void PARTICLEINTERACTION::SPHMomentum::InitArtificialViscosityHandler()
{
  // create artificial viscosity handler
  artificialviscosity_ = std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity>(
      new PARTICLEINTERACTION::SPHArtificialViscosity());

  // init artificial viscosity handler
  artificialviscosity_->Init();
}

void PARTICLEINTERACTION::SPHMomentum::SetupParticleInteractionWriter()
{
  // register specific runtime vtp writer
  if (writeparticlewallinteraction_)
    particleinteractionwriter_->RegisterSpecificRuntimeVtuWriter("particle-wall-momentum");
}

void PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleContribution");

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(allfluidtypes_, relindices);

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

    // get material for particle types
    const MAT::PAR::ParticleMaterialSPHFluid* material_i = fluidmaterial_[type_i];
    const MAT::PAR::ParticleMaterialSPHFluid* material_j = fluidmaterial_[type_j];

    // get pointer to particle states
    const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* press_i = container_i->GetPtrToState(PARTICLEENGINE::Pressure, particle_i);
    const double* vel_i = container_i->GetPtrToState(PARTICLEENGINE::Velocity, particle_i);

    double* acc_i = nullptr;
    if (intfluidtypes_.count(type_i))
      acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

    const double* mod_vel_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedVelocity, particle_i);
    double* mod_acc_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedAcceleration, particle_i);

    // get pointer to particle states
    const double* rad_j = container_j->GetPtrToState(PARTICLEENGINE::Radius, particle_j);
    const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j);
    const double* press_j = container_j->GetPtrToState(PARTICLEENGINE::Pressure, particle_j);
    const double* vel_j = container_j->GetPtrToState(PARTICLEENGINE::Velocity, particle_j);

    double* acc_j = nullptr;
    if (intfluidtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
      acc_j = container_j->GetPtrToState(PARTICLEENGINE::Acceleration, particle_j);

    const double* mod_vel_j =
        container_j->CondGetPtrToState(PARTICLEENGINE::ModifiedVelocity, particle_j);

    double* mod_acc_j = nullptr;
    if (status_j == PARTICLEENGINE::Owned)
      mod_acc_j = container_j->CondGetPtrToState(PARTICLEENGINE::ModifiedAcceleration, particle_j);

    // evaluate specific coefficient
    double speccoeff_ij(0.0);
    double speccoeff_ji(0.0);
    momentumformulation_->SpecificCoefficient(dens_i, dens_j, mass_i, mass_j, particlepair.dWdrij_,
        particlepair.dWdrji_, &speccoeff_ij, &speccoeff_ji);

    // evaluate pressure gradient
    momentumformulation_->PressureGradient(dens_i, dens_j, press_i, press_j, speccoeff_ij,
        speccoeff_ji, particlepair.e_ij_, acc_i, acc_j);

    // evaluate shear forces
    {
      // get factor from kernel space dimension
      int kernelfac = 0;
      kernel_->KernelSpaceDimension(kernelfac);
      kernelfac += 2;

      // evaluate shear forces
      momentumformulation_->ShearForces(dens_i, dens_j, vel_i, vel_j, kernelfac,
          material_i->dynamicViscosity_, material_j->dynamicViscosity_, material_i->bulkViscosity_,
          material_j->bulkViscosity_, particlepair.absdist_, speccoeff_ij, speccoeff_ji,
          particlepair.e_ij_, acc_i, acc_j);
    }

    // apply transport velocity formulation
    if (transportvelocityformulation_ ==
        INPAR::PARTICLE::TransportVelocityFormulation::StandardTransportVelocity)
    {
      // evaluate background pressure (standard formulation)
      momentumformulation_->StandardBackgroundPressure(dens_i, dens_j,
          material_i->backgroundPressure_, material_j->backgroundPressure_, speccoeff_ij,
          speccoeff_ji, particlepair.e_ij_, mod_acc_i, mod_acc_j);

      // evaluate convection of momentum with relative velocity
      momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j, mod_vel_i,
          mod_vel_j, speccoeff_ij, speccoeff_ji, particlepair.e_ij_, acc_i, acc_j);
    }
    else if (transportvelocityformulation_ ==
             INPAR::PARTICLE::TransportVelocityFormulation::GeneralizedTransportVelocity)
    {
      // modified first derivative of kernel
      const double mod_dWdrij =
          (mod_acc_i) ? kernel_->dWdrij(particlepair.absdist_, kernel_->SmoothingLength(rad_i[0]))
                      : 0.0;
      const double mod_dWdrji =
          (mod_acc_j) ? kernel_->dWdrij(particlepair.absdist_, kernel_->SmoothingLength(rad_j[0]))
                      : 0.0;

      // modified background pressure
      const double mod_bg_press_i =
          (mod_acc_i) ? std::min(std::abs(10.0 * press_i[0]), material_i->backgroundPressure_)
                      : 0.0;
      const double mod_bg_press_j =
          (mod_acc_j) ? std::min(std::abs(10.0 * press_j[0]), material_j->backgroundPressure_)
                      : 0.0;

      // evaluate background pressure (generalized formulation)
      momentumformulation_->GeneralizedBackgroundPressure(dens_i, dens_j, mass_i, mass_j,
          mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, particlepair.e_ij_, mod_acc_i,
          mod_acc_j);

      // evaluate convection of momentum with relative velocity
      momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j, mod_vel_i,
          mod_vel_j, speccoeff_ij, speccoeff_ji, particlepair.e_ij_, acc_i, acc_j);
    }

    // evaluate artificial viscosity
    if (material_i->artificialViscosity_ > 0.0 or material_j->artificialViscosity_ > 0.0)
    {
      // particle averaged smoothing length
      const double h_ij =
          0.5 * (kernel_->SmoothingLength(rad_i[0]) + kernel_->SmoothingLength(rad_j[0]));

      // get speed of sound
      const double c_i = material_i->SpeedOfSound();
      const double c_j = (type_i == type_j) ? c_i : material_j->SpeedOfSound();

      // particle averaged speed of sound
      const double c_ij = 0.5 * (c_i + c_j);

      // particle averaged density
      const double dens_ij = 0.5 * (dens_i[0] + dens_j[0]);

      // evaluate artificial viscosity
      artificialviscosity_->ArtificialViscosity(vel_i, vel_j, mass_i, mass_j,
          material_i->artificialViscosity_, material_j->artificialViscosity_, particlepair.dWdrij_,
          particlepair.dWdrji_, dens_ij, h_ij, c_ij, particlepair.absdist_, particlepair.e_ij_,
          acc_i, acc_j);
    }
  }
}

void PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleBoundaryContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleBoundaryContribution");

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForDisjointCombination(
      intfluidtypes_, boundarytypes_, relindices);

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
    UTILS::VecSet(e_ij, particlepair.e_ij_);
    if (swapparticles) UTILS::VecScale(e_ij, -1.0);

    // first derivative of kernel
    const double dWdrij = (swapparticles) ? particlepair.dWdrji_ : particlepair.dWdrij_;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get material for particle types
    const MAT::PAR::ParticleMaterialSPHFluid* material_i = fluidmaterial_[type_i];

    // get equation of state for particle types
    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_i =
        equationofstatebundle_->GetPtrToSpecificEquationOfState(type_i);

    // get pointer to particle states
    const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* press_i = container_i->GetPtrToState(PARTICLEENGINE::Pressure, particle_i);
    const double* vel_i = container_i->GetPtrToState(PARTICLEENGINE::Velocity, particle_i);

    double* acc_i = nullptr;
    if (status_i == PARTICLEENGINE::Owned)
      acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

    const double* mod_vel_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedVelocity, particle_i);

    double* mod_acc_i = nullptr;
    if (status_i == PARTICLEENGINE::Owned)
      mod_acc_i = container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedAcceleration, particle_i);

    // get pointer to boundary particle states
    const double* mass_j = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* press_j =
        container_j->GetPtrToState(PARTICLEENGINE::BoundaryPressure, particle_j);
    const double* vel_j = container_j->GetPtrToState(PARTICLEENGINE::BoundaryVelocity, particle_j);

    double temp_dens(0.0);
    temp_dens = equationofstate_i->PressureToDensity(press_j[0], material_i->initDensity_);
    const double* dens_j = &temp_dens;

    double* force_j = nullptr;
    if (status_j == PARTICLEENGINE::Owned)
      force_j = container_j->CondGetPtrToState(PARTICLEENGINE::Force, particle_j);

    // contribution from neighboring boundary particle j
    double acc_ij[3] = {0.0, 0.0, 0.0};
    double mod_acc_ij[3] = {0.0, 0.0, 0.0};

    // evaluate specific coefficient
    double speccoeff_ij(0.0);
    momentumformulation_->SpecificCoefficient(
        dens_i, dens_j, mass_i, mass_j, dWdrij, 0.0, &speccoeff_ij, nullptr);

    // evaluate pressure gradient
    momentumformulation_->PressureGradient(
        dens_i, dens_j, press_i, press_j, speccoeff_ij, 0.0, e_ij, acc_ij, nullptr);

    // evaluate shear forces
    if (boundaryparticleinteraction_ == INPAR::PARTICLE::NoSlipBoundaryParticle)
    {
      // get factor from kernel space dimension
      int kernelfac = 0;
      kernel_->KernelSpaceDimension(kernelfac);
      kernelfac += 2;

      // evaluate shear forces
      momentumformulation_->ShearForces(dens_i, dens_j, vel_i, vel_j, kernelfac,
          material_i->dynamicViscosity_, material_i->dynamicViscosity_, material_i->bulkViscosity_,
          material_i->bulkViscosity_, absdist, speccoeff_ij, 0.0, e_ij, acc_ij, nullptr);
    }

    // apply transport velocity formulation
    if (transportvelocityformulation_ ==
        INPAR::PARTICLE::TransportVelocityFormulation::StandardTransportVelocity)
    {
      // evaluate background pressure (standard formulation)
      momentumformulation_->StandardBackgroundPressure(dens_i, dens_j,
          material_i->backgroundPressure_, 0.0, speccoeff_ij, 0.0, e_ij, mod_acc_ij, nullptr);

      // evaluate convection of momentum with relative velocity
      momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j, mod_vel_i,
          nullptr, speccoeff_ij, 0.0, e_ij, acc_ij, nullptr);
    }
    else if (transportvelocityformulation_ ==
             INPAR::PARTICLE::TransportVelocityFormulation::GeneralizedTransportVelocity)
    {
      // modified first derivative of kernel
      const double mod_dWdrij = kernel_->dWdrij(absdist, kernel_->SmoothingLength(rad_i[0]));

      // modified background pressure
      const double mod_bg_press_i =
          std::min(std::abs(10.0 * press_i[0]), material_i->backgroundPressure_);

      // evaluate background pressure (generalized formulation)
      momentumformulation_->GeneralizedBackgroundPressure(dens_i, dens_j, mass_i, mass_j,
          mod_bg_press_i, 0.0, mod_dWdrij, 0.0, e_ij, mod_acc_ij, nullptr);

      // evaluate convection of momentum with relative velocity
      momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j, mod_vel_i,
          nullptr, speccoeff_ij, 0.0, e_ij, acc_ij, nullptr);
    }

    // evaluate artificial viscosity
    if (boundaryparticleinteraction_ == INPAR::PARTICLE::NoSlipBoundaryParticle and
        material_i->artificialViscosity_ > 0.0)
    {
      // get smoothing length
      const double h_i = kernel_->SmoothingLength(rad_i[0]);

      // get speed of sound
      const double c_i = material_i->SpeedOfSound();

      // particle averaged density
      const double dens_ij = 0.5 * (dens_i[0] + dens_j[0]);

      // evaluate artificial viscosity
      artificialviscosity_->ArtificialViscosity(vel_i, vel_j, mass_i, mass_j,
          material_i->artificialViscosity_, 0.0, dWdrij, 0.0, dens_ij, h_i, c_i, absdist, e_ij,
          acc_ij, nullptr);
    }

    // add contribution from neighboring boundary particle j
    if (acc_i) UTILS::VecAdd(acc_i, acc_ij);
    if (mod_acc_i) UTILS::VecAdd(mod_acc_i, mod_acc_ij);

    // add contribution to neighboring boundary particle j
    if (force_j) UTILS::VecAddScale(force_j, -mass_i[0], acc_ij);
  }
}

void PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleWallContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleWallContribution");

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // get reference to particle-wall pair data
  const SPHParticleWallPairData& particlewallpairdata =
      neighborpairs_->GetRefToParticleWallPairData();

  // get number of particle-wall pairs
  const int numparticlewallpairs = particlewallpairdata.size();

  // write interaction output
  const bool writeinteractionoutput =
      particleinteractionwriter_->GetCurrentWriteResultFlag() and writeparticlewallinteraction_;

  // init storage for interaction output
  std::vector<double> attackpoints;
  std::vector<double> contactforces;
  std::vector<double> normaldirection;

  // prepare storage for interaction output
  if (writeinteractionoutput)
  {
    attackpoints.reserve(3 * numparticlewallpairs);
    contactforces.reserve(3 * numparticlewallpairs);
    normaldirection.reserve(3 * numparticlewallpairs);
  }

  // get reference to weighted fluid particle pressure
  const std::vector<double>& weightedpressure = virtualwallparticle_->GetWeightedPressure();

  // get reference to weighted fluid particle pressure gradient
  const std::vector<std::vector<double>>& weightedpressuregradient =
      virtualwallparticle_->GetWeightedPressureGradient();

  // get reference to weighted fluid particle distance vector
  const std::vector<std::vector<double>>& weighteddistancevector =
      virtualwallparticle_->GetWeightedDistanceVector();

  // get reference to weighted fluid particle velocity
  const std::vector<std::vector<double>>& weightedvelocity =
      virtualwallparticle_->GetWeightedVelocity();

  // get relevant particle wall pair indices for specific particle types
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticleWallPairIndices(intfluidtypes_, relindices);

  // iterate over relevant particle-wall pairs
  for (const int particlewallpairindex : relindices)
  {
    const SPHParticleWallPair& particlewallpair = particlewallpairdata[particlewallpairindex];

    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get material for particle types
    const MAT::PAR::ParticleMaterialSPHFluid* material_i = fluidmaterial_[type_i];

    // get equation of state for particle types
    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_i =
        equationofstatebundle_->GetPtrToSpecificEquationOfState(type_i);

    // get pointer to particle states
    const double* pos_i = container_i->GetPtrToState(PARTICLEENGINE::Position, particle_i);
    const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* press_i = container_i->GetPtrToState(PARTICLEENGINE::Pressure, particle_i);
    const double* vel_i = container_i->GetPtrToState(PARTICLEENGINE::Velocity, particle_i);
    double* acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

    const double* mod_vel_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedVelocity, particle_i);
    double* mod_acc_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::ModifiedAcceleration, particle_i);

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // number of nodes of wall element
    const int numnodes = ele->NumNode();

    // shape functions and location vector of wall element
    CORE::LINALG::SerialDenseVector funct(numnodes);
    std::vector<int> lmele;

    if (walldatastate->GetVelCol() != Teuchos::null or
        walldatastate->GetForceCol() != Teuchos::null)
    {
      // evaluate shape functions of element at wall contact point
      CORE::DRT::UTILS::shape_function_2D(
          funct, particlewallpair.elecoords_[0], particlewallpair.elecoords_[1], ele->Shape());

      // get location vector of wall element
      lmele.reserve(numnodes * 3);
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      ele->LocationVector(
          *particlewallinterface_->GetWallDiscretization(), lmele, lmowner, lmstride);
    }

    // velocity of wall contact point j
    double vel_j[3] = {0.0, 0.0, 0.0};

    if (walldatastate->GetVelCol() != Teuchos::null)
    {
      // get nodal velocities
      std::vector<double> nodal_vel(numnodes * 3);
      DRT::UTILS::ExtractMyValues(*walldatastate->GetVelCol(), nodal_vel, lmele);

      // determine velocity of wall contact point j
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim) vel_j[dim] += funct[node] * nodal_vel[node * 3 + dim];
    }

    // sum contribution from neighboring virtual particle k
    double sumk_acc_ik[3] = {0.0, 0.0, 0.0};
    double sumk_mod_acc_ik[3] = {0.0, 0.0, 0.0};

    // compute vector from wall contact point j to particle i
    double r_ij[3];
    UTILS::VecSetScale(r_ij, particlewallpair.absdist_, particlewallpair.e_ij_);

    // vector from weighted fluid particle positions l to wall contact point j
    double r_jl_weighted[3];
    UTILS::VecSet(r_jl_weighted, weighteddistancevector[particlewallpairindex].data());

    // inverse normal distance from weighted fluid particle positions l to wall contact point j
    const double inv_norm_dist_jl_weighted =
        1.0 / UTILS::VecDot(r_jl_weighted, particlewallpair.e_ij_);

    // unit surface tangent vectors in wall contact point j
    double t_j_1[3];
    double t_j_2[3];
    UTILS::UnitSurfaceTangents(particlewallpair.e_ij_, t_j_1, t_j_2);

    // iterate over virtual particles
    for (const std::vector<double>& virtualparticle :
        virtualwallparticle_->GetRelativePositionsOfVirtualParticles())
    {
      // vector from virtual particle k to wall contact point j
      double r_jk[3];
      UTILS::VecSetScale(r_jk, virtualparticle[0], particlewallpair.e_ij_);
      UTILS::VecAddScale(r_jk, virtualparticle[1], t_j_1);
      UTILS::VecAddScale(r_jk, virtualparticle[2], t_j_2);

      // vector from virtual particle k to particle i
      double r_ik[3];
      UTILS::VecSet(r_ik, r_ij);
      UTILS::VecAdd(r_ik, r_jk);

      // absolute distance between virtual particle k and particle i
      const double absdist = UTILS::VecNormTwo(r_ik);

      // virtual particle within interaction distance
      if (absdist < rad_i[0])
      {
        // vector from weighted fluid particle positions l to virtual particle k
        double r_kl_weighted[3];
        UTILS::VecSet(r_kl_weighted, r_jl_weighted);
        UTILS::VecSub(r_kl_weighted, r_jk);

        // get pointer to virtual particle states
        const double* mass_k = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);

        const double temp_press_k =
            weightedpressure[particlewallpairindex] +
            UTILS::VecDot(r_kl_weighted, weightedpressuregradient[particlewallpairindex].data());
        const double* press_k = &temp_press_k;

        const double temp_dens_k =
            equationofstate_i->PressureToDensity(press_k[0], material_i->initDensity_);
        const double* dens_k = &temp_dens_k;

        double temp_vel_k[3];
        double fac = -virtualparticle[0] * inv_norm_dist_jl_weighted;
        UTILS::VecSetScale(temp_vel_k, 1 + fac, vel_j);
        UTILS::VecAddScale(temp_vel_k, -fac, weightedvelocity[particlewallpairindex].data());
        const double* vel_k = temp_vel_k;

        // versor from virtual particle k to particle i
        double e_ik[3];
        UTILS::VecSetScale(e_ik, 1.0 / absdist, r_ik);

        // evaluate first derivative of kernel
        const double dWdrik = kernel_->dWdrij(absdist, rad_i[0]);

        // evaluate specific coefficient
        double speccoeff_ik(0.0);
        momentumformulation_->SpecificCoefficient(
            dens_i, dens_k, mass_i, mass_k, dWdrik, 0.0, &speccoeff_ik, nullptr);

        // evaluate pressure gradient
        momentumformulation_->PressureGradient(
            dens_i, dens_k, press_i, press_k, speccoeff_ik, 0.0, e_ik, sumk_acc_ik, nullptr);

        // evaluate shear forces
        if (boundaryparticleinteraction_ == INPAR::PARTICLE::NoSlipBoundaryParticle)
        {
          // get factor from kernel space dimension
          int kernelfac = 0;
          kernel_->KernelSpaceDimension(kernelfac);
          kernelfac += 2;

          // evaluate shear forces
          momentumformulation_->ShearForces(dens_i, dens_k, vel_i, vel_k, kernelfac,
              material_i->dynamicViscosity_, material_i->dynamicViscosity_,
              material_i->bulkViscosity_, material_i->bulkViscosity_, absdist, speccoeff_ik, 0.0,
              e_ik, sumk_acc_ik, nullptr);
        }

        // apply transport velocity formulation
        if (transportvelocityformulation_ ==
            INPAR::PARTICLE::TransportVelocityFormulation::StandardTransportVelocity)
        {
          // evaluate background pressure (standard formulation)
          momentumformulation_->StandardBackgroundPressure(dens_i, dens_k,
              material_i->backgroundPressure_, 0.0, speccoeff_ik, 0.0, e_ik, sumk_mod_acc_ik,
              nullptr);

          // evaluate convection of momentum with relative velocity
          momentumformulation_->ModifiedVelocityContribution(dens_i, dens_k, vel_i, vel_k,
              mod_vel_i, nullptr, speccoeff_ik, 0.0, e_ik, sumk_acc_ik, nullptr);
        }
        else if (transportvelocityformulation_ ==
                 INPAR::PARTICLE::TransportVelocityFormulation::GeneralizedTransportVelocity)
        {
          // modified first derivative of kernel
          const double mod_dWdrij = kernel_->dWdrij(absdist, kernel_->SmoothingLength(rad_i[0]));

          // modified background pressure
          const double mod_bg_press_i =
              std::min(std::abs(10.0 * press_i[0]), material_i->backgroundPressure_);

          // evaluate background pressure (generalized formulation)
          momentumformulation_->GeneralizedBackgroundPressure(dens_i, dens_k, mass_i, mass_k,
              mod_bg_press_i, 0.0, mod_dWdrij, 0.0, e_ik, sumk_mod_acc_ik, nullptr);

          // evaluate convection of momentum with relative velocity
          momentumformulation_->ModifiedVelocityContribution(dens_i, dens_k, vel_i, vel_k,
              mod_vel_i, nullptr, speccoeff_ik, 0.0, e_ik, sumk_acc_ik, nullptr);
        }

        // evaluate artificial viscosity
        if (boundaryparticleinteraction_ == INPAR::PARTICLE::NoSlipBoundaryParticle and
            material_i->artificialViscosity_ > 0.0)
        {
          // get smoothing length
          const double h_i = kernel_->SmoothingLength(rad_i[0]);

          // get speed of sound
          const double c_i = material_i->SpeedOfSound();

          // particle averaged density
          const double dens_ij = 0.5 * (dens_i[0] + dens_k[0]);

          // evaluate artificial viscosity
          artificialviscosity_->ArtificialViscosity(vel_i, vel_k, mass_i, mass_k,
              material_i->artificialViscosity_, 0.0, dWdrik, 0.0, dens_ij, h_i, c_i, absdist, e_ik,
              sumk_acc_ik, nullptr);
        }
      }
    }

    // add contribution from neighboring virtual particle k
    UTILS::VecAdd(acc_i, sumk_acc_ik);
    if (mod_acc_i) UTILS::VecAdd(mod_acc_i, sumk_mod_acc_ik);

    // calculation of wall contact force
    double wallcontactforce[3] = {0.0, 0.0, 0.0};
    if (writeinteractionoutput or walldatastate->GetForceCol() != Teuchos::null)
      UTILS::VecSetScale(wallcontactforce, -mass_i[0], sumk_acc_ik);

    // write interaction output
    if (writeinteractionoutput)
    {
      // calculate wall contact point
      double wallcontactpoint[3];
      UTILS::VecSet(wallcontactpoint, pos_i);
      UTILS::VecSub(wallcontactpoint, r_ij);

      // set wall attack point and states
      for (int dim = 0; dim < 3; ++dim) attackpoints.push_back(wallcontactpoint[dim]);
      for (int dim = 0; dim < 3; ++dim) contactforces.push_back(wallcontactforce[dim]);
      for (int dim = 0; dim < 3; ++dim) normaldirection.push_back(particlewallpair.e_ij_[dim]);
    }

    // assemble contact force acting on wall element
    if (walldatastate->GetForceCol() != Teuchos::null)
    {
      // determine nodal forces
      std::vector<double> nodal_force(numnodes * 3);
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim)
          nodal_force[node * 3 + dim] = funct[node] * wallcontactforce[dim];

      // assemble nodal forces
      const int err = walldatastate->GetMutableForceCol()->SumIntoGlobalValues(
          numnodes * 3, nodal_force.data(), lmele.data());
      if (err < 0) dserror("sum into Epetra_Vector failed!");
    }
  }

  if (writeinteractionoutput)
  {
    // get specific runtime vtp writer
    IO::VisualizationManager* visualization_manager =
        particleinteractionwriter_->GetSpecificRuntimeVtuWriter("particle-wall-momentum");
    auto& visualization_data = visualization_manager->GetVisualizationDataMutable();

    // set wall attack points
    visualization_data.GetPointCoordinatesMutable() = attackpoints;

    // append states
    visualization_data.SetPointDataVector<double>("contact force", contactforces, 3);
    visualization_data.SetPointDataVector<double>("normal direction", normaldirection, 3);
  }
}
