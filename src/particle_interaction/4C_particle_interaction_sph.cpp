/*---------------------------------------------------------------------------*/
/*! \file
\brief smoothed particle hydrodynamics (SPH) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph.hpp"

#include "4C_io_pstream.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_sph_boundary_particle.hpp"
#include "4C_particle_interaction_sph_density.hpp"
#include "4C_particle_interaction_sph_equationofstate.hpp"
#include "4C_particle_interaction_sph_equationofstate_bundle.hpp"
#include "4C_particle_interaction_sph_kernel.hpp"
#include "4C_particle_interaction_sph_momentum.hpp"
#include "4C_particle_interaction_sph_neighbor_pairs.hpp"
#include "4C_particle_interaction_sph_open_boundary.hpp"
#include "4C_particle_interaction_sph_phase_change.hpp"
#include "4C_particle_interaction_sph_pressure.hpp"
#include "4C_particle_interaction_sph_rigid_particle_contact.hpp"
#include "4C_particle_interaction_sph_surface_tension.hpp"
#include "4C_particle_interaction_sph_temperature.hpp"
#include "4C_particle_interaction_sph_virtual_wall_particle.hpp"
#include "4C_particle_interaction_utils.hpp"
#include "4C_particle_wall_interface.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionSPH::ParticleInteractionSPH(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params), params_sph_(params.sublist("SPH"))
{
  // empty constructor
}

PARTICLEINTERACTION::ParticleInteractionSPH::~ParticleInteractionSPH() = default;

void PARTICLEINTERACTION::ParticleInteractionSPH::Init()
{
  // call base class init
  ParticleInteractionBase::Init();

  // init kernel handler
  init_kernel_handler();

  // init equation of state bundle
  init_equation_of_state_bundle();

  // init neighbor pair handler
  init_neighbor_pair_handler();

  // init density handler
  init_density_handler();

  // init pressure handler
  init_pressure_handler();

  // init temperature handler
  init_temperature_handler();

  // init momentum handler
  init_momentum_handler();

  // init surface tension handler
  init_surface_tension_handler();

  // init boundary particle handler
  init_boundary_particle_handler();

  // init dirichlet open boundary handler
  init_dirichlet_open_boundary_handler();

  // init neumann open boundary handler
  init_neumann_open_boundary_handler();

  // init virtual wall particle handler
  init_virtual_wall_particle_handler();

  // init phase change handler
  init_phase_change_handler();

  // init rigid particle contact handler
  init_rigid_particle_contact_handler();

  // safety check
  if (surfacetension_ and virtualwallparticle_)
    FOUR_C_THROW("surface tension formulation with wall interaction not implemented!");
}

void PARTICLEINTERACTION::ParticleInteractionSPH::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface, particlewallinterface);

  // setup kernel handler
  kernel_->Setup();

  // setup equation of state bundle
  equationofstatebundle_->Setup();

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, particlewallinterface, kernel_);

  // setup density handler
  density_->Setup(particleengineinterface, particlewallinterface, kernel_, particlematerial_,
      equationofstatebundle_, neighborpairs_, virtualwallparticle_);

  // setup pressure handler
  pressure_->Setup(particleengineinterface, particlematerial_, equationofstatebundle_);

  // setup temperature handler
  if (temperature_) temperature_->Setup(particleengineinterface, particlematerial_, neighborpairs_);

  // setup momentum handler
  momentum_->Setup(particleengineinterface, particlewallinterface, kernel_, particlematerial_,
      particleinteractionwriter_, equationofstatebundle_, neighborpairs_, virtualwallparticle_);

  // setup surface tension handler
  if (surfacetension_)
    surfacetension_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Setup(particleengineinterface, neighborpairs_);

  // setup dirichlet open boundary handler
  if (dirichletopenboundary_)
    dirichletopenboundary_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup neumann open boundary handler
  if (neumannopenboundary_)
    neumannopenboundary_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup virtual wall particle handler
  if (virtualwallparticle_)
    virtualwallparticle_->Setup(
        particleengineinterface, particlewallinterface, kernel_, neighborpairs_);

  // setup phase change handler
  if (phasechange_)
    phasechange_->Setup(particleengineinterface, particlematerial_, equationofstatebundle_);

  // setup rigid particle contact handler
  if (rigidparticlecontact_)
    rigidparticlecontact_->Setup(
        particleengineinterface, particlewallinterface, particleinteractionwriter_, neighborpairs_);

  // short screen output
  if ((dirichletopenboundary_ or neumannopenboundary_) and
      particleengineinterface_->have_periodic_boundary_conditions())
  {
    if (myrank_ == 0)
      CORE::IO::cout << "Warning: periodic boundary and open boundary conditions applied!"
                     << CORE::IO::endl;
  }
}

void PARTICLEINTERACTION::ParticleInteractionSPH::write_restart() const
{
  // call base class function
  ParticleInteractionBase::write_restart();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::read_restart(
    const std::shared_ptr<CORE::IO::DiscretizationReader> reader)
{
  // call base class function
  ParticleInteractionBase::read_restart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::insert_particle_states_of_particle_types(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase)
    {
      // insert states of boundary and rigid particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius,
          PARTICLEENGINE::BoundaryPressure, PARTICLEENGINE::BoundaryVelocity});
    }
    else if (type == PARTICLEENGINE::DirichletPhase or type == PARTICLEENGINE::NeumannPhase)
    {
      // insert states of open boundary particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius, PARTICLEENGINE::Density,
          PARTICLEENGINE::Pressure});
    }
    else
    {
      // insert states of regular phase particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius, PARTICLEENGINE::Density,
          PARTICLEENGINE::Pressure});
    }
  }

  // states for density evaluation scheme
  density_->insert_particle_states_of_particle_types(particlestatestotypes);

  // states for temperature evaluation scheme
  if (temperature_) temperature_->insert_particle_states_of_particle_types(particlestatestotypes);

  // insert momentum evaluation dependent states
  momentum_->insert_particle_states_of_particle_types(particlestatestotypes);

  // additional states for surface tension formulation
  if (surfacetension_)
    surfacetension_->insert_particle_states_of_particle_types(particlestatestotypes);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::SetInitialStates()
{
  // get kernel space dimension
  int kernelspacedim = 0;
  kernel_->kernel_space_dimension(kernelspacedim);

  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // safety check
  if (not(initialparticlespacing > 0.0)) FOUR_C_THROW("negative initial particle spacing!");

  // compute initial particle volume
  const double initialparticlevolume = std::pow(initialparticlespacing, kernelspacedim);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->get_ptr_to_particle_mat_parameter(type_i);

    // initial density of current phase
    std::vector<double> initdensity(1);
    initdensity[0] = material->initDensity_;

    // (initial) mass of current phase
    std::vector<double> initmass(1);
    initmass[0] = initdensity[0] * initialparticlevolume;

    // (initial) radius of current phase
    std::vector<double> initradius(1);
    initradius[0] = material->initRadius_;

    // set initial density for respective particles of current type
    if (container->HaveStoredState(PARTICLEENGINE::Density))
      container->set_state(initdensity, PARTICLEENGINE::Density);

    // set initial mass and radius for all particles of current type
    container->set_state(initmass, PARTICLEENGINE::Mass);
    container->set_state(initradius, PARTICLEENGINE::Radius);

    // evaluate initial inertia for respective particles of current type
    if (container->HaveStoredState(PARTICLEENGINE::Inertia))
    {
      // (initial) inertia of current phase
      std::vector<double> initinertia(1);

      if (kernelspacedim == 2)
      {
        // effective particle radius considering initial particle volume in disk shape
        const double effectiveradius = std::sqrt(M_1_PI * initialparticlevolume);

        // inertia for disk shape
        initinertia[0] = 0.5 * initmass[0] * UTILS::Pow<2>(effectiveradius);
      }
      else if (kernelspacedim == 3)
      {
        // effective particle radius considering initial particle volume in spherical shape
        const double effectiveradius = std::pow(0.75 * M_1_PI * initialparticlevolume, 1.0 / 3.0);

        // inertia for spherical shape
        initinertia[0] = 0.4 * initmass[0] * UTILS::Pow<2>(effectiveradius);
      }
      else
      {
        FOUR_C_THROW("inertia for particles only in two and three dimensional evaluation given!");
      }

      // set initial inertia for respective particles of current type
      container->set_state(initinertia, PARTICLEENGINE::Inertia);
    }

    // initial states for temperature evaluation
    if (temperature_)
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialThermo* material =
          dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
              particlematerial_->get_ptr_to_particle_mat_parameter(type_i));

      // initial temperature of current phase
      std::vector<double> inittemperature(1);
      inittemperature[0] = material->initTemperature_;

      // set initial temperature for all particles of current type
      container->set_state(inittemperature, PARTICLEENGINE::Temperature);
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionSPH::pre_evaluate_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::pre_evaluate_time_step");

  // prescribe open boundary states
  if (dirichletopenboundary_) dirichletopenboundary_->prescribe_open_boundary_states(time_);

  // prescribe open boundary states
  if (neumannopenboundary_) neumannopenboundary_->prescribe_open_boundary_states(time_);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::evaluate_interactions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::evaluate_interactions");

  // evaluate particle neighbor pairs
  neighborpairs_->evaluate_neighbor_pairs();

  // init relative positions of virtual particles
  if (virtualwallparticle_)
    virtualwallparticle_->init_relative_positions_of_virtual_particles(max_interaction_distance());

  // compute density field
  density_->ComputeDensity();

  // compute pressure using equation of state and density
  pressure_->ComputePressure();

  // compute interface quantities
  if (surfacetension_) surfacetension_->compute_interface_quantities();

  // compute temperature field
  if (temperature_) temperature_->ComputeTemperature();

  // interpolate open boundary states
  if (dirichletopenboundary_) dirichletopenboundary_->interpolate_open_boundary_states();

  // interpolate open boundary states
  if (neumannopenboundary_) neumannopenboundary_->interpolate_open_boundary_states();

  // init boundary particle states
  if (boundaryparticle_) boundaryparticle_->init_boundary_particle_states(gravity_);

  // init states at wall contact points
  if (virtualwallparticle_) virtualwallparticle_->init_states_at_wall_contact_points(gravity_);

  // add momentum contribution to acceleration field
  momentum_->add_acceleration_contribution();

  // add surface tension contribution to acceleration field
  if (surfacetension_) surfacetension_->add_acceleration_contribution();

  // add rigid particle contact contribution to force field
  if (rigidparticlecontact_) rigidparticlecontact_->add_force_contribution();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::post_evaluate_time_step(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::post_evaluate_time_step");

  // check open boundary phase change
  if (dirichletopenboundary_)
    dirichletopenboundary_->check_open_boundary_phase_change(max_interaction_distance());

  // check open boundary phase change
  if (neumannopenboundary_)
    neumannopenboundary_->check_open_boundary_phase_change(max_interaction_distance());

  // evaluate phase change
  if (phasechange_) phasechange_->EvaluatePhaseChange(particlesfromphasetophase);
}

double PARTICLEINTERACTION::ParticleInteractionSPH::max_interaction_distance() const
{
  return max_particle_radius();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::distribute_interaction_history() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionSPH::communicate_interaction_history() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionSPH::set_current_time(const double currenttime)
{
  // call base class method
  ParticleInteractionBase::set_current_time(currenttime);

  // set current time
  if (temperature_) temperature_->set_current_time(currenttime);

  // set current time
  if (surfacetension_) surfacetension_->set_current_time(currenttime);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::set_current_step_size(
    const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::set_current_step_size(currentstepsize);

  // set current step size
  density_->set_current_step_size(currentstepsize);

  // set current step size
  if (temperature_) temperature_->set_current_step_size(currentstepsize);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_kernel_handler()
{
  // get type of smoothed particle hydrodynamics kernel
  INPAR::PARTICLE::KernelType kerneltype =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::KernelType>(params_sph_, "KERNEL");

  // create kernel handler
  switch (kerneltype)
  {
    case INPAR::PARTICLE::CubicSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_sph_);
      break;
    }
    case INPAR::PARTICLE::QuinticSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_sph_);
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown kernel type!");
      break;
    }
  }

  // init kernel handler
  kernel_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_equation_of_state_bundle()
{
  // create equation of state bundle
  equationofstatebundle_ =
      std::make_shared<PARTICLEINTERACTION::SPHEquationOfStateBundle>(params_sph_);

  // init equation of state bundle
  equationofstatebundle_->Init(particlematerial_);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_neighbor_pair_handler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::SPHNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_density_handler()
{
  // get type of smoothed particle hydrodynamics density evaluation scheme
  INPAR::PARTICLE::DensityEvaluationScheme densityevaluationscheme =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::DensityEvaluationScheme>(
          params_sph_, "DENSITYEVALUATION");

  // create density handler
  switch (densityevaluationscheme)
  {
    case INPAR::PARTICLE::DensitySummation:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensitySummation>(
          new PARTICLEINTERACTION::SPHDensitySummation(params_sph_));
      break;
    }
    case INPAR::PARTICLE::DensityIntegration:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityIntegration>(
          new PARTICLEINTERACTION::SPHDensityIntegration(params_sph_));
      break;
    }
    case INPAR::PARTICLE::DensityPredictCorrect:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityPredictCorrect>(
          new PARTICLEINTERACTION::SPHDensityPredictCorrect(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown density evaluation scheme type!");
      break;
    }
  }

  // init density handler
  density_->Init();

  // safety check
  if (densityevaluationscheme != INPAR::PARTICLE::DensityPredictCorrect and
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::DensityCorrectionScheme>(
          params_sph_, "DENSITYCORRECTION") != INPAR::PARTICLE::NoCorrection)
    FOUR_C_THROW(
        "the density correction scheme set is not valid with the current density evaluation "
        "scheme!");
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_pressure_handler()
{
  // create pressure handler
  pressure_ =
      std::unique_ptr<PARTICLEINTERACTION::SPHPressure>(new PARTICLEINTERACTION::SPHPressure());

  // init pressure handler
  pressure_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_temperature_handler()
{
  // get type of smoothed particle hydrodynamics temperature evaluation scheme
  INPAR::PARTICLE::TemperatureEvaluationScheme temperatureevaluationscheme =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::TemperatureEvaluationScheme>(
          params_sph_, "TEMPERATUREEVALUATION");

  // create temperature handler
  switch (temperatureevaluationscheme)
  {
    case INPAR::PARTICLE::NoTemperatureEvaluation:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperature>(nullptr);
      break;
    }
    case INPAR::PARTICLE::TemperatureIntegration:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperature>(
          new PARTICLEINTERACTION::SPHTemperature(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown temperature evaluation scheme!");
      break;
    }
  }

  // init temperature handler
  if (temperature_) temperature_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_momentum_handler()
{
  // create momentum handler
  momentum_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentum>(
      new PARTICLEINTERACTION::SPHMomentum(params_sph_));

  // init momentum handler
  momentum_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_surface_tension_handler()
{
  // get type of smoothed particle hydrodynamics surface tension formulation
  INPAR::PARTICLE::SurfaceTensionFormulation surfacetensionformulation =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::SurfaceTensionFormulation>(
          params_sph_, "SURFACETENSIONFORMULATION");

  // create surface tension handler
  switch (surfacetensionformulation)
  {
    case INPAR::PARTICLE::NoSurfaceTension:
    {
      surfacetension_ = std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTension>(nullptr);
      break;
    }
    case INPAR::PARTICLE::ContinuumSurfaceForce:
    {
      surfacetension_ = std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTension>(
          new PARTICLEINTERACTION::SPHSurfaceTension(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown surface tension formulation type!");
      break;
    }
  }

  // init surface tension handler
  if (surfacetension_) surfacetension_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_boundary_particle_handler()
{
  // get type of boundary particle formulation
  INPAR::PARTICLE::BoundaryParticleFormulationType boundaryparticleformulation =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::BoundaryParticleFormulationType>(
          params_sph_, "BOUNDARYPARTICLEFORMULATION");

  // create boundary particle handler
  switch (boundaryparticleformulation)
  {
    case INPAR::PARTICLE::NoBoundaryFormulation:
    {
      boundaryparticle_ = std::unique_ptr<PARTICLEINTERACTION::SPHBoundaryParticleBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::AdamiBoundaryFormulation:
    {
      boundaryparticle_ = std::unique_ptr<PARTICLEINTERACTION::SPHBoundaryParticleAdami>(
          new PARTICLEINTERACTION::SPHBoundaryParticleAdami(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown boundary particle formulation type!");
      break;
    }
  }

  // init boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_dirichlet_open_boundary_handler()
{
  // get type of dirichlet open boundary
  INPAR::PARTICLE::DirichletOpenBoundaryType dirichletopenboundarytype =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::DirichletOpenBoundaryType>(
          params_sph_, "DIRICHLETBOUNDARYTYPE");

  // create open boundary handler
  switch (dirichletopenboundarytype)
  {
    case INPAR::PARTICLE::NoDirichletOpenBoundary:
    {
      dirichletopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::DirichletNormalToPlane:
    {
      dirichletopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryDirichlet>(
          new PARTICLEINTERACTION::SPHOpenBoundaryDirichlet(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown dirichlet open boundary type!");
      break;
    }
  }

  // init open boundary handler
  if (dirichletopenboundary_) dirichletopenboundary_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_neumann_open_boundary_handler()
{
  // get type of neumann open boundary
  INPAR::PARTICLE::NeumannOpenBoundaryType neumannopenboundarytype =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::NeumannOpenBoundaryType>(
          params_sph_, "NEUMANNBOUNDARYTYPE");

  // create open boundary handler
  switch (neumannopenboundarytype)
  {
    case INPAR::PARTICLE::NoNeumannOpenBoundary:
    {
      neumannopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::NeumannNormalToPlane:
    {
      neumannopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryNeumann>(
          new PARTICLEINTERACTION::SPHOpenBoundaryNeumann(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown neumann open boundary type!");
      break;
    }
  }

  // init open boundary handler
  if (neumannopenboundary_) neumannopenboundary_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_virtual_wall_particle_handler()
{
  // get type of wall formulation
  INPAR::PARTICLE::WallFormulationType wallformulation =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::WallFormulationType>(
          params_sph_, "WALLFORMULATION");

  // create virtual wall particle handler
  switch (wallformulation)
  {
    case INPAR::PARTICLE::NoWallFormulation:
    {
      virtualwallparticle_ = std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle>(nullptr);
      break;
    }
    case INPAR::PARTICLE::VirtualParticleWallFormulation:
    {
      virtualwallparticle_ =
          std::make_shared<PARTICLEINTERACTION::SPHVirtualWallParticle>(params_sph_);
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown wall formulation type!");
      break;
    }
  }

  // init virtual wall particle handler
  if (virtualwallparticle_) virtualwallparticle_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_phase_change_handler()
{
  // get type of phase change
  INPAR::PARTICLE::PhaseChangeType phasechangetype =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::PhaseChangeType>(params_sph_, "PHASECHANGETYPE");

  // create phase change handler
  switch (phasechangetype)
  {
    case INPAR::PARTICLE::NoPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::OneWayScalarBelowToAbovePhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove>(
          new PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove(params_sph_));
      break;
    }
    case INPAR::PARTICLE::OneWayScalarAboveToBelowPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow>(
          new PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow(params_sph_));
      break;
    }
    case INPAR::PARTICLE::TwoWayScalarPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar>(
          new PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown boundary particle formulation type!");
      break;
    }
  }

  // init phase change handler
  if (phasechange_) phasechange_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::init_rigid_particle_contact_handler()
{
  // get type of rigid particle contact
  INPAR::PARTICLE::RigidParticleContactType rigidparticlecontacttype =
      CORE::UTILS::IntegralValue<INPAR::PARTICLE::RigidParticleContactType>(
          params_sph_, "RIGIDPARTICLECONTACTTYPE");

  // create rigid particle contact handler
  switch (rigidparticlecontacttype)
  {
    case INPAR::PARTICLE::NoRigidParticleContact:
    {
      rigidparticlecontact_ =
          std::unique_ptr<PARTICLEINTERACTION::SPHRigidParticleContactBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::ElasticRigidParticleContact:
    {
      rigidparticlecontact_ = std::unique_ptr<PARTICLEINTERACTION::SPHRigidParticleContactElastic>(
          new PARTICLEINTERACTION::SPHRigidParticleContactElastic(params_sph_));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown rigid particle contact type!");
      break;
    }
  }

  // init rigid particle contact handler
  if (rigidparticlecontact_) rigidparticlecontact_->Init();
}

FOUR_C_NAMESPACE_CLOSE
