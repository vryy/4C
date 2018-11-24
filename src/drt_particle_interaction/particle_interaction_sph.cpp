/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph.cpp

\brief smoothed particle hydrodynamics (SPH) interaction handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph.H"

#include "particle_interaction_material_handler.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_density.H"
#include "particle_interaction_sph_pressure.H"
#include "particle_interaction_sph_temperature.H"
#include "particle_interaction_sph_momentum.H"
#include "particle_interaction_sph_surface_tension.H"
#include "particle_interaction_sph_boundary_particle.H"
#include "particle_interaction_sph_phase_change.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionSPH::ParticleInteractionSPH(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params), params_sph_(params.sublist("SPH"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionSPH::~ParticleInteractionSPH()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init particle interaction handler                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::Init()
{
  // call base class init
  ParticleInteractionBase::Init();

  // init kernel handler
  InitKernelHandler();

  // init equation of state bundle
  InitEquationOfStateBundle();

  // init neighbor pair handler
  InitNeighborPairHandler();

  // init density handler
  InitDensityHandler();

  // init pressure handler
  InitPressureHandler();

  // init temperature handler
  InitTemperatureHandler();

  // init momentum handler
  InitMomentumHandler();

  // init surface tension handler
  InitSurfaceTensionHandler();

  // init boundary particle handler
  InitBoundaryParticleHandler();

  // init phase change handler
  InitPhaseChangeHandler();
}

/*---------------------------------------------------------------------------*
 | setup particle interaction handler                         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface);

  // setup kernel handler
  kernel_->Setup();

  // setup equation of state bundle
  equationofstatebundle_->Setup();

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, kernel_);

  // setup density handler
  density_->Setup(
      particleengineinterface, kernel_, particlematerial_, equationofstatebundle_, neighborpairs_);

  // setup pressure handler
  pressure_->Setup(particleengineinterface, particlematerial_, equationofstatebundle_);

  // setup temperature handler
  if (temperature_) temperature_->Setup(particleengineinterface, particlematerial_, neighborpairs_);

  // setup momentum handler
  momentum_->Setup(
      particleengineinterface, kernel_, particlematerial_, equationofstatebundle_, neighborpairs_);

  // setup surface tension handler
  if (surfacetension_)
    surfacetension_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Setup(particleengineinterface, neighborpairs_);

  // setup phase change handler
  if (phasechange_) phasechange_->Setup(particleengineinterface);
}

/*---------------------------------------------------------------------------*
 | write restart of particle interaction handler              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::WriteRestart(
    const int step, const double time) const
{
  // call base class function
  ParticleInteractionBase::WriteRestart(step, time);

  // write restart of kernel handler
  kernel_->WriteRestart(step, time);

  // write restart of equation of state bundle
  equationofstatebundle_->WriteRestart(step, time);

  // write restart of neighbor pair handler
  neighborpairs_->WriteRestart(step, time);

  // write restart of density handler
  density_->WriteRestart(step, time);

  // write restart of pressure handler
  pressure_->WriteRestart(step, time);

  // write restart of temperature handler
  if (temperature_) temperature_->WriteRestart(step, time);

  // write restart of momentum handler
  momentum_->WriteRestart(step, time);

  // write restart of surface tension handler
  if (surfacetension_) surfacetension_->WriteRestart(step, time);

  // write restart of boundary particle handler
  if (boundaryparticle_) boundaryparticle_->WriteRestart(step, time);

  // write restart of phase change handler
  if (phasechange_) phasechange_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of particle interaction handler               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  ParticleInteractionBase::ReadRestart(reader);

  // read restart of kernel handler
  kernel_->ReadRestart(reader);

  // read restart of equation of state bundle
  equationofstatebundle_->ReadRestart(reader);

  // read restart of neighbor pair handler
  neighborpairs_->ReadRestart(reader);

  // read restart of density handler
  density_->ReadRestart(reader);

  // read restart of pressure handler
  pressure_->ReadRestart(reader);

  // read restart of temperature handler
  if (temperature_) temperature_->ReadRestart(reader);

  // read restart of momentum handler
  momentum_->ReadRestart(reader);

  // read restart of surface tension handler
  if (surfacetension_) surfacetension_->ReadRestart(reader);

  // read restart of boundary particle handler
  if (boundaryparticle_) boundaryparticle_->ReadRestart(reader);

  // read restart of phase change handler
  if (phasechange_) phasechange_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | insert interaction dependent states of all particle types  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InsertParticleStatesOfParticleTypes(
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
    else
    {
      // insert states of regular phase particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius, PARTICLEENGINE::Density,
          PARTICLEENGINE::Pressure});
    }
  }

  // states for density evaluation scheme
  density_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // states for temperature evaluation scheme
  if (temperature_) temperature_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // insert momentum evaluation dependent states
  momentum_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // additional states for surface tension formulation
  if (surfacetension_) surfacetension_->InsertParticleStatesOfParticleTypes(particlestatestotypes);
}

/*---------------------------------------------------------------------------*
 | set initial states                                         sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::SetInitialStates()
{
  // compute initial consistent particle volume
  double consistentparticlevolume = ComputeConsistentParticleVolume();

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(type, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type);

    // initial density of current phase
    std::vector<double> initdensity(1);
    initdensity[0] = material->initDensity_;

    // (initial) mass of current phase
    std::vector<double> initmass(1);
    initmass[0] = initdensity[0] * consistentparticlevolume;

    // (initial) radius of current phase
    std::vector<double> initradius(1);
    initradius[0] = material->initRadius_;

    // set initial density for all non-boundary and non-rigid particles
    if (type != PARTICLEENGINE::BoundaryPhase and type != PARTICLEENGINE::RigidPhase)
      particlecontainerbundle_->SetStateSpecificContainer(
          initdensity, PARTICLEENGINE::Density, type);

    // set initial mass and radius for all particles of current type
    particlecontainerbundle_->SetStateSpecificContainer(initmass, PARTICLEENGINE::Mass, type);
    particlecontainerbundle_->SetStateSpecificContainer(initradius, PARTICLEENGINE::Radius, type);

    // initial states for temperature evaluation
    if (temperature_)
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialThermo* material =
          dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
              particlematerial_->GetPtrToParticleMatParameter(type));

      // initial temperature of current phase
      std::vector<double> inittemperature(1);
      inittemperature[0] = material->initTemperature_;

      // set initial temperature for all particles of current type
      particlecontainerbundle_->SetStateSpecificContainer(
          inittemperature, PARTICLEENGINE::Temperature, type);
    }
  }
}

/*---------------------------------------------------------------------------*
 | evaluate particle interactions                             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::EvaluateInteractions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::EvaluateInteractions");

  // evaluate particle neighbor pairs
  neighborpairs_->EvaluateNeighborPairs();

  // compute density field
  density_->ComputeDensity();

  // compute pressure using equation of state and density
  pressure_->ComputePressure();

  // compute temperature field
  if (temperature_) temperature_->ComputeTemperature();

  // init boundary particles
  if (boundaryparticle_) boundaryparticle_->InitBoundaryParticles(gravity_);

  // add momentum contribution to acceleration field
  momentum_->AddAccelerationContribution();

  // add surface tension contribution to acceleration field
  if (surfacetension_) surfacetension_->AddAccelerationContribution();

  // evaluate phase change
  if (phasechange_) phasechange_->EvaluatePhaseChange();
}

/*---------------------------------------------------------------------------*
 | maximum interaction distance (on this processor)           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::ParticleInteractionSPH::MaxInteractionDistance() const
{
  return MaxParticleRadius();
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::SetCurrentTime(const double currenttime)
{
  // call base class method
  ParticleInteractionBase::SetCurrentTime(currenttime);

  // set current time
  if (surfacetension_) surfacetension_->SetCurrentTime(currenttime);
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::SetCurrentStepSize(const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::SetCurrentStepSize(currentstepsize);

  // set current step size
  density_->SetCurrentStepSize(currentstepsize);

  // set current step size
  if (temperature_) temperature_->SetCurrentStepSize(currentstepsize);
}

/*---------------------------------------------------------------------------*
 | init kernel handler                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitKernelHandler()
{
  // get type of smoothed particle hydrodynamics kernel
  INPAR::PARTICLE::KernelType kerneltype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::KernelType>(params_sph_, "KERNEL");

  // create kernel handler
  switch (kerneltype)
  {
    case INPAR::PARTICLE::CubicSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_);
      break;
    }
    case INPAR::PARTICLE::QuinticSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_);
      break;
    }
    default:
    {
      dserror("unknown kernel type!");
      break;
    }
  }

  // init kernel handler
  kernel_->Init();
}

/*---------------------------------------------------------------------------*
 | init equation of state bundle                              sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitEquationOfStateBundle()
{
  // create equation of state bundle
  equationofstatebundle_ = std::make_shared<PARTICLEINTERACTION::SPHEquationOfStateBundle>(params_);

  // init equation of state bundle
  equationofstatebundle_->Init(particlematerial_);
}

/*---------------------------------------------------------------------------*
 | init neighbor pair handler                                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitNeighborPairHandler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::SPHNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

/*---------------------------------------------------------------------------*
 | init density handler                                       sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitDensityHandler()
{
  // get type of smoothed particle hydrodynamics density evaluation scheme
  INPAR::PARTICLE::DensityEvaluationScheme densityevaluationscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityEvaluationScheme>(
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
      dserror("unknown density evaluation scheme type!");
      break;
    }
  }

  // init density handler
  density_->Init();

  // safety check
  if (densityevaluationscheme != INPAR::PARTICLE::DensityPredictCorrect and
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityCorrectionScheme>(
          params_sph_, "DENSITYCORRECTION") != INPAR::PARTICLE::NoCorrection)
    dserror(
        "the density correction scheme set is not valid with the current density evaluation "
        "scheme!");
}

/*---------------------------------------------------------------------------*
 | init pressure handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitPressureHandler()
{
  // create pressure handler
  pressure_ =
      std::unique_ptr<PARTICLEINTERACTION::SPHPressure>(new PARTICLEINTERACTION::SPHPressure());

  // init pressure handler
  pressure_->Init();
}

/*---------------------------------------------------------------------------*
 | init temperature handler                                    meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitTemperatureHandler()
{
  // get type of smoothed particle hydrodynamics temperature evaluation scheme
  INPAR::PARTICLE::TemperatureEvaluationScheme temperatureevaluationscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::TemperatureEvaluationScheme>(
          params_sph_, "TEMPERATUREEVALUATION");

  // create temperature handler
  switch (temperatureevaluationscheme)
  {
    case INPAR::PARTICLE::NoTemperatureEvaluation:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperatureBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::TemperatureIntegration:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperatureIntegration>(
          new PARTICLEINTERACTION::SPHTemperatureIntegration(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown surface tension formulation type!");
      break;
    }
  }

  // init temperature handler
  if (temperature_) temperature_->Init();
}

/*---------------------------------------------------------------------------*
 | init momentum handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitMomentumHandler()
{
  // create momentum handler
  momentum_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentum>(
      new PARTICLEINTERACTION::SPHMomentum(params_sph_));

  // init momentum handler
  momentum_->Init();
}

/*---------------------------------------------------------------------------*
 | init surface tension handler                               sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitSurfaceTensionHandler()
{
  // get type of smoothed particle hydrodynamics surface tension formulation
  INPAR::PARTICLE::SurfaceTensionFormulation surfacetensionformulation =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionFormulation>(
          params_sph_, "SURFACETENSIONFORMULATION");

  // create surface tension handler
  switch (surfacetensionformulation)
  {
    case INPAR::PARTICLE::NoSurfaceTension:
    {
      surfacetension_ = std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTensionBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::ContinuumSurfaceForce:
    {
      surfacetension_ =
          std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce>(
              new PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown surface tension formulation type!");
      break;
    }
  }

  // init surface tension handler
  if (surfacetension_) surfacetension_->Init();
}

/*---------------------------------------------------------------------------*
 | init boundary particle handler                             sfuchs 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitBoundaryParticleHandler()
{
  // get type of boundary particle formulation
  INPAR::PARTICLE::BoundaryParticleFormulationType boundaryparticleformulation =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::BoundaryParticleFormulationType>(
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
      dserror("unknown boundary particle formulation type!");
      break;
    }
  }

  // init boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Init();
}

/*---------------------------------------------------------------------------*
 | init phase change handler                                  sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionSPH::InitPhaseChangeHandler()
{
  // get type phase change
  INPAR::PARTICLE::PhaseChangeType phasechangetype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::PhaseChangeType>(params_sph_, "PHASECHANGETYPE");

  // create phase change handler
  switch (phasechangetype)
  {
    case INPAR::PARTICLE::NoPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeBase>(nullptr);
      break;
    }
    default:
    {
      dserror("unknown boundary particle formulation type!");
      break;
    }
  }

  // init phase change handler
  if (phasechange_) phasechange_->Init();
}

/*---------------------------------------------------------------------------*
 | compute initial consistent particle volume                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::ParticleInteractionSPH::ComputeConsistentParticleVolume() const
{
  const double consistentproblemvolume = params_sph_.get<double>("CONSISTENTPROBLEMVOLUME");

  if (not(consistentproblemvolume > 0.0))
    dserror("consistent problem volume needs to be defined in SPH simulations!");

  // get number of particles on this processor
  int numberofparticles = particleengineinterface_->GetNumberOfParticles();

  // get number of boundary and rigid particles on this processor
  int numberofboundaryparticles =
      particleengineinterface_->GetNumberOfParticlesOfSpecificType(PARTICLEENGINE::BoundaryPhase);
  int numberofrigidparticles =
      particleengineinterface_->GetNumberOfParticlesOfSpecificType(PARTICLEENGINE::RigidPhase);

  // number of non-boundary particles on this processor
  int numberofnonboundaryparticles =
      numberofparticles - numberofboundaryparticles - numberofrigidparticles;

  // get total number of non-boundary particles on all processors
  int totalnumberofnonboundaryparticles = 0;
  comm_.SumAll(&numberofnonboundaryparticles, &totalnumberofnonboundaryparticles, 1);

  // consistent volume of non-boundary particles
  return (consistentproblemvolume / totalnumberofnonboundaryparticles);
}
