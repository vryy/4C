/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_temperature.cpp

\brief temperature handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                     meier 09/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_temperature.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_heatsource.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                 meier 09/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHTemperature::SPHTemperature(const Teuchos::ParameterList& params)
    : params_sph_(params), time_(0.0), dt_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHTemperature::~SPHTemperature()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init temperature handler                                    meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::Init()
{
  // init heat source handler
  InitHeatSourceHandler();
}

/*---------------------------------------------------------------------------*
 | setup temperature handler                                   meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // setup temperature of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Temperature};

    for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
    {
      // no refreshing of density states for boundary or rigid particles
      if (typeEnum == PARTICLEENGINE::BoundaryPhase) continue;

      temperaturetorefresh_.push_back(std::make_pair(typeEnum, states));
    }
  }

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  thermomaterial_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    thermomaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));

    // safety check
    if (not(thermomaterial_[type_i]->thermalCapacity_ > 0.0))
      dserror("thermal conductivity for particles of type '%s' not positive!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
  }

  // setup heat source handler
  if (heatsource_) heatsource_->Setup(particleengineinterface, particlematerial);
}

/*---------------------------------------------------------------------------*
 | write restart of temperature handler                        meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::WriteRestart(const int step, const double time) const
{
  // write restart of heat source handler
  if (heatsource_) heatsource_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of temperature handler                         meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of heat source handler
  if (heatsource_) heatsource_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

/*---------------------------------------------------------------------------*
 | insert temperature evaluation dependent states              meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // states for temperature evaluation scheme
    particlestates.insert(PARTICLEENGINE::Temperature);
    particlestates.insert(PARTICLEENGINE::TemperatureDot);
  }
}

/*---------------------------------------------------------------------------*
 | compute temperature field using energy equation             meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::ComputeTemperature() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHTemperature::ComputeTemperature");

  // evaluate energy equation
  EnergyEquation();

  // evaluate heat source
  if (heatsource_) heatsource_->EvaluateHeatSource(time_);

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no temperature calculation for boundary particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase) continue;

    // update temperature of all particles
    particlecontainerbundle_->UpdateStateSpecificContainer(
        1.0, PARTICLEENGINE::Temperature, dt_, PARTICLEENGINE::TemperatureDot, typeEnum);
  }

  // refresh temperature of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(temperaturetorefresh_);
}

/*---------------------------------------------------------------------------*
 | init heat source handler                                   sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::InitHeatSourceHandler()
{
  // get type of heat source
  INPAR::PARTICLE::HeatSourceType heatsourcetype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::HeatSourceType>(params_sph_, "HEATSOURCETYPE");

  // create heat source handler
  switch (heatsourcetype)
  {
    case INPAR::PARTICLE::NoHeatSource:
    {
      heatsource_ = std::unique_ptr<PARTICLEINTERACTION::SPHHeatSourceBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::VolumeHeatSource:
    {
      heatsource_ = std::unique_ptr<PARTICLEINTERACTION::SPHHeatSourceVolume>(
          new PARTICLEINTERACTION::SPHHeatSourceVolume(params_sph_));
      break;
    }
    case INPAR::PARTICLE::SurfaceHeatSource:
    {
      heatsource_ = std::unique_ptr<PARTICLEINTERACTION::SPHHeatSourceSurface>(
          new PARTICLEINTERACTION::SPHHeatSourceSurface(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown type of heat source!");
      break;
    }
  }

  // init heat source handler
  if (heatsource_) heatsource_->Init();
}

/*---------------------------------------------------------------------------*
 | evaluate energy equation                                    meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHTemperature::EnergyEquation() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // no temperature integration for boundary particles
    if (type_i == PARTICLEENGINE::BoundaryPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear temperature dot state
    container_i->ClearState(PARTICLEENGINE::TemperatureDot);
  }

  // iterate over neighbor pairs
  for (auto& neighborpair : neighborpairs_->GetRefToNeighborPairData())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = neighborpair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = neighborpair.tuple_j_;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* basematerial_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);
    const MAT::PAR::ParticleMaterialBase* basematerial_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];
    const MAT::PAR::ParticleMaterialThermo* thermomaterial_j = thermomaterial_[type_j];

    // declare pointer variables for particle i and j
    const double *mass_i, *dens_i, *temp_i;
    double* tempdot_i;

    const double *mass_j, *dens_j, *temp_j;
    double* tempdot_j;

    // get pointer to particle states
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    temp_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Temperature, particle_i);
    tempdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureDot, particle_i);

    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase)
      dens_i = &(basematerial_i->initDensity_);
    else
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);

    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    temp_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Temperature, particle_j);
    tempdot_j = container_j->GetPtrToParticleState(PARTICLEENGINE::TemperatureDot, particle_j);

    if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
      dens_j = &(basematerial_j->initDensity_);
    else
      dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);

    // factor containing temperature difference
    const double fac =
        (thermomaterial_i->thermalConductivity_ + thermomaterial_j->thermalConductivity_) *
        (temp_i[0] - temp_j[0]) / (dens_i[0] * dens_j[0] * neighborpair.absdist_);

    // no temperature integration for boundary particles
    if (type_i != PARTICLEENGINE::BoundaryPhase)
      tempdot_i[0] +=
          mass_j[0] * fac * neighborpair.dWdrij_ * thermomaterial_i->invThermalCapacity_;

    // no temperature integration for boundary particles
    if (type_j != PARTICLEENGINE::BoundaryPhase and status_j == PARTICLEENGINE::Owned)
      tempdot_j[0] -=
          mass_i[0] * fac * neighborpair.dWdrji_ * thermomaterial_j->invThermalCapacity_;
  }
}
