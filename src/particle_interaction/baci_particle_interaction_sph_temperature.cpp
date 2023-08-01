/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_sph_temperature.H"

#include "baci_particle_interaction_sph_kernel.H"
#include "baci_particle_interaction_material_handler.H"
#include "baci_particle_interaction_sph_neighbor_pairs.H"
#include "baci_particle_interaction_sph_heatsource.H"
#include "baci_particle_interaction_sph_heatloss_evaporation.H"

#include "baci_particle_interaction_utils.H"

#include "baci_particle_engine_interface.H"
#include "baci_particle_engine_container.H"

#include "baci_utils_exceptions.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHTemperature::SPHTemperature(const Teuchos::ParameterList& params)
    : params_sph_(params),
      time_(0.0),
      dt_(0.0),
      temperaturegradient_(DRT::INPUT::IntegralValue<int>(params_sph_, "TEMPERATUREGRADIENT"))
{
  // empty constructor
}

PARTICLEINTERACTION::SPHTemperature::~SPHTemperature() = default;

void PARTICLEINTERACTION::SPHTemperature::Init()
{
  // init heat source handler
  InitHeatSourceHandler();

  // init evaporation induced heat loss handler
  InitHeatLossEvaporationHandler();

  // init with potential integrated thermo particle types
  intthermotypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::RigidPhase};
}

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

  // update with actual integrated thermo particle types
  const auto intthermotypes = intthermotypes_;
  for (const auto& type_i : intthermotypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      intthermotypes_.erase(type_i);

  // setup temperature of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Temperature};

    for (const auto& type_i : intthermotypes_)
      temptorefresh_.push_back(std::make_pair(type_i, states));
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
      dserror("thermal capacity for particles of type '%s' not positive!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
  }

  // setup heat source handler
  if (heatsource_) heatsource_->Setup(particleengineinterface, particlematerial, neighborpairs);

  // setup evaporation induced heat loss handler
  if (heatlossevaporation_) heatlossevaporation_->Setup(particleengineinterface, particlematerial);
}

void PARTICLEINTERACTION::SPHTemperature::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

void PARTICLEINTERACTION::SPHTemperature::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

void PARTICLEINTERACTION::SPHTemperature::InsertParticleStatesOfParticleTypes(
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

    // set temperature state
    particlestates.insert(PARTICLEENGINE::Temperature);

    // current particle type is not an integrated thermo particle type
    if (not intthermotypes_.count(type)) continue;

    // state for temperature evaluation scheme
    particlestates.insert(PARTICLEENGINE::TemperatureDot);

    // state for temperature gradient evaluation
    if (temperaturegradient_) particlestates.insert(PARTICLEENGINE::TemperatureGradient);
  }
}

void PARTICLEINTERACTION::SPHTemperature::ComputeTemperature() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHTemperature::ComputeTemperature");

  // evaluate energy equation
  EnergyEquation();

  // evaluate heat source
  if (heatsource_) heatsource_->EvaluateHeatSource(time_);

  // evaluate evaporation induced heat loss
  if (heatlossevaporation_) heatlossevaporation_->EvaluateEvaporationInducedHeatLoss();

  // iterate over integrated thermo particle types
  for (const auto& type_i : intthermotypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // update temperature of all particles
    container_i->UpdateState(1.0, PARTICLEENGINE::Temperature, dt_, PARTICLEENGINE::TemperatureDot);
  }

  // refresh temperature of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(temptorefresh_);

  // evaluate temperature gradient
  if (temperaturegradient_) TemperatureGradient();
}

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

void PARTICLEINTERACTION::SPHTemperature::InitHeatLossEvaporationHandler()
{
  if (DRT::INPUT::IntegralValue<int>(params_sph_, "VAPOR_HEATLOSS"))
    heatlossevaporation_ = std::unique_ptr<PARTICLEINTERACTION::SPHHeatLossEvaporation>(
        new PARTICLEINTERACTION::SPHHeatLossEvaporation(params_sph_));

  // init evaporation induced heat loss handler
  if (heatlossevaporation_) heatlossevaporation_->Init();
}

void PARTICLEINTERACTION::SPHTemperature::EnergyEquation() const
{
  // iterate over integrated thermo particle types
  for (const auto& type_i : intthermotypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear temperature dot state
    container_i->ClearState(PARTICLEENGINE::TemperatureDot);
  }

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
    const MAT::PAR::ParticleMaterialBase* basematerial_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);
    const MAT::PAR::ParticleMaterialBase* basematerial_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];
    const MAT::PAR::ParticleMaterialThermo* thermomaterial_j = thermomaterial_[type_j];

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);

    const double* dens_i = (container_i->HaveStoredState(PARTICLEENGINE::Density))
                               ? container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i)
                               : &(basematerial_i->initDensity_);

    const double* temp_i = container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i);
    double* tempdot_i = container_i->CondGetPtrToState(PARTICLEENGINE::TemperatureDot, particle_i);

    const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);

    const double* dens_j = (container_j->HaveStoredState(PARTICLEENGINE::Density))
                               ? container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j)
                               : &(basematerial_j->initDensity_);

    const double* temp_j = container_j->GetPtrToState(PARTICLEENGINE::Temperature, particle_j);
    double* tempdot_j = container_j->CondGetPtrToState(PARTICLEENGINE::TemperatureDot, particle_j);

    // thermal conductivities
    const double& k_i = thermomaterial_i->thermalConductivity_;
    const double& k_j = thermomaterial_j->thermalConductivity_;

    // factor containing effective conductivity and temperature difference
    const double fac = (4.0 * k_i * k_j) * (temp_i[0] - temp_j[0]) /
                       (dens_i[0] * dens_j[0] * (k_i + k_j) * particlepair.absdist_);

    // sum contribution of neighboring particle j
    if (tempdot_i)
      tempdot_i[0] +=
          mass_j[0] * fac * particlepair.dWdrij_ * thermomaterial_i->invThermalCapacity_;

    // sum contribution of neighboring particle i
    if (tempdot_j and status_j == PARTICLEENGINE::Owned)
      tempdot_j[0] -=
          mass_i[0] * fac * particlepair.dWdrji_ * thermomaterial_j->invThermalCapacity_;
  }
}

void PARTICLEINTERACTION::SPHTemperature::TemperatureGradient() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHTemperature::TemperatureGradient");

  // iterate over integrated thermo particle types
  for (const auto& type_i : intthermotypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear temperature gradient state
    container_i->ClearState(PARTICLEENGINE::TemperatureGradient);
  }

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
    const MAT::PAR::ParticleMaterialBase* basematerial_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);
    const MAT::PAR::ParticleMaterialBase* basematerial_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);

    const double* dens_i = (container_i->HaveStoredState(PARTICLEENGINE::Density))
                               ? container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i)
                               : &(basematerial_i->initDensity_);

    const double* temp_i = container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i);
    double* tempgrad_i =
        container_i->CondGetPtrToState(PARTICLEENGINE::TemperatureGradient, particle_i);

    const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);

    const double* dens_j = (container_j->HaveStoredState(PARTICLEENGINE::Density))
                               ? container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j)
                               : &(basematerial_j->initDensity_);

    const double* temp_j = container_j->GetPtrToState(PARTICLEENGINE::Temperature, particle_j);
    double* tempgrad_j =
        container_j->CondGetPtrToState(PARTICLEENGINE::TemperatureGradient, particle_j);

    const double temp_ji = temp_j[0] - temp_i[0];

    // sum contribution of neighboring particle j
    if (tempgrad_i)
      UTILS::VecAddScale(
          tempgrad_i, (mass_j[0] / dens_j[0]) * temp_ji * particlepair.dWdrij_, particlepair.e_ij_);

    // sum contribution of neighboring particle i
    if (tempgrad_j and status_j == PARTICLEENGINE::Owned)
      UTILS::VecAddScale(
          tempgrad_j, (mass_i[0] / dens_i[0]) * temp_ji * particlepair.dWdrji_, particlepair.e_ij_);
  }
}
