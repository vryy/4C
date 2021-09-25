/*---------------------------------------------------------------------------*/
/*! \file
\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_phase_change.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"
#include "../drt_particle_engine/particle_object.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHPhaseChangeBase::SPHPhaseChangeBase(const Teuchos::ParameterList& params)
    : params_sph_(params),
      belowphase_(PARTICLEENGINE::Phase1),
      abovephase_(PARTICLEENGINE::Phase2),
      transitionstate_(PARTICLEENGINE::Density),
      transitionvalue_(0.0),
      hysteresisgap_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::Init()
{
  // read from input file
  std::string word;
  std::istringstream phasechangedefinition(
      Teuchos::getNumericStringParameter(params_sph_, "PHASECHANGEDEFINITION"));

  // get phase below transition value
  if (phasechangedefinition >> word)
    belowphase_ = PARTICLEENGINE::EnumFromTypeName(word);
  else
    dserror("expecting particle type for phase below transition value!");

  // get phase above transition value
  if (phasechangedefinition >> word)
    abovephase_ = PARTICLEENGINE::EnumFromTypeName(word);
  else
    dserror("expecting particle type for phase above transition value!");

  // safety check
  if (belowphase_ == abovephase_)
    dserror("equal particle types for phase below and above transition value!");

  // get transition state of phase change
  if (phasechangedefinition >> word)
    transitionstate_ = PARTICLEENGINE::EnumFromStateName(word);
  else
    dserror("expecting particle state of phase change!");

  // safety check
  if (PARTICLEENGINE::EnumToStateDim(transitionstate_) != 1)
    dserror("expecting scalar particle state for phase change!");

  // get transition value of phase change
  if (phasechangedefinition >> word)
  {
    try
    {
      transitionvalue_ = std::stod(word);
    }
    catch (...)
    {
      dserror("invalid argument: expecting a double value for transition value of phase change!");
    };
  }
  else
    dserror("expecting transition value of phase change!");

  // probe for optional hysteresis gap at transition value
  if (phasechangedefinition >> word)
  {
    if (not(word == "hysteresis")) dserror("expecting optional keyword 'hysteresis'!");

    // get hysteresis gap at transition value
    if (phasechangedefinition >> word)
    {
      try
      {
        hysteresisgap_ = std::stod(word);
      }
      catch (...)
      {
        dserror("invalid argument: expecting a double value for hysteresis gap!");
      };
    }
    else
      dserror("expecting hysteresis gap at transition value!");
  }
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;

  // safety check
  for (const auto& type_i : {belowphase_, abovephase_})
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::EvaluatePhaseChangeFromBelowToAbovePhase(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase,
    std::vector<std::set<int>>& particlestoremove,
    std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>>& particlestoinsert)
    const
{
  // set source and target type of particles
  PARTICLEENGINE::TypeEnum type_source = belowphase_;
  PARTICLEENGINE::TypeEnum type_target = abovephase_;

  // check for boundary or rigid particles
  bool isboundaryrigid_source =
      (type_source == PARTICLEENGINE::BoundaryPhase or type_source == PARTICLEENGINE::RigidPhase);
  bool isboundaryrigid_target =
      (type_target == PARTICLEENGINE::BoundaryPhase or type_target == PARTICLEENGINE::RigidPhase);

  // get container of owned particles of source particle type
  PARTICLEENGINE::ParticleContainer* container =
      particlecontainerbundle_->GetSpecificContainer(type_source, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  int particlestored = container->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  // get pointer to particle state
  const double* state = container->GetPtrToState(transitionstate_, 0);

  // get material for particle types
  const MAT::PAR::ParticleMaterialBase* material_source =
      particlematerial_->GetPtrToParticleMatParameter(type_source);
  const MAT::PAR::ParticleMaterialBase* material_target =
      particlematerial_->GetPtrToParticleMatParameter(type_target);

  // get equation of state of target particle type
  const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_target;
  if (not isboundaryrigid_target)
    equationofstate_target = equationofstatebundle_->GetPtrToSpecificEquationOfState(type_target);

  // iterate over owned particles of current type
  for (int index = 0; index < particlestored; ++index)
  {
    // evaluate phase change condition of current particle
    if (state[index] > (transitionvalue_ + 0.5 * hysteresisgap_))
    {
      int globalid(0);
      PARTICLEENGINE::ParticleStates particlestates;
      container->GetParticle(index, globalid, particlestates);

      // add density and pressure state for boundary or rigid particles
      if (isboundaryrigid_source and (not isboundaryrigid_target))
      {
        particlestates[PARTICLEENGINE::Density].assign(1, material_source->initDensity_);

        const double press = equationofstate_target->DensityToPressure(
            material_source->initDensity_, material_target->initDensity_);

        particlestates[PARTICLEENGINE::Pressure].assign(1, press);
      }

      // clear velocity and acceleration state of boundary or rigid particles
      if (isboundaryrigid_target and (not isboundaryrigid_source))
      {
        particlestates[PARTICLEENGINE::Velocity].assign(3, 0.0);
        particlestates[PARTICLEENGINE::Acceleration].assign(3, 0.0);
      }

      PARTICLEENGINE::ParticleObjShrdPtr particleobject =
          std::make_shared<PARTICLEENGINE::ParticleObject>(type_target, globalid, particlestates);

      // append particle to be insert
      particlestoinsert[type_target].push_back(std::make_pair(-1, particleobject));

      // store index of particle to be removed from containers
      particlestoremove[type_source].insert(index);

      // append source and target type together with global id of particle
      particlesfromphasetophase.push_back(std::make_tuple(type_source, type_target, globalid));
    }
  }
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::EvaluatePhaseChangeFromAboveToBelowPhase(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase,
    std::vector<std::set<int>>& particlestoremove,
    std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>>& particlestoinsert)
    const
{
  // set source and target type of particles
  PARTICLEENGINE::TypeEnum type_source = abovephase_;
  PARTICLEENGINE::TypeEnum type_target = belowphase_;

  // check for boundary or rigid particles
  bool isboundaryrigid_source =
      (type_source == PARTICLEENGINE::BoundaryPhase or type_source == PARTICLEENGINE::RigidPhase);
  bool isboundaryrigid_target =
      (type_target == PARTICLEENGINE::BoundaryPhase or type_target == PARTICLEENGINE::RigidPhase);

  // get container of owned particles of source particle type
  PARTICLEENGINE::ParticleContainer* container =
      particlecontainerbundle_->GetSpecificContainer(type_source, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  int particlestored = container->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  // get pointer to particle state
  const double* state = container->GetPtrToState(transitionstate_, 0);

  // get material for particle types
  const MAT::PAR::ParticleMaterialBase* material_source =
      particlematerial_->GetPtrToParticleMatParameter(type_source);
  const MAT::PAR::ParticleMaterialBase* material_target =
      particlematerial_->GetPtrToParticleMatParameter(type_target);

  // get equation of state of target particle type
  const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_target;
  if (not isboundaryrigid_target)
    equationofstate_target = equationofstatebundle_->GetPtrToSpecificEquationOfState(type_target);

  // iterate over owned particles of current type
  for (int index = 0; index < particlestored; ++index)
  {
    // evaluate phase change condition of current particle
    if (state[index] < (transitionvalue_ - 0.5 * hysteresisgap_))
    {
      int globalid(0);
      PARTICLEENGINE::ParticleStates particlestates;
      container->GetParticle(index, globalid, particlestates);

      // add density and pressure state for boundary or rigid particles
      if (isboundaryrigid_source and (not isboundaryrigid_target))
      {
        particlestates[PARTICLEENGINE::Density].assign(1, material_source->initDensity_);

        const double press = equationofstate_target->DensityToPressure(
            material_source->initDensity_, material_target->initDensity_);

        particlestates[PARTICLEENGINE::Pressure].assign(1, press);
      }

      // clear velocity and acceleration state of boundary or rigid particles
      if (isboundaryrigid_target and (not isboundaryrigid_source))
      {
        particlestates[PARTICLEENGINE::Velocity].assign(3, 0.0);
        particlestates[PARTICLEENGINE::Acceleration].assign(3, 0.0);
      }

      PARTICLEENGINE::ParticleObjShrdPtr particleobject =
          std::make_shared<PARTICLEENGINE::ParticleObject>(type_target, globalid, particlestates);

      // append particle to be insert
      particlestoinsert[type_target].push_back(std::make_pair(-1, particleobject));

      // store index of particle to be removed from containers
      particlestoremove[type_source].insert(index);

      // append source and target type together with global id of particle
      particlesfromphasetophase.push_back(std::make_tuple(type_source, type_target, globalid));
    }
  }
}

PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove::SPHPhaseChangeOneWayScalarBelowToAbove(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from below to above phase
  EvaluatePhaseChangeFromBelowToAbovePhase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->HandOverParticlesToBeRemoved(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->HandOverParticlesToBeInserted(particlestoinsert);
}

PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow::SPHPhaseChangeOneWayScalarAboveToBelow(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from above to below phase
  EvaluatePhaseChangeFromAboveToBelowPhase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->HandOverParticlesToBeRemoved(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->HandOverParticlesToBeInserted(particlestoinsert);
}

PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::SPHPhaseChangeTwoWayScalar(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from below to above phase
  EvaluatePhaseChangeFromBelowToAbovePhase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // evaluate phase change from above to below phase
  EvaluatePhaseChangeFromAboveToBelowPhase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->HandOverParticlesToBeRemoved(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->HandOverParticlesToBeInserted(particlestoinsert);
}
