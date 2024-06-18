/*---------------------------------------------------------------------------*/
/*! \file
\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph_phase_change.hpp"

#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_engine_object.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_sph_equationofstate.hpp"
#include "4C_particle_interaction_sph_equationofstate_bundle.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHPhaseChangeBase::SPHPhaseChangeBase(const Teuchos::ParameterList& params)
    : params_sph_(params),
      belowphase_(PARTICLEENGINE::Phase1),
      abovephase_(PARTICLEENGINE::Phase2),
      transitionstate_(PARTICLEENGINE::Density),
      transitionvalue_(0.0),
      hysteresisgap_(0.0)
{
  // empty constructor
}

void ParticleInteraction::SPHPhaseChangeBase::Init()
{
  // read from input file
  std::string word;
  std::istringstream phasechangedefinition(
      Teuchos::getNumericStringParameter(params_sph_, "PHASECHANGEDEFINITION"));

  // get phase below transition value
  if (phasechangedefinition >> word)
    belowphase_ = PARTICLEENGINE::EnumFromTypeName(word);
  else
    FOUR_C_THROW("expecting particle type for phase below transition value!");

  // get phase above transition value
  if (phasechangedefinition >> word)
    abovephase_ = PARTICLEENGINE::EnumFromTypeName(word);
  else
    FOUR_C_THROW("expecting particle type for phase above transition value!");

  // safety check
  if (belowphase_ == abovephase_)
    FOUR_C_THROW("equal particle types for phase below and above transition value!");

  // get transition state of phase change
  if (phasechangedefinition >> word)
    transitionstate_ = PARTICLEENGINE::EnumFromStateName(word);
  else
    FOUR_C_THROW("expecting particle state of phase change!");

  // safety check
  if (PARTICLEENGINE::EnumToStateDim(transitionstate_) != 1)
    FOUR_C_THROW("expecting scalar particle state for phase change!");

  // get transition value of phase change
  if (phasechangedefinition >> word)
  {
    try
    {
      transitionvalue_ = std::stod(word);
    }
    catch (...)
    {
      FOUR_C_THROW(
          "invalid argument: expecting a double value for transition value of phase change!");
    };
  }
  else
    FOUR_C_THROW("expecting transition value of phase change!");

  // probe for optional hysteresis gap at transition value
  if (phasechangedefinition >> word)
  {
    if (not(word == "hysteresis")) FOUR_C_THROW("expecting optional keyword 'hysteresis'!");

    // get hysteresis gap at transition value
    if (phasechangedefinition >> word)
    {
      try
      {
        hysteresisgap_ = std::stod(word);
      }
      catch (...)
      {
        FOUR_C_THROW("invalid argument: expecting a double value for hysteresis gap!");
      };
    }
    else
      FOUR_C_THROW("expecting hysteresis gap at transition value!");
  }
}

void ParticleInteraction::SPHPhaseChangeBase::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<ParticleInteraction::MaterialHandler> particlematerial,
    const std::shared_ptr<ParticleInteraction::SPHEquationOfStateBundle> equationofstatebundle)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->get_particle_container_bundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;

  // safety check
  for (const auto& type_i : {belowphase_, abovephase_})
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      FOUR_C_THROW("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
}

void ParticleInteraction::SPHPhaseChangeBase::evaluate_phase_change_from_below_to_above_phase(
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
      particlecontainerbundle_->get_specific_container(type_source, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  int particlestored = container->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  // get pointer to particle state
  const double* state = container->GetPtrToState(transitionstate_, 0);

  // get material for particle types
  const Mat::PAR::ParticleMaterialBase* material_source =
      particlematerial_->get_ptr_to_particle_mat_parameter(type_source);
  const Mat::PAR::ParticleMaterialBase* material_target =
      particlematerial_->get_ptr_to_particle_mat_parameter(type_target);

  // get equation of state of target particle type
  const ParticleInteraction::SPHEquationOfStateBase* equationofstate_target;
  if (not isboundaryrigid_target)
    equationofstate_target =
        equationofstatebundle_->get_ptr_to_specific_equation_of_state(type_target);

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

void ParticleInteraction::SPHPhaseChangeBase::evaluate_phase_change_from_above_to_below_phase(
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
      particlecontainerbundle_->get_specific_container(type_source, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  int particlestored = container->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  // get pointer to particle state
  const double* state = container->GetPtrToState(transitionstate_, 0);

  // get material for particle types
  const Mat::PAR::ParticleMaterialBase* material_source =
      particlematerial_->get_ptr_to_particle_mat_parameter(type_source);
  const Mat::PAR::ParticleMaterialBase* material_target =
      particlematerial_->get_ptr_to_particle_mat_parameter(type_target);

  // get equation of state of target particle type
  const ParticleInteraction::SPHEquationOfStateBase* equationofstate_target;
  if (not isboundaryrigid_target)
    equationofstate_target =
        equationofstatebundle_->get_ptr_to_specific_equation_of_state(type_target);

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

ParticleInteraction::SPHPhaseChangeOneWayScalarBelowToAbove::SPHPhaseChangeOneWayScalarBelowToAbove(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void ParticleInteraction::SPHPhaseChangeOneWayScalarBelowToAbove::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from below to above phase
  evaluate_phase_change_from_below_to_above_phase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->hand_over_particles_to_be_removed(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->hand_over_particles_to_be_inserted(particlestoinsert);
}

ParticleInteraction::SPHPhaseChangeOneWayScalarAboveToBelow::SPHPhaseChangeOneWayScalarAboveToBelow(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void ParticleInteraction::SPHPhaseChangeOneWayScalarAboveToBelow::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from above to below phase
  evaluate_phase_change_from_above_to_below_phase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->hand_over_particles_to_be_removed(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->hand_over_particles_to_be_inserted(particlestoinsert);
}

ParticleInteraction::SPHPhaseChangeTwoWayScalar::SPHPhaseChangeTwoWayScalar(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params)
{
  // empty constructor
}

void ParticleInteraction::SPHPhaseChangeTwoWayScalar::EvaluatePhaseChange(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // evaluate phase change from below to above phase
  evaluate_phase_change_from_below_to_above_phase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // evaluate phase change from above to below phase
  evaluate_phase_change_from_above_to_below_phase(
      particlesfromphasetophase, particlestoremove, particlestoinsert);

  // hand over particles to be removed
  particleengineinterface_->hand_over_particles_to_be_removed(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->hand_over_particles_to_be_inserted(particlestoinsert);
}

FOUR_C_NAMESPACE_CLOSE
