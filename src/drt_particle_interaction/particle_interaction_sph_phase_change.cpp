/*---------------------------------------------------------------------------*/
/*! \file
\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
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
 | declarations                                                              |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHPhaseChangeBase::SPHPhaseChangeBase(const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::Init()
{
  // nothing to do
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
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHPhaseChangeBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::SPHPhaseChangeTwoWayScalar(
    const Teuchos::ParameterList& params)
    : SPHPhaseChangeBase::SPHPhaseChangeBase(params),
      belowphase_(PARTICLEENGINE::Phase1),
      abovephase_(PARTICLEENGINE::Phase2),
      transitionstate_(PARTICLEENGINE::Density),
      transitionvalue_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::Init()
{
  // call base class init
  SPHPhaseChangeBase::Init();

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
}

void PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle)
{
  // call base class setup
  SPHPhaseChangeBase::Setup(particleengineinterface, particlematerial, equationofstatebundle);

  // safety check
  for (const auto& typeEnum : {belowphase_, abovephase_})
    if (not particlecontainerbundle_->GetParticleTypes().count(typeEnum))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
}

void PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar::EvaluatePhaseChange() const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // iterate over source type of particles
  for (const auto& type_source : {belowphase_, abovephase_})
  {
    // determine target type of particles
    PARTICLEENGINE::TypeEnum type_target = (type_source == belowphase_) ? abovephase_ : belowphase_;

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
    if (particlestored <= 0) continue;

    // get pointer to particle state
    const double* state = container->GetPtrToParticleState(transitionstate_, 0);

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
      // evaluate transition condition for phase change
      bool havephasechange = (type_source == belowphase_) ? (state[index] > transitionvalue_)
                                                          : (state[index] < transitionvalue_);

      // phase change of current particle
      if (havephasechange)
      {
        int globalid(0);
        PARTICLEENGINE::ParticleStates particleStates;
        container->GetParticle(index, globalid, particleStates);

        // add density and pressure state for boundary or rigid particles
        if (isboundaryrigid_source and (not isboundaryrigid_target))
        {
          particleStates[PARTICLEENGINE::Density].resize(1, material_source->initDensity_);

          const double press = equationofstate_target->DensityToPressure(
              material_source->initDensity_, material_target->initDensity_);

          particleStates[PARTICLEENGINE::Pressure].resize(1, press);
        }

        PARTICLEENGINE::ParticleObjShrdPtr particleobject =
            std::make_shared<PARTICLEENGINE::ParticleObject>(type_target, globalid, particleStates);

        // append particle to be insert
        particlestoinsert[type_target].push_back(std::make_pair(-1, particleobject));

        // store index of particle to be removed from containers
        particlestoremove[type_source].insert(index);
      }
    }
  }

  // change type of particles
  particleengineinterface_->TypeChangeParticles(particlestoremove, particlestoinsert);
}
