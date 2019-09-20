/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature boundary condition handler for particle simulations

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_temperature_bc.H"

#include "particle_algorithm_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::TemperatureBoundaryConditionHandler(
    const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::Init()
{
  // get control parameters for conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_conditions, "TEMPERATURE_BOUNDARY_CONDITION", temperaturebctypetofunctid_);

  // iterate over particle types and insert into set
  for (auto& typeIt : temperaturebctypetofunctid_)
    typessubjectedtotemperaturebc_.insert(typeIt.first);
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types subjected to temperature boundary conditions
  for (auto& particleType : typessubjectedtotemperaturebc_)
  {
    // insert states for types subjected to temperature boundary conditions
    particlestatestotypes[particleType].insert(PARTICLEENGINE::ReferencePosition);
  }
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::SetParticleReferencePosition() const
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types subjected to temperature boundary conditions
  for (auto& particleType : typessubjectedtotemperaturebc_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // set particle reference position
    container->UpdateState(0.0, PARTICLEENGINE::ReferencePosition, 1.0, PARTICLEENGINE::Position);
  }
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::EvaluateTemperatureBoundaryCondition(
    const double& evaltime) const
{
  // init vector containing evaluated function
  std::vector<double> funct(1);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types subjected to temperature boundary conditions
  for (auto& typeIt : temperaturebctypetofunctid_)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get id of function
    const int functid = typeIt.second;

    // get reference to function
    DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(functid - 1);

    // declare pointer variables
    const double* refpos;
    double* temp;

    // get pointer to particle states
    refpos = container->GetPtrToParticleState(PARTICLEENGINE::ReferencePosition, 0);
    temp = container->GetPtrToParticleState(PARTICLEENGINE::Temperature, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // safety check
    if (function.NumberComponents() != 1)
      dserror("dimension of function defining temperature boundary condition is not one!");

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // evaluate function
      funct = function.EvaluateTimeDerivative(0, &(refpos[statedim * i]), evaltime, 0);
      temp[i] = funct[0];
    }
  }
}
