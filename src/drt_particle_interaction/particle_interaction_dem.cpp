/*---------------------------------------------------------------------------*/
/*! \file
\brief discrete element method (DEM) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_history_pairs.H"
#include "particle_interaction_dem_contact.H"
#include "particle_interaction_dem_adhesion.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_particle_wall/particle_wall_interface.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionDEM::ParticleInteractionDEM(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params), params_dem_(params.sublist("DEM"))
{
  // empty constructor
}

PARTICLEINTERACTION::ParticleInteractionDEM::~ParticleInteractionDEM() = default;

void PARTICLEINTERACTION::ParticleInteractionDEM::Init()
{
  // call base class init
  ParticleInteractionBase::Init();

  // init neighbor pair handler
  InitNeighborPairHandler();

  // init history pair handler
  InitHistoryPairHandler();

  // init contact handler
  InitContactHandler();

  // init adhesion handler
  InitAdhesionHandler();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface, particlewallinterface);

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, particlewallinterface);

  // setup history pair handler
  historypairs_->Setup(particleengineinterface);

  // setup contact handler
  contact_->Setup(particleengineinterface, particlewallinterface, particlematerial_,
      particleinteractionwriter_, neighborpairs_, historypairs_);

  // setup adhesion handler
  if (adhesion_)
    adhesion_->Setup(particleengineinterface, particlewallinterface, particleinteractionwriter_,
        neighborpairs_, historypairs_, contact_->GetNormalContactStiffness());
}

void PARTICLEINTERACTION::ParticleInteractionDEM::WriteRestart(
    const int step, const double time) const
{
  // call base class function
  ParticleInteractionBase::WriteRestart(step, time);

  // write restart of neighbor pair handler
  neighborpairs_->WriteRestart(step, time);

  // write restart of history pair handler
  historypairs_->WriteRestart(step, time);

  // write restart of contact handler
  contact_->WriteRestart(step, time);

  // write restart of adhesion handler
  if (adhesion_) adhesion_->WriteRestart(step, time);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  ParticleInteractionBase::ReadRestart(reader);

  // read restart of neighbor pair handler
  neighborpairs_->ReadRestart(reader);

  // read restart of history pair handler
  historypairs_->ReadRestart(reader);

  // read restart of contact handler
  contact_->ReadRestart(reader);

  // read restart of adhesion handler
  if (adhesion_) adhesion_->ReadRestart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // insert states of regular phase particles
    particlestates.insert({PARTICLEENGINE::Force, PARTICLEENGINE::Mass, PARTICLEENGINE::Radius});
  }

  // states for contact evaluation scheme
  contact_->InsertParticleStatesOfParticleTypes(particlestatestotypes);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialStates()
{
  // set initial radius
  SetInitialRadius();

  // set initial mass
  SetInitialMass();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::PrepareTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::PrepareTimeStep");
}

void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions");

  // clear force and moment states of particles
  ClearForceAndMomentStates();

  // evaluate neighbor pairs
  neighborpairs_->EvaluateNeighborPairs();

  // evaluate adhesion neighbor pairs
  if (adhesion_) neighborpairs_->EvaluateNeighborPairsAdhesion(adhesion_->GetAdhesionDistance());

  // check critical time step
  contact_->CheckCriticalTimeStep();

  // add contact contribution to force and moment field
  contact_->AddForceAndMomentContribution();

  // add adhesion contribution to force field
  if (adhesion_) adhesion_->AddForceContribution();

  // compute acceleration from force and moment
  ComputeAcceleration();

  // update history pairs
  historypairs_->UpdateHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::PostEvaluateTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::PostEvaluateTimeStep");
}

double PARTICLEINTERACTION::ParticleInteractionDEM::MaxInteractionDistance() const
{
  // particle contact interaction distance
  double interactiondistance = 2.0 * MaxParticleRadius();

  // add adhesion distance
  if (adhesion_) interactiondistance += adhesion_->GetAdhesionDistance();

  return interactiondistance;
}

void PARTICLEINTERACTION::ParticleInteractionDEM::DistributeInteractionHistory() const
{
  // distribute history pairs
  historypairs_->DistributeHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::CommunicateInteractionHistory() const
{
  // communicate history pairs
  historypairs_->CommunicateHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetCurrentStepSize(const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::SetCurrentStepSize(currentstepsize);

  // set current step size
  contact_->SetCurrentStepSize(currentstepsize);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitNeighborPairHandler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::DEMNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitHistoryPairHandler()
{
  // create history pair handler
  historypairs_ = std::make_shared<PARTICLEINTERACTION::DEMHistoryPairs>(comm_);

  // init history pair handler
  historypairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitContactHandler()
{
  // create contact handler
  contact_ = std::unique_ptr<PARTICLEINTERACTION::DEMContact>(
      new PARTICLEINTERACTION::DEMContact(params_dem_));

  // init contact handler
  contact_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitAdhesionHandler()
{
  // get type of adhesion law
  INPAR::PARTICLE::AdhesionLaw adhesionlaw =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::AdhesionLaw>(params_dem_, "ADHESIONLAW");

  // create adhesion handler
  if (adhesionlaw != INPAR::PARTICLE::NoAdhesion)
    adhesion_ = std::unique_ptr<PARTICLEINTERACTION::DEMAdhesion>(
        new PARTICLEINTERACTION::DEMAdhesion(params_dem_));

  // init adhesion handler
  if (adhesion_) adhesion_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialRadius()
{
  // get allowed bounds for particle radius
  double r_min = params_dem_.get<double>("MIN_RADIUS");
  double r_max = params_dem_.get<double>("MAX_RADIUS");

  // safety checks
  if (r_min < 0.0) dserror("negative minimum allowed particle radius!");
  if (not(r_max > 0.0)) dserror("non-positive maximum allowed particle radius!");
  if (r_min > r_max)
    dserror("minimum allowed particle radius larger than maximum allowed particle radius!");

  // get type of initial particle radius assignment
  INPAR::PARTICLE::InitialRadiusAssignment radiusdistributiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::InitialRadiusAssignment>(
          params_dem_, "INITIAL_RADIUS");

  switch (radiusdistributiontype)
  {
    // particle radius from particle material
    case INPAR::PARTICLE::RadiusFromParticleMaterial:
    {
      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(type_i);

        // safety checks
        if (material->initRadius_ < r_min)
          dserror("material particle radius smaller than minimum allowed particle radius!");

        if (material->initRadius_ > r_max)
          dserror("material particle radius larger than maximum allowed particle radius!");

        // (initial) radius of current phase
        std::vector<double> initradius(1);
        initradius[0] = material->initRadius_;

        // set initial radius for all particles of current type
        container->SetState(initradius, PARTICLEENGINE::Radius);
      }

      break;
    }
    // particle radius from particle input
    case INPAR::PARTICLE::RadiusFromParticleInput:
    {
      // note: particle radius set as read in from input file, only safety checks here

      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // safety checks
        if (container->GetMinValueOfState(PARTICLEENGINE::Radius) < r_min)
          dserror("minimum particle radius smaller than minimum allowed particle radius!");

        if (container->GetMaxValueOfState(PARTICLEENGINE::Radius) > r_max)
          dserror("maximum particle radius larger than maximum allowed particle radius!");
      }

      break;
    }
    // normal or log-normal random particle radius distribution
    case INPAR::PARTICLE::NormalRadiusDistribution:
    case INPAR::PARTICLE::LogNormalRadiusDistribution:
    {
      // get sigma of random particle radius distribution
      double sigma = params_dem_.get<double>("RADIUSDISTRIBUTION_SIGMA");

      // safety check
      if (not(sigma > 0.0)) dserror("non-positive sigma of random particle radius distribution!");

      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(type_i);

        // get pointer to particle state
        double* radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);

        // determine mu of random particle radius distribution
        const double mu = (radiusdistributiontype == INPAR::PARTICLE::NormalRadiusDistribution)
                              ? material->initRadius_
                              : std::log(material->initRadius_);

        // initialize random number generator
        DRT::Problem::Instance()->Random()->SetMeanVariance(mu, sigma);

        // iterate over particles stored in container
        for (int i = 0; i < particlestored; ++i)
        {
          // generate random value
          const double randomvalue = DRT::Problem::Instance()->Random()->Normal();

          // set normal or log-normal distributed random value for particle radius
          radius[i] = (radiusdistributiontype == INPAR::PARTICLE::NormalRadiusDistribution)
                          ? randomvalue
                          : std::exp(randomvalue);

          // adjust radius to allowed bounds
          if (radius[i] > r_max)
            radius[i] = r_max;
          else if (radius[i] < r_min)
            radius[i] = r_min;
        }
      }
      break;
    }
    default:
    {
      dserror("invalid type of (random) particle radius distribution!");
      break;
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialMass()
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // get pointer to particle states
    const double* radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);
    double* mass = container->GetPtrToParticleState(PARTICLEENGINE::Mass, 0);

    // compute mass via particle volume and initial density
    const double fac = material->initDensity_ * 4.0 / 3.0 * M_PI;
    for (int i = 0; i < particlestored; ++i) mass[i] = fac * UTILS::pow<3>(radius[i]);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ClearForceAndMomentStates() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear force of all particles
    container->ClearState(PARTICLEENGINE::Force);

    // clear moment of all particles
    if (container->HaveStoredState(PARTICLEENGINE::Moment))
      container->ClearState(PARTICLEENGINE::Moment);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get particle state dimension
    const int statedim = container->GetParticleStateDim(PARTICLEENGINE::Acceleration);

    // get pointer to particle states
    const double* radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);
    const double* mass = container->GetPtrToParticleState(PARTICLEENGINE::Mass, 0);
    const double* force = container->GetPtrToParticleState(PARTICLEENGINE::Force, 0);
    double* acc = container->GetPtrToParticleState(PARTICLEENGINE::Acceleration, 0);

    const double* moment = container->HaveStoredState(PARTICLEENGINE::Moment)
                               ? container->GetPtrToParticleState(PARTICLEENGINE::Moment, 0)
                               : nullptr;

    double* angacc = container->HaveStoredState(PARTICLEENGINE::AngularAcceleration)
                         ? container->GetPtrToParticleState(PARTICLEENGINE::AngularAcceleration, 0)
                         : nullptr;

    // compute acceleration
    for (int i = 0; i < particlestored; ++i)
      UTILS::vec_addscale(&acc[statedim * i], (1.0 / mass[i]), &force[statedim * i]);

    // compute angular acceleration
    if (angacc and moment)
      for (int i = 0; i < particlestored; ++i)
        UTILS::vec_addscale(&angacc[statedim * i],
            (5.0 / (2.0 * mass[i] * UTILS::pow<2>(radius[i]))), &moment[statedim * i]);
  }
}
