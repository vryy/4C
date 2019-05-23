/*---------------------------------------------------------------------------*/
/*!
\brief discrete element method (DEM) interaction handler

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_history_pairs.H"
#include "particle_interaction_dem_contact.H"

#include "../drt_particle_algorithm/particle_wall_interface.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionDEM::ParticleInteractionDEM(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params), params_dem_(params.sublist("DEM"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionDEM::~ParticleInteractionDEM()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init particle interaction handler                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
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
}

/*---------------------------------------------------------------------------*
 | setup particle interaction handler                         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface, particlewallinterface);

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, particlewallinterface);

  // setup history pair handler
  historypairs_->Setup(particleengineinterface);

  // setup contact handler
  contact_->Setup(particleengineinterface, particlewallinterface, particlematerial_, neighborpairs_,
      historypairs_);
}

/*---------------------------------------------------------------------------*
 | write restart of particle interaction handler              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
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
}

/*---------------------------------------------------------------------------*
 | read restart of particle interaction handler               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
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
}

/*---------------------------------------------------------------------------*
 | insert interaction dependent states of all particle types  sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | set initial states                                         sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialStates()
{
  // set initial radius
  SetInitialRadius();

  // set initial mass
  SetInitialMass();
}

/*---------------------------------------------------------------------------*
 | evaluate particle interactions                             sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions");

  // clear force and moment states of particles
  ClearForceAndMomentStates();

  // evaluate neighbor pairs
  neighborpairs_->EvaluateNeighborPairs();

  // check critical time step
  contact_->CheckCriticalTimeStep();

  // add contact contribution to force and moment field
  contact_->AddForceAndMomentContribution();

  // compute acceleration from force and moment
  ComputeAcceleration();

  // update history pairs
  historypairs_->UpdateHistoryPairs();
}

/*---------------------------------------------------------------------------*
 | maximum interaction distance (on this processor)           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::ParticleInteractionDEM::MaxInteractionDistance() const
{
  return (2.0 * MaxParticleRadius());
}

/*---------------------------------------------------------------------------*
 | distribute interaction history                             sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::DistributeInteractionHistory() const
{
  // distribute history pairs
  historypairs_->DistributeHistoryPairs();
}

/*---------------------------------------------------------------------------*
 | communicate interaction history                            sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::CommunicateInteractionHistory() const
{
  // communicate history pairs
  historypairs_->CommunicateHistoryPairs();
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::SetCurrentStepSize(const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::SetCurrentStepSize(currentstepsize);

  // set current step size
  contact_->SetCurrentStepSize(currentstepsize);
}

/*---------------------------------------------------------------------------*
 | init neighbor pair handler                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::InitNeighborPairHandler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::DEMNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

/*---------------------------------------------------------------------------*
 | init history pair handler                                  sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::InitHistoryPairHandler()
{
  // create history pair handler
  historypairs_ = std::make_shared<PARTICLEINTERACTION::DEMHistoryPairs>(comm_);

  // init history pair handler
  historypairs_->Init();
}

/*---------------------------------------------------------------------------*
 | init contact handler                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::InitContactHandler()
{
  // create contact handler
  contact_ = std::unique_ptr<PARTICLEINTERACTION::DEMContact>(
      new PARTICLEINTERACTION::DEMContact(params_dem_));

  // init contact handler
  contact_->Init();
}

/*---------------------------------------------------------------------------*
 | set initial radius                                         sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialRadius()
{
  // get type of (random) particle radius distribution
  INPAR::PARTICLE::RadiusDistribution radiusdistributiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::RadiusDistribution>(
          params_dem_, "RADIUSDISTRIBUTION");

  switch (radiusdistributiontype)
  {
    // no (random) particle radius distribution
    case INPAR::PARTICLE::NoRadiusDistribution:
    {
      // iterate over particle types
      for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(typeEnum);

        // (initial) radius of current phase
        std::vector<double> initradius(1);
        initradius[0] = material->initRadius_;

        // set initial radius for all particles of current type
        particlecontainerbundle_->SetStateSpecificContainer(
            initradius, PARTICLEENGINE::Radius, typeEnum);
      }

      break;
    }
    // normal or log-normal random particle radius distribution
    case INPAR::PARTICLE::NormalRadiusDistribution:
    case INPAR::PARTICLE::LogNormalRadiusDistribution:
    {
      // get allowed bounds for particle radius
      double r_min = params_dem_.get<double>("MIN_RADIUS");
      double r_max = params_dem_.get<double>("MAX_RADIUS");

      // safety checks
      if (not(r_min > 0.0)) dserror("non-positive minimum particle radius!");
      if (not(r_max > 0.0)) dserror("non-positive maximum particle radius!");
      if (r_min > r_max) dserror("minimum particle radius larger than maximum particle radius!");

      // get variance of random particle radius distribution
      double variance = params_dem_.get<double>("RADIUSDISTRIBUTION_VAR");

      // safety check
      if (not(variance > 0.0))
        dserror("non-positive variance of random particle radius distribution!");

      // iterate over particle types
      for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(typeEnum);

        // get pointer to particle state
        double* radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);

        // determine mean of normal random number generator
        const double mean = (radiusdistributiontype == INPAR::PARTICLE::NormalRadiusDistribution)
                                ? material->initRadius_
                                : std::log(material->initRadius_);

        // initialize random number generator
        DRT::Problem::Instance()->Random()->SetMeanVariance(mean, variance);

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

/*---------------------------------------------------------------------------*
 | set initial mass                                           sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialMass()
{
  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(typeEnum);

    // declare pointer variables
    const double* radius;
    double* mass;

    // get pointer to particle states
    radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);
    mass = container->GetPtrToParticleState(PARTICLEENGINE::Mass, 0);

    // compute mass via particle volume and initial density
    const double fac = material->initDensity_ * 4.0 / 3.0 * M_PI;
    for (int i = 0; i < particlestored; ++i) mass[i] = fac * UTILS::pow<3>(radius[i]);
  }
}

/*---------------------------------------------------------------------------*
 | clear force and moment states of particles                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::ClearForceAndMomentStates() const
{
  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // clear force of all particles
    container->ClearState(PARTICLEENGINE::Force);

    // clear moment of all particles
    if (contact_->HaveTangentialContact()) container->ClearState({PARTICLEENGINE::Moment});
  }
}

/*---------------------------------------------------------------------------*
 | compute acceleration from force and moment                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration");

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get particle state dimension
    const int statedim = container->GetParticleStateDim(PARTICLEENGINE::Acceleration);

    // declare pointer variables
    const double *mass, *radius, *force, *moment;
    double *acc, *angacc;

    // get pointer to particle states
    mass = container->GetPtrToParticleState(PARTICLEENGINE::Mass, 0);
    force = container->GetPtrToParticleState(PARTICLEENGINE::Force, 0);
    acc = container->GetPtrToParticleState(PARTICLEENGINE::Acceleration, 0);

    if (contact_->HaveTangentialContact())
    {
      radius = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);
      moment = container->GetPtrToParticleState(PARTICLEENGINE::Moment, 0);
      angacc = container->GetPtrToParticleState(PARTICLEENGINE::AngularAcceleration, 0);
    }

    // compute acceleration
    for (int i = 0; i < particlestored; ++i)
      UTILS::vec_addscale(&acc[statedim * i], (1.0 / mass[i]), &force[statedim * i]);

    // compute angular acceleration
    if (contact_->HaveTangentialContact())
      for (int i = 0; i < particlestored; ++i)
        UTILS::vec_addscale(&angacc[statedim * i],
            (5.0 / (2.0 * mass[i] * UTILS::pow<2>(radius[i]))), &moment[statedim * i]);
  }
}
