/*---------------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

\brief algorithm to control particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_algorithm.H"

#include "particle_algorithm_utils.H"

#include "particle_timint.H"
#include "particle_input_generator.H"
#include "particle_gravity.H"
#include "particle_viscous_damping.H"
#include "particle_wall.H"
#include "particle_initial_field.H"
#include "particle_result_test.H"

#include "../drt_particle_interaction/particle_interaction_base.H"
#include "../drt_particle_interaction/particle_interaction_sph.H"
#include "../drt_particle_interaction/particle_interaction_dem.H"

#include "../drt_particle_engine/particle_engine.H"
#include "../drt_particle_engine/particle_communication_utils.H"
#include "../drt_particle_engine/particle_object.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_resulttest.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ParticleAlgorithm::ParticleAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : AlgorithmBase(comm, params),
      myrank_(comm.MyPID()),
      params_(params),
      numparticlesafterlastloadbalance_(0),
      transferevery_(DRT::INPUT::IntegralValue<int>(params_, "TRANSFER_EVERY")),
      writeresultsevery_(params.get<int>("RESULTSEVRY")),
      writerestartevery_(params.get<int>("RESTARTEVRY")),
      writeresultsthisstep_(false),
      writerestartthisstep_(false),
      isrestarted_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ParticleAlgorithm::~ParticleAlgorithm()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init particle algorithm                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::Init(
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles)
{
  // init particle engine
  InitParticleEngine();

  // init particle time integration
  InitParticleTimeIntegration();

  // init particle interaction handler
  InitParticleInteraction();

  // init particle gravity handler
  InitParticleGravity();

  // init viscous damping handler
  InitViscousDamping();

  // init particle wall handler
  InitParticleWall();

  // set initial particles to vector of particles to be distributed
  particlestodistribute_ = initialparticles;

  // clear vector of initial particles in global problem
  initialparticles.clear();
}

/*---------------------------------------------------------------------------*
 | setup particle algorithm                                   sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::Setup()
{
  // generate initial particles
  if (not isrestarted_) GenerateInitialParticles();

  // determine all particle types
  DetermineParticleTypes();

  // determine particle states of all particle types
  DetermineParticleStatesOfParticleTypes();

  // setup particle engine
  particleengine_->Setup(particlestatestotypes_);

  // check particle types for modified states
  const bool havemodifiedstates = HaveModifiedStates();

  // setup particle time integration
  particletimint_->Setup(particleengine_, havemodifiedstates);

  // setup particle interaction handler
  if (particleinteraction_) particleinteraction_->Setup(particleengine_);

  // setup gravity handler
  if (particlegravity_) particlegravity_->Setup();

  // setup viscous damping handler
  if (viscousdamping_) viscousdamping_->Setup(particleengine_);

  // setup wall handler
  if (particlewall_) particlewall_->Setup(particleengine_, particleengine_->GetBinningStrategy());

  // setup initial particles
  SetupInitialParticles();

  // setup initial wall distribution
  if (particlewall_) SetupInitalWallDistribution();

  // setup initial states
  if (not isrestarted_) SetupInitialStates();
}

/*---------------------------------------------------------------------------*
 | read restart information for given time step               sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::ReadRestart(const int restartstep)
{
  // clear vector of particles to be distributed
  particlestodistribute_.clear();

  // create discretization reader
  const std::shared_ptr<IO::DiscretizationReader> reader =
      particleengine_->BinDisReader(restartstep);

  // safety check
  if (restartstep != reader->ReadInt("step")) dserror("time step on file not equal to given step!");

  // get restart time
  double restarttime = reader->ReadDouble("time");

  // read restart of particle engine
  particleengine_->ReadRestart(reader, particlestodistribute_);

  // read restart of particle time integration
  particletimint_->ReadRestart(reader);

  // read restart of particle interaction handler
  if (particleinteraction_) particleinteraction_->ReadRestart(reader);

  // read restart of gravity handler
  if (particlegravity_) particlegravity_->ReadRestart(reader);

  // read restart of viscous damping handler
  if (viscousdamping_) viscousdamping_->ReadRestart(reader);

  // read restart of wall handler
  if (particlewall_) particlewall_->ReadRestart(restartstep);

  // set time and step after restart
  SetTimeStep(restarttime, restartstep);

  // set flag indicating restart to true
  isrestarted_ = true;

  // short screen output
  if (myrank_ == 0)
    IO::cout << "====== restart of the particle simulation from step " << restartstep << IO::endl;
}

/*---------------------------------------------------------------------------*
 | time loop for particle problem                             sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    PrepareTimeStep();

    // integrate particle time step
    Integrate();

    // output particle time step
    Output();
  }
}

/*---------------------------------------------------------------------------*
 | prepare time step                                          sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::PrepareTimeStep(bool print_header)
{
  // increment time and step
  IncrementTimeAndStep();

  // set current time
  SetCurrentTime();

  // set current step size
  SetCurrentStepSize();

  // print header
  if (print_header) PrintHeader();

  // update result and restart control flags
  writeresultsthisstep_ = (writeresultsevery_ and (Step() % writeresultsevery_ == 0));
  writerestartthisstep_ = (writerestartevery_ and (Step() % writerestartevery_ == 0));
}

/*---------------------------------------------------------------------------*
 | integrate particle time step                               sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::Integrate()
{
  // time integration scheme specific pre-interaction routine
  particletimint_->PreInteractionRoutine();

  // update connectivity
  UpdateConnectivity();

  // set gravity acceleration
  if (particlegravity_) SetGravityAcceleration();

  // evaluate particle interactions
  if (particleinteraction_) particleinteraction_->EvaluateInteractions();

  // apply viscous damping contribution
  if (viscousdamping_) viscousdamping_->ApplyViscousDamping();

  // time integration scheme specific post-interaction routine
  particletimint_->PostInteractionRoutine();
}

/*---------------------------------------------------------------------------*
 | output particle time step                                  sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::Output() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::Output");

  // write result step
  if (writeresultsthisstep_)
  {
    // write particle runtime vtp output
    particleengine_->WriteParticleRuntimeVtpOutput(Step(), Time());

    // write binning discretization output (debug feature)
    particleengine_->WriteBinDisOutput(Step(), Time());

    // write wall runtime vtu output
    if (particlewall_) particlewall_->WriteWallRuntimeVtuOutput(Step(), Time());
  }

  // write restart step
  if (writerestartthisstep_)
  {
    // write restart of particle engine
    particleengine_->WriteRestart(Step(), Time());

    // write restart of particle time integration
    particletimint_->WriteRestart(Step(), Time());

    // write restart of particle interaction handler
    if (particleinteraction_) particleinteraction_->WriteRestart(Step(), Time());

    // write restart of gravity handler
    if (particlegravity_) particlegravity_->WriteRestart(Step(), Time());

    // write restart of viscous damping handler
    if (viscousdamping_) viscousdamping_->WriteRestart(Step(), Time());

    // write restart of wall handler
    if (particlewall_) particlewall_->WriteRestart(Step(), Time());

    // short screen output
    if (myrank_ == 0)
      IO::cout(IO::verbose) << "====== restart of the particle simulation written in step "
                            << Step() << IO::endl;
  }
}

/*---------------------------------------------------------------------------*
 | create result test                                         sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
std::shared_ptr<DRT::ResultTest> PARTICLEALGORITHM::ParticleAlgorithm::CreateResultTest()
{
  // create and init particle result test
  std::shared_ptr<PARTICLEALGORITHM::ParticleResultTest> resulttest =
      std::make_shared<PARTICLEALGORITHM::ParticleResultTest>();
  resulttest->Init();

  // setup particle result test
  resulttest->Setup(particleengine_);

  return resulttest;
}

/*---------------------------------------------------------------------------*
 | init particle engine                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleEngine()
{
  // create and init particle engine
  particleengine_ = std::make_shared<PARTICLEENGINE::ParticleEngine>(Comm(), params_);
  particleengine_->Init();
}

/*---------------------------------------------------------------------------*
 | init particle time integration                             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleTimeIntegration()
{
  // get particle time integration scheme
  INPAR::PARTICLE::DynamicType timinttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(params_, "DYNAMICTYP");

  // create particle time integration
  switch (timinttype)
  {
    case INPAR::PARTICLE::dyna_semiimpliciteuler:
    {
      particletimint_ = std::unique_ptr<PARTICLEALGORITHM::TimIntSemiImplicitEuler>(
          new PARTICLEALGORITHM::TimIntSemiImplicitEuler(params_));
      break;
    }
    case INPAR::PARTICLE::dyna_velocityverlet:
    {
      particletimint_ = std::unique_ptr<PARTICLEALGORITHM::TimIntVelocityVerlet>(
          new PARTICLEALGORITHM::TimIntVelocityVerlet(params_));
      break;
    }
    default:
    {
      dserror("unknown particle time integration scheme!");
      break;
    }
  }

  // init particle time integration
  particletimint_->Init();
}

/*---------------------------------------------------------------------------*
 | init particle interaction handler                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleInteraction()
{
  // get particle interaction type
  INPAR::PARTICLE::InteractionType interactiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::InteractionType>(params_, "INTERACTION");

  // create particle interaction handler
  switch (interactiontype)
  {
    case INPAR::PARTICLE::interaction_none:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::interaction_sph:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionSPH>(
          new PARTICLEINTERACTION::ParticleInteractionSPH(Comm(), params_));
      break;
    }
    case INPAR::PARTICLE::interaction_dem:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionDEM>(
          new PARTICLEINTERACTION::ParticleInteractionDEM(Comm(), params_));
      break;
    }
    default:
    {
      dserror("unknown particle interaction type!");
      break;
    }
  }

  // init particle interaction handler
  if (particleinteraction_) particleinteraction_->Init();
}

/*---------------------------------------------------------------------------*
 | init particle gravity handler                              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleGravity()
{
  // init gravity acceleration vector
  std::vector<double> gravity;
  double value;
  std::istringstream gravitystream(
      Teuchos::getNumericStringParameter(params_, "GRAVITY_ACCELERATION"));

  while (gravitystream >> value) gravity.push_back(value);

  // safety check
  if ((int)gravity.size() != 3)
    dserror("dimension (dim = %d) of gravity acceleration vector is wrong!", (int)gravity.size());

  // get magnitude of gravity
  double temp = 0.0;
  for (int i = 0; i < (int)gravity.size(); ++i) temp += gravity[i] * gravity[i];
  const double gravity_norm = std::sqrt(temp);

  // no gravity defined
  if (gravity_norm == 0.0)
    particlegravity_ = std::unique_ptr<PARTICLEALGORITHM::GravityHandler>(nullptr);
  // create particle gravity handler
  else
    particlegravity_ = std::unique_ptr<PARTICLEALGORITHM::GravityHandler>(
        new PARTICLEALGORITHM::GravityHandler(params_));

  // init particle gravity handler
  if (particlegravity_) particlegravity_->Init(gravity);
}

/*---------------------------------------------------------------------------*
 | init viscous damping handler                               sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitViscousDamping()
{
  // get viscous damping factor
  const double viscdampfac = params_.get<double>("VISCOUS_DAMPING");

  // create viscous damping handler
  if (viscdampfac > 0.0)
    viscousdamping_ = std::unique_ptr<PARTICLEALGORITHM::ViscousDampingHandler>(
        new PARTICLEALGORITHM::ViscousDampingHandler(viscdampfac));
  else
    viscousdamping_ = std::unique_ptr<PARTICLEALGORITHM::ViscousDampingHandler>(nullptr);

  // init viscous damping handler
  if (viscousdamping_) viscousdamping_->Init();
}

/*---------------------------------------------------------------------------*
 | init particle wall handler                                 sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleWall()
{
  // get type of particle wall source
  INPAR::PARTICLE::ParticleWallSource particlewallsource =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleWallSource>(
          params_, "PARTICLE_WALL_SOURCE");

  // create particle wall handler
  switch (particlewallsource)
  {
    case INPAR::PARTICLE::NoParticleWall:
    {
      particlewall_ = std::shared_ptr<PARTICLEALGORITHM::WallHandlerBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::DiscretCondition:
    {
      particlewall_ =
          std::make_shared<PARTICLEALGORITHM::WallHandlerDiscretCondition>(Comm(), params_);
      break;
    }
    case INPAR::PARTICLE::BoundingBox:
    {
      particlewall_ = std::make_shared<PARTICLEALGORITHM::WallHandlerBoundingBox>(Comm(), params_);
      break;
    }
    default:
    {
      dserror("unknown type of particle wall source!");
      break;
    }
  }

  // init particle wall handler
  if (particlewall_) particlewall_->Init();
}

/*---------------------------------------------------------------------------*
 | generate initial particles                                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::GenerateInitialParticles()
{
  // create and init particle input generator
  std::unique_ptr<PARTICLEALGORITHM::InputGenerator> particleinputgenerator =
      std::unique_ptr<PARTICLEALGORITHM::InputGenerator>(
          new PARTICLEALGORITHM::InputGenerator(Comm(), params_));
  particleinputgenerator->Init();

  // generate particles
  particleinputgenerator->GenerateParticles(particlestodistribute_);
}

/*---------------------------------------------------------------------------*
 | determine all particle types                               sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::DetermineParticleTypes()
{
  // init map relating particle types to dynamic load balance factor
  std::map<PARTICLEENGINE::TypeEnum, double> typetodynloadbal;

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_, "PHASE_TO_DYNLOADBALFAC", typetodynloadbal);

  // insert into map of particle types and corresponding states with empty set
  for (auto& typeIt : typetodynloadbal)
    particlestatestotypes_.insert(
        std::make_pair(typeIt.first, std::set<PARTICLEENGINE::StateEnum>()));

  // safety check
  for (auto& particle : particlestodistribute_)
    if (not particlestatestotypes_.count(particle->ReturnParticleType()))
      dserror("particle type '%s' of initial particle not defined!",
          PARTICLEENGINE::EnumToTypeName(particle->ReturnParticleType()).c_str());
}

/*---------------------------------------------------------------------------*
 | determine particle states of all particle types            sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::DetermineParticleStatesOfParticleTypes()
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes_)
  {
    // insert default particle states
    (typeIt.second)
        .insert({PARTICLEENGINE::Position, PARTICLEENGINE::Velocity, PARTICLEENGINE::Acceleration,
            PARTICLEENGINE::LastTransferPosition});
  }

  // insert integration dependent states of all particle types
  particletimint_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);

  // insert interaction dependent states of all particle types
  if (particleinteraction_)
    particleinteraction_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);
}

/*---------------------------------------------------------------------------*
 | setup initial particles                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitialParticles()
{
  // erase particles outside bounding box
  particleengine_->EraseParticlesOutsideBoundingBox(particlestodistribute_);

  // distribute particles to owning processor
  particleengine_->DistributeParticles(particlestodistribute_);

  // dynamic load balancing
  particleengine_->DynamicLoadBalancing();

  // get number of particles on this processor
  numparticlesafterlastloadbalance_ = particleengine_->GetNumberOfParticles();

  // ghost particles on other processors
  particleengine_->GhostParticles();

  // build particle to particle neighbors
  if (particleinteraction_) particleengine_->BuildParticleToParticleNeighbors();

  // distribute interaction history
  if (particleinteraction_) particleinteraction_->DistributeInteractionHistory();

  // build global id to local index map
  particleengine_->BuildGlobalIDToLocalIndexMap();
}

/*---------------------------------------------------------------------------*
 | setup initial wall distribution                            sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitalWallDistribution()
{
  // update bin row and column map
  particlewall_->UpdateBinRowAndColMap(
      particleengine_->GetBinRowMap(), particleengine_->GetBinColMap());

  // distribute wall elements and nodes
  particlewall_->DistributeWallElementsAndNodes();

  // relate bins to column wall elements
  particlewall_->RelateBinsToColWallEles();

  // build particle to wall neighbors
  if (particleinteraction_)
    particlewall_->BuildParticleToWallNeighbors(particleengine_->GetParticlesToBins());
}

/*---------------------------------------------------------------------------*
 | setup initial states                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitialStates()
{
  // set initial states
  if (particleinteraction_) particleinteraction_->SetInitialStates();

  // set initial conditions
  SetInitialConditions();

  // time integration scheme specific initialization routine
  particletimint_->SetInitialStates();

  // update connectivity
  UpdateConnectivity();

  // set gravity acceleration
  if (particlegravity_) SetGravityAcceleration();

  // evaluate particle interactions
  if (particleinteraction_) particleinteraction_->EvaluateInteractions();
}

/*---------------------------------------------------------------------------*
 | check particle types for modified states                   sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEALGORITHM::ParticleAlgorithm::HaveModifiedStates()
{
  bool havemodifiedstates = false;

  int numtypes = 0;
  int numtypeswithmodifiedstates = 0;

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes_)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // no modified velocity and acceleration for boundary or rigid particles
    if (particleType == PARTICLEENGINE::BoundaryPhase or particleType == PARTICLEENGINE::RigidPhase)
      continue;

    // get reference to set of particle states of current type
    std::set<PARTICLEENGINE::StateEnum>& stateEnumSet = typeIt.second;

    // check for modified velocity and acceleration of current type
    bool modifiedvelocity = stateEnumSet.count(PARTICLEENGINE::ModifiedVelocity);
    bool modifiedacceleration = stateEnumSet.count(PARTICLEENGINE::ModifiedAcceleration);

    // safety check
    if (modifiedvelocity != modifiedacceleration)
      dserror("modified states of both velocity and acceleration need to be set!");

    // increase counter of type
    ++numtypes;

    // modified states for current type
    if (modifiedvelocity) ++numtypeswithmodifiedstates;
  }

  if (numtypeswithmodifiedstates > 0)
  {
    // safety check
    if (numtypes != numtypeswithmodifiedstates)
      dserror(
          "modified states need to be set for all particle types except boundary phase particles!");

    havemodifiedstates = true;
  }

  return havemodifiedstates;
}

/*---------------------------------------------------------------------------*
 | update connectivity                                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::UpdateConnectivity()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::UpdateConnectivity");

  // check particle transfer including bin size safety checks
  bool transferneeded = CheckParticleTransfer();

  // check for valid particle connectivity
  bool invalidparticleconnectivity = (not particleengine_->HaveValidParticleConnectivity());

  if (transferevery_ or transferneeded or invalidparticleconnectivity or writeresultsthisstep_ or
      writerestartthisstep_)
  {
    // transfer particles to new bins and processors
    particleengine_->TransferParticles();

    // communicate interaction history
    if (particleinteraction_) particleinteraction_->CommunicateInteractionHistory();

    // check load balancing
    bool loadbalanceneeded = CheckLoadBalancing();

    if (loadbalanceneeded or writerestartthisstep_)
    {
      // dynamic load balancing
      particleengine_->DynamicLoadBalancing();

      // communicate interaction history
      if (particleinteraction_) particleinteraction_->CommunicateInteractionHistory();

      // get number of particles on this processor
      numparticlesafterlastloadbalance_ = particleengine_->GetNumberOfParticles();

      if (particlewall_)
      {
        // update bin row and column map
        particlewall_->UpdateBinRowAndColMap(
            particleengine_->GetBinRowMap(), particleengine_->GetBinColMap());

        // distribute wall elements and nodes
        particlewall_->DistributeWallElementsAndNodes();

        // relate bins to column wall elements
        particlewall_->RelateBinsToColWallEles();
      }
    }

    // ghost particles on other processors
    particleengine_->GhostParticles();

    // build particle to particle neighbors
    if (particleinteraction_) particleengine_->BuildParticleToParticleNeighbors();

    // build particle to wall neighbors
    if (particleinteraction_ and particlewall_)
      particlewall_->BuildParticleToWallNeighbors(particleengine_->GetParticlesToBins());

    // build global id to local index map
    particleengine_->BuildGlobalIDToLocalIndexMap();

    // short screen output
    if (myrank_ == 0)
    {
      if (loadbalanceneeded or writerestartthisstep_)
        IO::cout(IO::verbose) << "dynamic load balancing and particle transfer in step " << Step()
                              << IO::endl;
      else
        IO::cout(IO::verbose) << "particle transfer in step " << Step() << IO::endl;
    }
  }
  else
  {
    // refresh particles being ghosted on other processors
    particleengine_->RefreshParticles();
  }
}

/*---------------------------------------------------------------------------*
 | check particle transfer including bin size safety checks   sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEALGORITHM::ParticleAlgorithm::CheckParticleTransfer()
{
  // get maximum particle interaction distance
  double allprocinteractiondistance = 0.0;
  if (particleinteraction_)
  {
    double interactiondistance = particleinteraction_->MaxInteractionDistance();
    Comm().MaxAll(&interactiondistance, &allprocinteractiondistance, 1);
  }

  // get minimum relevant bin size
  const double minbinsize = particleengine_->MinBinSize();

  // bin size safety check
  if (allprocinteractiondistance > minbinsize)
    dserror("the particle interaction distance is larger than the minimal bin size (%f > %f)!",
        allprocinteractiondistance, minbinsize);

  // loop over all spatial directions
  for (int dim = 0; dim < 3; ++dim)
  {
    // check for periodic boundary condition in current spatial direction
    if (particleengine_->HavePBC(dim))
    {
      // check periodic length in current spatial direction
      if ((2.0 * allprocinteractiondistance) > particleengine_->PBCDelta(dim))
        dserror("particles are not allowed to interact directly and across the periodic boundary!");
    }
  }

  // maximum position increment since last particle transfer
  double maxpositionincrement = 0.0;

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengine_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored == 0) continue;

    // get pointer to particle states
    const double* pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);
    const double* lasttransferpos =
        container->GetPtrToParticleState(PARTICLEENGINE::LastTransferPosition, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // iterate over coordinate values of owned particles of current type
    for (int i = 0; i < (statedim * particlestored); ++i)
    {
      // get position increment of particle in current spatial dimension since last transfer
      double absolutpositionincrement = std::abs(pos[i] - lasttransferpos[i]);

      // compare to current maximum
      maxpositionincrement = std::max(maxpositionincrement, absolutpositionincrement);
    }
  }

  // get maximum position increment on all processors
  double allprocmaxpositionincrement = 0.0;
  Comm().MaxAll(&maxpositionincrement, &allprocmaxpositionincrement, 1);

  // check if a particle transfer is needed: it is assumed that (in a worst case scenario)
  // two particles approach each other with maximum position increment in one spatial dimension
  bool transferneeded =
      ((allprocinteractiondistance + 2.0 * allprocmaxpositionincrement) > minbinsize);

  // safety check
  if (transferneeded and maxpositionincrement > minbinsize)
    dserror("a particle traveled more than one bin on this processor!");

  return transferneeded;
}

/*---------------------------------------------------------------------------*
 | check load balancing                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEALGORITHM::ParticleAlgorithm::CheckLoadBalancing()
{
  // percentage limit
  const double percentagelimit = 0.1;

  // get number of particles on this processor
  int numberofparticles = particleengine_->GetNumberOfParticles();

  // percentage change of particles on this processor
  double percentagechange = 0.0;
  if (numparticlesafterlastloadbalance_ > 0)
    percentagechange =
        std::abs(static_cast<double>(numberofparticles - numparticlesafterlastloadbalance_) /
                 numparticlesafterlastloadbalance_);

  // get maximum percentage change of particles
  double maxpercentagechange = 0.0;
  Comm().MaxAll(&percentagechange, &maxpercentagechange, 1);

  // check if dynamic load balance is needed based on the maximum percentage change of the number of
  // particles on one processor
  bool loadbalanceneeded = (maxpercentagechange > percentagelimit);

  return loadbalanceneeded;
}

/*---------------------------------------------------------------------------*
 | set initial conditions                                     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetInitialConditions()
{
  // create and init particle initial field handler
  std::unique_ptr<PARTICLEALGORITHM::InitialFieldHandler> initialfield =
      std::unique_ptr<PARTICLEALGORITHM::InitialFieldHandler>(
          new PARTICLEALGORITHM::InitialFieldHandler(params_));
  initialfield->Init();

  // setup particle initial field handler
  initialfield->Setup(particleengine_);

  // set initial fields
  initialfield->SetInitialFields();
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetCurrentTime()
{
  // set current time in particle time integration
  particletimint_->SetCurrentTime(Time());

  // set current time in particle interaction
  if (particleinteraction_) particleinteraction_->SetCurrentTime(Time());
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetCurrentStepSize()
{
  // set current step size in particle interaction
  if (particleinteraction_) particleinteraction_->SetCurrentStepSize(Dt());
}

/*---------------------------------------------------------------------------*
 | set gravity acceleration                                   sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ParticleAlgorithm::SetGravityAcceleration()
{
  std::vector<double> scaled_gravity(3);

  // get gravity acceleration
  particlegravity_->GetGravityAcceleration(Time(), scaled_gravity);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengine_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
  {
    // gravity is not set for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // set gravity acceleration for all particles of current type
    particlecontainerbundle->SetStateSpecificContainer(
        scaled_gravity, PARTICLEENGINE::Acceleration, typeEnum);
  }

  // set scaled gravity in particle interaction handler
  if (particleinteraction_) particleinteraction_->SetGravity(scaled_gravity);
}
