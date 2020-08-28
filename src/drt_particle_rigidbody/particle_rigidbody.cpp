/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body handler for particle problem
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_rigidbody.H"

#include "particle_rigidbody_datastate.H"
#include "particle_rigidbody_affiliation_pairs.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_unique_global_id.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLERIGIDBODY::RigidBodyHandler::RigidBodyHandler(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), myrank_(comm.MyPID()), params_(params)
{
  // empty constructor
}

PARTICLERIGIDBODY::RigidBodyHandler::~RigidBodyHandler() = default;

void PARTICLERIGIDBODY::RigidBodyHandler::Init()
{
  // init rigid body unique global identifier handler
  InitRigidBodyUniqueGlobalIdHandler();

  // init rigid body data state container
  InitRigidBodyDataState();

  // init affiliation pair handler
  InitAffiliationPairHandler();
}

void PARTICLERIGIDBODY::RigidBodyHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // setup unique global identifier handler
  rigidbodyuniqueglobalidhandler_->Setup();

  // setup rigid body data state container
  rigidbodydatastate_->Setup();

  // setup affiliation pair handler
  affiliationpairs_->Setup(particleengineinterface);

  // safety check
  {
    // get particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
        particleengineinterface_->GetParticleContainerBundle();

    if (not particlecontainerbundle->GetParticleTypes().count(PARTICLEENGINE::RigidPhase))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(PARTICLEENGINE::RigidPhase).c_str());
  }

  // short screen output
  if (particleengineinterface_->HavePeriodicBoundaryConditions() and myrank_ == 0)
    IO::cout << "Warning: rigid bodies not transferred over periodic boundary!" << IO::endl;
}

void PARTICLERIGIDBODY::RigidBodyHandler::WriteRestart() const
{
  // get bin discretization writer
  std::shared_ptr<IO::DiscretizationWriter> binwriter =
      particleengineinterface_->GetBinDiscretizationWriter();

  // write restart of unique global identifier handler
  rigidbodyuniqueglobalidhandler_->WriteRestart(binwriter);

  // write restart of affiliation pair handler
  affiliationpairs_->WriteRestart();
}

void PARTICLERIGIDBODY::RigidBodyHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of unique global identifier handler
  rigidbodyuniqueglobalidhandler_->ReadRestart(reader);

  // read restart of affiliation pair handler
  affiliationpairs_->ReadRestart(reader);
}

void PARTICLERIGIDBODY::RigidBodyHandler::DistributeRigidBody()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::DistributeRigidBody");

  // distribute affiliation pairs
  affiliationpairs_->DistributeAffiliationPairs();
}

void PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody");

  // communicate affiliation pairs
  affiliationpairs_->CommunicateAffiliationPairs();
}

void PARTICLERIGIDBODY::RigidBodyHandler::InitRigidBodyUniqueGlobalIdHandler()
{
  // create and init unique global identifier handler
  rigidbodyuniqueglobalidhandler_ = std::unique_ptr<PARTICLEENGINE::UniqueGlobalIdHandler>(
      new PARTICLEENGINE::UniqueGlobalIdHandler(comm_, "rigidbody"));
  rigidbodyuniqueglobalidhandler_->Init();
}

void PARTICLERIGIDBODY::RigidBodyHandler::InitRigidBodyDataState()
{
  // create rigid body data state container
  rigidbodydatastate_ = std::make_shared<PARTICLERIGIDBODY::RigidBodyDataState>();

  // init rigid body data state container
  rigidbodydatastate_->Init();
}

void PARTICLERIGIDBODY::RigidBodyHandler::InitAffiliationPairHandler()
{
  // create affiliation pair handler
  affiliationpairs_ = std::unique_ptr<PARTICLERIGIDBODY::RigidBodyAffiliationPairs>(
      new PARTICLERIGIDBODY::RigidBodyAffiliationPairs(comm_));

  // init affiliation pair handler
  affiliationpairs_->Init();
}
