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

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_unique_global_id.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

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
}

void PARTICLERIGIDBODY::RigidBodyHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // setup unique global identifier handler
  rigidbodyuniqueglobalidhandler_->Setup();

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
}

void PARTICLERIGIDBODY::RigidBodyHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of unique global identifier handler
  rigidbodyuniqueglobalidhandler_->ReadRestart(reader);
}

void PARTICLERIGIDBODY::RigidBodyHandler::InitRigidBodyUniqueGlobalIdHandler()
{
  // create and init unique global identifier handler
  rigidbodyuniqueglobalidhandler_ = std::unique_ptr<PARTICLEENGINE::UniqueGlobalIdHandler>(
      new PARTICLEENGINE::UniqueGlobalIdHandler(comm_, "rigidbody"));
  rigidbodyuniqueglobalidhandler_->Init();
}
