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
#include "particle_rigidbody_runtime_vtp_writer.H"
#include "particle_rigidbody_affiliation_pairs.H"
#include "particle_rigidbody_utils.H"

#include "../drt_particle_interaction/particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_communication_utils.H"
#include "../drt_particle_engine/particle_unique_global_id.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_pack_buffer.H"
#include "../drt_lib/drt_parobject.H"

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

  // init rigid body runtime vtp writer
  InitRigidBodyVtpWriter();

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

  // setup rigid body runtime vtp writer
  SetupRigidBodyVtpWriter();

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

  // get packed rigid body state data
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);
  GetPackedRigidBodyStates(*buffer);

  // write rigid body state data
  binwriter->WriteCharVector("RigidBodyStateData", buffer);
}

void PARTICLERIGIDBODY::RigidBodyHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of unique global identifier handler
  rigidbodyuniqueglobalidhandler_->ReadRestart(reader);

  // read restart of runtime vtp writer
  rigidbodyvtpwriter_->ReadRestart(reader);

  // read restart of affiliation pair handler
  affiliationpairs_->ReadRestart(reader);

  // allocate rigid body states
  AllocateRigidBodyStates();

  // read rigid body state data
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);
  reader->ReadCharVector(buffer, "RigidBodyStateData");

  // extract packed rigid body state data
  ExtractPackedRigidBodyStates(*buffer);
}

void PARTICLERIGIDBODY::RigidBodyHandler::InsertParticleStatesOfParticleTypes(
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

    if (type == PARTICLEENGINE::RigidPhase)
    {
      // insert states of rigid particles
      particlestates.insert(
          {PARTICLEENGINE::RigidBodyColor, PARTICLEENGINE::ReferenceRelativePosition,
              PARTICLEENGINE::RelativePosition, PARTICLEENGINE::Inertia});
    }
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::WriteRigidBodyRuntimeOutput(
    const int step, const double time) const
{
  rigidbodyvtpwriter_->ResetTimeAndTimeStep(time, step);
  rigidbodyvtpwriter_->SetRigidBodyPositionsAndStates(ownedrigidbodies_);
  rigidbodyvtpwriter_->WriteFiles();
  rigidbodyvtpwriter_->WriteCollectionFileOfAllWrittenFiles();
}

void PARTICLERIGIDBODY::RigidBodyHandler::SetUniqueGlobalIdsForAllRigidBodies()
{
  // get reference to affiliation pair data
  std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

  // maximum global id of rigid bodies on this processor
  int maxglobalid = -1;

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    // get pointer to particle states
    const double* rigidbodycolor_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::RigidBodyColor, particle_i);

    // get global id of affiliated rigid body k
    const int rigidbody_k = std::round(rigidbodycolor_i[0]);

    // insert affiliation pair
    affiliationpairdata.insert(std::make_pair(globalid_i[0], rigidbody_k));

    // get maximum global id of rigid bodies on this processor
    maxglobalid = std::max(maxglobalid, rigidbody_k);
  }

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // get maximum global id of rigid bodies on all processors
  int allprocmaxglobalid = -1;
  comm_.MaxAll(&maxglobalid, &allprocmaxglobalid, 1);

  // number of global ids on all processors
  const int numglobalids = allprocmaxglobalid + 1;

#ifdef DEBUG
  if (not(rigidbodyuniqueglobalidhandler_->GetMaxGlobalId() < 0))
    dserror("maximum global id of rigid body unique global identifier handler already touched!");
#endif

  // request number of global ids of all rigid bodies on processor 0
  std::vector<int> requesteduniqueglobalids;
  if (myrank_ == 0) requesteduniqueglobalids.reserve(numglobalids);

  // draw requested number of global ids
  rigidbodyuniqueglobalidhandler_->DrawRequestedNumberOfGlobalIds(requesteduniqueglobalids);

#ifdef DEBUG
  if (myrank_ == 0)
    for (int i = 0; i < numglobalids; ++i)
      if (requesteduniqueglobalids[i] != i) dserror("drawn requested global ids not consecutive!");
#endif

  // used global ids on all processors
  std::vector<int> usedglobalids(numglobalids, 0);

  // get used global ids on this processor
  for (const auto& it : affiliationpairdata) usedglobalids[it.second] = 1;

  // mpi communicator
  const Epetra_MpiComm* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm_);
  if (!mpicomm) dserror("dynamic cast to Epetra_MpiComm failed!");

  // get used global ids on all processors
  MPI_Allreduce(MPI_IN_PLACE, &usedglobalids[0], numglobalids, MPI_INT, MPI_MAX, mpicomm->Comm());

  // free unused global ids on processor 0
  if (myrank_ == 0)
    for (int i = 0; i < numglobalids; ++i)
      if (usedglobalids[i] == 0)
        rigidbodyuniqueglobalidhandler_->InsertFreedGlobalId(requesteduniqueglobalids[i]);
}

void PARTICLERIGIDBODY::RigidBodyHandler::AllocateRigidBodyStates()
{
  // number of global ids
  const int numglobalids = rigidbodyuniqueglobalidhandler_->GetMaxGlobalId() + 1;

  // allocate stored states
  rigidbodydatastate_->AllocateStoredStates(numglobalids);
}

void PARTICLERIGIDBODY::RigidBodyHandler::DistributeRigidBody()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::DistributeRigidBody");

  // distribute affiliation pairs
  affiliationpairs_->DistributeAffiliationPairs();

  // store rigid bodies previously owned by this processor
  std::vector<int> previouslyownedrigidbodies = ownedrigidbodies_;

  // update rigid body ownership
  UpdateRigidBodyOwnership();

  // relate owned rigid bodies to all hosting processors
  RelateOwnedRigidBodiesToHostingProcs();

  // communicate rigid body states
  CommunicateRigidBodyStates(previouslyownedrigidbodies);
}

void PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody");

  // communicate affiliation pairs
  affiliationpairs_->CommunicateAffiliationPairs();

  // store rigid bodies previously owned by this processor
  std::vector<int> previouslyownedrigidbodies = ownedrigidbodies_;

  // update rigid body ownership
  UpdateRigidBodyOwnership();

  // relate owned rigid bodies to all hosting processors
  RelateOwnedRigidBodiesToHostingProcs();

  // communicate rigid body states
  CommunicateRigidBodyStates(previouslyownedrigidbodies);
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputeMassQuantities()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::ComputeMassQuantities");

  // clear partial mass quantities of rigid bodies
  ClearPartialMassQuantities();

  // compute partial mass quantities of rigid bodies
  ComputePartialMassQuantities();

  // gathered partial mass quantities of rigid bodies from all corresponding processors
  std::unordered_map<int, std::vector<double>> gatheredpartialmass;
  std::unordered_map<int, std::vector<std::vector<double>>> gatheredpartialinertia;
  std::unordered_map<int, std::vector<std::vector<double>>> gatheredpartialposition;

  // gather partial mass quantities of rigid bodies
  GatherPartialMassQuantities(gatheredpartialmass, gatheredpartialinertia, gatheredpartialposition);

  // compute full mass quantities of rigid bodies
  ComputeFullMassQuantities(gatheredpartialmass, gatheredpartialinertia, gatheredpartialposition);
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputeAccelerations()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::ComputeAccelerations");

  // clear force and torque acting on rigid bodies
  ClearForceAndTorque();

  // compute partial force and torque acting on rigid bodies
  ComputePartialForceAndTorque();

  // gather partial and compute full force and torque acting on rigid bodies
  GatherPartialAndComputeFullForceAndTorque();

  // compute accelerations of rigid bodies from force and torque
  ComputeAccelerationsFromForceAndTorque();
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdatePositions(const double timeincrement)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::UpdatePositions");

  // update positions of rigid bodies with given time increment
  UpdateRigidBodyPositions(timeincrement);

  // broadcast positions of rigid bodies
  BroadcastPositions();

  // update relative position of rigid particles
  UpdateRigidParticleRelativePosition();

  // set position of rigid particles
  SetRigidParticlePosition();
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdateVelocities(const double timeincrement)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::UpdateVelocities");

  // update velocities of rigid bodies with given time increment
  UpdateRigidBodyVelocities(timeincrement);

  // broadcast velocities of rigid bodies
  BroadcastVelocities();

  // set velocities of rigid particles
  SetRigidParticleVelocities();
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

void PARTICLERIGIDBODY::RigidBodyHandler::InitRigidBodyVtpWriter()
{
  // construct and init rigid body runtime vtp writer
  rigidbodyvtpwriter_ = std::unique_ptr<PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter>(
      new PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter(comm_));
  rigidbodyvtpwriter_->Init(rigidbodydatastate_);
}

void PARTICLERIGIDBODY::RigidBodyHandler::InitAffiliationPairHandler()
{
  // create affiliation pair handler
  affiliationpairs_ = std::unique_ptr<PARTICLERIGIDBODY::RigidBodyAffiliationPairs>(
      new PARTICLERIGIDBODY::RigidBodyAffiliationPairs(comm_));

  // init affiliation pair handler
  affiliationpairs_->Init();
}

void PARTICLERIGIDBODY::RigidBodyHandler::SetupRigidBodyVtpWriter()
{
  // get data format for written numeric data via vtp
  bool write_binary_output = (DRT::INPUT::IntegralValue<INPAR::PARTICLE::OutputDataFormat>(
                                  params_, "OUTPUT_DATA_FORMAT") == INPAR::PARTICLE::binary);

  // setup rigid body runtime vtp writer
  rigidbodyvtpwriter_->Setup(write_binary_output);
}

void PARTICLERIGIDBODY::RigidBodyHandler::GetPackedRigidBodyStates(std::vector<char>& buffer) const
{
  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get reference to rigid body states
    const double& mass_k = rigidbodydatastate_->GetRefMass()[rigidbody_k];
    const std::vector<double>& inertia_k = rigidbodydatastate_->GetRefInertia()[rigidbody_k];
    const std::vector<double>& pos_k = rigidbodydatastate_->GetRefPosition()[rigidbody_k];
    const std::vector<double>& rot_k = rigidbodydatastate_->GetRefRotation()[rigidbody_k];
    const std::vector<double>& vel_k = rigidbodydatastate_->GetRefVelocity()[rigidbody_k];
    const std::vector<double>& angvel_k = rigidbodydatastate_->GetRefAngularVelocity()[rigidbody_k];
    const std::vector<double>& acc_k = rigidbodydatastate_->GetRefAcceleration()[rigidbody_k];
    const std::vector<double>& angacc_k =
        rigidbodydatastate_->GetRefAngularAcceleration()[rigidbody_k];

    // pack data for sending
    DRT::PackBuffer data;
    data.StartPacking();

    data.AddtoPack(rigidbody_k);
    data.AddtoPack(mass_k);
    for (int i = 0; i < 6; ++i) data.AddtoPack(inertia_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(pos_k[i]);
    for (int i = 0; i < 4; ++i) data.AddtoPack(rot_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(vel_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(angvel_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(acc_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(angacc_k[i]);

    buffer.insert(buffer.end(), data().begin(), data().end());
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ExtractPackedRigidBodyStates(std::vector<char>& buffer)
{
  std::vector<char>::size_type position = 0;

  while (position < buffer.size())
  {
    const int rigidbody_k = DRT::ParObject::ExtractInt(position, buffer);

    // get global ids of rigid bodies owned by this processor
    ownedrigidbodies_.push_back(rigidbody_k);

    // get reference to rigid body states
    double& mass_k = rigidbodydatastate_->GetRefMutableMass()[rigidbody_k];
    std::vector<double>& inertia_k = rigidbodydatastate_->GetRefMutableInertia()[rigidbody_k];
    std::vector<double>& pos_k = rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k];
    std::vector<double>& rot_k = rigidbodydatastate_->GetRefMutableRotation()[rigidbody_k];
    std::vector<double>& vel_k = rigidbodydatastate_->GetRefMutableVelocity()[rigidbody_k];
    std::vector<double>& angvel_k =
        rigidbodydatastate_->GetRefMutableAngularVelocity()[rigidbody_k];
    std::vector<double>& acc_k = rigidbodydatastate_->GetRefMutableAcceleration()[rigidbody_k];
    std::vector<double>& angacc_k =
        rigidbodydatastate_->GetRefMutableAngularAcceleration()[rigidbody_k];

    DRT::ParObject::ExtractfromPack(position, buffer, mass_k);
    for (int i = 0; i < 6; ++i) DRT::ParObject::ExtractfromPack(position, buffer, inertia_k[i]);
    for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, buffer, pos_k[i]);
    for (int i = 0; i < 4; ++i) DRT::ParObject::ExtractfromPack(position, buffer, rot_k[i]);
    for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, buffer, vel_k[i]);
    for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, buffer, angvel_k[i]);
    for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, buffer, acc_k[i]);
    for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, buffer, angacc_k[i]);
  }

  if (position != buffer.size())
    dserror("mismatch in size of data %d <-> %d", static_cast<int>(buffer.size()), position);
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdateRigidBodyOwnership()
{
  ownedrigidbodies_.clear();
  hostedrigidbodies_.clear();
  ownerofrigidbodies_.clear();

  // number of global ids
  const int numglobalids = rigidbodyuniqueglobalidhandler_->GetMaxGlobalId() + 1;

  // maximum number of particles per rigid body over all processors
  std::vector<std::pair<int, int>> maxnumberofparticlesperrigidbodyonproc(
      numglobalids, std::make_pair(0, myrank_));

  // get number of particle per rigid body on this processor
  for (const auto& it : affiliationpairs_->GetRefToAffiliationPairData())
    maxnumberofparticlesperrigidbodyonproc[it.second].first++;

  // get global ids of rigid bodies hosted (owned and non-owned) by this processor
  for (int rigidbody_k = 0; rigidbody_k < numglobalids; ++rigidbody_k)
    if (maxnumberofparticlesperrigidbodyonproc[rigidbody_k].first > 0)
      hostedrigidbodies_.push_back(rigidbody_k);

  // mpi communicator
  const Epetra_MpiComm* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm_);
  if (!mpicomm) dserror("dynamic cast to Epetra_MpiComm failed!");

  // get maximum number of particles per rigid body over all processors
  MPI_Allreduce(MPI_IN_PLACE, &maxnumberofparticlesperrigidbodyonproc[0], numglobalids, MPI_2INT,
      MPI_MAXLOC, mpicomm->Comm());

  // get owner of all rigid bodies
  ownerofrigidbodies_.reserve(numglobalids);
  for (const auto& it : maxnumberofparticlesperrigidbodyonproc)
    ownerofrigidbodies_.push_back(it.second);

  // get global ids of rigid bodies owned by this processor
  for (const int rigidbody_k : hostedrigidbodies_)
    if (ownerofrigidbodies_[rigidbody_k] == myrank_) ownedrigidbodies_.push_back(rigidbody_k);
}

void PARTICLERIGIDBODY::RigidBodyHandler::RelateOwnedRigidBodiesToHostingProcs()
{
  // number of global ids
  const int numglobalids = rigidbodyuniqueglobalidhandler_->GetMaxGlobalId() + 1;

  // allocate memory
  ownedrigidbodiestohostingprocs_.assign(numglobalids, std::vector<int>(0));

  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // owner of rigid body k
    const int owner_k = ownerofrigidbodies_[rigidbody_k];

    // communicate global id of rigid body to owning processor
    if (owner_k != myrank_)
    {
      // pack data for sending
      DRT::PackBuffer data;
      data.StartPacking();

      data.AddtoPack(rigidbody_k);

      sdata[owner_k].insert(sdata[owner_k].end(), data().begin(), data().end());
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    int msgsource = p.first;
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);

      // insert processor id the gathered global id of rigid body is received from
      ownedrigidbodiestohostingprocs_[rigidbody_k].push_back(msgsource);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBodyStates(
    std::vector<int>& previouslyownedrigidbodies)
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over previously owned rigid bodies
  for (const int rigidbody_k : previouslyownedrigidbodies)
  {
    // owner of rigid body k
    const int owner_k = ownerofrigidbodies_[rigidbody_k];

    // get reference to rigid body states
    const double& mass_k = rigidbodydatastate_->GetRefMass()[rigidbody_k];
    const std::vector<double>& inertia_k = rigidbodydatastate_->GetRefInertia()[rigidbody_k];
    const std::vector<double>& pos_k = rigidbodydatastate_->GetRefPosition()[rigidbody_k];
    const std::vector<double>& rot_k = rigidbodydatastate_->GetRefRotation()[rigidbody_k];
    const std::vector<double>& vel_k = rigidbodydatastate_->GetRefVelocity()[rigidbody_k];
    const std::vector<double>& angvel_k = rigidbodydatastate_->GetRefAngularVelocity()[rigidbody_k];
    const std::vector<double>& acc_k = rigidbodydatastate_->GetRefAcceleration()[rigidbody_k];
    const std::vector<double>& angacc_k =
        rigidbodydatastate_->GetRefAngularAcceleration()[rigidbody_k];

    // communicate states to owning processor
    if (owner_k != myrank_)
    {
      // pack data for sending
      DRT::PackBuffer data;
      data.StartPacking();

      data.AddtoPack(rigidbody_k);
      data.AddtoPack(mass_k);
      for (int i = 0; i < 6; ++i) data.AddtoPack(inertia_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(pos_k[i]);
      for (int i = 0; i < 4; ++i) data.AddtoPack(rot_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(vel_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(angvel_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(acc_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(angacc_k[i]);

      sdata[owner_k].insert(sdata[owner_k].end(), data().begin(), data().end());
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);

      // get reference to rigid body states
      double& mass_k = rigidbodydatastate_->GetRefMutableMass()[rigidbody_k];
      std::vector<double>& inertia_k = rigidbodydatastate_->GetRefMutableInertia()[rigidbody_k];
      std::vector<double>& pos_k = rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k];
      std::vector<double>& rot_k = rigidbodydatastate_->GetRefMutableRotation()[rigidbody_k];
      std::vector<double>& vel_k = rigidbodydatastate_->GetRefMutableVelocity()[rigidbody_k];
      std::vector<double>& angvel_k =
          rigidbodydatastate_->GetRefMutableAngularVelocity()[rigidbody_k];
      std::vector<double>& acc_k = rigidbodydatastate_->GetRefMutableAcceleration()[rigidbody_k];
      std::vector<double>& angacc_k =
          rigidbodydatastate_->GetRefMutableAngularAcceleration()[rigidbody_k];

      DRT::ParObject::ExtractfromPack(position, rmsg, mass_k);
      for (int i = 0; i < 6; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, inertia_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, pos_k[i]);
      for (int i = 0; i < 4; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, rot_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, vel_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, angvel_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, acc_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, angacc_k[i]);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ClearPartialMassQuantities()
{
  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // get pointer to rigid body states
    double* mass_k = &rigidbodydatastate_->GetRefMutableMass()[rigidbody_k];
    double* inertia_k = &rigidbodydatastate_->GetRefMutableInertia()[rigidbody_k][0];
    double* pos_k = &rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k][0];

    // clear mass quantities
    mass_k[0] = 0.0;
    for (int i = 0; i < 6; ++i) inertia_k[i] = 0.0;
    PARTICLEINTERACTION::UTILS::vec_clear(pos_k);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputePartialMassQuantities()
{
  // get reference to affiliation pair data
  const std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    double* mass_k = &rigidbodydatastate_->GetRefMutableMass()[rigidbody_k];
    double* pos_k = &rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k][0];

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);

    // sum contribution of particle i
    mass_k[0] += mass_i[0];
    PARTICLEINTERACTION::UTILS::vec_addscale(pos_k, mass_i[0], pos_i);
  }

  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // get pointer to rigid body states
    const double* mass_k = &rigidbodydatastate_->GetRefMass()[rigidbody_k];
    double* pos_k = &rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k][0];

#ifdef DEBUG
    if (not(mass_k[0] > 0.0)) dserror("partial mass of rigid body %d is zero!", rigidbody_k);
#endif

    // determine center of gravity of (partial) rigid body k
    PARTICLEINTERACTION::UTILS::vec_scale(pos_k, 1.0 / mass_k[0]);
  }

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    const double* pos_k = &rigidbodydatastate_->GetRefPosition()[rigidbody_k][0];
    double* inertia_k = &rigidbodydatastate_->GetRefMutableInertia()[rigidbody_k][0];

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    const double* inertia_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::Inertia, particle_i);

    double r_ki[3];
    PARTICLEINTERACTION::UTILS::vec_set(r_ki, pos_k);
    PARTICLEINTERACTION::UTILS::vec_sub(r_ki, pos_i);

    // sum contribution of particle i
    inertia_k[0] += inertia_i[0] + (r_ki[1] * r_ki[1] + r_ki[2] * r_ki[2]) * mass_i[0];
    inertia_k[1] += inertia_i[0] + (r_ki[0] * r_ki[0] + r_ki[2] * r_ki[2]) * mass_i[0];
    inertia_k[2] += inertia_i[0] + (r_ki[0] * r_ki[0] + r_ki[1] * r_ki[1]) * mass_i[0];
    inertia_k[3] -= r_ki[0] * r_ki[1] * mass_i[0];
    inertia_k[4] -= r_ki[0] * r_ki[2] * mass_i[0];
    inertia_k[5] -= r_ki[1] * r_ki[2] * mass_i[0];
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::GatherPartialMassQuantities(
    std::unordered_map<int, std::vector<double>>& gatheredpartialmass,
    std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialinertia,
    std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialposition)
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // owner of rigid body k
    const int owner_k = ownerofrigidbodies_[rigidbody_k];

    // get reference to rigid body states
    const double& mass_k = rigidbodydatastate_->GetRefMass()[rigidbody_k];
    const std::vector<double>& inertia_k = rigidbodydatastate_->GetRefInertia()[rigidbody_k];
    const std::vector<double>& pos_k = rigidbodydatastate_->GetRefPosition()[rigidbody_k];

    // rigid body k owned by this processor
    if (owner_k == myrank_)
    {
      // append to gathered partial mass quantities
      gatheredpartialmass[rigidbody_k].push_back(mass_k);
      gatheredpartialinertia[rigidbody_k].push_back(inertia_k);
      gatheredpartialposition[rigidbody_k].push_back(pos_k);
    }
    // communicate partial mass quantities to owning processor
    else
    {
      // pack data for sending
      DRT::PackBuffer data;
      data.StartPacking();

      data.AddtoPack(rigidbody_k);
      data.AddtoPack(mass_k);
      for (int i = 0; i < 6; ++i) data.AddtoPack(inertia_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(pos_k[i]);

      sdata[owner_k].insert(sdata[owner_k].end(), data().begin(), data().end());
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);
      double mass_k = DRT::ParObject::ExtractDouble(position, rmsg);

      std::vector<double> inertia_k(6);
      for (int i = 0; i < 6; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, inertia_k[i]);

      std::vector<double> pos_k(3);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, pos_k[i]);

      // append to gathered partial mass quantities
      gatheredpartialmass[rigidbody_k].push_back(mass_k);
      gatheredpartialinertia[rigidbody_k].push_back(inertia_k);
      gatheredpartialposition[rigidbody_k].push_back(pos_k);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputeFullMassQuantities(
    std::unordered_map<int, std::vector<double>>& gatheredpartialmass,
    std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialinertia,
    std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialposition)
{
  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    std::vector<double>& partialmass_k = gatheredpartialmass[rigidbody_k];
    std::vector<std::vector<double>>& partialpos_k = gatheredpartialposition[rigidbody_k];

    // number of partial mass quantities of rigid body k including this processor
    const int numpartial_k = ownedrigidbodiestohostingprocs_[rigidbody_k].size() + 1;

#ifdef DEBUG
    if (static_cast<int>(partialmass_k.size()) != numpartial_k or
        static_cast<int>(partialpos_k.size()) != numpartial_k)
      dserror("the number of partial mass quantities of rigid body %d do not match!", rigidbody_k);
#endif

    // get pointer to rigid body states
    double* mass_k = &rigidbodydatastate_->GetRefMutableMass()[rigidbody_k];
    double* pos_k = &rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k][0];

    // clear mass and position
    mass_k[0] = 0.0;
    PARTICLEINTERACTION::UTILS::vec_clear(pos_k);

    // iterate over partial quantities
    for (int p = 0; p < numpartial_k; ++p)
    {
      // sum contribution of partial quantity
      mass_k[0] += partialmass_k[p];
      PARTICLEINTERACTION::UTILS::vec_addscale(pos_k, partialmass_k[p], &partialpos_k[p][0]);
    }

    // determine center of gravity of rigid body k
    PARTICLEINTERACTION::UTILS::vec_scale(pos_k, 1.0 / mass_k[0]);
  }

  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    std::vector<double>& partialmass_k = gatheredpartialmass[rigidbody_k];
    std::vector<std::vector<double>>& partialpos_k = gatheredpartialposition[rigidbody_k];
    std::vector<std::vector<double>>& partialinertia_k = gatheredpartialinertia[rigidbody_k];

    // number of partial mass quantities of rigid body k including this processor
    const int numpartial_k = ownedrigidbodiestohostingprocs_[rigidbody_k].size() + 1;

#ifdef DEBUG
    if (static_cast<int>(partialmass_k.size()) != numpartial_k or
        static_cast<int>(partialinertia_k.size()) != numpartial_k)
      dserror("the number of partial mass quantities of rigid body %d do not match!", rigidbody_k);
#endif

    // get pointer to rigid body states
    const double* pos_k = &rigidbodydatastate_->GetRefPosition()[rigidbody_k][0];
    double* inertia_k = &rigidbodydatastate_->GetRefMutableInertia()[rigidbody_k][0];

    // clear inertia
    for (int i = 0; i < 6; ++i) inertia_k[i] = 0.0;

    // iterate over partial quantities
    for (int p = 0; p < numpartial_k; ++p)
    {
      double r_kp[3];
      PARTICLEINTERACTION::UTILS::vec_set(r_kp, pos_k);
      PARTICLEINTERACTION::UTILS::vec_sub(r_kp, &partialpos_k[p][0]);

      // sum contribution of partial quantity
      for (int i = 0; i < 6; ++i) inertia_k[i] += partialinertia_k[p][i];

      inertia_k[0] += (r_kp[1] * r_kp[1] + r_kp[2] * r_kp[2]) * partialmass_k[p];
      inertia_k[1] += (r_kp[0] * r_kp[0] + r_kp[2] * r_kp[2]) * partialmass_k[p];
      inertia_k[2] += (r_kp[0] * r_kp[0] + r_kp[1] * r_kp[1]) * partialmass_k[p];
      inertia_k[3] -= r_kp[0] * r_kp[1] * partialmass_k[p];
      inertia_k[4] -= r_kp[0] * r_kp[2] * partialmass_k[p];
      inertia_k[5] -= r_kp[1] * r_kp[2] * partialmass_k[p];
    }
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ClearForceAndTorque()
{
  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // get pointer to rigid body states
    double* force_k = &rigidbodydatastate_->GetRefMutableForce()[rigidbody_k][0];
    double* torque_k = &rigidbodydatastate_->GetRefMutableTorque()[rigidbody_k][0];

    // clear force and torque of rigid body k
    PARTICLEINTERACTION::UTILS::vec_clear(force_k);
    PARTICLEINTERACTION::UTILS::vec_clear(torque_k);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputePartialForceAndTorque()
{
  // get reference to affiliation pair data
  const std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    double* force_k = &rigidbodydatastate_->GetRefMutableForce()[rigidbody_k][0];
    double* torque_k = &rigidbodydatastate_->GetRefMutableTorque()[rigidbody_k][0];

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* relpos_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::RelativePosition, particle_i);
    const double* acc_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

    // compute force of particle i
    double force_i[3];
    PARTICLEINTERACTION::UTILS::vec_setscale(force_i, mass_i[0], acc_i);

    // sum contribution of particle i
    PARTICLEINTERACTION::UTILS::vec_add(force_k, force_i);
    PARTICLEINTERACTION::UTILS::vec_addcross(torque_k, relpos_i, force_i);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::GatherPartialAndComputeFullForceAndTorque()
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over hosted rigid bodies
  for (const int rigidbody_k : hostedrigidbodies_)
  {
    // owner of rigid body k
    const int owner_k = ownerofrigidbodies_[rigidbody_k];

    // get reference to rigid body states
    const std::vector<double>& force_k = rigidbodydatastate_->GetRefForce()[rigidbody_k];
    const std::vector<double>& torque_k = rigidbodydatastate_->GetRefTorque()[rigidbody_k];

    // communicate partial force and torque to owning processor
    if (owner_k != myrank_)
    {
      // pack data for sending
      DRT::PackBuffer data;
      data.StartPacking();

      data.AddtoPack(rigidbody_k);
      for (int i = 0; i < 3; ++i) data.AddtoPack(force_k[i]);
      for (int i = 0; i < 3; ++i) data.AddtoPack(torque_k[i]);

      sdata[owner_k].insert(sdata[owner_k].end(), data().begin(), data().end());
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);

      std::vector<double> tmp_force_k(3);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, tmp_force_k[i]);

      std::vector<double> tmp_torque_k(3);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, tmp_torque_k[i]);

      // get pointer to rigid body states
      double* force_k = &rigidbodydatastate_->GetRefMutableForce()[rigidbody_k][0];
      double* torque_k = &rigidbodydatastate_->GetRefMutableTorque()[rigidbody_k][0];

      // sum gathered contribution to full force and torque
      PARTICLEINTERACTION::UTILS::vec_add(force_k, &tmp_force_k[0]);
      PARTICLEINTERACTION::UTILS::vec_add(torque_k, &tmp_torque_k[0]);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::ComputeAccelerationsFromForceAndTorque()
{
  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get pointer to rigid body states
    const double* mass_k = &rigidbodydatastate_->GetRefMass()[rigidbody_k];
    const double* inertia_k = &rigidbodydatastate_->GetRefInertia()[rigidbody_k][0];
    const double* rot_k = &rigidbodydatastate_->GetRefRotation()[rigidbody_k][0];
    const double* force_k = &rigidbodydatastate_->GetRefForce()[rigidbody_k][0];
    const double* torque_k = &rigidbodydatastate_->GetRefTorque()[rigidbody_k][0];
    double* acc_k = &rigidbodydatastate_->GetRefMutableAcceleration()[rigidbody_k][0];
    double* angacc_k = &rigidbodydatastate_->GetRefMutableAngularAcceleration()[rigidbody_k][0];

    // compute acceleration of rigid body k
    PARTICLEINTERACTION::UTILS::vec_setscale(acc_k, 1 / mass_k[0], force_k);

    // compute inverse of rotation
    double invrot_k[4];
    UTILS::quaternion_invert(invrot_k, rot_k);

    // get torque in the reference frame
    double reftorque_k[3];
    UTILS::quaternion_rotate_vector(reftorque_k, invrot_k, torque_k);

    // determinant of mass moment of inertia
    const double det_inertia_k =
        inertia_k[0] * inertia_k[1] * inertia_k[2] + inertia_k[3] * inertia_k[4] * inertia_k[5] +
        inertia_k[3] * inertia_k[4] * inertia_k[5] - inertia_k[1] * inertia_k[4] * inertia_k[4] -
        inertia_k[2] * inertia_k[3] * inertia_k[3] - inertia_k[0] * inertia_k[5] * inertia_k[5];

    // no mass moment of inertia
    if (std::abs(det_inertia_k) < 1E-14) continue;

    // evaluate angular acceleration of rigid body k in the reference frame
    double refangacc_k[3];
    refangacc_k[0] = reftorque_k[0] * (inertia_k[1] * inertia_k[2] - inertia_k[5] * inertia_k[5]) +
                     reftorque_k[1] * (inertia_k[4] * inertia_k[5] - inertia_k[2] * inertia_k[3]) +
                     reftorque_k[2] * (inertia_k[3] * inertia_k[5] - inertia_k[1] * inertia_k[4]);

    refangacc_k[1] = reftorque_k[0] * (inertia_k[4] * inertia_k[5] - inertia_k[2] * inertia_k[3]) +
                     reftorque_k[1] * (inertia_k[0] * inertia_k[2] - inertia_k[4] * inertia_k[4]) +
                     reftorque_k[2] * (inertia_k[3] * inertia_k[4] - inertia_k[0] * inertia_k[5]);

    refangacc_k[2] = reftorque_k[0] * (inertia_k[3] * inertia_k[5] - inertia_k[1] * inertia_k[4]) +
                     reftorque_k[1] * (inertia_k[3] * inertia_k[4] - inertia_k[0] * inertia_k[5]) +
                     reftorque_k[2] * (inertia_k[0] * inertia_k[1] - inertia_k[3] * inertia_k[3]);

    PARTICLEINTERACTION::UTILS::vec_scale(refangacc_k, 1.0 / det_inertia_k);

    // compute angular acceleration of rigid body k in the rotating frame
    UTILS::quaternion_rotate_vector(angacc_k, rot_k, refangacc_k);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdateRigidBodyPositions(const double timeincrement)
{
  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get pointer to rigid body states
    double* pos_k = &rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k][0];
    double* rot_k = &rigidbodydatastate_->GetRefMutableRotation()[rigidbody_k][0];
    const double* vel_k = &rigidbodydatastate_->GetRefVelocity()[rigidbody_k][0];
    const double* angvel_k = &rigidbodydatastate_->GetRefAngularVelocity()[rigidbody_k][0];

    // update position
    PARTICLEINTERACTION::UTILS::vec_addscale(pos_k, timeincrement, vel_k);

    // save current rotation
    double curr_rot_k[4];
    UTILS::quaternion_set(curr_rot_k, rot_k);

    // get rotation increment
    double phi_k[3];
    PARTICLEINTERACTION::UTILS::vec_setscale(phi_k, timeincrement, angvel_k);

    double incr_rot_k[4];
    UTILS::quaternion_from_angle(incr_rot_k, phi_k);

    // update rotation
    UTILS::quaternion_product(rot_k, incr_rot_k, curr_rot_k);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdateRigidBodyVelocities(const double timeincrement)
{
  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get pointer to rigid body states
    double* vel_k = &rigidbodydatastate_->GetRefMutableVelocity()[rigidbody_k][0];
    double* angvel_k = &rigidbodydatastate_->GetRefMutableAngularVelocity()[rigidbody_k][0];
    const double* acc_k = &rigidbodydatastate_->GetRefAcceleration()[rigidbody_k][0];
    const double* angacc_k = &rigidbodydatastate_->GetRefAngularAcceleration()[rigidbody_k][0];

    // update velocities
    PARTICLEINTERACTION::UTILS::vec_addscale(vel_k, timeincrement, acc_k);
    PARTICLEINTERACTION::UTILS::vec_addscale(angvel_k, timeincrement, angacc_k);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::BroadcastPositions()
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get reference to hosting processors of rigid body k
    std::vector<int>& hostingprocs_k = ownedrigidbodiestohostingprocs_[rigidbody_k];

    // get reference to rigid body states
    const std::vector<double>& pos_k = rigidbodydatastate_->GetRefPosition()[rigidbody_k];
    const std::vector<double>& rot_k = rigidbodydatastate_->GetRefRotation()[rigidbody_k];

    // pack data for sending
    DRT::PackBuffer data;
    data.StartPacking();

    data.AddtoPack(rigidbody_k);

    for (int i = 0; i < 3; ++i) data.AddtoPack(pos_k[i]);
    for (int i = 0; i < 4; ++i) data.AddtoPack(rot_k[i]);

    for (int torank : hostingprocs_k)
      sdata[torank].insert(sdata[torank].end(), data().begin(), data().end());
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);

      // get reference to rigid body states
      std::vector<double>& pos_k = rigidbodydatastate_->GetRefMutablePosition()[rigidbody_k];
      std::vector<double>& rot_k = rigidbodydatastate_->GetRefMutableRotation()[rigidbody_k];

      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, pos_k[i]);
      for (int i = 0; i < 4; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, rot_k[i]);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::BroadcastVelocities()
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // iterate over owned rigid bodies
  for (const int rigidbody_k : ownedrigidbodies_)
  {
    // get reference to hosting processors of rigid body k
    std::vector<int>& hostingprocs_k = ownedrigidbodiestohostingprocs_[rigidbody_k];

    // get reference to rigid body states
    const std::vector<double>& vel_k = rigidbodydatastate_->GetRefVelocity()[rigidbody_k];
    const std::vector<double>& angvel_k = rigidbodydatastate_->GetRefAngularVelocity()[rigidbody_k];

    // pack data for sending
    DRT::PackBuffer data;
    data.StartPacking();

    data.AddtoPack(rigidbody_k);

    for (int i = 0; i < 3; ++i) data.AddtoPack(vel_k[i]);
    for (int i = 0; i < 3; ++i) data.AddtoPack(angvel_k[i]);

    for (int torank : hostingprocs_k)
      sdata[torank].insert(sdata[torank].end(), data().begin(), data().end());
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      const int rigidbody_k = DRT::ParObject::ExtractInt(position, rmsg);

      // get reference to rigid body states
      std::vector<double>& vel_k = rigidbodydatastate_->GetRefMutableVelocity()[rigidbody_k];
      std::vector<double>& angvel_k =
          rigidbodydatastate_->GetRefMutableAngularVelocity()[rigidbody_k];

      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, vel_k[i]);
      for (int i = 0; i < 3; ++i) DRT::ParObject::ExtractfromPack(position, rmsg, angvel_k[i]);
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::UpdateRigidParticleRelativePosition()
{
  // get reference to affiliation pair data
  const std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    const double* rot_k = &rigidbodydatastate_->GetRefRotation()[rigidbody_k][0];

    // get pointer to particle states
    const double* refrelpos_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::ReferenceRelativePosition, particle_i);
    double* relpos_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::RelativePosition, particle_i);

    // update relative position of particle i
    UTILS::quaternion_rotate_vector(relpos_i, rot_k, refrelpos_i);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::SetRigidParticlePosition()
{
  // get reference to affiliation pair data
  const std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    const double* pos_k = &rigidbodydatastate_->GetRefPosition()[rigidbody_k][0];

    // get pointer to particle states
    const double* relpos_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::RelativePosition, particle_i);
    double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);

    // set position of particle i
    PARTICLEINTERACTION::UTILS::vec_set(pos_i, pos_k);
    PARTICLEINTERACTION::UTILS::vec_add(pos_i, relpos_i);
  }
}

void PARTICLERIGIDBODY::RigidBodyHandler::SetRigidParticleVelocities()
{
  // get reference to affiliation pair data
  const std::unordered_map<int, int>& affiliationpairdata =
      affiliationpairs_->GetRefToAffiliationPairData();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get container of owned particles of rigid phase
  PARTICLEENGINE::ParticleContainer* container_i = particlecontainerbundle->GetSpecificContainer(
      PARTICLEENGINE::RigidPhase, PARTICLEENGINE::Owned);

#ifdef DEBUG
  if (static_cast<int>(affiliationpairdata.size()) != container_i->ParticlesStored())
    dserror("number of affiliation pairs and rigid particles not equal!");
#endif

  // loop over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    auto it = affiliationpairdata.find(globalid_i[0]);

#ifdef DEBUG
    // no affiliation pair for current global id
    if (it == affiliationpairdata.end())
      dserror("no affiliated rigid body found for particle with global id %d", globalid_i[0]);
#endif

    // get global id of affiliated rigid body k
    const int rigidbody_k = it->second;

    // get pointer to rigid body states
    const double* vel_k = &rigidbodydatastate_->GetRefVelocity()[rigidbody_k][0];
    const double* angvel_k = &rigidbodydatastate_->GetRefAngularVelocity()[rigidbody_k][0];

    // get pointer to particle states
    const double* relpos_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::RelativePosition, particle_i);
    double* vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    double* angvel_i = nullptr;
    if (container_i->HaveStoredState(PARTICLEENGINE::AngularVelocity))
      angvel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::AngularVelocity, particle_i);

    // set velocities of particle i
    PARTICLEINTERACTION::UTILS::vec_set(vel_i, vel_k);
    PARTICLEINTERACTION::UTILS::vec_addcross(vel_i, angvel_k, relpos_i);
    if (angvel_i) PARTICLEINTERACTION::UTILS::vec_set(angvel_i, angvel_k);
  }
}
