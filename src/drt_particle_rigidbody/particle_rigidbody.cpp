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
#include "../drt_particle_engine/particle_communication_utils.H"
#include "../drt_particle_engine/particle_unique_global_id.H"

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

  // update rigid body ownership
  UpdateRigidBodyOwnership();

  // relate owned rigid bodies to all hosting processors
  RelateOwnedRigidBodiesToHostingProcs();
}

void PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLERIGIDBODY::RigidBodyHandler::CommunicateRigidBody");

  // communicate affiliation pairs
  affiliationpairs_->CommunicateAffiliationPairs();

  // update rigid body ownership
  UpdateRigidBodyOwnership();

  // relate owned rigid bodies to all hosting processors
  RelateOwnedRigidBodiesToHostingProcs();
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
