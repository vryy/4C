/*---------------------------------------------------------------------------*/
/*! \file
\brief particle engine to control particle problem
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_engine.H"

#include "particle_communication_utils.H"

#include "particle_container.H"
#include "particle_container_bundle.H"
#include "particle_object.H"
#include "particle_unique_global_id.H"
#include "particle_runtime_vtp_writer.H"

#include "../drt_particle_algorithm/particle_algorithm_utils.H"

#include "../drt_binstrategy/binning_strategy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleEngine::ParticleEngine(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm),
      myrank_(comm.MyPID()),
      params_(params),
      minbinsize_(0.0),
      typevectorsize_(0),
      validownedparticles_(false),
      validghostedparticles_(false),
      validparticleneighbors_(false),
      validglobalidtolocalindex_(false),
      validdirectghosting_(false),
      validhalfneighboringbins_(false)
{
  // empty constructor
}

PARTICLEENGINE::ParticleEngine::~ParticleEngine() = default;

void PARTICLEENGINE::ParticleEngine::Init()
{
  // init binning strategy
  InitBinningStrategy();

  // init particle container bundle
  InitParticleContainerBundle();

  // init particle unique global identifier handler
  InitParticleUniqueGlobalIdHandler();

  // init particle runtime vtp writer
  InitParticleVtpWriter();
}

void PARTICLEENGINE::ParticleEngine::Setup(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes)
{
  // setup binning strategy
  SetupBinningStrategy();

  // setup particle container bundle
  SetupParticleContainerBundle(particlestatestotypes);

  // setup particle unique global identifier handler
  SetupParticleUniqueGlobalIdHandler();

  // setup data storage
  SetupDataStorage(particlestatestotypes);

  // setup particle runtime vtp writer
  SetupParticleVtpWriter();

  // setup particle type weights for dynamic load balancing
  SetupTypeWeights();
}

void PARTICLEENGINE::ParticleEngine::WriteRestart(const int step, const double time) const
{
  // get bin discretization writer
  std::shared_ptr<IO::DiscretizationWriter> binwriter =
      Teuchos::get_shared_ptr(binstrategy_->BinDiscret()->Writer());

  binwriter->NewStep(step, time);

  // pack particles of all containers
  std::shared_ptr<std::vector<char>> particlebuffer = std::make_shared<std::vector<char>>();
  particlecontainerbundle_->GetPackedParticleObjectsOfAllContainers(particlebuffer);

  // write particle data
  binwriter->WriteCharVector("ParticleData", Teuchos::rcp(particlebuffer));

  // write restart of unique global identifier handler
  particleuniqueglobalidhandler_->WriteRestart(binwriter);
}

void PARTICLEENGINE::ParticleEngine::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader,
    std::vector<ParticleObjShrdPtr>& particlestoread) const
{
  // read particle data
  Teuchos::RCP<std::vector<char>> particledata = Teuchos::rcp(new std::vector<char>);
  reader->ReadCharVector(particledata, "ParticleData");

  std::vector<char>::size_type position = 0;

  while (position < particledata->size())
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(position, *particledata, data);

    // this std::shared_ptr holds the memory
    std::shared_ptr<DRT::ParObject> object(DRT::UTILS::Factory(data));
    ParticleObjShrdPtr particleobject =
        std::dynamic_pointer_cast<PARTICLEENGINE::ParticleObject>(object);
    if (particleobject == nullptr) dserror("received object is not a particle object!");

    // store read particle
    particlestoread.push_back(particleobject);
  }

  if (position != particledata->size())
    dserror("mismatch in size of data %d <-> %d", static_cast<int>(particledata->size()), position);

  // read restart of unique global identifier handler
  particleuniqueglobalidhandler_->ReadRestart(reader);

  // read restart of runtime vtp writer
  particlevtpwriter_->ReadRestart(reader);
}

void PARTICLEENGINE::ParticleEngine::WriteParticleRuntimeOutput(
    const int step, const double time) const
{
  particlevtpwriter_->ResetTimeAndTimeStep(time, step);
  particlevtpwriter_->SetParticlePositionsAndStates();
  particlevtpwriter_->WriteFiles();
  particlevtpwriter_->WriteCollectionFileOfAllWrittenFiles();
}

void PARTICLEENGINE::ParticleEngine::FreeUniqueGlobalIds(std::vector<int>& freeuniquegids)
{
  // insert freed global id
  for (const int currfreegid : freeuniquegids)
    particleuniqueglobalidhandler_->InsertFreedGlobalId(currfreegid);

  // clear after all unique global ids are freed
  freeuniquegids.clear();
}

void PARTICLEENGINE::ParticleEngine::GetUniqueGlobalIdsForAllParticles(
    std::vector<ParticleObjShrdPtr>& particlestogetuniquegids)
{
  std::vector<int> requesteduniqueglobalids;
  requesteduniqueglobalids.reserve(particlestogetuniquegids.size());

  // draw requested number of global ids
  particleuniqueglobalidhandler_->DrawRequestedNumberOfGlobalIds(requesteduniqueglobalids);

  // number of particles to get unique ids
  int numparticles = particlestogetuniquegids.size();

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
    particlestogetuniquegids[i]->SetParticleGlobalID(requesteduniqueglobalids[i]);
}

void PARTICLEENGINE::ParticleEngine::CheckNumberOfUniqueGlobalIds()
{
  // mpi communicator
  const auto* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm_);
  if (!mpicomm) dserror("dynamic cast to Epetra_MpiComm failed!");

  // get number of particles on all processors
  int numberofparticles = GetNumberOfParticles();
  MPI_Allreduce(MPI_IN_PLACE, &numberofparticles, 1, MPI_INT, MPI_SUM, mpicomm->Comm());

  // get number of reusable global ids on all processors
  int numberofreusableglobalids = particleuniqueglobalidhandler_->GetNumberOfReusableGlobalIds();
  MPI_Allreduce(MPI_IN_PLACE, &numberofreusableglobalids, 1, MPI_INT, MPI_SUM, mpicomm->Comm());

  // maximum global id
  const int maxglobalid = particleuniqueglobalidhandler_->GetMaxGlobalId();

  // safety check
  if (numberofparticles + numberofreusableglobalids != (maxglobalid + 1))
    dserror("sum of particles and reusable global ids unequal total global ids: %d + %d != %d",
        numberofparticles, numberofreusableglobalids, (maxglobalid + 1));
}

void PARTICLEENGINE::ParticleEngine::EraseParticlesOutsideBoundingBox(
    std::vector<ParticleObjShrdPtr>& particlestocheck)
{
  // get bounding box dimensions
  LINALG::Matrix<3, 2> boundingbox = binstrategy_->DomainBoundingBoxCornerPositions();

  // set of particles located outside bounding box
  std::set<int> particlesoutsideboundingbox;

  // number of particles to check
  int numparticles = particlestocheck.size();

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
  {
    // get states of particle
    const ParticleStates& particleStates = particlestocheck[i]->ReturnParticleStates();

    // get position of particle
    const std::vector<double>& pos = particleStates[PARTICLEENGINE::Position];

#ifdef DEBUG
    // get type of particles
    TypeEnum typeEnum = particlestocheck[i]->ReturnParticleType();

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    if (static_cast<int>(pos.size()) != container->GetParticleStateDim(PARTICLEENGINE::Position))
      dserror("dimension of particle state '%s' not valid!",
          PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());
#endif

    // check particle location with respect to bounding box in each spatial directions
    for (int dim = 0; dim < 3; ++dim)
    {
      // particle located outside bounding box
      if ((pos[dim] < boundingbox(dim, 0)) or (pos[dim] > boundingbox(dim, 1)))
      {
        // insert particle into set
        particlesoutsideboundingbox.insert(i);

#ifdef DEBUG
        if (particlestocheck[i]->ReturnParticleGlobalID() < 0)
          dserror("no global id assigned to particle!");
#endif

        // insert freed global id
        particleuniqueglobalidhandler_->InsertFreedGlobalId(
            particlestocheck[i]->ReturnParticleGlobalID());

        break;
      }
    }
  }

  // number of particles located outside bounding box
  const int numparticlesoutside = particlesoutsideboundingbox.size();

  // no particles located outside of bounding box
  if (numparticlesoutside == 0) return;

  // put particles to be erased at the end of the vector
  int swapposition = numparticles - 1;

  // iterate in reversed order over particles to be erased
  std::set<int>::reverse_iterator rit;
  for (rit = particlesoutsideboundingbox.rbegin(); rit != particlesoutsideboundingbox.rend(); ++rit)
  {
    // swap position of particle to be erase
    std::swap(particlestocheck[*rit], particlestocheck[swapposition]);

    // set new swap position
    --swapposition;
  }

  // erase particles located outside bounding box
  particlestocheck.resize(numparticles - numparticlesoutside);

  // short screen output
  if (numparticlesoutside)
    IO::cout << "on processor " << myrank_ << " removed " << numparticlesoutside
             << " particle(s) being outside the computational domain!" << IO::endl;
}

void PARTICLEENGINE::ParticleEngine::DistributeParticles(
    std::vector<ParticleObjShrdPtr>& particlestodistribute)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::DistributeParticles");

  std::vector<std::vector<ParticleObjShrdPtr>> particlestosend(comm_.NumProc());
  std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert(typevectorsize_);

  // determine particles that need to be distributed
  DetermineParticlesToBeDistributed(particlestodistribute, particlestosend, particlestoinsert);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert owned particles received from other processors
  InsertOwnedParticles(particlestoinsert);

  // store particle positions after transfer of particles
  StorePositionsAfterParticleTransfer();

  // relate owned particles to bins
  RelateOwnedParticlesToBins();
}

void PARTICLEENGINE::ParticleEngine::TransferParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::TransferParticles");

  std::vector<std::set<int>> particlestoremove(typevectorsize_);
  std::vector<std::vector<ParticleObjShrdPtr>> particlestosend(comm_.NumProc());
  std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert(typevectorsize_);

  // relate owned particles to bins
  if (not validownedparticles_) RelateOwnedParticlesToBins();

  // check particles for periodic boundaries/leaving domain
  CheckParticlesAtBoundaries(particlestoremove);

  // determine particles that need to be transfered
  DetermineParticlesToBeTransfered(particlestoremove, particlestosend);

  // remove particles from containers
  RemoveParticlesFromContainers(particlestoremove);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert owned particles received from other processors
  InsertOwnedParticles(particlestoinsert);

  // store particle positions after transfer of particles
  StorePositionsAfterParticleTransfer();

  // relate owned particles to bins
  RelateOwnedParticlesToBins();

  // check and decrease the size of all containers of owned particles
  particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
      PARTICLEENGINE::Owned);
}

void PARTICLEENGINE::ParticleEngine::GhostParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::GhostParticles");

  std::vector<std::vector<ParticleObjShrdPtr>> particlestosend(comm_.NumProc());
  std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert(typevectorsize_);
  std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>> directghosting;

  // clear all containers of ghosted particles
  particlecontainerbundle_->ClearAllContainersOfSpecificStatus(PARTICLEENGINE::Ghosted);

  // determine particles that need to be ghosted
  DetermineParticlesToBeGhosted(particlestosend);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert ghosted particles received from other processors
  InsertGhostedParticles(particlestoinsert, directghosting);

  // communicate and build map for direct ghosting
  CommunicateDirectGhostingMap(directghosting);

  // check and decrease the size of all containers of ghosted particles
  particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
      PARTICLEENGINE::Ghosted);
}

void PARTICLEENGINE::ParticleEngine::RefreshParticles() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::RefreshParticles");

  std::vector<std::vector<ParticleObjShrdPtr>> particlestosend(comm_.NumProc());
  std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert(typevectorsize_);

  // determine particles that need to be refreshed
  DetermineParticlesToBeRefreshed(particlestosend);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert refreshed particles received from other processors
  InsertRefreshedParticles(particlestoinsert);
}

void PARTICLEENGINE::ParticleEngine::RefreshParticlesOfSpecificStatesAndTypes(
    const StatesOfTypesToRefresh& particlestatestotypes) const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEENGINE::ParticleEngine::RefreshParticlesOfSpecificStatesAndTypes");

  std::vector<std::vector<ParticleObjShrdPtr>> particlestosend(comm_.NumProc());
  std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert(typevectorsize_);

  // determine particles that need to be refreshed
  DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
      particlestatestotypes, particlestosend);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert refreshed particles received from other processors
  InsertRefreshedParticles(particlestoinsert);
}

void PARTICLEENGINE::ParticleEngine::DynamicLoadBalancing()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::DynamicLoadBalancing");

  // determine bin weights needed for repartitioning
  DetermineBinWeights();

  // distribute bins via recursive coordinate bisection
  binstrategy_->DistributeBinsRecursCoordBisection(binrowmap_, bincenters_, binweights_);

  // export elements to new layout
  binstrategy_->BinDiscret()->ExportRowElements(*binrowmap_);

  // setup ghosting of bins
  SetupBinGhosting();

  // determine bin distribution dependent maps/sets
  DetermineBinDisDependentMapsAndSets();

  // determine ghosting dependent maps/sets for communication
  DetermineGhostingDependentMapsAndSets();

  // prepare storage for particle objects
  std::vector<ParticleObjShrdPtr> particlestodistribute;
  particlestodistribute.reserve(GetNumberOfParticles());

  // get vector of particle objects of all containers
  particlecontainerbundle_->GetVectorOfParticleObjectsOfAllContainers(particlestodistribute);

  // clear all containers of owned particles
  particlecontainerbundle_->ClearAllContainersOfSpecificStatus(PARTICLEENGINE::Owned);

  // invalidate particle safety flags
  InvalidateParticleSafetyFlags();

  // invalidate flag denoting valid relation of half surrounding neighboring bins to owned bins
  validhalfneighboringbins_ = false;

  // distribute particles to owning processor
  DistributeParticles(particlestodistribute);

  // check and decrease the size of all containers of owned particles
  particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
      PARTICLEENGINE::Owned);
}

void PARTICLEENGINE::ParticleEngine::HandOverParticlesToBeRemoved(
    std::vector<std::set<int>>& particlestoremove)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::HandOverParticlesToBeRemoved");

  // number of particles to remove
  int numparticlestoremove = 0;
  for (auto typeIt : particlestoremove) numparticlestoremove += typeIt.size();

  if (numparticlestoremove)
  {
    // remove particles from containers
    RemoveParticlesFromContainers(particlestoremove);

    // check and decrease the size of all containers of owned particles
    particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
        PARTICLEENGINE::Owned);
  }
}

void PARTICLEENGINE::ParticleEngine::HandOverParticlesToBeInserted(
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::HandOverParticlesToBeInserted");

  // number of particles to insert
  int numparticlestoinsert = 0;
  for (auto typeIt : particlestoinsert) numparticlestoinsert += typeIt.size();

  if (numparticlestoinsert)
  {
    // insert owned particles
    InsertOwnedParticles(particlestoinsert);
  }
}

void PARTICLEENGINE::ParticleEngine::BuildParticleToParticleNeighbors()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::BuildParticleToParticleNeighbors");

  // safety check
  if ((not validownedparticles_) or (not validghostedparticles_))
    dserror("invalid relation of particles to bins!");

  // relate half neighboring bins to owned bins
  if (not validhalfneighboringbins_) RelateHalfNeighboringBinsToOwnedBins();

  // clear potential particle neighbors
  potentialparticleneighbors_.clear();

  // invalidate flag denoting validity of particle neighbors map
  validparticleneighbors_ = false;

  // loop over row bins
  for (int rowlidofbin = 0; rowlidofbin < binrowmap_->NumMyElements(); ++rowlidofbin)
  {
    // get global id of bin
    const int gidofbin = binrowmap_->GID(rowlidofbin);

    // get local id of bin
    const int collidofbin = bincolmap_->LID(gidofbin);

    // check if current bin contains owned particles
    if (particlestobins_[collidofbin].empty()) continue;

    // iterate over owned particles in current bin
    for (auto& particleIt : particlestobins_[collidofbin])
    {
      // get type of particle
      TypeEnum typeEnum = particleIt.first;

      // get local index of owned particle
      const int ownedindex = particleIt.second;

      // get container of owned particles of current particle type
      ParticleContainer* container =
          particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

      // get global id of particle
      const int* currglobalid = container->GetPtrToParticleGlobalID(ownedindex);

      // get position of particle
      const double* currpos =
          container->GetPtrToParticleState(PARTICLEENGINE::Position, ownedindex);

      // iterate over neighboring bins (including current bin)
      for (int gidofneighborbin : halfneighboringbinstobins_[rowlidofbin])
      {
        // get local id of neighboring bin
        const int collidofneighboringbin = bincolmap_->LID(gidofneighborbin);

        // check if current neighboring bin contains particles
        if (particlestobins_[collidofneighboringbin].empty()) continue;

        // get status of neighboring particles
        StatusEnum neighborStatusEnum = (binrowmap_->LID(gidofneighborbin) < 0)
                                            ? PARTICLEENGINE::Ghosted
                                            : PARTICLEENGINE::Owned;

        // iterate over particles in current neighboring bin
        for (auto& neighborParticleIt : particlestobins_[collidofneighboringbin])
        {
          // get type of neighboring particle
          TypeEnum neighborTypeEnum = neighborParticleIt.first;

          // get local index of neighboring particle
          const int neighborindex = neighborParticleIt.second;

          // get container of neighboring particle of current particle type
          ParticleContainer* neighborcontainer =
              particlecontainerbundle_->GetSpecificContainer(neighborTypeEnum, neighborStatusEnum);

          // get global id of neighboring particle
          const int* neighborglobalid = neighborcontainer->GetPtrToParticleGlobalID(neighborindex);

          // avoid duplicate neighbor pairs and self-neighboring
          if (gidofbin == gidofneighborbin and neighborglobalid[0] <= currglobalid[0]) continue;

          // get position of neighboring particle
          const double* neighborpos =
              neighborcontainer->GetPtrToParticleState(PARTICLEENGINE::Position, neighborindex);

          // distance vector from owned particle to neighboring particle
          double dist[3];

          // distance between particles considering periodic boundaries
          DistanceBetweenParticles(currpos, neighborpos, dist);

          // distance between particles larger than minimum bin size
          if (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2] >
              (minbinsize_ * minbinsize_))
            continue;

          // append potential particle neighbor pair
          potentialparticleneighbors_.push_back(
              std::make_pair(std::make_tuple(typeEnum, PARTICLEENGINE::Owned, ownedindex),
                  std::make_tuple(neighborTypeEnum, neighborStatusEnum, neighborindex)));
        }
      }
    }
  }

  // validate flag denoting validity of particle neighbors map
  validparticleneighbors_ = true;
}

void PARTICLEENGINE::ParticleEngine::BuildGlobalIDToLocalIndexMap()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::BuildGlobalIDToLocalIndexMap");

  // clear map
  globalidtolocalindex_.clear();

  // invalidate flag denoting validity of map relating particle global ids to local index
  validglobalidtolocalindex_ = false;

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (auto& statusEnum : {PARTICLEENGINE::Owned, PARTICLEENGINE::Ghosted})
    {
      // get container of current particle type and current status
      ParticleContainer* container =
          particlecontainerbundle_->GetSpecificContainer(typeEnum, statusEnum);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // no particles of current type and current status
      if (particlestored <= 0) continue;

      // get pointer to global id of particles
      int* globalids = container->GetPtrToParticleGlobalID(0);

      // loop over particles in container
      for (int index = 0; index < particlestored; ++index)
      {
        // get global id of particle
        int globalid = globalids[index];

        // add entry to map
        globalidtolocalindex_[globalid] =
            std::make_shared<LocalIndexTuple>(typeEnum, statusEnum, index);
      }
    }
  }

  // validate flag denoting validity of map relating particle global ids to local index
  validglobalidtolocalindex_ = true;
}

bool PARTICLEENGINE::ParticleEngine::HaveValidParticleConnectivity() const
{
  int localcheck = ((validownedparticles_ and validghostedparticles_ and
                     validglobalidtolocalindex_ and validdirectghosting_));

  // check among all processors
  int globalcheck = 0;
  comm_.MinAll(&localcheck, &globalcheck, 1);

  return globalcheck;
}

bool PARTICLEENGINE::ParticleEngine::HaveValidParticleNeighbors() const
{
  int localcheck = validparticleneighbors_;

  // check among all processors
  int globalcheck = 0;
  comm_.MinAll(&localcheck, &globalcheck, 1);

  return globalcheck;
}

const PARTICLEENGINE::ParticlesToBins& PARTICLEENGINE::ParticleEngine::GetParticlesToBins() const
{
  // safety check
  if ((not validownedparticles_) or (not validghostedparticles_))
    dserror("invalid map relating particles to bins!");

  return particlestobins_;
}

const PARTICLEENGINE::PotentialParticleNeighbors&
PARTICLEENGINE::ParticleEngine::GetPotentialParticleNeighbors() const
{
  // safety check
  if (not validparticleneighbors_) dserror("invalid particle neighbors!");

  return potentialparticleneighbors_;
}

const PARTICLEENGINE::LocalIndexTupleShrdPtr
PARTICLEENGINE::ParticleEngine::GetLocalIndexInSpecificContainer(int globalid) const
{
  // safety check
  if (not validglobalidtolocalindex_) dserror("invalid global id to local index map!");

  // get local index of particle in specific container
  auto globalidIt = globalidtolocalindex_.find(globalid);
  if (globalidIt == globalidtolocalindex_.end()) return nullptr;

  return globalidIt->second;
}

const std::shared_ptr<IO::DiscretizationWriter>
PARTICLEENGINE::ParticleEngine::GetBinDiscretizationWriter() const
{
  return Teuchos::get_shared_ptr(binstrategy_->BinDiscret()->Writer());
}

void PARTICLEENGINE::ParticleEngine::RelateAllParticlesToAllProcs(
    std::vector<int>& particlestoproc) const
{
  // global ids on this processor
  std::vector<int> thisprocglobalids;

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no particles of current type and current status
    if (particlestored <= 0) continue;

    // get pointer to global id of particles
    int* globalids = container->GetPtrToParticleGlobalID(0);

    // insert global id of particles
    thisprocglobalids.insert(thisprocglobalids.end(), globalids, globalids + particlestored);
  }

  // get maximum global id on this processor
  int thisprocmaxglobalid = 0;
  if (not thisprocglobalids.empty())
    thisprocmaxglobalid = *std::max_element(thisprocglobalids.begin(), thisprocglobalids.end());

  // get maximum global id on all processors
  int allprocmaxglobalid(0);
  comm_.MaxAll(&thisprocmaxglobalid, &allprocmaxglobalid, 1);

  // resize to hold all particles
  const int vecsize = allprocmaxglobalid + 1;
  particlestoproc.resize(vecsize, -1);

  // relate this processor id to its global ids
  for (int globalid : thisprocglobalids) particlestoproc[globalid] = myrank_;

  // mpi communicator
  const auto* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm_);
  if (!mpicomm) dserror("dynamic cast to Epetra_MpiComm failed!");

  // communicate global ids between all processors
  MPI_Allreduce(MPI_IN_PLACE, &particlestoproc[0], vecsize, MPI_INT, MPI_MAX, mpicomm->Comm());
}

void PARTICLEENGINE::ParticleEngine::GetParticlesWithinRadius(const double* position,
    const double radius, std::vector<LocalIndexTuple>& neighboringparticles) const
{
  // safety check
  if ((not validownedparticles_) or (not validghostedparticles_))
    dserror("invalid relation of particles to bins!");

  // get global id of bin
  const int gidofbin = binstrategy_->ConvertPosToGid(position);

#ifdef DEBUG
  // position outside computational domain
  if (gidofbin == -1) dserror("position outside of computational domain!");

  // bin not owned or ghosted by this processor
  if (bincolmap_->LID(gidofbin) < 0)
    dserror("position not in owned or ghosted bin on this processor!");
#endif

  // get neighboring bins to current bin
  std::vector<int> binvec;
  binstrategy_->GetNeighborAndOwnBinIds(gidofbin, binvec);

  // iterate over neighboring bins
  for (int gidofneighborbin : binvec)
  {
    // get local id of neighboring bin
    const int collidofneighboringbin = bincolmap_->LID(gidofneighborbin);

    // neighboring bin not ghosted by this processor
    if (collidofneighboringbin < 0) continue;

    // check if current neighboring bin contains particles
    if (particlestobins_[collidofneighboringbin].empty()) continue;

    // get status of neighboring particles
    StatusEnum neighborStatusEnum =
        (binrowmap_->LID(gidofneighborbin) < 0) ? PARTICLEENGINE::Ghosted : PARTICLEENGINE::Owned;

    // iterate over particles in current neighboring bin
    for (auto& neighborParticleIt : particlestobins_[collidofneighboringbin])
    {
      // get type of neighboring particle
      TypeEnum neighborTypeEnum = neighborParticleIt.first;

      // get local index of neighboring particle
      const int neighborindex = neighborParticleIt.second;

      // get container of neighboring particle of current particle type
      ParticleContainer* neighborcontainer =
          particlecontainerbundle_->GetSpecificContainer(neighborTypeEnum, neighborStatusEnum);

      // get position of neighboring particle
      const double* neighborpos =
          neighborcontainer->GetPtrToParticleState(PARTICLEENGINE::Position, neighborindex);

      // distance vector from position to neighboring particle
      double dist[3];

      // distance between position and neighboring particle considering periodic boundaries
      DistanceBetweenParticles(position, neighborpos, dist);

      // distance between position and neighboring particle larger than radius
      if (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2] > (radius * radius)) continue;

      // append neighboring particle
      neighboringparticles.push_back(
          std::make_tuple(neighborTypeEnum, neighborStatusEnum, neighborindex));
    }
  }
}

const double* PARTICLEENGINE::ParticleEngine::BinSize() const { return binstrategy_->BinSize(); }

bool PARTICLEENGINE::ParticleEngine::HavePeriodicBoundaryConditions() const
{
  return binstrategy_->HavePeriodicBoundaryConditionsApplied();
}

bool PARTICLEENGINE::ParticleEngine::HavePeriodicBoundaryConditionsInSpatialDirection(
    const int dim) const
{
  return binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim);
}

double PARTICLEENGINE::ParticleEngine::LengthOfBinningDomainInASpatialDirection(const int dim) const
{
  return binstrategy_->LengthOfBinningDomainInASpatialDirection(dim);
}

LINALG::Matrix<3, 2> const& PARTICLEENGINE::ParticleEngine::DomainBoundingBoxCornerPositions() const
{
  return binstrategy_->DomainBoundingBoxCornerPositions();
}

void PARTICLEENGINE::ParticleEngine::DistanceBetweenParticles(
    const double* pos_i, const double* pos_j, double* r_ji) const
{
  // loop over all spatial directions
  for (int dim = 0; dim < 3; ++dim)
  {
    // vector from particle i to j
    r_ji[dim] = pos_j[dim] - pos_i[dim];

    // check for periodic boundary condition in current spatial direction
    if (binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim))
    {
      // binning domain length in current spatial direction
      double binningdomainlength = binstrategy_->LengthOfBinningDomainInASpatialDirection(dim);

      // shift by periodic length if particles are closer over periodic boundary
      if (std::abs(r_ji[dim]) > (0.5 * binningdomainlength))
      {
        if (pos_i[dim] < pos_j[dim])
          r_ji[dim] -= binningdomainlength;
        else
          r_ji[dim] += binningdomainlength;
      }
    }
  }
}

const std::shared_ptr<IO::DiscretizationReader> PARTICLEENGINE::ParticleEngine::BinDisReader(
    int restartstep) const
{
  return std::make_shared<IO::DiscretizationReader>(binstrategy_->BinDiscret(), restartstep);
}

int PARTICLEENGINE::ParticleEngine::GetNumberOfParticles() const
{
  int numberofparticles = 0;

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // add number of particles stored in container
    numberofparticles += container->ParticlesStored();
  }

  return numberofparticles;
}

int PARTICLEENGINE::ParticleEngine::GetNumberOfParticlesOfSpecificType(
    const TypeEnum typeEnum) const
{
  if (not particlecontainerbundle_->GetParticleTypes().count(typeEnum)) return 0;

  // get container of owned particles of specific particle type
  ParticleContainer* container =
      particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

  return container->ParticlesStored();
}

void PARTICLEENGINE::ParticleEngine::WriteBinDisOutput(const int step, const double time) const
{
  // write bins to output file
  binstrategy_->WriteBinOutput(step, time);
}

void PARTICLEENGINE::ParticleEngine::InitBinningStrategy()
{
  // create and init binning strategy and create bins
  binstrategy_ = std::make_shared<BINSTRATEGY::BinningStrategy>();
  binstrategy_->Init();
}

void PARTICLEENGINE::ParticleEngine::SetupBinningStrategy()
{
  // determine minimum relevant bin size
  DetermineMinRelevantBinSize();

  // create an initial linear distribution of row bins
  binrowmap_ = binstrategy_->CreateLinearMapForNumbin(comm_);

  // initialize vector for storage of bin center coordinates and bin weights
  bincenters_ = Teuchos::rcp(new Epetra_MultiVector(*binrowmap_, 3));
  binweights_ = Teuchos::rcp(new Epetra_MultiVector(*binrowmap_, 1));

  // get all bin centers needed for repartitioning
  binstrategy_->GetAllBinCenters(binrowmap_, bincenters_);

  // initialize weights of all bins
  binweights_->PutScalar(1.0e-05);

  // distribute bins via recursive coordinate bisection
  binstrategy_->DistributeBinsRecursCoordBisection(binrowmap_, bincenters_, binweights_);

  // create bins and fill bins into binning discretization
  binstrategy_->FillBinsIntoBinDiscretization(binrowmap_);

  // setup ghosting of bins
  SetupBinGhosting();

  // determine bin distribution dependent maps/sets
  DetermineBinDisDependentMapsAndSets();

  // determine ghosting dependent maps/sets for communication
  DetermineGhostingDependentMapsAndSets();
}

void PARTICLEENGINE::ParticleEngine::SetupBinGhosting()
{
  // gather bins of row map and all its neighbors (row + ghost)
  std::set<int> bins;
  for (int lid = 0; lid < binrowmap_->NumMyElements(); ++lid)
  {
    int gidofbin = binrowmap_->GID(lid);
    std::vector<int> binvec;
    // get neighboring bins
    binstrategy_->GetNeighborAndOwnBinIds(gidofbin, binvec);
    bins.insert(binvec.begin(), binvec.end());
  }

  // remove non-existing ghost bins from original bin set
  {
    // create copy of column bins
    std::set<int> ghostbins(bins);
    // find ghost bins and check for existence
    for (int lid = 0; lid < binrowmap_->NumMyElements(); ++lid)
    {
      const int gid = binrowmap_->GID(lid);
      std::set<int>::iterator iter = ghostbins.find(gid);
      if (iter != ghostbins.end()) ghostbins.erase(iter);
    }
    // only ghost bins remain
    std::vector<int> ghostbins_vec(ghostbins.begin(), ghostbins.end());
    const int size = static_cast<int>(ghostbins.size());
    std::vector<int> pidlist(size);
    const int err = binrowmap_->RemoteIDList(size, ghostbins_vec.data(), pidlist.data(), NULL);
    if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d", err);

    for (int i = 0; i < size; ++i)
    {
      if (pidlist[i] == -1)
      {
        std::set<int>::iterator iter = bins.find(ghostbins_vec[i]);
        if (iter == bins.end()) dserror("bin id is missing in bin set");
        // erase non-existing id
        bins.erase(iter);
      }
    }
  }

  // copy bin gids to a vector and create bincolmap
  std::vector<int> bincolmapvec(bins.begin(), bins.end());
  bincolmap_ = Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(bincolmapvec.size()), &bincolmapvec[0], 0, comm_));

  if (bincolmap_->NumGlobalElements() == 1 && comm_.NumProc() > 1)
    dserror("one bin cannot be run in parallel -> reduce BIN_SIZE_LOWER_BOUND");

  // make sure that all processors are either filled or unfilled
  binstrategy_->BinDiscret()->CheckFilledGlobally();

  // create ghosting for bins
  binstrategy_->BinDiscret()->ExtendedGhosting(*bincolmap_, true, false, true, false);
}

void PARTICLEENGINE::ParticleEngine::InitParticleContainerBundle()
{
  // create and init particle container bundle
  particlecontainerbundle_ = std::make_shared<PARTICLEENGINE::ParticleContainerBundle>();
  particlecontainerbundle_->Init();
}

void PARTICLEENGINE::ParticleEngine::SetupParticleContainerBundle(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes) const
{
  // setup particle container bundle
  particlecontainerbundle_->Setup(particlestatestotypes);
}

void PARTICLEENGINE::ParticleEngine::InitParticleUniqueGlobalIdHandler()
{
  // create and init unique global identifier handler
  particleuniqueglobalidhandler_ = std::unique_ptr<PARTICLEENGINE::UniqueGlobalIdHandler>(
      new PARTICLEENGINE::UniqueGlobalIdHandler(comm_, "particle"));
  particleuniqueglobalidhandler_->Init();
}

void PARTICLEENGINE::ParticleEngine::SetupParticleUniqueGlobalIdHandler() const
{
  // setup unique global identifier handler
  particleuniqueglobalidhandler_->Setup();
}

void PARTICLEENGINE::ParticleEngine::SetupDataStorage(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes)
{
  // determine size of vectors indexed by particle types
  typevectorsize_ = ((--particlestatestotypes.end())->first) + 1;

  // allocate memory to hold particle types
  directghostingtargets_.resize(typevectorsize_);

  // allocate memory for particles being communicated to target processors
  communicatedparticletargets_.assign(comm_.NumProc(), std::vector<int>(0));
}

void PARTICLEENGINE::ParticleEngine::InitParticleVtpWriter()
{
  // construct and init particle runtime vtp writer
  particlevtpwriter_ = std::unique_ptr<PARTICLEENGINE::ParticleRuntimeVtpWriter>(
      new PARTICLEENGINE::ParticleRuntimeVtpWriter(comm_));
  particlevtpwriter_->Init(particlecontainerbundle_);
}

void PARTICLEENGINE::ParticleEngine::SetupParticleVtpWriter() const
{
  // get data format for written numeric data via vtp
  bool write_binary_output = (DRT::INPUT::IntegralValue<INPAR::PARTICLE::OutputDataFormat>(
                                  params_, "OUTPUT_DATA_FORMAT") == INPAR::PARTICLE::binary);

  // get flag to determine output of ghosted particles (debug feature)
  bool write_ghosted_particles = DRT::INPUT::IntegralValue<int>(params_, "WRITE_GHOSTED_PARTICLES");

  // setup particle runtime vtp writer
  particlevtpwriter_->Setup(write_binary_output, write_ghosted_particles);
}

void PARTICLEENGINE::ParticleEngine::SetupTypeWeights()
{
  // allocate memory to hold particle types
  typeweights_.resize(typevectorsize_);

  // init map relating particle types to dynamic load balance factor
  std::map<PARTICLEENGINE::TypeEnum, double> typetodynloadbal;

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_, "PHASE_TO_DYNLOADBALFAC", typetodynloadbal);

  // insert weight of particle type
  for (auto& typeIt : typetodynloadbal) typeweights_[typeIt.first] = typeIt.second;
}

void PARTICLEENGINE::ParticleEngine::DetermineBinDisDependentMapsAndSets()
{
  // clear sets and maps
  boundarybins_.clear();
  touchedbins_.clear();
  firstlayerbinsownedby_.clear();

  // check for finalized construction of binning discretization
  if (binstrategy_->BinDiscret()->Filled() == false)
    dserror("construction of binning discretization not finalized!");

  // loop over row bins
  for (int rowlidofbin = 0; rowlidofbin < binrowmap_->NumMyElements(); ++rowlidofbin)
  {
    int currbin = binrowmap_->GID(rowlidofbin);

    // first insert all owned bins
    boundarybins_.insert(currbin);

    // get neighboring bins
    std::vector<int> binvec;
    binstrategy_->GetNeighborBinIds(currbin, binvec);

    // iterate over neighboring bins
    for (int neighbin : binvec)
    {
      // neighboring bin not owned by this processor
      if (binrowmap_->LID(neighbin) < 0)
      {
        // insert owned bin
        touchedbins_.insert(currbin);

        // insert owner of neighbouring bin
        int neighbinowner = binstrategy_->BinDiscret()->gElement(neighbin)->Owner();
        firstlayerbinsownedby_.insert(std::make_pair(neighbin, neighbinowner));
      }
    }
  }

  // determine all non-boundary bins
  std::set<int> innerbinids;

  // get number of bins in all spatial directions
  const int* binperdir = binstrategy_->BinPerDir();

  // safety check
  for (int dim = 0; dim < 3; ++dim)
    if (binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim) and
        binperdir[dim] < 3)
      dserror("at least 3 bins in direction with periodic boundary conditions necessary!");

  // determine range of all inner bins
  std::vector<int> ijk_min(3);
  std::vector<int> ijk_max(3);
  for (int dim = 0; dim < 3; ++dim)
  {
    ijk_min[dim] = (binperdir[dim] > 2) ? 1 : 0;
    ijk_max[dim] = (binperdir[dim] > 2) ? (binperdir[dim] - 2) : (binperdir[dim] - 1);
  }

  // ijk_range of inner bins (contains: i_min i_max j_min j_max k_min k_max)
  int ijk_range[] = {ijk_min[0], ijk_max[0], ijk_min[1], ijk_max[1], ijk_min[2], ijk_max[2]};

  // get corresponding owned bin ids in ijk range
  binstrategy_->GidsInijkRange(&ijk_range[0], innerbinids, true);

  // substract non-boundary bins from all owned bins to obtain boundary bins
  for (int currbin : innerbinids) boundarybins_.erase(currbin);
}

void PARTICLEENGINE::ParticleEngine::DetermineGhostingDependentMapsAndSets()
{
  // clear sets and maps
  ghostedbins_.clear();
  thisbinsghostedby_.clear();

  // check for finalized construction of binning discretization
  if (binstrategy_->BinDiscret()->Filled() == false)
    dserror("construction of binning discretization not finalized!");

  // loop over col bins
  for (int collidofbin = 0; collidofbin < bincolmap_->NumMyElements(); ++collidofbin)
  {
    int currbin = bincolmap_->GID(collidofbin);

    // current bin not owned by this processor
    if (binrowmap_->LID(currbin) < 0) ghostedbins_.insert(currbin);
  }

  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack data for sending
  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data, ghostedbins_);
  data.StartPacking();
  DRT::ParObject::AddtoPack(data, ghostedbins_);

  // communicate ghosted bins between all processors
  for (int torank = 0; torank < comm_.NumProc(); ++torank)
  {
    if (torank == myrank_) continue;

    sdata[torank].insert(sdata[torank].end(), data().begin(), data().end());
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // init receiving vector
  std::vector<int> receivedbins;

  // unpack and store received data
  for (auto& p : rdata)
  {
    int msgsource = p.first;
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      DRT::ParObject::ExtractfromPack(position, rmsg, receivedbins);

      // iterate over received bins
      for (int receivedbin : receivedbins)
      {
        // received bin is owned by this processor
        if (binrowmap_->LID(receivedbin) >= 0) (thisbinsghostedby_[receivedbin]).insert(msgsource);
      }
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLEENGINE::ParticleEngine::RelateHalfNeighboringBinsToOwnedBins()
{
  // allocate memory for neighbors of owned bins
  halfneighboringbinstobins_.assign(binrowmap_->NumMyElements(), std::set<int>());

  // loop over row bins
  for (int rowlidofbin = 0; rowlidofbin < binrowmap_->NumMyElements(); ++rowlidofbin)
  {
    // get global id of bin
    const int gidofbin = binrowmap_->GID(rowlidofbin);

    // get ijk of current bin
    int ijk[3];
    binstrategy_->ConvertGidToijk(gidofbin, ijk);

    // get reference to neighboring bins (including current bin) of current bin
    std::set<int>& neighboringbins = halfneighboringbinstobins_[rowlidofbin];

    // insert current bin id
    neighboringbins.insert(gidofbin);

    // insert half of the surrounding bins following a specific stencil
    int ijk_range_9bin[] = {ijk[0] - 1, ijk[0] + 1, ijk[1] - 1, ijk[1] + 1, ijk[2] + 1, ijk[2] + 1};
    binstrategy_->GidsInijkRange(&ijk_range_9bin[0], neighboringbins, false);

    int ijk_range_3bin[] = {ijk[0] + 1, ijk[0] + 1, ijk[1] - 1, ijk[1] + 1, ijk[2], ijk[2]};
    binstrategy_->GidsInijkRange(&ijk_range_3bin[0], neighboringbins, false);

    int ijk_range_1bin[] = {ijk[0], ijk[0], ijk[1] + 1, ijk[1] + 1, ijk[2], ijk[2]};
    binstrategy_->GidsInijkRange(&ijk_range_1bin[0], neighboringbins, false);
  }

  // iterate over bins being ghosted on this processor
  for (int gidofbin : ghostedbins_)
  {
    // get neighboring bins
    std::vector<int> binvec;
    binstrategy_->GetNeighborBinIds(gidofbin, binvec);

    // iterate over neighboring bins
    for (int neighbin : binvec)
    {
      // get local id of bin
      const int rowlidofbin = binrowmap_->LID(neighbin);

      // neighboring bin not owned by this processor
      if (rowlidofbin < 0) continue;

      // insert neighboring bins being ghosted on this processor
      halfneighboringbinstobins_[rowlidofbin].insert(gidofbin);
    }
  }

  // validate flag denoting valid relation of half surrounding neighboring bins to owned bins
  validhalfneighboringbins_ = true;
}

void PARTICLEENGINE::ParticleEngine::CheckParticlesAtBoundaries(
    std::vector<std::set<int>>& particlestoremove) const
{
  // safety check
  if (not validownedparticles_) dserror("invalid relation of owned particles to bins!");

  // get bounding box dimensions
  LINALG::Matrix<3, 2> boundingbox = binstrategy_->DomainBoundingBoxCornerPositions();

  // count particles that left the computational domain
  int numparticlesoutside = 0;

  // iterate over owned bins at the boundary
  for (int bdrybin : boundarybins_)
  {
    // get local id of bin
    const int collidofbin = bincolmap_->LID(bdrybin);

    // check if current bin contains owned particles
    if (particlestobins_[collidofbin].empty()) continue;

    // iterate over owned particles in current bin
    for (auto& particleIt : particlestobins_[collidofbin])
    {
      // get type of particle
      TypeEnum typeEnum = particleIt.first;

      // get local index of owned particle
      const int ownedindex = particleIt.second;

      // get container of owned particle of current particle type
      ParticleContainer* container =
          particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

      // get position of particle
      double* currpos = container->GetPtrToParticleState(PARTICLEENGINE::Position, ownedindex);

      // get global id of bin
      const int gidofbin = binstrategy_->ConvertPosToGid(currpos);

      // particle left computational domain
      if (gidofbin == -1)
      {
        (particlestoremove[typeEnum]).insert(ownedindex);

        // get global id of particle
        const int* currglobalid = container->GetPtrToParticleGlobalID(ownedindex);

#ifdef DEBUG
        if (currglobalid[0] < 0) dserror("no global id assigned to particle!");
#endif

        // insert freed global id
        particleuniqueglobalidhandler_->InsertFreedGlobalId(currglobalid[0]);

        ++numparticlesoutside;

        continue;
      }

      // no periodic boundary conditions
      if (not binstrategy_->HavePeriodicBoundaryConditionsApplied()) continue;

      // check for periodic boundary in each spatial directions
      for (int dim = 0; dim < 3; ++dim)
      {
        if (binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim))
        {
          // binning domain length in current spatial direction
          double binningdomainlength = binstrategy_->LengthOfBinningDomainInASpatialDirection(dim);

          // shift position by periodic length
          if (currpos[dim] < boundingbox(dim, 0))
            currpos[dim] += binningdomainlength;
          else if (currpos[dim] > boundingbox(dim, 1))
            currpos[dim] -= binningdomainlength;
        }
      }
    }
  }

  // short screen output
  if (numparticlesoutside)
    IO::cout << "on processor " << myrank_ << " removed " << numparticlesoutside
             << " particle(s) being outside the computational domain!" << IO::endl;
}

void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeDistributed(
    std::vector<ParticleObjShrdPtr>& particlestodistribute,
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestokeep)
{
  // clear particles being communicated to target processors
  communicatedparticletargets_.assign(comm_.NumProc(), std::vector<int>(0));

  // number of particles to distribute
  int numparticles = particlestodistribute.size();

  // global ids of bins containing particles
  std::vector<int> bingidlist(numparticles);

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
  {
    // get states of particle
    const ParticleStates& particleStates = particlestodistribute[i]->ReturnParticleStates();

    // get position of particle
    const std::vector<double>& pos = particleStates[PARTICLEENGINE::Position];

#ifdef DEBUG
    // get type of particles
    TypeEnum typeEnum = particlestodistribute[i]->ReturnParticleType();

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    if (static_cast<int>(pos.size()) != container->GetParticleStateDim(PARTICLEENGINE::Position))
      dserror("dimension of particle state '%s' not valid!",
          PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());
#endif

    // get global id of bin
    bingidlist[i] = binstrategy_->ConvertPosToGid(pos.data());
  }

  // get corresponding processor id
  std::vector<int> pidlist(numparticles);
  {
    // only unique id lists are accepted in RemoteIDList
    // 1) make gid list unique
    std::set<int> unique_bingidlist(bingidlist.begin(), bingidlist.end());
    std::vector<int> uniquevec_bingidlist(unique_bingidlist.begin(), unique_bingidlist.end());
    const int uniquesize = static_cast<int>(uniquevec_bingidlist.size());

    // 2) communication
    std::vector<int> unique_pidlist(uniquesize);
    int err = binrowmap_->RemoteIDList(
        uniquesize, uniquevec_bingidlist.data(), unique_pidlist.data(), NULL);
    if (err < 0) dserror("RemoteIDList returned err=%d", err);

    // 3) build full pid list via lookup table
    std::map<int, int> lookuptable;
    for (int s = 0; s < uniquesize; ++s)
      lookuptable.insert(
          lookuptable.end(), std::pair<int, int>(uniquevec_bingidlist[s], unique_pidlist[s]));
    for (int s = 0; s < numparticles; ++s) pidlist[s] = lookuptable[bingidlist[s]];
  }

  // count particles that are outside of the computational domain
  int numparticlesoutside = 0;

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
  {
    // get type of particle
    TypeEnum typeEnum = particlestodistribute[i]->ReturnParticleType();

    // get owner of particle
    int ownerofparticle = pidlist[i];

    // particle outside of computational domain
    if (ownerofparticle == -1)
    {
      ++numparticlesoutside;

#ifdef DEBUG
      if (particlestodistribute[i]->ReturnParticleGlobalID() < 0)
        dserror("no global id assigned to particle!");
#endif

      // insert freed global id
      particleuniqueglobalidhandler_->InsertFreedGlobalId(
          particlestodistribute[i]->ReturnParticleGlobalID());
    }
    // particle is owned by this processor
    else if (myrank_ == ownerofparticle)
      particlestokeep[typeEnum].push_back(
          std::make_pair(ownerofparticle, particlestodistribute[i]));
    // particle is owned by another processor
    else
    {
      // append particle to be send
      particlestosend[ownerofparticle].push_back(particlestodistribute[i]);

      // append global id of particle to be send
      communicatedparticletargets_[ownerofparticle].emplace_back(
          particlestodistribute[i]->ReturnParticleGlobalID());
    }
  }

  // short screen output
  if (numparticlesoutside)
    IO::cout << "on processor " << myrank_ << " removed " << numparticlesoutside
             << " particle(s) being outside the computational domain!" << IO::endl;

  // clear after all particles are prepared for distribution
  particlestodistribute.clear();
}

void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeTransfered(
    std::vector<std::set<int>>& particlestoremove,
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend)
{
  // safety check
  if (not validownedparticles_) dserror("invalid relation of owned particles to bins!");

  // clear particles being communicated to target processors
  communicatedparticletargets_.assign(comm_.NumProc(), std::vector<int>(0));

  // iterate over this processors bins being touched by other processors
  for (int touchedbin : touchedbins_)
  {
    // get local id of bin
    const int collidofbin = bincolmap_->LID(touchedbin);

    // check if current bin contains owned particles
    if (particlestobins_[collidofbin].empty()) continue;

    // iterate over owned particles in current bin
    for (auto& particleIt : particlestobins_[collidofbin])
    {
      // get type of particle
      TypeEnum typeEnum = particleIt.first;

      // get local index of owned particle
      const int ownedindex = particleIt.second;

      // get container of owned particle of current particle type
      ParticleContainer* container =
          particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

      // get position of particle
      const double* currpos =
          container->GetPtrToParticleState(PARTICLEENGINE::Position, ownedindex);

      // get global id of bin
      const int gidofbin = binstrategy_->ConvertPosToGid(currpos);

      // particle left computational domain
      if (gidofbin == -1)
      {
#ifdef DEBUG
        if (not particlestoremove[typeEnum].count(ownedindex))
          dserror(
              "on processor %d a particle left the computational domain without being detected!",
              myrank_);
#endif
        continue;
      }

      // particle remains owned on this processor
      if (binrowmap_->LID(gidofbin) >= 0) continue;

      // get owning processor
      auto targetIt = firstlayerbinsownedby_.find(gidofbin);
      if (targetIt == firstlayerbinsownedby_.end())
        dserror("particle not owned on this proc but target processor is unknown!");
      int sendtoproc = targetIt->second;

      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(ownedindex, globalid, particleStates);

      // append particle to be send
      particlestosend[sendtoproc].emplace_back(std::make_shared<PARTICLEENGINE::ParticleObject>(
          typeEnum, globalid, particleStates, gidofbin));

      // store index of particle to be removed from containers after particle transfer
      (particlestoremove[typeEnum]).insert(ownedindex);

      // append global id of particle to be send
      communicatedparticletargets_[sendtoproc].emplace_back(globalid);
    }
  }
}

void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeGhosted(
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // safety check
  if (not validownedparticles_) dserror("invalid relation of owned particles to bins!");

  // iterate over this processors bins being ghosted by other processors
  for (auto& targetIt : thisbinsghostedby_)
  {
    // bin being ghosted on other processors
    const int ghostedbin = targetIt.first;

    // get local id of bin
    const int collidofbin = bincolmap_->LID(ghostedbin);

    // check if current bin contains owned particles
    if (particlestobins_[collidofbin].empty()) continue;

    // iterate over owned particles in current bin
    for (auto& particleIt : particlestobins_[collidofbin])
    {
      // get type of particle
      TypeEnum typeEnum = particleIt.first;

      // get local index of owned particle
      const int ownedindex = particleIt.second;

      // get container of owned particle of current particle type
      ParticleContainer* container =
          particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(ownedindex, globalid, particleStates);

      // iterate over target processors
      for (int sendtoproc : targetIt.second)
      {
        // append particle to be send
        particlestosend[sendtoproc].emplace_back(std::make_shared<PARTICLEENGINE::ParticleObject>(
            typeEnum, globalid, particleStates, ghostedbin, ownedindex));
      }
    }
  }
}

void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeRefreshed(
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // safety check
  if (not validdirectghosting_) dserror("invalid direct ghosting!");

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // check for particles of current type to be sent
    if (directghostingtargets_[typeEnum].empty()) continue;

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // iterate over owned particles of current type
    for (auto& indexIt : directghostingtargets_[typeEnum])
    {
      int ownedindex = indexIt.first;

      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(ownedindex, globalid, particleStates);

      // iterate over target processors
      for (auto& targetIt : indexIt.second)
      {
        int sendtoproc = targetIt.first;
        int ghostedindex = targetIt.second;

        // append particle to be send
        particlestosend[sendtoproc].emplace_back(std::make_shared<PARTICLEENGINE::ParticleObject>(
            typeEnum, -1, particleStates, -1, ghostedindex));
      }
    }
  }
}

void PARTICLEENGINE::ParticleEngine::DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
    const StatesOfTypesToRefresh& particlestatestotypes,
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // safety check
  if (not validdirectghosting_) dserror("invalid direct ghosting!");

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    TypeEnum typeEnum = typeIt.first;

    // check for particles of current type to be sent
    if (directghostingtargets_[typeEnum].empty()) continue;

    // determine necessary size of vector for states
    int statesvectorsize = *std::max_element(typeIt.second.begin(), typeIt.second.end()) + 1;

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // iterate over owned particles of current type
    for (auto& indexIt : directghostingtargets_[typeEnum])
    {
      int ownedindex = indexIt.first;

      // allocate memory to hold particle states
      ParticleStates particleStates;
      particleStates.assign(statesvectorsize, std::vector<double>(0));

      // iterate over states to be sent
      for (auto& stateEnum : typeIt.second)
      {
        // get particle state dimension
        int statedim = container->GetParticleStateDim(stateEnum);

        // get pointer to particle state
        const double* state_ptr = container->GetPtrToParticleState(stateEnum, ownedindex);

        // fill particle state
        particleStates[stateEnum].assign(state_ptr, state_ptr + statedim);
      }

      // iterate over target processors
      for (auto& targetIt : indexIt.second)
      {
        int sendtoproc = targetIt.first;
        int ghostedindex = targetIt.second;

        // append particle to be send
        particlestosend[sendtoproc].emplace_back(std::make_shared<PARTICLEENGINE::ParticleObject>(
            typeEnum, -1, particleStates, -1, ghostedindex));
      }
    }
  }
}

void PARTICLEENGINE::ParticleEngine::CommunicateParticles(
    std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoreceive) const
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack data for sending
  for (int torank = 0; torank < comm_.NumProc(); ++torank)
  {
    if (particlestosend[torank].empty()) continue;

    for (auto& iter : particlestosend[torank])
    {
      DRT::PackBuffer data;
      iter->Pack(data);
      data.StartPacking();
      iter->Pack(data);
      sdata[torank].insert(sdata[torank].end(), data().begin(), data().end());
    }
  }

  // clear after all particles are packed
  particlestosend.clear();

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
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack(position, rmsg, data);

      // this std::shared_ptr holds the memory
      std::shared_ptr<DRT::ParObject> object(DRT::UTILS::Factory(data));
      ParticleObjShrdPtr particleobject =
          std::dynamic_pointer_cast<PARTICLEENGINE::ParticleObject>(object);
      if (particleobject == nullptr) dserror("received object is not a particle object!");

      // store received particle
      particlestoreceive[particleobject->ReturnParticleType()].push_back(
          std::make_pair(msgsource, particleobject));
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }
}

void PARTICLEENGINE::ParticleEngine::CommunicateDirectGhostingMap(
    std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting)
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
    directghostingtargets_[typeEnum].clear();

  // invalidate flags denoting validity of direct ghosting
  validdirectghosting_ = false;

  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack data for sending
  for (auto& p : directghosting)
  {
    DRT::PackBuffer data;
    DRT::ParObject::AddtoPack(data, p.second);
    data.StartPacking();
    DRT::ParObject::AddtoPack(data, p.second);
    std::swap(sdata[p.first], data());
  }

  // clear after all ghosting information is packed
  directghosting.clear();

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // init receiving map
  std::map<TypeEnum, std::map<int, std::pair<int, int>>> receiveddirectghosting;

  // unpack and store received data
  for (auto& p : rdata)
  {
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      DRT::ParObject::ExtractfromPack(position, rmsg, receiveddirectghosting);

      // iterate over particle types
      for (auto& typeIt : receiveddirectghosting)
      {
        // get type of particles
        TypeEnum typeEnum = typeIt.first;

        // iterate over this processors local indices of owned particles
        for (auto& indexIt : typeIt.second)
        {
          // get index of owned particle
          int ownedindex = indexIt.first;

          (directghostingtargets_[typeEnum])[ownedindex].push_back(indexIt.second);
        }
      }
    }

    if (position != rmsg.size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>(rmsg.size()), position);
  }

  // validate flags denoting validity of direct ghosting
  validdirectghosting_ = true;
}

void PARTICLEENGINE::ParticleEngine::InsertOwnedParticles(
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert)
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // check for particles of current type
    if (particlestoinsert[typeEnum].empty()) continue;

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // iterate over particle objects pairs
    for (auto& objectpair : particlestoinsert[typeEnum])
    {
      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get global id of particle
      int globalid = particleobject->ReturnParticleGlobalID();

      // get states of particle
      const ParticleStates& particleStates = particleobject->ReturnParticleStates();

#ifdef DEBUG
      if (globalid < 0) dserror("no global id assigned to particle!");

      // get bin of particle
      int gidofbin = particleobject->ReturnBinGid();

      // bin particle
      if (gidofbin < 0)
      {
        // get position of particle
        const std::vector<double>& pos = particleStates[PARTICLEENGINE::Position];

        // get type of particles
        TypeEnum typeEnum = particleobject->ReturnParticleType();

        // get container of owned particles of current particle type
        ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

        if (static_cast<int>(pos.size()) !=
            container->GetParticleStateDim(PARTICLEENGINE::Position))
          dserror("dimension of particle state '%s' not valid!",
              PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());

        // get global id of bin
        gidofbin = binstrategy_->ConvertPosToGid(pos.data());
      }

      // particle not owned by this processor
      if (binrowmap_->LID(gidofbin) < 0) dserror("particle received not owned on this proc!");
#endif

      // add particle to container of owned particles
      int index(0);
      container->AddParticle(index, globalid, particleStates);
    }
  }

  // clear after all particles are inserted
  particlestoinsert.clear();

  // invalidate particle safety flags
  InvalidateParticleSafetyFlags();
}

void PARTICLEENGINE::ParticleEngine::InsertGhostedParticles(
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert,
    std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting)
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // check for particles of current type
    if (particlestoinsert[typeEnum].empty()) continue;

    // get container of ghosted particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Ghosted);

    // iterate over particle objects pairs
    for (auto& objectpair : particlestoinsert[typeEnum])
    {
      // get owner of sending processor
      int sendingproc = objectpair.first;

      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get global id of particle
      int globalid = particleobject->ReturnParticleGlobalID();

      // get states of particle
      const ParticleStates& particleStates = particleobject->ReturnParticleStates();

      // get bin of particle
      const int gidofbin = particleobject->ReturnBinGid();
      if (gidofbin < 0)
        dserror("received ghosted particle contains no information about its bin gid!");

      // add particle to container of ghosted particles
      int ghostedindex(0);
      container->AddParticle(ghostedindex, globalid, particleStates);

      // add index relating (owned and ghosted) particles to col bins
      particlestobins_[bincolmap_->LID(gidofbin)].push_back(std::make_pair(typeEnum, ghostedindex));

      // get local index of particle in container of owned particles of sending processor
      int ownedindex = particleobject->ReturnContainerIndex();

      // insert necessary information being communicated to other processors for direct ghosting
      (((directghosting[sendingproc])[typeEnum])[ownedindex]) =
          std::make_pair(myrank_, ghostedindex);
    }
  }

  // clear after all particles are inserted
  particlestoinsert.clear();

  // validate flag denoting valid relation of ghosted particles to bins
  validghostedparticles_ = true;

  // invalidate safety flags dependent on ghosting
  validparticleneighbors_ = false;
  validglobalidtolocalindex_ = false;
  validdirectghosting_ = false;
}

void PARTICLEENGINE::ParticleEngine::InsertRefreshedParticles(
    std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) const
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // check for particles of current type
    if (particlestoinsert[typeEnum].empty()) continue;

    // get container of ghosted particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Ghosted);

    // iterate over particle objects pairs
    for (auto& objectpair : particlestoinsert[typeEnum])
    {
      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get states of particle
      const ParticleStates& particleStates = particleobject->ReturnParticleStates();

      // get local index of particle in container of ghosted particles on this processor
      int ghostedindex = particleobject->ReturnContainerIndex();

      // replace particle in container of ghosted particles
      container->ReplaceParticle(ghostedindex, -1, particleStates);
    }
  }

  // clear after all particles are inserted
  particlestoinsert.clear();
}

void PARTICLEENGINE::ParticleEngine::RemoveParticlesFromContainers(
    std::vector<std::set<int>>& particlestoremove)
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // check for particles of current type
    if (particlestoremove[typeEnum].empty()) continue;

    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // iterate in reversed order over particles to be removed
    std::set<int>::reverse_iterator rit;
    for (rit = particlestoremove[typeEnum].rbegin(); rit != particlestoremove[typeEnum].rend();
         ++rit)
      container->RemoveParticle(*rit);
  }

  // clear after all particles are removed
  particlestoremove.clear();

  // invalidate particle safety flags
  InvalidateParticleSafetyFlags();
}

void PARTICLEENGINE::ParticleEngine::StorePositionsAfterParticleTransfer()
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored == 0) continue;

    // get pointer to particle states
    const double* pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);
    double* lasttransferpos =
        container->GetPtrToParticleState(PARTICLEENGINE::LastTransferPosition, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // copy particle position data
    for (int i = 0; i < (statedim * particlestored); ++i) lasttransferpos[i] = pos[i];
  }
}

void PARTICLEENGINE::ParticleEngine::RelateOwnedParticlesToBins()
{
  // clear vector relating (owned and ghosted) particles to col bins
  particlestobins_.resize(bincolmap_->NumMyElements());
  for (auto& binIt : particlestobins_) binIt.clear();

  // invalidate particle safety flags
  InvalidateParticleSafetyFlags();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to position of particle after last transfer
    const double* lasttransferpos =
        container->GetPtrToParticleState(PARTICLEENGINE::LastTransferPosition, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // loop over particles in container
    for (int index = 0; index < particlestored; ++index)
    {
      // get global id of bin
      const int gidofbin = binstrategy_->ConvertPosToGid(&(lasttransferpos[statedim * index]));

#ifdef DEBUG
      if (gidofbin < 0) dserror("particle out of bounding box but not removed from container!");

      if (binrowmap_->LID(gidofbin) < 0)
        dserror("particle not owned by this proc but not removed from container!");
#endif

      // add index relating (owned and ghosted) particles to col bins
      particlestobins_[bincolmap_->LID(gidofbin)].push_back(std::make_pair(typeEnum, index));
    }
  }

  // validate flag denoting valid relation of owned particles to bins
  validownedparticles_ = true;
}

void PARTICLEENGINE::ParticleEngine::DetermineMinRelevantBinSize()
{
  // get number of bins in all spatial directions
  const int* binperdir = binstrategy_->BinPerDir();

  // get bin size
  const double* binsize = binstrategy_->BinSize();

  // initialize minimum bin size to maximum bin size
  minbinsize_ = binstrategy_->GetMaxBinSize();

  // check for minimum bin size in spatial directions with more than one bin layer
  for (int i = 0; i < 3; ++i)
    if (binperdir[i] > 1) minbinsize_ = std::min(minbinsize_, binsize[i]);
}

void PARTICLEENGINE::ParticleEngine::DetermineBinWeights()
{
  // safety check
  if (not validownedparticles_) dserror("invalid relation of owned particles to bins!");

  // initialize weights of all bins
  binweights_->PutScalar(1.0e-05);

  // loop over row bins
  for (int rowlidofbin = 0; rowlidofbin < binrowmap_->NumMyElements(); ++rowlidofbin)
  {
    // get global id of bin
    const int gidofbin = binrowmap_->GID(rowlidofbin);

    // iterate over owned particles in current bin
    for (auto& particleIt : particlestobins_[bincolmap_->LID(gidofbin)])
    {
      // add weight of particle of specific type
      (*binweights_)[0][rowlidofbin] += typeweights_[particleIt.first];
    }
  }
}

void PARTICLEENGINE::ParticleEngine::InvalidateParticleSafetyFlags()
{
  validownedparticles_ = false;
  validghostedparticles_ = false;
  validparticleneighbors_ = false;
  validglobalidtolocalindex_ = false;
  validdirectghosting_ = false;
}
