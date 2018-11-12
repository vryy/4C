/*---------------------------------------------------------------------------*/
/*!
\file particle_engine.cpp

\brief particle engine to control particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_engine.H"

#include "particle_communication_utils.H"

#include "particle_container.H"
#include "particle_container_bundle.H"
#include "particle_object.H"
#include "particle_runtime_vtp_writer.H"

#include "../drt_binstrategy/binning_strategy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleEngine::ParticleEngine(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm),
      myrank_(comm.MyPID()),
      params_(params),
      minbinsize_(0.0),
      validparticleconnectivity_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleEngine::~ParticleEngine()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init particle engine                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::Init()
{
  // init binning strategy
  InitBinningStrategy();

  // init particle container bundle
  InitParticleContainerBundle();

  // init particle runtime vtp writer
  InitParticleVtpWriter();
}

/*---------------------------------------------------------------------------*
 | setup particle engine                                      sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::Setup(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes)
{
  // setup binning strategy
  SetupBinningStrategy();

  // setup particle container bundle
  SetupParticleContainerBundle(particlestatestotypes);

  // setup particle runtime vtp writer
  SetupParticleVtpWriter();
}

/*---------------------------------------------------------------------------*
 | write restart of particle engine                           sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::WriteRestart(const int step, const double time) const
{
  // pack particles of all containers
  Teuchos::RCP<std::vector<char>> particlebuffer = Teuchos::rcp(new std::vector<char>);
  particlecontainerbundle_->PackParticleContainerBundle(particlebuffer);

  // get bin discretization writer
  Teuchos::RCP<IO::DiscretizationWriter> binwriter = binstrategy_->BinDiscret()->Writer();

  binwriter->NewStep(step, time);

  // write particle data
  binwriter->WriteCharVector("ParticleData", particlebuffer);

  // write restart of runtime vtp writer
  particlevtpwriter_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of particle engine                            sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
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

  // read restart of runtime vtp writer
  particlevtpwriter_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | write particle runtime vtp output                          sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::WriteParticleRuntimeVtpOutput(
    const int step, const double time) const
{
  particlevtpwriter_->ResetTimeAndTimeStep(time, step);
  particlevtpwriter_->SetParticlePositionsAndStates();
  particlevtpwriter_->WriteFiles();
  particlevtpwriter_->WriteCollectionFileOfAllWrittenFiles();
}

/*---------------------------------------------------------------------------*
 | erase particles outside bounding box                       sfuchs 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::EraseParticlesOutsideBoundingBox(
    std::vector<ParticleObjShrdPtr>& particlestocheck)
{
  // get bounding box dimensions
  LINALG::Matrix<3, 2> xaabb = binstrategy_->XAABB();

  // set of particles located outside bounding box
  std::set<int> particlesoutsideboundingbox;

  // number of particles to check
  int numparticles = particlestocheck.size();

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
  {
    // get particle object
    ParticleObjShrdPtr particleobject = particlestocheck[i];

    // get states of particle
    ParticleStates particleStates = particleobject->ReturnParticleStates();

    // get position of particle
    auto pos = particleStates.find(PARTICLEENGINE::Position);
    if (pos == particleStates.end())
      dserror("particle state '%s' not found!",
          PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());
    double* currpos = (pos->second).data();

    // check particle location with respect to bounding box in each spatial directions
    for (int dim = 0; dim < 3; ++dim)
    {
      // particle located outside bounding box
      if ((currpos[dim] < xaabb(dim, 0)) or (currpos[dim] > xaabb(dim, 1)))
      {
        // insert particle into set
        particlesoutsideboundingbox.insert(i);

        continue;
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
    std::cout << "on processor " << myrank_ << " a total of " << numparticlesoutside
              << " particles are outside of the computational domain and therefore removed!"
              << std::endl;
}

/*---------------------------------------------------------------------------*
 | distribute particles to owning processor                   sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DistributeParticles(
    std::vector<ParticleObjShrdPtr>& particlestodistribute)
{
  // init maps
  std::map<int, std::vector<ParticleObjShrdPtr>> particlestosend;
  std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert;

  // determine particles that need to be distributed
  DetermineParticlesToBeDistributed(particlestodistribute, particlestosend, particlestoinsert);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert owned particles received from other processors
  InsertOwnedParticles(particlestoinsert);

  // rebuild index of owned particles in bin content map
  RebuildIndexOfOwnedParticlesInBinContentMap();
}

/*---------------------------------------------------------------------------*
 | transfer particles to new bins and processors              sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::TransferParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::TransferParticles");

  // init maps
  std::map<TypeEnum, std::set<int>> particlestoremove;
  std::map<int, std::vector<ParticleObjShrdPtr>> particlestosend;
  std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert;

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

  // rebuild index of owned particles in bin content map
  RebuildIndexOfOwnedParticlesInBinContentMap();
}

/*---------------------------------------------------------------------------*
 | ghost particles on other processors                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::GhostParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::GhostParticles");

  // init maps
  std::map<int, std::vector<ParticleObjShrdPtr>> particlestosend;
  std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert;
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
}

/*---------------------------------------------------------------------------*
 | refresh particles being ghosted on other processors        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::RefreshParticles() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::RefreshParticles");

  // init maps
  std::map<int, std::vector<ParticleObjShrdPtr>> particlestosend;
  std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert;

  // determine particles that need to be refreshed
  DetermineParticlesToBeRefreshed(particlestosend);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert refreshed particles received from other processors
  InsertRefreshedParticles(particlestoinsert);
}

/*---------------------------------------------------------------------------*
 | refresh specific states of particles of specific types     sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::RefreshSpecificStatesOfParticlesOfSpecificTypes(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes) const
{
  // init maps
  std::map<int, std::vector<ParticleObjShrdPtr>> particlestosend;
  std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>> particlestoinsert;

  // determine particles that need to be refreshed
  DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
      particlestatestotypes, particlestosend);

  // communicate particles
  CommunicateParticles(particlestosend, particlestoinsert);

  // insert refreshed particles received from other processors
  InsertRefreshedParticles(particlestoinsert);
}

/*---------------------------------------------------------------------------*
 | dynamic load balancing                                     sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DynamicLoadBalancing()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::DynamicLoadBalancing");

  // init maps
  std::vector<ParticleObjShrdPtr> particlestodistribute;

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

  // get vector of particle objects of all containers
  particlecontainerbundle_->GetVectorOfParticleObjectsOfAllContainers(particlestodistribute);

  // clear all containers of owned particles
  particlecontainerbundle_->ClearAllContainersOfSpecificStatus(PARTICLEENGINE::Owned);

  // distribute particles to owning processor
  DistributeParticles(particlestodistribute);
}

/*---------------------------------------------------------------------------*
 | change type of particles                                   sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::TypeChangeParticles(
    std::map<TypeEnum, std::set<int>>& particlestoremove,
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::TypeChangeParticles");

  // remove particles from containers
  RemoveParticlesFromContainers(particlestoremove);

  // insert owned particles undergoing a type change
  InsertOwnedParticles(particlestoinsert);

  // invalidate particle connectivity flag
  validparticleconnectivity_ = false;
}

/*---------------------------------------------------------------------------*
 | build overlapping particle to particle neighbor map        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::BuildParticleToParticleNeighbors()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::BuildParticleToParticleNeighbors");

  // clear map
  particletoparticleneighbors_.clear();

  // iterate over bins containing particles
  for (auto& binIt : bincontent_)
  {
    // get global id of bin
    const int gidofbin = binIt.first;

    // get neighboring bins
    std::vector<int> binvec;
    binstrategy_->GetNeighborAndOwnBinIds(gidofbin, binvec);

    // iterate over particle types
    for (auto& typeIt : binIt.second)
    {
      auto ownedIt = (typeIt.second).find(PARTICLEENGINE::Owned);
      // check if current bin contains owned particles of current type
      if (ownedIt == (typeIt.second).end()) continue;

      // get type of particles
      TypeEnum particleType = typeIt.first;

      // get reference to sub-map
      std::map<int, TypeStatusIndexMap>& currentTypeMap =
          particletoparticleneighbors_[particleType];

      // get set of owned particles of current type
      std::set<int>& setofownedparticles = ownedIt->second;

      // get container of owned particles of current type
      ParticleContainerShrdPtr ownedcontainer =
          particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // iterate over owned particles of current type
      for (int ownedparticle : setofownedparticles)
      {
        // get position of owned particle
        const double* ownedcurrpos =
            ownedcontainer->GetPtrToParticleState(PARTICLEENGINE::Position, ownedparticle);

        // get reference to sub-map
        TypeStatusIndexMap& currentTypeCurrentParticleMap = currentTypeMap[ownedparticle];

        // iterate over neighboring bins (including current bin)
        for (int neighboringbin : binvec)
        {
          auto neighborBinIt = bincontent_.find(neighboringbin);
          // check if neighboring bin contains particles
          if (neighborBinIt == bincontent_.end()) continue;

          // iterate over particle types in neighboring bin
          for (auto& neighborTypeIt : neighborBinIt->second)
          {
            // get type of neighboring particles
            TypeEnum neighborParticleType = neighborTypeIt.first;

            // get reference to sub-map
            std::map<StatusEnum, std::set<int>>& neighborTypeMap =
                currentTypeCurrentParticleMap[neighborParticleType];

            // iterate over particle statuses
            for (auto& neighborStatusIt : neighborTypeIt.second)
            {
              // get status of neighboring particles of current type
              StatusEnum neighborParticleStatus = neighborStatusIt.first;

              // get reference to sub-set
              std::set<int>& neighborTypeStatusSet = neighborTypeMap[neighborParticleStatus];

              // get set of neighboring particles of current type and status
              const std::set<int>& currentTypeCurrentStatusNeighbors = neighborStatusIt.second;

              // get container of neighboring particles of current type and status
              ParticleContainerShrdPtr neighborcontainer =
                  particlecontainerbundle_->GetSpecificContainer(
                      neighborParticleType, neighborParticleStatus);

              // iterate over neighboring particles of current type and status
              for (int neighborparticle : currentTypeCurrentStatusNeighbors)
              {
                // no self-neighboring
                if (particleType == neighborParticleType and gidofbin == neighboringbin and
                    neighborParticleStatus == PARTICLEENGINE::Owned and
                    ownedparticle == neighborparticle)
                  continue;

                // get position of neighboring particle
                const double* neighborcurrpos = neighborcontainer->GetPtrToParticleState(
                    PARTICLEENGINE::Position, neighborparticle);

                // distance vector from owned particle to neighboring particle
                double dist[3];

                // distance between particles considering periodic boundaries
                DistanceBetweenParticles(ownedcurrpos, neighborcurrpos, dist);

                // distance between particles larger than minimum bin size
                if (std::sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]) >
                    minbinsize_)
                  continue;

                // insert neighboring particle of current type and status
                neighborTypeStatusSet.insert(neighborparticle);
              }
            }
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | build global id to local index map                         sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::BuildGlobalIDToLocalIndexMap()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEENGINE::ParticleEngine::BuildGlobalIDToLocalIndexMap");

  // clear map
  globalidtolocalindex_.clear();

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // get status of neighboring particles of current type
      StatusEnum particleStatus = statusIt.first;

      // get container of current type and current status
      ParticleContainerShrdPtr container = statusIt.second;

      // get number of particles stored in container
      int particlestored = container->ParticlesStored();

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
            std::make_shared<LocalIndexTuple>(particleType, particleStatus, index);
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | get reference to particle neighbors map                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
const PARTICLEENGINE::ParticleNeighborsMap&
PARTICLEENGINE::ParticleEngine::GetParticleNeighborsMap() const
{
  // safety check
  if (not validparticleconnectivity_) dserror("invalid particle connectivity!");

  return particletoparticleneighbors_;
}

/*---------------------------------------------------------------------------*
 | get local index in specific particle container             sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
const PARTICLEENGINE::LocalIndexTupleShrdPtr
PARTICLEENGINE::ParticleEngine::GetLocalIndexInSpecificContainer(int globalid) const
{
  // safety check
  if (not validparticleconnectivity_) dserror("invalid particle connectivity!");

  // get local index of particle in specific container
  auto globalidIt = globalidtolocalindex_.find(globalid);
  if (globalidIt == globalidtolocalindex_.end()) return nullptr;

  return globalidIt->second;
}

/*---------------------------------------------------------------------------*
 | return bin size                                            sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
const double* PARTICLEENGINE::ParticleEngine::BinSize() const { return binstrategy_->BinSize(); }

/*---------------------------------------------------------------------------*
 | return flag whether pbc are applied                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEENGINE::ParticleEngine::HavePBC(const int dim) const
{
  return binstrategy_->HavePBC(dim);
}

/*---------------------------------------------------------------------------*
 | return delta for pbc in x, y, or z direction               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEENGINE::ParticleEngine::PBCDelta(const int dim) const
{
  return binstrategy_->PBCDelta(dim);
}

/*---------------------------------------------------------------------------*
 | get bounding box dimensions                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
LINALG::Matrix<3, 2>& PARTICLEENGINE::ParticleEngine::XAABB() const
{
  return binstrategy_->XAABB();
}

/*---------------------------------------------------------------------------*
 | distance between particles considering periodic boundaries sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DistanceBetweenParticles(
    const double* pos_i, const double* pos_j, double* r_ji) const
{
  // loop over all spatial directions
  for (int dim = 0; dim < 3; ++dim)
  {
    // vector from particle i to j
    r_ji[dim] = pos_j[dim] - pos_i[dim];

    // check for periodic boundary condition in current spatial direction
    if (binstrategy_->HavePBC(dim))
    {
      // periodic length in current spatial direction
      const double pbcdelta = binstrategy_->PBCDelta(dim);

      // shift by periodic length if particles are closer over periodic boundary
      if (std::abs(r_ji[dim]) > (0.5 * pbcdelta))
      {
        if (pos_i[dim] < pos_j[dim])
          r_ji[dim] -= pbcdelta;
        else
          r_ji[dim] += pbcdelta;
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | create binning discretization reader                       sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
const std::shared_ptr<IO::DiscretizationReader> PARTICLEENGINE::ParticleEngine::BinDisReader(
    int restartstep) const
{
  return std::make_shared<IO::DiscretizationReader>(binstrategy_->BinDiscret(), restartstep);
}

/*---------------------------------------------------------------------------*
 | get number of particles on this processors                 sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
int PARTICLEENGINE::ParticleEngine::GetNumberOfParticles() const
{
  int numberofparticles = 0;

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // add number of particles stored in container
    numberofparticles += container->ParticlesStored();
  }

  return numberofparticles;
}

/*---------------------------------------------------------------------------*
 | get number of particles on this processor of specific type sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
int PARTICLEENGINE::ParticleEngine::GetNumberOfParticlesOfSpecificType(
    const TypeEnum particleType) const
{
  // get container of owned particles of specific particle type
  auto typeIt = particlecontainerbundle_->GetRefToAllContainersMap().find(particleType);
  if (typeIt == particlecontainerbundle_->GetRefToAllContainersMap().end()) return 0;

  auto statusIt = (typeIt->second).find(PARTICLEENGINE::Owned);
  if (statusIt == (typeIt->second).end()) return 0;

  ParticleContainerShrdPtr container = statusIt->second;

  return container->ParticlesStored();
}

/*---------------------------------------------------------------------------*
 | write binning discretization output (debug feature)        sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::WriteBinDisOutput(const int step, const double time) const
{
  // write bins to output file
  binstrategy_->WriteBinOutput(step, time);
}

/*---------------------------------------------------------------------------*
 | init binning strategy                                      sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InitBinningStrategy()
{
  // create and init binning strategy
  binstrategy_ = std::unique_ptr<BINSTRATEGY::BinningStrategy>(new BINSTRATEGY::BinningStrategy());
  binstrategy_->Init(comm_);
}

/*---------------------------------------------------------------------------*
 | setup binning strategy                                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::SetupBinningStrategy()
{
  // create bins
  binstrategy_->CreateBinsBasedOnCutoffAndXAABB();

  // determine minimum relevant bin size
  DetermineMinRelevantBinSize();

  // build periodic boundary condition
  binstrategy_->BuildPeriodicBC();

  // create an initial linear distribution of row bins
  binrowmap_ = binstrategy_->CreateLinearMapForNumbin(comm_);

  // initialize vector for storage of bin center coordinates and bin weights
  bincenters_ = Teuchos::rcp(new Epetra_MultiVector(*binrowmap_, 3));
  binweights_ = Teuchos::rcp(new Epetra_MultiVector(*binrowmap_, 1));

  // get all bin centers needed for repartitioning
  binstrategy_->GetAllBinCenters(binrowmap_, bincenters_);

  // determine bin weights needed for repartitioning
  DetermineBinWeights();

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

/*---------------------------------------------------------------------------*
 | setup ghosting of bins                                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::SetupBinGhosting() const
{
  // gather bins of rowmap and all its neighbors (row + ghost)
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

  // copy bingids to a vector and create bincolmap
  std::vector<int> bincolmapvec(bins.begin(), bins.end());
  Teuchos::RCP<Epetra_Map> bincolmap = Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(bincolmapvec.size()), &bincolmapvec[0], 0, comm_));

  if (bincolmap->NumGlobalElements() == 1 && comm_.NumProc() > 1)
    dserror("one bin cannot be run in parallel -> reduce CUTOFF_RADIUS");

  // make sure that all processors are either filled or unfilled
  binstrategy_->BinDiscret()->CheckFilledGlobally();

  // create ghosting for bins
  binstrategy_->BinDiscret()->ExtendedGhosting(*bincolmap, true, false, true, false);
}

/*---------------------------------------------------------------------------*
 | init particle container bundle                             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InitParticleContainerBundle()
{
  // create and init particle container bundle
  particlecontainerbundle_ = std::make_shared<PARTICLEENGINE::ParticleContainerBundle>(myrank_);
  particlecontainerbundle_->Init();
}

/*---------------------------------------------------------------------------*
 | setup particle container bundle                            sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::SetupParticleContainerBundle(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes) const
{
  // setup particle container bundle
  particlecontainerbundle_->Setup(particlestatestotypes);
}

/*---------------------------------------------------------------------------*
 | init particle runtime vtp writer                           sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InitParticleVtpWriter()
{
  // construct and init particle runtime vtp writer
  particlevtpwriter_ = std::unique_ptr<PARTICLEENGINE::ParticleRuntimeVtpWriter>(
      new PARTICLEENGINE::ParticleRuntimeVtpWriter(comm_));
  particlevtpwriter_->Init(particlecontainerbundle_);
}

/*---------------------------------------------------------------------------*
 | setup particle runtime vtp writer                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | determine bin distribution dependent maps/sets             sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
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
  for (int lid = 0; lid < binrowmap_->NumMyElements(); ++lid)
  {
    int currbin = binrowmap_->GID(lid);

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
    if (binstrategy_->HavePBC(dim) and binperdir[dim] < 3)
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

/*---------------------------------------------------------------------------*
 | determine ghosting dependent maps/sets for communication   sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineGhostingDependentMapsAndSets()
{
  // clear sets and maps
  ghostedbins_.clear();
  thisbinsghostedby_.clear();

  // check for finalized construction of binning discretization
  if (binstrategy_->BinDiscret()->Filled() == false)
    dserror("construction of binning discretization not finalized!");

  const Epetra_Map* bincolmap = binstrategy_->BinDiscret()->ElementColMap();

  // -----------------------------------------------------------------------
  // determine set ghostedbins_
  // -----------------------------------------------------------------------

  // loop over col bins
  for (int lid = 0; lid < bincolmap->NumMyElements(); ++lid)
  {
    int colbinid = bincolmap->GID(lid);

    // current bin not owned by this processor
    if (binrowmap_->LID(colbinid) < 0) ghostedbins_.insert(colbinid);
  }

  // -----------------------------------------------------------------------
  // determine map thisbinsghostedby_
  // -----------------------------------------------------------------------

  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // number of processors
  int const numproc = comm_.NumProc();

  // pack data for sending
  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data, ghostedbins_);
  data.StartPacking();
  DRT::ParObject::AddtoPack(data, ghostedbins_);

  // communicate ghosted bins between all processors
  for (int torank = 0; torank < numproc; ++torank)
  {
    if (torank == myrank_) continue;

    sdata[torank].insert(sdata[torank].end(), data().begin(), data().end());
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // init receiving vector
  std::vector<int> receivedbins;

  // insert received bins
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

    if (position != (rdata[msgsource]).size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>((rdata[msgsource]).size()),
          position);
  }
}

/*---------------------------------------------------------------------------*
 | check particles for periodic boundaries/leaving domain     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::CheckParticlesAtBoundaries(
    std::map<TypeEnum, std::set<int>>& particlestoremove) const
{
  // get bounding box dimensions
  LINALG::Matrix<3, 2> xaabb = binstrategy_->XAABB();

  // count particles that left the computational domain
  int numparticlesoutside = 0;

  // iterate over owned bins at the boundary
  for (int bdrybin : boundarybins_)
  {
    auto binIt = bincontent_.find(bdrybin);
    // check if current bin contains particles
    if (binIt == bincontent_.end()) continue;

    // iterate over particle types
    for (auto& typeIt : binIt->second)
    {
      auto ownedIt = (typeIt.second).find(PARTICLEENGINE::Owned);
      // check if current bin contains owned particles of current type
      if (ownedIt == (typeIt.second).end()) continue;

      // get type of particles
      TypeEnum particleType = typeIt.first;

      // get container of owned particles of current particle type
      ParticleContainerShrdPtr container =
          particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // iterate over owned particles of current type in current bin
      for (int ownedindex : ownedIt->second)
      {
        // get position of particle
        double* currpos = container->GetPtrToParticleState(PARTICLEENGINE::Position, ownedindex);

        // get global id of bin
        const int gidofbin = binstrategy_->ConvertPosToGid(currpos);

        // particle left computational domain
        if (gidofbin == -1)
        {
          (particlestoremove[particleType]).insert(ownedindex);

          ++numparticlesoutside;

          continue;
        }

        // no periodic boundary conditions
        if (not binstrategy_->HavePBC()) continue;

        // check for periodic boundary in each spatial directions
        for (int dim = 0; dim < 3; ++dim)
        {
          if (binstrategy_->HavePBC(dim))
          {
            // periodic length in current spatial direction
            double pbc_length = binstrategy_->PBCDelta(dim);

            // shift position by periodic length
            if (currpos[dim] < xaabb(dim, 0))
              currpos[dim] += pbc_length;
            else if (currpos[dim] > xaabb(dim, 1))
              currpos[dim] -= pbc_length;
          }
        }
      }
    }
  }

  // short screen output
  if (numparticlesoutside)
    std::cout << "on processor " << myrank_ << " a total of " << numparticlesoutside
              << " particles left the computational domain and therefore removed!" << std::endl;
}

/*---------------------------------------------------------------------------*
 | determine particles that need to be distributed            sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeDistributed(
    std::vector<ParticleObjShrdPtr>& particlestodistribute,
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend,
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestokeep) const
{
  // number of particles to distribute
  int numparticles = particlestodistribute.size();

  // global ids of bins containing particles
  std::vector<int> bingidlist(numparticles);

  // iterate over particles objects
  for (int i = 0; i < numparticles; ++i)
  {
    // get particle object
    ParticleObjShrdPtr particleobject = particlestodistribute[i];

    // get states of particle
    ParticleStates particleStates = particleobject->ReturnParticleStates();

    // get position of particle
    auto pos = particleStates.find(PARTICLEENGINE::Position);
    if (pos == particleStates.end())
      dserror("particle state '%s' not found!",
          PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());
    double* currpos = (pos->second).data();

    // get global id of bin
    bingidlist[i] = binstrategy_->ConvertPosToGid(currpos);
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
    // get particle object
    ParticleObjShrdPtr particleobject = particlestodistribute[i];

    // get type of particle
    TypeEnum particleType = particleobject->ReturnParticleType();

    // get owner of particle
    int ownerofparticle = pidlist[i];

    // particle outside of computational domain
    if (ownerofparticle == -1) ++numparticlesoutside;
    // particle is owned by this processor
    else if (myrank_ == ownerofparticle)
      particlestokeep[particleType].push_back(std::make_pair(ownerofparticle, particleobject));
    // particle is owned by another processor
    else
      particlestosend[ownerofparticle].push_back(particleobject);
  }

  // short screen output
  if (numparticlesoutside)
    std::cout << "on processor " << myrank_ << " a total of " << numparticlesoutside
              << " particles are outside of the computational domain and therefore removed!"
              << std::endl;

  // clear map after all particles are prepared for distribution
  particlestodistribute.clear();
}

/*---------------------------------------------------------------------------*
 | determine particles that need to be transfered             sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeTransfered(
    std::map<TypeEnum, std::set<int>>& particlestoremove,
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // iterate over this processors bins being touched by other processors
  for (int touchedbin : touchedbins_)
  {
    auto binIt = bincontent_.find(touchedbin);
    // check if current bin contains particles
    if (binIt == bincontent_.end()) continue;

    // iterate over particle types
    for (auto& typeIt : binIt->second)
    {
      auto ownedIt = (typeIt.second).find(PARTICLEENGINE::Owned);
      // check if current bin contains owned particles of current type
      if (ownedIt == (typeIt.second).end()) continue;

      // get type of particles
      TypeEnum particleType = typeIt.first;

      // get container of owned particles of current particle type
      ParticleContainerShrdPtr container =
          particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // iterate over owned particles of current type in current bin
      for (int ownedindex : ownedIt->second)
      {
        // get position of particle
        double* currpos = container->GetPtrToParticleState(PARTICLEENGINE::Position, ownedindex);

        // get global id of bin
        const int gidofbin = binstrategy_->ConvertPosToGid(currpos);

        // particle left computational domain
        if (gidofbin == -1) continue;

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

        ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
        particleobject->Init(particleType, globalid, particleStates, gidofbin);

        // append particle to be send
        particlestosend[sendtoproc].push_back(particleobject);

        // store index of particle to be removed from containers after particle transfer
        (particlestoremove[particleType]).insert(ownedindex);
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | determine particles that need to be ghosted                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeGhosted(
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // iterate over this processors bins being ghosted by other processors
  for (auto& targetIt : thisbinsghostedby_)
  {
    // bin being ghosted on other processors
    int ghostedbin = targetIt.first;

    auto binIt = bincontent_.find(ghostedbin);
    // check if current bin contains particles
    if (binIt == bincontent_.end()) continue;

    // iterate over particle types
    for (auto& typeIt : binIt->second)
    {
      auto ownedIt = (typeIt.second).find(PARTICLEENGINE::Owned);
      // check if current bin contains owned particles of current type
      if (ownedIt == (typeIt.second).end()) continue;

      // get type of particles
      TypeEnum particleType = typeIt.first;

      // get container of owned particles of current particle type
      ParticleContainerShrdPtr container =
          particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // iterate over owned particles of current type in current bin
      for (int ownedindex : ownedIt->second)
      {
        int globalid(0);
        ParticleStates particleStates;
        container->GetParticle(ownedindex, globalid, particleStates);

        ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
        particleobject->Init(particleType, globalid, particleStates, ghostedbin, ownedindex);

        // iterate over target processors
        for (int sendtoproc : targetIt.second)
        {
          // append particle to be send
          particlestosend[sendtoproc].push_back(particleobject);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | determine particles that need to be refreshed              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineParticlesToBeRefreshed(
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // safety check
  if (not validparticleconnectivity_) dserror("invalid particle connectivity!");

  // iterate over particle types
  for (auto& typeIt : directghostingmap_)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // iterate over owned particles of current type to be sent
    for (auto& indexIt : typeIt.second)
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

        ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
        particleobject->Init(particleType, -1, particleStates, -1, ghostedindex);

        // append particle to be send
        particlestosend[sendtoproc].push_back(particleobject);
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | determine particles that need to be refreshed              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes,
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend) const
{
  // safety check
  if (not validparticleconnectivity_) dserror("invalid particle connectivity!");

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get state enum set
    const std::set<StateEnum>& stateEnumSet = typeIt.second;

    // get iterator to current particle type
    auto ghostingTypeIt = directghostingmap_.find(particleType);
    // check if owned particles of current type need to be refreshed
    if (ghostingTypeIt == directghostingmap_.end()) continue;

    // get container of owned particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // iterate over owned particles of current type to be sent
    for (auto& indexIt : ghostingTypeIt->second)
    {
      int ownedindex = indexIt.first;

      ParticleStates particleStates;

      // iterate over states to be sent
      for (auto& stateEnum : stateEnumSet)
        particleStates[stateEnum] = container->GetParticleState(stateEnum, ownedindex);

      // iterate over target processors
      for (auto& targetIt : indexIt.second)
      {
        int sendtoproc = targetIt.first;
        int ghostedindex = targetIt.second;

        ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
        particleobject->Init(particleType, -1, particleStates, -1, ghostedindex);

        // append particle to be send
        particlestosend[sendtoproc].push_back(particleobject);
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | communicate particles                                      sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::CommunicateParticles(
    std::map<int, std::vector<ParticleObjShrdPtr>>& particlestosend,
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoreceive) const
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // ---- pack data for sending ----
  for (auto& p : particlestosend)
  {
    for (auto& iter : p.second)
    {
      DRT::PackBuffer data;
      iter->Pack(data);
      data.StartPacking();
      iter->Pack(data);
      sdata[p.first].insert(sdata[p.first].end(), data().begin(), data().end());
    }
  }

  // clear map after all particles are packed
  particlestosend.clear();

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // ---- unpack and store received data ----
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

/*---------------------------------------------------------------------------*
 | communicate and build map for direct ghosting              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::CommunicateDirectGhostingMap(
    std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting)
{
  // clear map
  directghostingmap_.clear();

  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // ---- pack data for sending ----
  for (auto& p : directghosting)
  {
    DRT::PackBuffer data;
    DRT::ParObject::AddtoPack(data, p.second);
    data.StartPacking();
    DRT::ParObject::AddtoPack(data, p.second);
    std::swap(sdata[p.first], data());
  }

  // clear map after all ghosting information is packed
  directghosting.clear();

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // init receiving map
  std::map<TypeEnum, std::map<int, std::pair<int, int>>> receiveddirectghosting;

  // ---- unpack and store received data ----
  for (auto& p : rdata)
  {
    int msgsource = p.first;
    std::vector<char>& rmsg = p.second;

    std::vector<char>::size_type position = 0;

    while (position < rmsg.size())
    {
      DRT::ParObject::ExtractfromPack(position, rmsg, receiveddirectghosting);

      // iterate over particle types
      for (auto& typeIt : receiveddirectghosting)
      {
        // get type of particles
        TypeEnum particleType = typeIt.first;

        // iterate over this processors local indices of owned particles
        for (auto& indexIt : typeIt.second)
        {
          // get index of owned particle
          int ownedindex = indexIt.first;

          (directghostingmap_[particleType])[ownedindex].push_back(indexIt.second);
        }
      }
    }

    if (position != (rdata[msgsource]).size())
      dserror("mismatch in size of data %d <-> %d", static_cast<int>((rdata[msgsource]).size()),
          position);
  }
}

/*---------------------------------------------------------------------------*
 | insert owned particles received from other processors      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InsertOwnedParticles(
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) const
{
  // iterate over particle types
  for (auto& typeIt : particlestoinsert)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // iterate over particle objects pairs
    for (auto& objectpair : typeIt.second)
    {
      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get global id of particle
      int globalid = particleobject->ReturnParticleGlobalID();

      // get states of particle
      ParticleStates particleStates = particleobject->ReturnParticleStates();

      // get bin of particle
      int gidofbin = particleobject->ReturnBinGid();

      // bin particle
      if (gidofbin < 0)
      {
        // get position of particle
        auto pos = particleStates.find(PARTICLEENGINE::Position);
        if (pos == particleStates.end())
          dserror("particle state '%s' not found!",
              PARTICLEENGINE::EnumToStateName(PARTICLEENGINE::Position).c_str());
        double* currpos = (pos->second).data();

        // get global id of bin
        gidofbin = binstrategy_->ConvertPosToGid(currpos);
      }

      // particle not owned by this processor
      if (binrowmap_->LID(gidofbin) < 0) dserror("particle received not owned on this proc!");

      // add particle to container of owned particles
      int index(0);
      container->AddParticle(index, globalid, particleStates);
    }
  }

  // clear map after all particles are inserted
  particlestoinsert.clear();
}

/*---------------------------------------------------------------------------*
 | insert ghosted particles received from other processors    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InsertGhostedParticles(
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert,
    std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting)
{
  // iterate over particle types
  for (auto& typeIt : particlestoinsert)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of ghosted particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Ghosted);

    // iterate over particle objects pairs
    for (auto& objectpair : typeIt.second)
    {
      // get owner of sending processor
      int sendingproc = objectpair.first;

      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get global id of particle
      int globalid = particleobject->ReturnParticleGlobalID();

      // get states of particle
      ParticleStates particleStates = particleobject->ReturnParticleStates();

      // get bin of particle
      const int gidofbin = particleobject->ReturnBinGid();
      if (gidofbin < 0)
        dserror("received ghosted particle contains no information about its bin gid!");

      // add particle to container of ghosted particles
      int ghostedindex(0);
      container->AddParticle(ghostedindex, globalid, particleStates);

      // add index to bin content map
      (((bincontent_[gidofbin])[particleType])[PARTICLEENGINE::Ghosted]).insert(ghostedindex);

      // get local index of particle in container of owned particles of sending processor
      int ownedindex = particleobject->ReturnContainerIndex();

      // insert necessary information into map being communicated to other processors needed for
      // direct ghosting
      (((directghosting[sendingproc])[particleType])[ownedindex]) =
          std::make_pair(myrank_, ghostedindex);
    }
  }

  // clear map after all particles are inserted
  particlestoinsert.clear();
}

/*---------------------------------------------------------------------------*
 | insert refreshed particles received from other processors  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::InsertRefreshedParticles(
    std::map<TypeEnum, std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) const
{
  // iterate over particle types
  for (auto& typeIt : particlestoinsert)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of ghosted particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Ghosted);

    // iterate over particle objects pairs
    for (auto& objectpair : typeIt.second)
    {
      // get particle object
      ParticleObjShrdPtr particleobject = objectpair.second;

      // get states of particle
      ParticleStates particleStates = particleobject->ReturnParticleStates();

      // get local index of particle in container of ghosted particles on this processor
      int ghostedindex = particleobject->ReturnContainerIndex();

      // replace particle in container of ghosted particles
      container->ReplaceParticle(ghostedindex, -1, particleStates);
    }
  }

  // clear map after all particles are inserted
  particlestoinsert.clear();
}

/*---------------------------------------------------------------------------*
 | remove particles from containers                           sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::RemoveParticlesFromContainers(
    std::map<TypeEnum, std::set<int>>& particlestoremove) const
{
  // iterate over particle types
  for (auto& typeIt : particlestoremove)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // iterate in reversed order over particles to be removed
    std::set<int>::reverse_iterator rit;
    for (rit = (typeIt.second).rbegin(); rit != (typeIt.second).rend(); ++rit)
      container->RemoveParticle(*rit);
  }

  // clear map after all particles are removed
  particlestoremove.clear();
}

/*---------------------------------------------------------------------------*
 | rebuild index of owned particles in bin content map        sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::RebuildIndexOfOwnedParticlesInBinContentMap()
{
  // clear map
  bincontent_.clear();

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    auto ownedIt = (typeIt.second).find(PARTICLEENGINE::Owned);
    // check for container of owned particles of current type
    if (ownedIt == (typeIt.second).end()) continue;

    // get container of owned particles of current type
    ParticleContainerShrdPtr container = ownedIt->second;

    // get number of particles stored in container
    int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to position of particle
    double* pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);

    // get dimension of particle position
    int statedim = PARTICLEENGINE::EnumToStateDim(PARTICLEENGINE::Position);

    // loop over particles in container
    for (int index = 0; index < particlestored; ++index)
    {
      // get global id of bin
      const int gidofbin = binstrategy_->ConvertPosToGid(&(pos[statedim * index]));

      // safety checks
      {
        if (gidofbin == -1) dserror("particle out of bounding box but not removed from container!");

        if (binrowmap_->LID(gidofbin) < 0)
          dserror("particle not owned by this proc but not removed from container!");
      }

      // add index to bin content map
      (((bincontent_[gidofbin])[particleType])[PARTICLEENGINE::Owned]).insert(index);
    }
  }
}

/*---------------------------------------------------------------------------*
 | determine minimum relevant bin size                        sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | determine bin weights needed for repartitioning            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleEngine::DetermineBinWeights()
{
  // initialize weights of all bins
  binweights_->PutScalar(1.0e-05);

  // iterate over bins containing particles
  for (auto& binIt : bincontent_)
  {
    // get global id of bin
    const int gidofbin = binIt.first;

    // number of particles in current bin
    int particlecounter = 0;

    // iterate over particle types
    for (auto& typeIt : binIt.second)
    {
      // iterate over particle statuses
      for (auto& statusIt : typeIt.second)
      {
        // only insert owned particles as weights
        if (statusIt.first == PARTICLEENGINE::Owned) particlecounter += (statusIt.second).size();
      }
    }

    // add number of particles in current bin to weights
    binweights_->SumIntoGlobalValue(gidofbin, 0, particlecounter);
  }
}
