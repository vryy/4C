/*---------------------------------------------------------------------------*/
/*!

\brief particle wall handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_wall.H"

#include "particle_wall_datastate.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_binstrategy/binning_strategy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_dofset_transparent.H"

#include "../drt_io/io.H"
#include "../drt_io/discretization_runtime_vtu_writer.H"

#include "../linalg/linalg_utils.H"

#include "../drt_geometry/searchtree_geometry_service.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::WallHandlerBase::WallHandlerBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : setuptime_(0.0),
      comm_(comm),
      myrank_(comm.MyPID()),
      params_(params),
      validwallelements_(false),
      validwallneighbors_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::WallHandlerBase::~WallHandlerBase()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init wall handler                                          sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::Init()
{
  // init wall discretization
  InitWallDiscretization();

  // init wall data state container
  InitWallDataState();

  // init wall discretization runtime vtu writer
  InitWallVtuWriter();
}

/*---------------------------------------------------------------------------*
 | setup wall handler                                         sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set interface to binning strategy
  binstrategy_ = binstrategy;

  // setup wall discretization
  SetupWallDiscretization();

  // setup wall data state container
  SetupWallDataState();

  // setup wall discretization runtime vtu writer
  SetupWallVtuWriter();
}

/*---------------------------------------------------------------------------*
 | write restart of wall handler                              sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::WriteRestart(const int step, const double time) const
{
  // get wall discretization writer
  Teuchos::RCP<IO::DiscretizationWriter> walldiscretizationwriter = walldiscretization_->Writer();

  walldiscretizationwriter->NewStep(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of wall handler                               sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::ReadRestart(const int restartstep)
{
  // create wall discretization reader
  const Teuchos::RCP<IO::DiscretizationReader> wallreader =
      Teuchos::rcp(new IO::DiscretizationReader(walldiscretization_, restartstep));

  // safety check
  if (restartstep != wallreader->ReadInt("step"))
    dserror("time step on file not equal to given step!");

  // get restart time
  setuptime_ = wallreader->ReadDouble("time");
}

/*---------------------------------------------------------------------------*
 | insert wall handler dependent states of all particle types sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // get flags defining considered states of particle wall
  bool ismoving = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_MOVING");
  bool isloaded = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_LOADED");

  if (not(ismoving and isloaded)) return;

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // insert states needed for iteration in particle structure interaction
    particlestates.insert({
        PARTICLEENGINE::LastIterPosition,
        PARTICLEENGINE::LastIterVelocity,
        PARTICLEENGINE::LastIterAcceleration,
    });

    if (particlestates.count(PARTICLEENGINE::AngularVelocity))
      particlestates.insert(PARTICLEENGINE::LastIterAngularVelocity);

    if (particlestates.count(PARTICLEENGINE::AngularAcceleration))
      particlestates.insert(PARTICLEENGINE::LastIterAngularAcceleration);

    if (particlestates.count(PARTICLEENGINE::ModifiedAcceleration))
      particlestates.insert(PARTICLEENGINE::LastIterModifiedAcceleration);
  }
}

/*---------------------------------------------------------------------------*
 | write wall runtime vtu output                              sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::WriteWallRuntimeVtuOutput(
    const int step, const double time) const
{
  // reset time and time step of the writer object
  wallvtuwriter_->ResetTimeAndTimeStep(time, step);


  // node displacements
  if (walldatastate_->GetDispCol() != Teuchos::null)
    wallvtuwriter_->AppendDofBasedResultDataVector(
        walldatastate_->GetRefMutableDispCol(), 3, 0, "disp");

  // node owner
  Teuchos::RCP<Epetra_Vector> nodeowner =
      Teuchos::rcp(new Epetra_Vector(*walldiscretization_->NodeColMap(), true));
  for (int inode = 0; inode < walldiscretization_->NumMyColNodes(); ++inode)
  {
    const DRT::Node* node = walldiscretization_->lColNode(inode);
    (*nodeowner)[inode] = node->Owner();
  }
  wallvtuwriter_->AppendNodeBasedResultDataVector(nodeowner, 1, "owner");

  // element owner
  Teuchos::RCP<Epetra_Vector> eleowner =
      Teuchos::rcp(new Epetra_Vector(*walldiscretization_->ElementColMap(), true));
  for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
  {
    const DRT::Element* ele = walldiscretization_->lColElement(iele);
    (*eleowner)[iele] = ele->Owner();
  }
  wallvtuwriter_->AppendElementBasedResultDataVector(eleowner, 1, "owner");

  // element id
  Teuchos::RCP<Epetra_Vector> eleid =
      Teuchos::rcp(new Epetra_Vector(*walldiscretization_->ElementColMap(), true));
  for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
  {
    const DRT::Element* ele = walldiscretization_->lColElement(iele);
    (*eleid)[iele] = ele->Id();
  }
  wallvtuwriter_->AppendElementBasedResultDataVector(eleid, 1, "id");


  // finalize everything and write all required files to filesystem
  wallvtuwriter_->WriteFiles();
  wallvtuwriter_->WriteCollectionFileOfAllWrittenFiles();
}

/*---------------------------------------------------------------------------*
 | update bin row and column map                              sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::UpdateBinRowAndColMap(
    const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap)
{
  binrowmap_ = binrowmap;
  bincolmap_ = bincolmap;
}

/*---------------------------------------------------------------------------*
 | get max wall position increment since last transfer        sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::GetMaxWallPositionIncrement(
    double& allprocmaxpositionincrement)
{
  if (walldatastate_->GetDispRow() != Teuchos::null)
  {
#ifdef DEBUG
    // safety checks
    if (walldatastate_->GetDispRowLastTransfer() == Teuchos::null)
      dserror("vector of wall displacements after last transfer not set!");

    if (not walldatastate_->GetDispRow()->Map().SameAs(
            walldatastate_->GetDispRowLastTransfer()->Map()))
      dserror("maps are not equal as expected!");
#endif

    // maximum position increment since last particle transfer
    double maxpositionincrement = 0.0;

    // iterate over coordinate values of wall displacements
    for (int i = 0; i < walldatastate_->GetDispRow()->MyLength(); ++i)
    {
      // get position increment of wall node in current spatial dimension since last transfer
      double absolutpositionincrement =
          std::abs(walldatastate_->GetDispRow()->operator[](i) -
                   walldatastate_->GetDispRowLastTransfer()->operator[](i));

      // compare to current maximum
      maxpositionincrement = std::max(maxpositionincrement, absolutpositionincrement);
    }

    // get maximum position increment on all processors
    walldiscretization_->Comm().MaxAll(&maxpositionincrement, &allprocmaxpositionincrement, 1);
  }
}

/*---------------------------------------------------------------------------*
 | relate bins to column wall elements                        sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::RelateBinsToColWallEles()
{
  // valid flag denoting validity of map relating bins to column wall elements
  if (validwallelements_) return;

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::WallHandlerBase::RelateBinsToColWallEles");

#ifdef DEBUG
  walldatastate_->CheckForCorrectMaps(walldiscretization_);
#endif

  // clear vector relating column wall elements to bins
  binstocolwalleles_.assign(walldiscretization_->NumMyColElements(), std::vector<int>(0));

  // invalidate flag denoting validity of wall neighbors map
  validwallneighbors_ = false;

  // iterate over column wall elements
  for (int collidofele = 0; collidofele < walldiscretization_->NumMyColElements(); ++collidofele)
  {
    // get pointer to current column wall elements
    DRT::Element* ele = walldiscretization_->lColElement(collidofele);

    // get corresponding bin ids for element
    std::vector<int> binids;
    binstrategy_->DistributeElementToBinsUsingEleXAABB(
        walldiscretization_, ele, binids, walldatastate_->GetRefMutableDispCol());

    // relate ids of owned bins to column wall elements
    for (int gidofbin : binids)
      if (not(bincolmap_->LID(gidofbin) < 0)) binstocolwalleles_[collidofele].push_back(gidofbin);
  }

  // validate flag denoting validity of map relating bins to column wall elements
  validwallelements_ = true;
}

/*---------------------------------------------------------------------------*
 | build particle to wall neighbors                           sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::BuildParticleToWallNeighbors(
    const PARTICLEENGINE::ParticlesToBins& particlestobins)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::WallHandlerBase::BuildParticleToWallNeighbors");

#ifdef DEBUG
  walldatastate_->CheckForCorrectMaps(walldiscretization_);
#endif

  // safety check
  if ((not validwallelements_)) dserror("invalid relation of bins to column wall elements!");

  // clear potential neighboring column wall elements
  potentialwallneighbors_.clear();

  // invalidate flag denoting validity of wall neighbors map
  validwallneighbors_ = false;

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // get minimum relevant bin size
  const double minbinsize = particleengineinterface_->MinBinSize();

  // iterate over column wall elements
  for (int collidofele = 0; collidofele < walldiscretization_->NumMyColElements(); ++collidofele)
  {
    // set of neighboring bins of current column wall element
    std::set<int> neighborbins;

    // iterate over bins related to current column wall element
    for (int gidofbin : binstocolwalleles_[collidofele])
    {
      // get neighboring bins to current bin
      std::vector<int> binvec;
      binstrategy_->GetNeighborAndOwnBinIds(gidofbin, binvec);

      // insert into set of neighboring bins of current column wall element
      neighborbins.insert(binvec.begin(), binvec.end());
    }

    // get pointer to current column wall elements
    DRT::Element* ele = walldiscretization_->lColElement(collidofele);

    // determine nodal positions of column wall element
    std::map<int, LINALG::Matrix<3, 1>> colelenodalpos;
    DetermineColWallEleNodalPos(ele, colelenodalpos);

    // iterate over neighboring bins
    for (int gidofneighborbin : neighborbins)
    {
      // consider only owned neighboring bins
      if (binrowmap_->LID(gidofneighborbin) < 0) continue;

      // get local id of neighboring bin
      const int collidofneighboringbin = bincolmap_->LID(gidofneighborbin);

      // check if current neighboring bin contains particles
      if (particlestobins[collidofneighboringbin].empty()) continue;

      // iterate over particles in current neighboring bin
      for (auto& neighborParticleIt : particlestobins[collidofneighboringbin])
      {
        // get type of neighboring particle
        PARTICLEENGINE::TypeEnum neighborTypeEnum = neighborParticleIt.first;

        // get local index of neighboring particle
        const int neighborindex = neighborParticleIt.second;

        // get container of neighboring particle of current particle type
        PARTICLEENGINE::ParticleContainer* neighborcontainer =
            particlecontainerbundle->GetSpecificContainer(neighborTypeEnum, PARTICLEENGINE::Owned);

        // get position of neighboring particle
        const LINALG::Matrix<3, 1> currpos(
            neighborcontainer->GetPtrToParticleState(PARTICLEENGINE::Position, neighborindex));

        // get coordinates of closest point on current column wall element to particle
        LINALG::Matrix<3, 1> closestpos;
        GEO::nearest3DObjectOnElement(ele, colelenodalpos, currpos, closestpos);

        // distance vector from particle to closest point on current column wall element
        double dist[3];
        for (int i = 0; i < 3; i++) dist[i] = closestpos(i) - currpos(i);

        // distance between particle and wall element larger than minimum bin size
        if (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2] > (minbinsize * minbinsize))
          continue;

        // append potential wall neighbor pair
        potentialwallneighbors_.push_back(std::make_pair(
            std::make_tuple(neighborTypeEnum, PARTICLEENGINE::Owned, neighborindex), ele));
      }
    }
  }

  // validate flag denoting validity of wall neighbors map
  validwallneighbors_ = true;
}

/*---------------------------------------------------------------------------*
 | get reference to potential wall neighbors                  sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
const PARTICLEENGINE::PotentialWallNeighbors&
PARTICLEALGORITHM::WallHandlerBase::GetPotentialWallNeighbors() const
{
  // safety check
  if (not validwallneighbors_) dserror("invalid wall neighbors!");

  return potentialwallneighbors_;
}

/*---------------------------------------------------------------------------*
 | determine nodal positions of column wall element           sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::DetermineColWallEleNodalPos(
    DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) const
{
#ifdef DEBUG
  if (walldiscretization_->ElementColMap()->LID(ele->Id()) < 0)
    dserror("element gid=%d not in element column map!", ele->Id());
#endif

  // get pointer to nodes of current column wall element
  DRT::Node** nodes = ele->Nodes();
  const int numnodes = ele->NumNode();

#ifdef DEBUG
  for (int i = 0; i < numnodes; ++i)
    if (walldiscretization_->NodeColMap()->LID(nodes[i]->Id()) < 0)
      dserror(
          "node gid=%d of column element gid=%d not in node column map", nodes[i]->Id(), ele->Id());
#endif

  // determine nodal displacements
  std::vector<double> nodal_disp(numnodes * 3, 0.0);
  if (walldatastate_->GetDispCol() != Teuchos::null)
  {
    std::vector<int> lm_wall;
    lm_wall.reserve(numnodes * 3);
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*walldiscretization_, lm_wall, lmowner, lmstride);

#ifdef DEBUG
    for (int i = 0; i < numnodes * 3; ++i)
      if (walldiscretization_->DofColMap()->LID(lm_wall[i]) < 0)
        dserror("dof gid=%d not in dof column map!", lm_wall[i]);
#endif

    DRT::UTILS::ExtractMyValues(*walldatastate_->GetDispCol(), nodal_disp, lm_wall);
  }

  // iterate over nodes of current column wall element
  for (int k = 0; k < numnodes; ++k)
  {
    // get reference to current nodal position
    LINALG::Matrix<3, 1>& currpos = colelenodalpos[nodes[k]->Id()];

    // determine nodal position
    for (int dim = 0; dim < 3; ++dim) currpos(dim) = nodes[k]->X()[dim] + nodal_disp[k * 3 + dim];
  }
}

/*---------------------------------------------------------------------------*
 | init wall discretization                                   sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::InitWallDiscretization()
{
  // create wall discretization
  walldiscretization_ =
      Teuchos::rcp(new DRT::Discretization("particlewalls", Teuchos::rcp(comm_.Clone())));

  // create wall discretization writer
  walldiscretization_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(walldiscretization_)));
}

/*---------------------------------------------------------------------------*
 | init wall data state container                             sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::InitWallDataState()
{
  // create and init wall data state container
  walldatastate_ = std::make_shared<PARTICLEALGORITHM::WallDataState>();
  walldatastate_->Init();
}

/*---------------------------------------------------------------------------*
 | setup wall data state container                            sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::SetupWallDataState()
{
  // get flags defining considered states of particle wall
  bool ismoving = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_MOVING");
  bool isloaded = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_LOADED");

  // setup wall data state container
  walldatastate_->Setup(walldiscretization_, ismoving, isloaded);
}

/*---------------------------------------------------------------------------*
 | init wall discretization runtime vtu writer                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::InitWallVtuWriter()
{
  // construct wall discretization runtime vtu writer
  wallvtuwriter_ =
      std::unique_ptr<DiscretizationRuntimeVtuWriter>(new DiscretizationRuntimeVtuWriter());
}

/*---------------------------------------------------------------------------*
 | setup wall discretization runtime vtu writer               sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::SetupWallVtuWriter()
{
  // get data format for written numeric data via vtu
  bool write_binary_output = (DRT::INPUT::IntegralValue<INPAR::PARTICLE::OutputDataFormat>(
                                  params_, "OUTPUT_DATA_FORMAT") == INPAR::PARTICLE::binary);

  // we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int max_number_timesteps_to_be_written = 1.0e+6;

  // initialize the writer object
  wallvtuwriter_->Initialize(
      walldiscretization_, max_number_timesteps_to_be_written, setuptime_, write_binary_output);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::WallHandlerDiscretCondition::WallHandlerDiscretCondition(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::WallHandlerBase(comm, params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | distribute wall elements and nodes                         sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerDiscretCondition::DistributeWallElementsAndNodes()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEALGORITHM::WallHandlerDiscretCondition::DistributeWallElementsAndNodes");

  // invalidate flags
  validwallelements_ = false;
  validwallneighbors_ = false;

  // distribute wall elements to bins with standard ghosting
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmapdummy;
  binstrategy_->StandardDiscretizationGhosting(walldiscretization_, binrowmap_,
      walldatastate_->GetRefMutableDispRow(), stdelecolmap, stdnodecolmapdummy);

  // export displacement vector
  Teuchos::RCP<Epetra_Vector> disn_col = Teuchos::null;
  if (walldatastate_->GetDispRow() != Teuchos::null)
  {
    disn_col = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap()));
    LINALG::Export(*walldatastate_->GetDispRow(), *disn_col);
  }

  // determine bin to row wall element distribution
  std::map<int, std::set<int>> bintorowelemap;
  binstrategy_->DistributeRowElementsToBinsUsingEleXAABB(
      walldiscretization_, bintorowelemap, disn_col);

  // extend wall element ghosting
  ExtendWallElementGhosting(bintorowelemap);

  // update maps of state vectors
  walldatastate_->UpdateMapsOfStateVectors(walldiscretization_);
}

/*---------------------------------------------------------------------------*
 | transfer wall elements and nodes                           sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerDiscretCondition::TransferWallElementsAndNodes()
{
  // transfer wall elements and nodes only if wall displacements are set
  if (walldatastate_->GetDispRow() == Teuchos::null) return;

  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEALGORITHM::WallHandlerDiscretCondition::TransferWallElementsAndNodes");

  // invalidate flags
  validwallelements_ = false;
  validwallneighbors_ = false;

  // transfer wall elements and nodes
  std::map<int, std::set<int>> bintorowelemap;
  binstrategy_->TransferNodesAndElements(
      walldiscretization_, walldatastate_->GetMutableDispCol(), bintorowelemap);

  // extend wall element ghosting
  ExtendWallElementGhosting(bintorowelemap);

  // update maps of state vectors
  walldatastate_->UpdateMapsOfStateVectors(walldiscretization_);
}

/*---------------------------------------------------------------------------*
 | extend wall element ghosting                               sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerDiscretCondition::ExtendWallElementGhosting(
    std::map<int, std::set<int>>& bintorowelemap)
{
  std::map<int, std::set<int>> colbintoelemap;
  Teuchos::RCP<Epetra_Map> extendedelecolmap =
      binstrategy_->ExtendGhosting(bintorowelemap, colbintoelemap, bincolmap_);

  BINSTRATEGY::UTILS::ExtendDiscretizationGhosting(
      walldiscretization_, extendedelecolmap, true, false, false);
}

/*---------------------------------------------------------------------------*
 | setup wall discretization                                  sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerDiscretCondition::SetupWallDiscretization() const
{
  // safety check
  if (binstrategy_->HavePBC())
    dserror("periodic boundary conditions not supported for particle wall from discretization!");

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structurediscretization =
      DRT::Problem::Instance()->GetDis("structure");

  // finalize structure discretization construction
  if (not structurediscretization->Filled()) structurediscretization->FillComplete();

  // declare structure objects in wall condition
  std::map<int, std::map<int, DRT::Node*>> nodes;
  std::map<int, std::map<int, DRT::Node*>> colnodes;
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> colelements;

  // get structure objects in wall condition
  DRT::UTILS::FindConditionObjects(
      *structurediscretization, nodes, colnodes, colelements, "ParticleWall");

  for (auto& condit : colelements)
  {
    // iterate over column wall nodes
    for (auto& nodeit : colnodes[condit.first])
    {
      // get current node
      DRT::Node* currnode = nodeit.second;

      // add current node to wall discretization
      walldiscretization_->AddNode(
          Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
    }

    // iterate over column wall elements
    for (auto& eleit : condit.second)
    {
      // get current element
      Teuchos::RCP<DRT::Element> currele = eleit.second;

      // create wall element
      Teuchos::RCP<DRT::Element> wallele =
          DRT::UTILS::Factory("BELE3_3", "Polynomial", currele->Id(), currele->Owner());

      // set node ids to element
      wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());

      // add wall element to discretization
      walldiscretization_->AddElement(wallele);
    }
  }

  // reuse dofs of structural discretization for wall discretization
  bool parallel = (comm_.NumProc() == 1) ? false : true;
  Teuchos::RCP<DRT::DofSet> newdofset =
      Teuchos::rcp(new DRT::TransparentDofSet(structurediscretization, parallel));
  walldiscretization_->ReplaceDofSet(newdofset);

  // finalize wall discretization construction
  walldiscretization_->FillComplete(true, false, false);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::WallHandlerBoundingBox::WallHandlerBoundingBox(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::WallHandlerBase(comm, params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | distribute wall elements and nodes                         sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBoundingBox::DistributeWallElementsAndNodes()
{
  // no need to distribute wall elements and nodes
}

/*---------------------------------------------------------------------------*
 | transfer wall elements and nodes                           sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBoundingBox::TransferWallElementsAndNodes()
{
  // no need to transfer wall elements and nodes
}

/*---------------------------------------------------------------------------*
 | setup wall discretization                                  sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBoundingBox::SetupWallDiscretization() const
{
  // init vector of node and element ids
  std::vector<int> nodeids(0);
  std::vector<int> eleids(0);

  // generate wall discretization from bounding box on first processor
  if (myrank_ == 0)
  {
    // prepare vector of node and element ids
    nodeids.reserve(8);
    eleids.reserve(6);

    // get bounding box dimension
    LINALG::Matrix<3, 2> xaabb = binstrategy_->XAABB();

    // reduce bounding box size to account for round-off errors
    for (int dim = 0; dim < 3; ++dim)
    {
      // periodic boundary conditions in current spatial direction
      if (binstrategy_->HavePBC(dim)) continue;

      xaabb(dim, 0) += 1.0e-12;
      xaabb(dim, 1) -= 1.0e-12;
    }

    // init vector of corner node positions
    std::vector<std::vector<double>> nodepositions;
    nodepositions.reserve(8);

    // determine corner node positions from bounding box dimension
    nodepositions.push_back({xaabb(0, 0), xaabb(1, 0), xaabb(2, 0)});
    nodepositions.push_back({xaabb(0, 0), xaabb(1, 1), xaabb(2, 0)});
    nodepositions.push_back({xaabb(0, 0), xaabb(1, 1), xaabb(2, 1)});
    nodepositions.push_back({xaabb(0, 0), xaabb(1, 0), xaabb(2, 1)});
    nodepositions.push_back({xaabb(0, 1), xaabb(1, 0), xaabb(2, 0)});
    nodepositions.push_back({xaabb(0, 1), xaabb(1, 1), xaabb(2, 0)});
    nodepositions.push_back({xaabb(0, 1), xaabb(1, 1), xaabb(2, 1)});
    nodepositions.push_back({xaabb(0, 1), xaabb(1, 0), xaabb(2, 1)});

    int nodeid = 0;
    for (auto& nodepos : nodepositions)
    {
      // add corner node to wall discretization
      walldiscretization_->AddNode(Teuchos::rcp(new DRT::Node(nodeid, &nodepos[0], myrank_)));

      // add node id
      nodeids.push_back(nodeid++);
    }

    // init vector of node ids to corresponding wall elements
    std::vector<std::vector<int>> nodeidsofelements;
    nodeidsofelements.reserve(6);

    // set corner node ids for each wall element
    nodeidsofelements.push_back({0, 3, 2, 1});
    nodeidsofelements.push_back({4, 5, 6, 7});
    nodeidsofelements.push_back({0, 4, 7, 3});
    nodeidsofelements.push_back({1, 2, 6, 5});
    nodeidsofelements.push_back({0, 1, 5, 4});
    nodeidsofelements.push_back({2, 3, 7, 6});

    int eleid = 0;
    for (int dim = 0; dim < 3; ++dim)
    {
      // periodic boundary conditions in current spatial direction
      if (binstrategy_->HavePBC(dim)) continue;

      // positive and negative end of bounding box in current spatial direction
      for (int sign = 0; sign < 2; ++sign)
      {
        // create wall element
        Teuchos::RCP<DRT::Element> wallele =
            DRT::UTILS::Factory("BELE3_3", "Polynomial", eleid, myrank_);

        // set node ids to element
        wallele->SetNodeIds(4, &(nodeidsofelements[dim * 2 + sign])[0]);

        // add wall element to discretization
        walldiscretization_->AddElement(wallele);

        // add element id
        eleids.push_back(eleid++);
      }
    }

    // no wall elements added
    if (eleids.size() == 0) dserror("no wall elements added, check periodic boundary conditions!");
  }

  // node row map of wall elements
  std::shared_ptr<Epetra_Map> noderowmap =
      std::make_shared<Epetra_Map>(-1, nodeids.size(), &nodeids[0], 0, walldiscretization_->Comm());

  // fully overlapping node column map
  Teuchos::RCP<Epetra_Map> nodecolmap = LINALG::AllreduceEMap(*noderowmap);

  // element row map of wall elements
  std::shared_ptr<Epetra_Map> elerowmap =
      std::make_shared<Epetra_Map>(-1, eleids.size(), &eleids[0], 0, walldiscretization_->Comm());

  // fully overlapping element column map
  Teuchos::RCP<Epetra_Map> elecolmap = LINALG::AllreduceEMap(*elerowmap);

  // fully overlapping ghosting of the wall elements to have everything redundant
  walldiscretization_->ExportColumnNodes(*nodecolmap);
  walldiscretization_->ExportColumnElements(*elecolmap);

  // finalize wall discretization construction
  walldiscretization_->FillComplete(true, false, false);
}
