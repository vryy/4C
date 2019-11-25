/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall handler for particle simulations

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_wall.H"

#include "particle_wall_datastate.H"
#include "particle_wall_discretization_runtime_vtu_writer.H"

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
#include "../drt_io/io_pstream.H"

#include "../linalg/linalg_utils.H"

#include "../drt_geometry/searchtree_geometry_service.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallHandlerBase::WallHandlerBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm),
      myrank_(comm.MyPID()),
      params_(params),
      validwallelements_(false),
      validwallneighbors_(false)
{
  // empty constructor
}

PARTICLEWALL::WallHandlerBase::~WallHandlerBase()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

void PARTICLEWALL::WallHandlerBase::Init(
    const std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy)
{
  // set interface to binning strategy
  binstrategy_ = binstrategy;

  // init wall discretization
  InitWallDiscretization();

  // init wall data state container
  InitWallDataState();

  // init wall discretization runtime vtu writer
  InitWallDiscretizationRuntimeVtuWriter();
}

void PARTICLEWALL::WallHandlerBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // setup wall discretization
  SetupWallDiscretization();

  // setup wall data state container
  walldatastate_->Setup();

  // setup wall runtime vtu/vtp writers
  {
    // get data format for written numeric data
    bool write_binary_output = (DRT::INPUT::IntegralValue<INPAR::PARTICLE::OutputDataFormat>(
                                    params_, "OUTPUT_DATA_FORMAT") == INPAR::PARTICLE::binary);

    // setup wall discretization runtime vtu writer
    walldiscretizationruntimevtuwriter_->Setup(write_binary_output);
  }
}

void PARTICLEWALL::WallHandlerBase::WriteRestart(const int step, const double time) const
{
  // get wall discretization writer
  Teuchos::RCP<IO::DiscretizationWriter> walldiscretizationwriter = walldiscretization_->Writer();

  walldiscretizationwriter->NewStep(step, time);

  // write restart of wall data state container
  walldatastate_->WriteRestart(step, time);

  // write restart of wall discretization runtime vtu writer
  walldiscretizationruntimevtuwriter_->WriteRestart(step, time);
}

void PARTICLEWALL::WallHandlerBase::ReadRestart(const int restartstep)
{
  // create discretization reader
  const std::shared_ptr<IO::DiscretizationReader> reader =
      std::make_shared<IO::DiscretizationReader>(walldiscretization_, restartstep);

  // safety check
  if (restartstep != reader->ReadInt("step")) dserror("time step on file not equal to given step!");

  // read restart of wall data state container
  walldatastate_->ReadRestart(reader);

  // read restart of wall discretization runtime vtu writer
  walldiscretizationruntimevtuwriter_->ReadRestart(reader);
}

void PARTICLEWALL::WallHandlerBase::InsertParticleStatesOfParticleTypes(
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

void PARTICLEWALL::WallHandlerBase::WriteWallRuntimeOutput(const int step, const double time) const
{
  // write wall discretization runtime output
  walldiscretizationruntimevtuwriter_->WriteWallDiscretizationRuntimeOutput(step, time);
}

void PARTICLEWALL::WallHandlerBase::UpdateBinRowAndColMap(
    const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap)
{
  binrowmap_ = binrowmap;
  bincolmap_ = bincolmap;
}

void PARTICLEWALL::WallHandlerBase::CheckWallNodesLocatedInBoundingBox() const
{
  // get bounding box dimension
  LINALG::Matrix<3, 2> boundingbox = binstrategy_->DomainBoundingBoxCornerPositions();

  // iterate over row wall nodes
  for (int rowlidofnode = 0; rowlidofnode < walldiscretization_->NumMyRowNodes(); ++rowlidofnode)
  {
    // get pointer to current row wall node
    DRT::Node* node = walldiscretization_->lRowNode(rowlidofnode);

    // init current position of node
    LINALG::Matrix<3, 1> currpos;
    for (int dim = 0; dim < 3; ++dim) currpos(dim) = node->X()[dim];

    if (walldatastate_->GetDispRow() != Teuchos::null)
    {
      // get nodal dofs
      std::vector<int> lm;
      walldiscretization_->Dof(static_cast<unsigned int>(0), node, lm);

      // iterate over spatial directions
      for (int dim = 0; dim < 3; ++dim)
      {
        // local id of nodal dof in current spatial direction
        const int lid = walldiscretization_->DofRowMap()->LID(lm[dim]);

#ifdef DEBUG
        // safety check
        if (lid < 0) dserror("dof gid=%d not in dof row map!", lm[dim]);
#endif

        currpos(dim) += walldatastate_->GetDispRow()->operator[](lid);
      }
    }

    // safety check
    for (int dim = 0; dim < 3; ++dim)
      if (currpos(dim) < boundingbox(dim, 0) or boundingbox(dim, 1) < currpos(dim))
        dserror("node gid=%d resides outside of bounding box!", node->Id());
  }
}

void PARTICLEWALL::WallHandlerBase::GetMaxWallPositionIncrement(
    double& allprocmaxpositionincrement) const
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

    // bin size safety check
    if (maxpositionincrement > particleengineinterface_->MinBinSize())
      dserror("a wall node traveled more than one bin on this processor!");

    // get maximum position increment on all processors
    walldiscretization_->Comm().MaxAll(&maxpositionincrement, &allprocmaxpositionincrement, 1);
  }
}

void PARTICLEWALL::WallHandlerBase::RelateBinsToColWallEles()
{
  // valid flag denoting validity of map relating bins to column wall elements
  if (validwallelements_) return;

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEWALL::WallHandlerBase::RelateBinsToColWallEles");

#ifdef DEBUG
  walldatastate_->CheckForCorrectMaps();
#endif

  // clear vector relating column wall elements to bins
  binstocolwalleles_.assign(walldiscretization_->NumMyColElements(), std::vector<int>(0));

  // invalidate flag denoting validity of wall neighbors map
  validwallneighbors_ = false;

  // iterate over column wall elements
  for (int collidofele = 0; collidofele < walldiscretization_->NumMyColElements(); ++collidofele)
  {
    // get pointer to current column wall element
    DRT::Element* ele = walldiscretization_->lColElement(collidofele);

    // get corresponding bin ids for element
    std::vector<int> binids;
    binstrategy_->DistributeSingleElementToBinsUsingEleAABB(
        walldiscretization_, ele, binids, walldatastate_->GetRefMutableDispCol());

    // relate ids of owned bins to column wall elements
    for (int gidofbin : binids)
      if (not(bincolmap_->LID(gidofbin) < 0)) binstocolwalleles_[collidofele].push_back(gidofbin);
  }

  // validate flag denoting validity of map relating bins to column wall elements
  validwallelements_ = true;
}

void PARTICLEWALL::WallHandlerBase::BuildParticleToWallNeighbors(
    const PARTICLEENGINE::ParticlesToBins& particlestobins)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEWALL::WallHandlerBase::BuildParticleToWallNeighbors");

#ifdef DEBUG
  walldatastate_->CheckForCorrectMaps();
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

    // get pointer to current column wall element
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

const PARTICLEENGINE::PotentialWallNeighbors&
PARTICLEWALL::WallHandlerBase::GetPotentialWallNeighbors() const
{
  // safety check
  if (not validwallneighbors_) dserror("invalid wall neighbors!");

  return potentialwallneighbors_;
}

void PARTICLEWALL::WallHandlerBase::DetermineColWallEleNodalPos(
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

void PARTICLEWALL::WallHandlerBase::InitWallDataState()
{
  // create wall data state container
  walldatastate_ = std::make_shared<PARTICLEWALL::WallDataState>(params_);

  // init wall data state container
  walldatastate_->Init(walldiscretization_);
}

void PARTICLEWALL::WallHandlerBase::InitWallDiscretizationRuntimeVtuWriter()
{
  // create wall discretization runtime vtu writer
  walldiscretizationruntimevtuwriter_ =
      std::unique_ptr<PARTICLEWALL::WallDiscretizationRuntimeVtuWriter>(
          new PARTICLEWALL::WallDiscretizationRuntimeVtuWriter());

  // init wall discretization runtime vtu writer
  walldiscretizationruntimevtuwriter_->Init(walldiscretization_, walldatastate_);
}

void PARTICLEWALL::WallHandlerBase::CreateWallDiscretization()
{
  // create wall discretization
  walldiscretization_ =
      Teuchos::rcp(new DRT::Discretization("particlewalls", Teuchos::rcp(comm_.Clone())));

  // create wall discretization writer
  walldiscretization_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(walldiscretization_)));
}

PARTICLEWALL::WallHandlerDiscretCondition::WallHandlerDiscretCondition(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEWALL::WallHandlerBase(comm, params)
{
  // empty constructor
}

void PARTICLEWALL::WallHandlerDiscretCondition::DistributeWallElementsAndNodes()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEWALL::WallHandlerDiscretCondition::DistributeWallElementsAndNodes");

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
  binstrategy_->DistributeRowElementsToBinsUsingEleAABB(
      walldiscretization_, bintorowelemap, disn_col);

  // extend wall element ghosting
  ExtendWallElementGhosting(bintorowelemap);

  // update maps of state vectors
  walldatastate_->UpdateMapsOfStateVectors();
}

void PARTICLEWALL::WallHandlerDiscretCondition::TransferWallElementsAndNodes()
{
  // transfer wall elements and nodes only if wall displacements are set
  if (walldatastate_->GetDispRow() == Teuchos::null) return;

  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEWALL::WallHandlerDiscretCondition::TransferWallElementsAndNodes");

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
  walldatastate_->UpdateMapsOfStateVectors();
}

void PARTICLEWALL::WallHandlerDiscretCondition::ExtendWallElementGhosting(
    std::map<int, std::set<int>>& bintorowelemap)
{
  std::map<int, std::set<int>> colbintoelemap;
  Teuchos::RCP<Epetra_Map> extendedelecolmap =
      binstrategy_->ExtendElementColMap(bintorowelemap, bintorowelemap, colbintoelemap, bincolmap_);

  BINSTRATEGY::UTILS::ExtendDiscretizationGhosting(
      walldiscretization_, extendedelecolmap, true, false, false);
}

void PARTICLEWALL::WallHandlerDiscretCondition::InitWallDiscretization()
{
  // create wall discretization
  CreateWallDiscretization();

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structurediscretization =
      DRT::Problem::Instance()->GetDis("structure");

  // finalize structure discretization construction
  if (not structurediscretization->Filled()) structurediscretization->FillComplete();

  // get all particle wall conditions
  std::vector<DRT::Condition*> conditions;
  structurediscretization->GetCondition("ParticleWall", conditions);

  // iterate over particle wall conditions
  for (int i = 0; i < (int)conditions.size(); ++i)
  {
    // set current particle wall condition
    std::vector<DRT::Condition*> currcondition(0);
    currcondition.push_back(conditions[i]);

    // get material id for current particle wall condition
    const int mat = currcondition[0]->GetInt("MAT");

    // initialize maps for particle wall conditions
    std::map<int, DRT::Node*> nodes;
    std::map<int, DRT::Node*> colnodes;
    std::map<int, Teuchos::RCP<DRT::Element>> colelements;

    // get structure objects in wall condition
    DRT::UTILS::FindConditionObjects(
        *structurediscretization, nodes, colnodes, colelements, currcondition);

    // iterate over column wall nodes
    for (auto& nodeit : colnodes)
    {
      // get current node
      DRT::Node* currnode = nodeit.second;

      // add current node to wall discretization
      walldiscretization_->AddNode(
          Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
    }

    // iterate over column wall elements
    for (auto& eleit : colelements)
    {
      // get current element
      Teuchos::RCP<DRT::Element> currele = eleit.second;

      // create wall element
      Teuchos::RCP<DRT::Element> wallele =
          DRT::UTILS::Factory("BELE3_3", "Polynomial", currele->Id(), currele->Owner());

      // set node ids to element
      wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());

      // create material for current wall element
      if (not(mat < 0)) wallele->SetMaterial(mat);

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

void PARTICLEWALL::WallHandlerDiscretCondition::SetupWallDiscretization() const
{
  // short screen output
  if (binstrategy_->HavePeriodicBoundaryConditionsApplied() and myrank_ == 0)
    IO::cout << "Warning: particle wall not transferred over periodic boundary!" << IO::endl;
}

PARTICLEWALL::WallHandlerBoundingBox::WallHandlerBoundingBox(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEWALL::WallHandlerBase(comm, params)
{
  // empty constructor
}

void PARTICLEWALL::WallHandlerBoundingBox::DistributeWallElementsAndNodes()
{
  // no need to distribute wall elements and nodes
}

void PARTICLEWALL::WallHandlerBoundingBox::TransferWallElementsAndNodes()
{
  // no need to transfer wall elements and nodes
}

void PARTICLEWALL::WallHandlerBoundingBox::InitWallDiscretization()
{
  // create wall discretization
  CreateWallDiscretization();

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
    LINALG::Matrix<3, 2> boundingbox = binstrategy_->DomainBoundingBoxCornerPositions();

    // reduce bounding box size to account for round-off errors
    for (int dim = 0; dim < 3; ++dim)
    {
      // periodic boundary conditions in current spatial direction
      if (binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim)) continue;

      boundingbox(dim, 0) += 1.0e-12;
      boundingbox(dim, 1) -= 1.0e-12;
    }

    // init vector of corner node positions
    std::vector<std::vector<double>> nodepositions;
    nodepositions.reserve(8);

    // determine corner node positions from bounding box dimension
    nodepositions.push_back({boundingbox(0, 0), boundingbox(1, 0), boundingbox(2, 0)});
    nodepositions.push_back({boundingbox(0, 0), boundingbox(1, 1), boundingbox(2, 0)});
    nodepositions.push_back({boundingbox(0, 0), boundingbox(1, 1), boundingbox(2, 1)});
    nodepositions.push_back({boundingbox(0, 0), boundingbox(1, 0), boundingbox(2, 1)});
    nodepositions.push_back({boundingbox(0, 1), boundingbox(1, 0), boundingbox(2, 0)});
    nodepositions.push_back({boundingbox(0, 1), boundingbox(1, 1), boundingbox(2, 0)});
    nodepositions.push_back({boundingbox(0, 1), boundingbox(1, 1), boundingbox(2, 1)});
    nodepositions.push_back({boundingbox(0, 1), boundingbox(1, 0), boundingbox(2, 1)});

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

    // get material id for particle wall
    const int mat = params_.get<int>("PARTICLE_WALL_MAT");

    int eleid = 0;
    for (int dim = 0; dim < 3; ++dim)
    {
      // periodic boundary conditions in current spatial direction
      if (binstrategy_->HavePeriodicBoundaryConditionsAppliedInSpatialDirection(dim)) continue;

      // positive and negative end of bounding box in current spatial direction
      for (int sign = 0; sign < 2; ++sign)
      {
        // create wall element
        Teuchos::RCP<DRT::Element> wallele =
            DRT::UTILS::Factory("BELE3_3", "Polynomial", eleid, myrank_);

        // set node ids to element
        wallele->SetNodeIds(4, &(nodeidsofelements[dim * 2 + sign])[0]);

        // create material for current wall element
        if (not(mat < 0)) wallele->SetMaterial(mat);

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

void PARTICLEWALL::WallHandlerBoundingBox::SetupWallDiscretization() const
{
  // nothing to do
}
