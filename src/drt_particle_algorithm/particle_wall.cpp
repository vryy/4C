/*---------------------------------------------------------------------------*/
/*!
\file particle_wall.cpp

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

#include "../drt_particle_engine/particle_engine_interface.H"

#include "../drt_binstrategy/binning_strategy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_factory.H"

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
    : setuptime_(0.0), comm_(comm), myrank_(comm.MyPID()), params_(params)
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
 | write wall runtime vtu output                              sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBase::WriteWallRuntimeVtuOutput(
    const int step, const double time) const
{
  // reset time and time step of the writer object
  wallvtuwriter_->ResetTimeAndTimeStep(time, step);


  // node owner
  Teuchos::RCP<Epetra_Vector> nodeowner =
      LINALG::CreateVector(*walldiscretization_->NodeColMap(), true);
  for (int inode = 0; inode < walldiscretization_->NumMyColNodes(); ++inode)
  {
    const DRT::Node* node = walldiscretization_->lColNode(inode);
    (*nodeowner)[inode] = node->Owner();
  }
  wallvtuwriter_->AppendNodeBasedResultDataVector(nodeowner, 1, "nodeowner");

  // element owner
  Teuchos::RCP<Epetra_Vector> eleowner =
      LINALG::CreateVector(*walldiscretization_->ElementColMap(), true);
  for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
  {
    const DRT::Element* ele = walldiscretization_->lColElement(iele);
    (*eleowner)[iele] = ele->Owner();
  }
  wallvtuwriter_->AppendElementBasedResultDataVector(eleowner, 1, "eleowner");

  // element id
  Teuchos::RCP<Epetra_Vector> eleid =
      LINALG::CreateVector(*walldiscretization_->ElementColMap(), true);
  for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
  {
    const DRT::Element* ele = walldiscretization_->lColElement(iele);
    (*eleid)[iele] = ele->Id();
  }
  wallvtuwriter_->AppendElementBasedResultDataVector(eleid, 1, "eleid");


  // finalize everything and write all required files to filesystem
  wallvtuwriter_->WriteFiles();
  wallvtuwriter_->WriteCollectionFileOfAllWrittenFiles();
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
void PARTICLEALGORITHM::WallHandlerDiscretCondition::DistributeWallElementsAndNodes(
    const Teuchos::RCP<Epetra_Map> binrowmap)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::WallHandlerBase::DistributeWallElesAndNodes");

  // wall discretization (currently) non-moving
  Teuchos::RCP<Epetra_Vector> disnp = Teuchos::null;

  // standard ghosting
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmapdummy;
  binstrategy_->StandardDiscretizationGhosting(
      walldiscretization_, binrowmap, disnp, stdelecolmap, stdnodecolmapdummy);

  // extended ghosting
  std::map<int, std::set<int>> bintorowelemap;
  binstrategy_->DistributeRowElementsToBinsUsingEleXAABB(
      walldiscretization_, bintorowelemap, disnp);

  std::map<int, std::set<int>> ext_bintoele_ghosting;
  Teuchos::RCP<Epetra_Map> extendedelecolmap = binstrategy_->ExtendGhosting(
      walldiscretization_->ElementColMap(), bintorowelemap, ext_bintoele_ghosting);

  // export column elements
  walldiscretization_->ExportColumnElements(*extendedelecolmap);

  // create node column map
  std::set<int> nodes;
  for (int lid = 0; lid < extendedelecolmap->NumMyElements(); ++lid)
  {
    DRT::Element* ele = walldiscretization_->gElement(extendedelecolmap->GID(lid));
    const int* nodeids = ele->NodeIds();
    for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
  }

  std::vector<int> colnodes(nodes.begin(), nodes.end());
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(
      new Epetra_Map(-1, (int)colnodes.size(), &colnodes[0], 0, walldiscretization_->Comm()));

  // export column nodes
  walldiscretization_->ExportColumnNodes(*nodecolmap);

  // finalize wall discretization construction with extended ghosting
  walldiscretization_->FillComplete();
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
void PARTICLEALGORITHM::WallHandlerBoundingBox::DistributeWallElementsAndNodes(
    const Teuchos::RCP<Epetra_Map> binrowmap)
{
  // no need to distribute wall elements and nodes
}

/*---------------------------------------------------------------------------*
 | setup wall discretization                                  sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallHandlerBoundingBox::SetupWallDiscretization() const
{
  // prepare vector of node and element ids
  std::vector<int> nodeids(8);
  std::vector<int> eleids;
  eleids.reserve(6);

  // generate wall discretization from bounding box on first processor
  if (myrank_ == 0)
  {
    // get bounding box dimension
    LINALG::Matrix<3, 2> xaabb = binstrategy_->XAABB();

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
