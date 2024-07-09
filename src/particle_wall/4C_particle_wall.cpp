/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_wall.hpp"

#include "4C_binstrategy.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_transparent.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_searchtree_service.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_particle.hpp"
#include "4C_io.hpp"
#include "4C_io_discretization_visualization_writer_mesh.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_wall_datastate.hpp"
#include "4C_particle_wall_discretization_runtime_vtu_writer.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

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

PARTICLEWALL::WallHandlerBase::~WallHandlerBase() = default;

void PARTICLEWALL::WallHandlerBase::init(
    const std::shared_ptr<Core::Binstrategy::BinningStrategy> binstrategy)
{
  // set interface to binning strategy
  binstrategy_ = binstrategy;

  // init wall discretization
  init_wall_discretization();

  // init wall data state container
  init_wall_data_state();
}

void PARTICLEWALL::WallHandlerBase::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const double restart_time)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // setup wall discretization
  setup_wall_discretization();

  // create wall discretization runtime vtu writer
  create_wall_discretization_runtime_vtu_writer(restart_time);

  // setup wall data state container
  walldatastate_->setup();
}

void PARTICLEWALL::WallHandlerBase::write_restart(const int step, const double time) const
{
  // get wall discretization writer
  Teuchos::RCP<Core::IO::DiscretizationWriter> walldiscretizationwriter =
      walldiscretization_->writer();

  walldiscretizationwriter->new_step(step, time);
}

void PARTICLEWALL::WallHandlerBase::read_restart(const int restartstep) {}

void PARTICLEWALL::WallHandlerBase::insert_particle_states_of_particle_types(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // get flags defining considered states of particle wall
  bool ismoving = Core::UTILS::IntegralValue<int>(params_, "PARTICLE_WALL_MOVING");
  bool isloaded = Core::UTILS::IntegralValue<int>(params_, "PARTICLE_WALL_LOADED");

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

    if (particlestates.count(PARTICLEENGINE::DensityDot))
      particlestates.insert(PARTICLEENGINE::LastIterDensity);

    if (particlestates.count(PARTICLEENGINE::TemperatureDot))
      particlestates.insert(PARTICLEENGINE::LastIterTemperature);
  }
}

void PARTICLEWALL::WallHandlerBase::write_wall_runtime_output(
    const int step, const double time) const
{
  // write wall discretization runtime output
  walldiscretizationruntimevtuwriter_->write_wall_discretization_runtime_output(step, time);
}

void PARTICLEWALL::WallHandlerBase::update_bin_row_and_col_map(
    const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap)
{
  binrowmap_ = binrowmap;
  bincolmap_ = bincolmap;
}

void PARTICLEWALL::WallHandlerBase::check_wall_nodes_located_in_bounding_box() const
{
  // get bounding box dimension
  Core::LinAlg::Matrix<3, 2> boundingbox = binstrategy_->domain_bounding_box_corner_positions();

  // iterate over row wall nodes
  for (int rowlidofnode = 0; rowlidofnode < walldiscretization_->num_my_row_nodes(); ++rowlidofnode)
  {
    // get pointer to current row wall node
    Core::Nodes::Node* node = walldiscretization_->l_row_node(rowlidofnode);

    // init current position of node
    Core::LinAlg::Matrix<3, 1> currpos;
    for (int dim = 0; dim < 3; ++dim) currpos(dim) = node->x()[dim];

    if (walldatastate_->get_disp_row() != Teuchos::null)
    {
      // get nodal dofs
      std::vector<int> lm;
      walldiscretization_->dof(static_cast<unsigned int>(0), node, lm);

      // iterate over spatial directions
      for (int dim = 0; dim < 3; ++dim)
      {
        // local id of nodal dof in current spatial direction
        const int lid = walldiscretization_->dof_row_map()->LID(lm[dim]);

#ifdef FOUR_C_ENABLE_ASSERTIONS
        // safety check
        if (lid < 0) FOUR_C_THROW("dof gid=%d not in dof row map!", lm[dim]);
#endif

        currpos(dim) += walldatastate_->get_disp_row()->operator[](lid);
      }
    }

    // safety check
    for (int dim = 0; dim < 3; ++dim)
      if (currpos(dim) < boundingbox(dim, 0) or boundingbox(dim, 1) < currpos(dim))
        FOUR_C_THROW("node gid=%d resides outside of bounding box!", node->id());
  }
}

void PARTICLEWALL::WallHandlerBase::get_max_wall_position_increment(
    double& allprocmaxpositionincrement) const
{
  if (walldatastate_->get_disp_row() != Teuchos::null)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    // safety checks
    if (walldatastate_->get_disp_row_last_transfer() == Teuchos::null)
      FOUR_C_THROW("vector of wall displacements after last transfer not set!");

    if (not walldatastate_->get_disp_row()->Map().SameAs(
            walldatastate_->get_disp_row_last_transfer()->Map()))
      FOUR_C_THROW("maps are not equal as expected!");
#endif

    // maximum position increment since last particle transfer
    double maxpositionincrement = 0.0;

    // iterate over coordinate values of wall displacements
    for (int i = 0; i < walldatastate_->get_disp_row()->MyLength(); ++i)
    {
      // get position increment of wall node in current spatial dimension since last transfer
      double absolutpositionincrement =
          std::abs(walldatastate_->get_disp_row()->operator[](i) -
                   walldatastate_->get_disp_row_last_transfer()->operator[](i));

      // compare to current maximum
      maxpositionincrement = std::max(maxpositionincrement, absolutpositionincrement);
    }

    // bin size safety check
    if (maxpositionincrement > particleengineinterface_->min_bin_size())
      FOUR_C_THROW("a wall node traveled more than one bin on this processor!");

    // get maximum position increment on all processors
    walldiscretization_->get_comm().MaxAll(&maxpositionincrement, &allprocmaxpositionincrement, 1);
  }
}

void PARTICLEWALL::WallHandlerBase::relate_bins_to_col_wall_eles()
{
  // valid flag denoting validity of map relating bins to column wall elements
  if (validwallelements_) return;

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEWALL::WallHandlerBase::relate_bins_to_col_wall_eles");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  walldatastate_->check_for_correct_maps();
#endif

  // clear vector relating column wall elements to bins
  binstocolwalleles_.assign(walldiscretization_->num_my_col_elements(), std::vector<int>(0));

  // invalidate flag denoting validity of wall neighbors map
  validwallneighbors_ = false;

  // iterate over column wall elements
  for (int collidofele = 0; collidofele < walldiscretization_->num_my_col_elements(); ++collidofele)
  {
    // get pointer to current column wall element
    Core::Elements::Element* ele = walldiscretization_->l_col_element(collidofele);

    // get corresponding bin ids for element
    std::vector<int> binids;
    binstrategy_->distribute_single_element_to_bins_using_ele_aabb(
        walldiscretization_, ele, binids, walldatastate_->get_disp_col());

    // relate ids of owned bins to column wall elements
    for (int gidofbin : binids)
      if (not(bincolmap_->LID(gidofbin) < 0)) binstocolwalleles_[collidofele].push_back(gidofbin);
  }

  // validate flag denoting validity of map relating bins to column wall elements
  validwallelements_ = true;
}

void PARTICLEWALL::WallHandlerBase::build_particle_to_wall_neighbors(
    const PARTICLEENGINE::ParticlesToBins& particlestobins)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEWALL::WallHandlerBase::build_particle_to_wall_neighbors");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  walldatastate_->check_for_correct_maps();
#endif

  // safety check
  if ((not validwallelements_)) FOUR_C_THROW("invalid relation of bins to column wall elements!");

  // clear potential neighboring column wall elements
  potentialwallneighbors_.clear();

  // invalidate flag denoting validity of wall neighbors map
  validwallneighbors_ = false;

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // get minimum relevant bin size
  const double minbinsize = particleengineinterface_->min_bin_size();

  // iterate over column wall elements
  for (int collidofele = 0; collidofele < walldiscretization_->num_my_col_elements(); ++collidofele)
  {
    // set of neighboring bins of current column wall element
    std::set<int> neighborbins;

    // iterate over bins related to current column wall element
    for (int gidofbin : binstocolwalleles_[collidofele])
    {
      // get neighboring bins to current bin
      std::vector<int> binvec;
      binstrategy_->get_neighbor_and_own_bin_ids(gidofbin, binvec);

      // insert into set of neighboring bins of current column wall element
      neighborbins.insert(binvec.begin(), binvec.end());
    }

    // get pointer to current column wall element
    Core::Elements::Element* ele = walldiscretization_->l_col_element(collidofele);

    // determine nodal positions of column wall element
    std::map<int, Core::LinAlg::Matrix<3, 1>> colelenodalpos;
    determine_col_wall_ele_nodal_pos(ele, colelenodalpos);

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
            particlecontainerbundle->get_specific_container(
                neighborTypeEnum, PARTICLEENGINE::Owned);

        // get position of neighboring particle
        const Core::LinAlg::Matrix<3, 1> currpos(
            neighborcontainer->get_ptr_to_state(PARTICLEENGINE::Position, neighborindex));

        // get coordinates of closest point on current column wall element to particle
        Core::LinAlg::Matrix<3, 1> closestpos;
        Core::Geo::nearest3DObjectOnElement(ele, colelenodalpos, currpos, closestpos);

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
PARTICLEWALL::WallHandlerBase::get_potential_wall_neighbors() const
{
  // safety check
  if (not validwallneighbors_) FOUR_C_THROW("invalid wall neighbors!");

  return potentialwallneighbors_;
}

void PARTICLEWALL::WallHandlerBase::determine_col_wall_ele_nodal_pos(
    Core::Elements::Element* ele, std::map<int, Core::LinAlg::Matrix<3, 1>>& colelenodalpos) const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (walldiscretization_->element_col_map()->LID(ele->id()) < 0)
    FOUR_C_THROW("element gid=%d not in element column map!", ele->id());
#endif

  // get pointer to nodes of current column wall element
  Core::Nodes::Node** nodes = ele->nodes();
  const int numnodes = ele->num_node();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  for (int i = 0; i < numnodes; ++i)
    if (walldiscretization_->node_col_map()->LID(nodes[i]->id()) < 0)
      FOUR_C_THROW(
          "node gid=%d of column element gid=%d not in node column map", nodes[i]->id(), ele->id());
#endif

  // determine nodal displacements
  std::vector<double> nodal_disp(numnodes * 3, 0.0);
  if (walldatastate_->get_disp_col() != Teuchos::null)
  {
    std::vector<int> lm_wall;
    lm_wall.reserve(numnodes * 3);
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->location_vector(*walldiscretization_, lm_wall, lmowner, lmstride);

#ifdef FOUR_C_ENABLE_ASSERTIONS
    for (int i = 0; i < numnodes * 3; ++i)
      if (walldiscretization_->dof_col_map()->LID(lm_wall[i]) < 0)
        FOUR_C_THROW("dof gid=%d not in dof column map!", lm_wall[i]);
#endif

    Core::FE::ExtractMyValues(*walldatastate_->get_disp_col(), nodal_disp, lm_wall);
  }

  // iterate over nodes of current column wall element
  for (int k = 0; k < numnodes; ++k)
  {
    // get reference to current nodal position
    Core::LinAlg::Matrix<3, 1>& currpos = colelenodalpos[nodes[k]->id()];

    // determine nodal position
    for (int dim = 0; dim < 3; ++dim) currpos(dim) = nodes[k]->x()[dim] + nodal_disp[k * 3 + dim];
  }
}

void PARTICLEWALL::WallHandlerBase::init_wall_data_state()
{
  // create wall data state container
  walldatastate_ = std::make_shared<PARTICLEWALL::WallDataState>(params_);

  // init wall data state container
  walldatastate_->init(walldiscretization_);
}

void PARTICLEWALL::WallHandlerBase::create_wall_discretization_runtime_vtu_writer(
    const double restart_time)
{
  // create wall discretization runtime vtu writer
  walldiscretizationruntimevtuwriter_ =
      std::make_unique<PARTICLEWALL::WallDiscretizationRuntimeVtuWriter>(
          walldiscretization_, walldatastate_, restart_time);
}

void PARTICLEWALL::WallHandlerBase::create_wall_discretization()
{
  // create wall discretization
  walldiscretization_ = Teuchos::rcp(new Core::FE::Discretization(
      "particlewalls", Teuchos::rcp(comm_.Clone()), Global::Problem::instance()->n_dim()));

  // create wall discretization writer
  walldiscretization_->set_writer(Teuchos::rcp(new Core::IO::DiscretizationWriter(
      walldiscretization_, Global::Problem::instance()->output_control_file(),
      Global::Problem::instance()->spatial_approximation_type())));
}

PARTICLEWALL::WallHandlerDiscretCondition::WallHandlerDiscretCondition(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEWALL::WallHandlerBase(comm, params)
{
  // empty constructor
}

void PARTICLEWALL::WallHandlerDiscretCondition::distribute_wall_elements_and_nodes()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEWALL::WallHandlerDiscretCondition::distribute_wall_elements_and_nodes");

  // invalidate flags
  validwallelements_ = false;
  validwallneighbors_ = false;

  // distribute wall elements to bins with standard ghosting
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmapdummy;
  binstrategy_->standard_discretization_ghosting(walldiscretization_, binrowmap_,
      walldatastate_->get_ref_disp_row(), stdelecolmap, stdnodecolmapdummy);

  // export displacement vector
  Teuchos::RCP<Epetra_Vector> disn_col = Teuchos::null;
  if (walldatastate_->get_disp_row() != Teuchos::null)
  {
    disn_col = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->dof_col_map()));
    Core::LinAlg::Export(*walldatastate_->get_disp_row(), *disn_col);
  }

  // determine bin to row wall element distribution
  std::map<int, std::set<int>> bintorowelemap;
  binstrategy_->distribute_row_elements_to_bins_using_ele_aabb(
      walldiscretization_, bintorowelemap, disn_col);

  // extend wall element ghosting
  extend_wall_element_ghosting(bintorowelemap);

  // update maps of state vectors
  walldatastate_->update_maps_of_state_vectors();
}

void PARTICLEWALL::WallHandlerDiscretCondition::transfer_wall_elements_and_nodes()
{
  // transfer wall elements and nodes only if wall displacements are set
  if (walldatastate_->get_disp_row() == Teuchos::null) return;

  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEWALL::WallHandlerDiscretCondition::transfer_wall_elements_and_nodes");

  // invalidate flags
  validwallelements_ = false;
  validwallneighbors_ = false;

  // transfer wall elements and nodes
  std::map<int, std::set<int>> bintorowelemap;
  binstrategy_->transfer_nodes_and_elements(
      walldiscretization_, walldatastate_->get_disp_col(), bintorowelemap);

  // extend wall element ghosting
  extend_wall_element_ghosting(bintorowelemap);

  // update maps of state vectors
  walldatastate_->update_maps_of_state_vectors();
}

void PARTICLEWALL::WallHandlerDiscretCondition::extend_wall_element_ghosting(
    std::map<int, std::set<int>>& bintorowelemap)
{
  std::map<int, std::set<int>> colbintoelemap;
  Teuchos::RCP<Epetra_Map> extendedelecolmap = binstrategy_->extend_element_col_map(
      bintorowelemap, bintorowelemap, colbintoelemap, bincolmap_);

  Core::Binstrategy::Utils::ExtendDiscretizationGhosting(
      walldiscretization_, extendedelecolmap, true, false, false);
}

void PARTICLEWALL::WallHandlerDiscretCondition::init_wall_discretization()
{
  // create wall discretization
  create_wall_discretization();

  // access the structural discretization
  Teuchos::RCP<Core::FE::Discretization> structurediscretization =
      Global::Problem::instance()->get_dis("structure");

  // finalize structure discretization construction
  if (not structurediscretization->filled()) structurediscretization->fill_complete();

  // get all particle wall conditions
  std::vector<Core::Conditions::Condition*> conditions;
  structurediscretization->get_condition("ParticleWall", conditions);

  // iterate over particle wall conditions
  for (int i = 0; i < static_cast<int>(conditions.size()); ++i)
  {
    // set current particle wall condition
    std::vector<Core::Conditions::Condition*> currcondition(0);
    currcondition.push_back(conditions[i]);

    // get material id for current particle wall condition
    const int mat = currcondition[0]->parameters().get<int>("MAT");

    // initialize maps for particle wall conditions
    std::map<int, Core::Nodes::Node*> nodes;
    std::map<int, Core::Nodes::Node*> colnodes;
    std::map<int, Teuchos::RCP<Core::Elements::Element>> colelements;

    // get structure objects in wall condition
    Core::Conditions::FindConditionObjects(
        *structurediscretization, nodes, colnodes, colelements, currcondition);

    // iterate over column wall nodes
    for (auto& nodeit : colnodes)
    {
      // get current node
      Core::Nodes::Node* currnode = nodeit.second;

      // add current node to wall discretization
      walldiscretization_->add_node(
          Teuchos::rcp(new Core::Nodes::Node(currnode->id(), currnode->x(), currnode->owner())));
    }

    // iterate over column wall elements
    for (auto& eleit : colelements)
    {
      // get current element
      Teuchos::RCP<Core::Elements::Element> currele = eleit.second;

      // create wall element
      Teuchos::RCP<Core::Elements::Element> wallele =
          Core::Communication::Factory("BELE3_3", "Polynomial", currele->id(), currele->owner());

      // set node ids to element
      wallele->set_node_ids(currele->num_node(), currele->node_ids());

      // create material for current wall element
      if (not(mat < 0)) wallele->set_material(0, Mat::Factory(mat));

      // add wall element to discretization
      walldiscretization_->add_element(wallele);
    }
  }

  // reuse dofs of structural discretization for wall discretization
  bool parallel = (comm_.NumProc() == 1) ? false : true;
  Teuchos::RCP<Core::DOFSets::DofSet> newdofset =
      Teuchos::rcp(new Core::DOFSets::TransparentDofSet(structurediscretization, parallel));
  walldiscretization_->replace_dof_set(newdofset);

  // finalize wall discretization construction
  walldiscretization_->fill_complete(true, false, false);
}

void PARTICLEWALL::WallHandlerDiscretCondition::setup_wall_discretization() const
{
  // short screen output
  if (binstrategy_->have_periodic_boundary_conditions_applied() and myrank_ == 0)
    Core::IO::cout << "Warning: particle wall not transferred over periodic boundary!"
                   << Core::IO::endl;
}

PARTICLEWALL::WallHandlerBoundingBox::WallHandlerBoundingBox(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEWALL::WallHandlerBase(comm, params)
{
  // empty constructor
}

void PARTICLEWALL::WallHandlerBoundingBox::distribute_wall_elements_and_nodes()
{
  // no need to distribute wall elements and nodes
}

void PARTICLEWALL::WallHandlerBoundingBox::transfer_wall_elements_and_nodes()
{
  // no need to transfer wall elements and nodes
}

void PARTICLEWALL::WallHandlerBoundingBox::init_wall_discretization()
{
  // create wall discretization
  create_wall_discretization();

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
    Core::LinAlg::Matrix<3, 2> boundingbox = binstrategy_->domain_bounding_box_corner_positions();

    // reduce bounding box size to account for round-off errors
    for (int dim = 0; dim < 3; ++dim)
    {
      // periodic boundary conditions in current spatial direction
      if (binstrategy_->have_periodic_boundary_conditions_applied_in_spatial_direction(dim))
        continue;

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
      walldiscretization_->add_node(Teuchos::rcp(new Core::Nodes::Node(nodeid, nodepos, myrank_)));

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
      if (binstrategy_->have_periodic_boundary_conditions_applied_in_spatial_direction(dim))
        continue;

      // positive and negative end of bounding box in current spatial direction
      for (int sign = 0; sign < 2; ++sign)
      {
        // create wall element
        Teuchos::RCP<Core::Elements::Element> wallele =
            Core::Communication::Factory("BELE3_3", "Polynomial", eleid, myrank_);

        // set node ids to element
        wallele->set_node_ids(4, nodeidsofelements[dim * 2 + sign].data());

        // create material for current wall element
        if (not(mat < 0)) wallele->set_material(0, Mat::Factory(mat));

        // add wall element to discretization
        walldiscretization_->add_element(wallele);

        // add element id
        eleids.push_back(eleid++);
      }
    }

    // no wall elements added
    if (eleids.empty()) FOUR_C_THROW("no wall elements added, check periodic boundary conditions!");
  }

  // node row map of wall elements
  std::shared_ptr<Epetra_Map> noderowmap = std::make_shared<Epetra_Map>(
      -1, nodeids.size(), nodeids.data(), 0, walldiscretization_->get_comm());

  // fully overlapping node column map
  Teuchos::RCP<Epetra_Map> nodecolmap = Core::LinAlg::AllreduceEMap(*noderowmap);

  // element row map of wall elements
  std::shared_ptr<Epetra_Map> elerowmap = std::make_shared<Epetra_Map>(
      -1, eleids.size(), eleids.data(), 0, walldiscretization_->get_comm());

  // fully overlapping element column map
  Teuchos::RCP<Epetra_Map> elecolmap = Core::LinAlg::AllreduceEMap(*elerowmap);

  // fully overlapping ghosting of the wall elements to have everything redundant
  walldiscretization_->export_column_nodes(*nodecolmap);
  walldiscretization_->export_column_elements(*elecolmap);

  // finalize wall discretization construction
  walldiscretization_->fill_complete(true, false, false);
}

void PARTICLEWALL::WallHandlerBoundingBox::setup_wall_discretization() const
{
  // nothing to do
}

FOUR_C_NAMESPACE_CLOSE
