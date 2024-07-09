/*----------------------------------------------------------------------*/
/*! \file

\brief provides the specific functionality for cutting a mesh with other meshes

\level 3
 *------------------------------------------------------------------------------------------------*/

#include "4C_cut_meshintersection.hpp"

#include "4C_cut_selfcut.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------*
 * add this background element if it falls within the bounding box of cut mesh
 * If it is not within BB, this element is never cut
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::ElementHandle* Core::Geo::Cut::MeshIntersection::add_element(int eid,
    const std::vector<int>& nids, const Core::LinAlg::SerialDenseMatrix& xyz,
    Core::FE::CellType distype, const double* lsv)
{
  for (std::vector<Teuchos::RCP<MeshHandle>>::iterator i = cut_mesh_.begin(); i != cut_mesh_.end();
       ++i)
  {
    MeshHandle& cut_mesh_handle = **i;
    Mesh& cut_mesh = cut_mesh_handle.linear_mesh();

    // element is to be added only when it falls within bounding box
    // generated over cut mesh. otherwise this element is never cut
    if (cut_mesh.within_bb(xyz))
    {
      int numnode = nids.size();
      if (numnode != xyz.numCols())
      {
        FOUR_C_THROW("node coordiante number mismatch");
      }

      // make sure all nodes are there
      for (int i = 0; i < numnode; ++i)
      {
        if (lsv != nullptr)
          normal_mesh().get_node(nids[i], &xyz(0, i), lsv[i]);
        else
          normal_mesh().get_node(nids[i], &xyz(0, i));
      }

      // create element
      return mesh_.create_element(eid, nids, distype);
    }
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic
 * sidehandle for quadratic sides)
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::add_cut_side(
    int sid, const std::vector<int>& nids, Core::FE::CellType distype, int mi)
{
  // create side
  return cut_mesh_[mi]->create_side(sid, nids, distype, options_);
}

/*----------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic
 * sidehandle for quadratic sides)
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::add_cut_side(int sid,
    const std::vector<int>& nids, const Core::LinAlg::SerialDenseMatrix& xyz,
    Core::FE::CellType distype, int mi)
{
  Mesh& cut_mesh = MeshIntersection::cut_mesh(mi);

  int numnode = nids.size();
  if (numnode != xyz.numCols())
  {
    FOUR_C_THROW("node coordiante number mismatch");
  }

  //   PointSet nodalpoints;

  // make sure all nodes are there
  for (int i = 0; i < numnode; ++i)
  {
    cut_mesh.get_node(nids[i], &xyz(0, i));
    //     nodalpoints.insert( n->point() );
    //     if ( n==nullptr )
    //     {
    //       // if there is no node with that id but a node at the given location,
    //       // the side is illegal and cannot be created
    //       return;
    //     }
  }

  //   // do not create degenerated cut sides
  //   if ( nodalpoints.size() < nids.size() )
  //     return;

  // create side
  return cut_mesh_[mi]->create_side(sid, nids, distype, options_);
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection in the self cut           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::build_self_cut_tree()
{
  Mesh& cm = cut_mesh();

  cm.build_self_cut_tree();
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::build_static_search_tree()
{
  Mesh& m = normal_mesh();

  m.build_static_search_tree();
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion via COMBUST-code                        *
 * where dofsets and node positions                                                               *
 * have not to be computed, standard cut for cut_test                                             *
 * !!!IS JUST USED FOR CUT TESTS                                                                  *
 *                                                                                   schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::cut_test_cut(bool include_inner,
    VCellGaussPts VCellgausstype, Core::Geo::Cut::BCellGaussPts BCellgausstype, bool tetcellsonly,
    bool screenoutput, bool do_Cut_Positions_Dofsets)
{
  int mypid = 0;
  int mpi_is_running = 0;
  MPI_Initialized(&mpi_is_running);
  if (mpi_is_running)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
  }

  // build the static search tree for the collision detection in the self cut
  if (mypid == 0)
  {
    std::cout << "Build self cut tree..............";
    fflush(stdout);
  }

  build_self_cut_tree();


  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // build the static search tree for the collision detection
  if (mypid == 0)
  {
    std::cout << "Build static search tree.........";
    fflush(stdout);
  }

  build_static_search_tree();

  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // handles cut sides which cut each other
  if (mypid == 0)
  {
    std::cout << "Perform self cut.................";
    fflush(stdout);
  }

  cut_self_cut(include_inner, screenoutput);

  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // detects if a side of the cut mesh possibly collides with an element of the background mesh
  if (mypid == 0)
  {
    std::cout << "Collision detection..............";
    fflush(stdout);
  }

  cut_collision_detection(include_inner, screenoutput);

  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // cut the mesh and create cutlines, facets, volumecells
  if (mypid == 0)
  {
    std::cout << "Mesh intersection................";
    fflush(stdout);
  }

  cut_mesh_intersection(screenoutput);

  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // determine inside-outside position and dofset-data, parallel communication if required
  // if false, position of the volumecell will be evaluated in direct divergence algorithm as it is
  // just required there for a cut test!
  if (do_Cut_Positions_Dofsets)  // for most of the cuttest this will not work (still no error is
                                 // thrown!), as there is no enclosed cut boundary!!!
  {
    if (mypid == 0)
    {
      std::cout << "Positions dof-sets...............";
      fflush(stdout);
    }

    cut_positions_dofsets(include_inner, screenoutput);

    if (mypid == 0)
    {
      std::cout << "success\n";
      fflush(stdout);
    }
  }

  // create integration points and/or subtetrahedralization
  if (mypid == 0)
  {
    std::cout << "Finalize.........................";
    fflush(stdout);
  }

  cut_finalize(include_inner, VCellgausstype, BCellgausstype, tetcellsonly, screenoutput);

  if (mypid == 0)
  {
    std::cout << "success\n";
    fflush(stdout);
  }

  // dump_gmsh_volume_cells("CUT_vc", true);
  // dump_gmsh_integration_cells("CUT_intcells");
}

/*------------------------------------------------------------------------------------------------*
 * handles cut sides which cut each other                                                         *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::cut_self_cut(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Geo::CUT --- 2/6 --- Cut_SelfCut");

  Teuchos::RCP<PointPool> point_pool = cut_mesh().points();
  if (cut_mesh().get_options().do_self_cut())
  {
    if (myrank_ == 0 and screenoutput) Core::IO::cout << "\t * 2/6 Cut_SelfCut ...      ";

    point_pool->set_merge_strategy(PointpoolMergeStrategy::SelfCutLoad);

    SelfCut selfcut(*cut_mesh_[0], myrank_);

    selfcut.perform_self_cut();
  }
  else if (myrank_ == 0 and screenoutput)
    Core::IO::cout << "\t *2/6 (Skip Cut_SelfCut) ...";

  point_pool->set_merge_strategy(PointpoolMergeStrategy::NormalCutLoad);
}

/*------------------------------------------------------------------------------------------------*
 * detects if a side of the cut mesh possibly collides with an element of the background mesh     *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::cut_collision_detection(
    bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Geo::CUT --- 3/6 --- cut_collision_detection");

  if (myrank_ == 0 and screenoutput) Core::IO::cout << "\t * 3/6 cut_collision_detection ...";

  Mesh& m = normal_mesh();

  Mesh& cutmesh = cut_mesh();

  m.search_collisions(cutmesh);
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI and XFLUIDFLUID where dofsets and node positions        *
 * have to be parallelized                                                           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::cut_mesh_intersection(bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Geo::CUT --- 4/6 --- cut_mesh_intersection");

  if (myrank_ == 0 and screenoutput) Core::IO::cout << "\t * 4/6 cut_mesh_intersection ...";

  //----------------------------------------------------------

  Mesh& m = normal_mesh();

  m.find_cut_points();
  m.make_cut_lines();
  m.make_facets();
  m.make_volume_cells();

}  // Core::Geo::Cut::MeshIntersection::cut_mesh_intersection

/*------------------------------------------------------------------------------------------------*
 * Routine for deciding the inside-outside position. This creates the dofset data,                *
 * also in parallel                                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::cut_positions_dofsets(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Geo::CUT --- 5/6 --- cut_positions_dofsets (serial)");

  if (myrank_ == 0 and screenoutput) Core::IO::cout << "\t * 5/6 cut_positions_dofsets ...";

  //----------------------------------------------------------

  Mesh& m = normal_mesh();

  if (options_.find_positions())
  {
    // find inside and outside positions of nodes, corresponding facets and volumecells using a
    // DFS-algorithm
    m.find_node_positions();

    // find the positions for all remaining facets ( and points, volumecells)
    m.find_facet_positions();

    // find number and connection of dofsets at nodes from cut volumes, also in parallel
    m.find_nodal_dof_sets(include_inner);
  }

}  // Core::Geo::Cut::MeshIntersection::cut_positions_dofsets


/*--------------------------------------------------------------------------------------*
 * get the cut mesh's side based on side id
 *-------------------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::get_cut_side(int sid, int mi) const
{
  return cut_mesh_[mi]->get_side(sid);
}

FOUR_C_NAMESPACE_CLOSE
