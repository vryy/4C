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
    Mesh& cut_mesh = cut_mesh_handle.LinearMesh();

    // element is to be added only when it falls within bounding box
    // generated over cut mesh. otherwise this element is never cut
    if (cut_mesh.WithinBB(xyz))
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
          NormalMesh().GetNode(nids[i], &xyz(0, i), lsv[i]);
        else
          NormalMesh().GetNode(nids[i], &xyz(0, i));
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
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::AddCutSide(
    int sid, const std::vector<int>& nids, Core::FE::CellType distype, int mi)
{
  // create side
  return cut_mesh_[mi]->create_side(sid, nids, distype, options_);
}

/*----------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic
 * sidehandle for quadratic sides)
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::AddCutSide(int sid,
    const std::vector<int>& nids, const Core::LinAlg::SerialDenseMatrix& xyz,
    Core::FE::CellType distype, int mi)
{
  Mesh& cut_mesh = CutMesh(mi);

  int numnode = nids.size();
  if (numnode != xyz.numCols())
  {
    FOUR_C_THROW("node coordiante number mismatch");
  }

  //   PointSet nodalpoints;

  // make sure all nodes are there
  for (int i = 0; i < numnode; ++i)
  {
    cut_mesh.GetNode(nids[i], &xyz(0, i));
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
void Core::Geo::Cut::MeshIntersection::BuildSelfCutTree()
{
  Mesh& cm = CutMesh();

  cm.BuildSelfCutTree();
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::build_static_search_tree()
{
  Mesh& m = NormalMesh();

  m.build_static_search_tree();
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion via COMBUST-code                        *
 * where dofsets and node positions                                                               *
 * have not to be computed, standard cut for cut_test                                             *
 * !!!IS JUST USED FOR CUT TESTS                                                                  *
 *                                                                                   schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void Core::Geo::Cut::MeshIntersection::CutTest_Cut(bool include_inner,
    Inpar::Cut::VCellGaussPts VCellgausstype, Inpar::Cut::BCellGaussPts BCellgausstype,
    bool tetcellsonly, bool screenoutput, bool do_Cut_Positions_Dofsets)
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

  BuildSelfCutTree();


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

  Cut_SelfCut(include_inner, screenoutput);

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

  Cut_Finalize(include_inner, VCellgausstype, BCellgausstype, tetcellsonly, screenoutput);

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
void Core::Geo::Cut::MeshIntersection::Cut_SelfCut(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::Geo::CUT --- 2/6 --- Cut_SelfCut");

  Teuchos::RCP<PointPool> point_pool = CutMesh().Points();
  if (CutMesh().GetOptions().Do_SelfCut())
  {
    if (myrank_ == 0 and screenoutput) Core::IO::cout << "\t * 2/6 Cut_SelfCut ...      ";

    point_pool->SetMergeStrategy(PointpoolMergeStrategy::SelfCutLoad);

    SelfCut selfcut(*cut_mesh_[0], myrank_);

    selfcut.PerformSelfCut();
  }
  else if (myrank_ == 0 and screenoutput)
    Core::IO::cout << "\t *2/6 (Skip Cut_SelfCut) ...";

  point_pool->SetMergeStrategy(PointpoolMergeStrategy::NormalCutLoad);
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

  Mesh& m = NormalMesh();

  Mesh& cutmesh = CutMesh();

  m.SearchCollisions(cutmesh);
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

  Mesh& m = NormalMesh();

  m.find_cut_points();
  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

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

  Mesh& m = NormalMesh();

  if (options_.FindPositions())
  {
    // find inside and outside positions of nodes, corresponding facets and volumecells using a
    // DFS-algorithm
    m.FindNodePositions();

    // find the positions for all remaining facets ( and points, volumecells)
    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes, also in parallel
    m.FindNodalDOFSets(include_inner);
  }

}  // Core::Geo::Cut::MeshIntersection::cut_positions_dofsets


/*--------------------------------------------------------------------------------------*
 * get the cut mesh's side based on side id
 *-------------------------------------------------------------------------------------*/
Core::Geo::Cut::SideHandle* Core::Geo::Cut::MeshIntersection::GetCutSide(int sid, int mi) const
{
  return cut_mesh_[mi]->get_side(sid);
}

FOUR_C_NAMESPACE_CLOSE
