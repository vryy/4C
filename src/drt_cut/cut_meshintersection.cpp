/*!-----------------------------------------------------------------------------------------------*
\file cut_meshintersection.cpp

\brief provides the specific functionality for cutting a mesh with other meshes

\level 3
<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
 *------------------------------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include <mpi.h>

#include "cut_selfcut.H"
#include "cut_meshintersection.H"

/*-----------------------------------------------------------------------------------------*
 * add this background element if it falls within the bounding box of cut mesh
 * If it is not within BB, this element is never cut
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle* GEO::CUT::MeshIntersection::AddElement(int eid,
    const std::vector<int>& nids, const Epetra_SerialDenseMatrix& xyz,
    DRT::Element::DiscretizationType distype, const double* lsv)
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
      if (numnode != xyz.N())
      {
        throw std::runtime_error("node coordiante number mismatch");
      }

      // make sure all nodes are there
      for (int i = 0; i < numnode; ++i)
      {
        if (lsv != NULL)
          NormalMesh().GetNode(nids[i], &xyz(0, i), lsv[i]);
        else
          NormalMesh().GetNode(nids[i], &xyz(0, i));
      }

      // create element
      return mesh_.CreateElement(eid, nids, distype);
    }
  }
  return NULL;
}

/*----------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic
 * sidehandle for quadratic sides)
 *----------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshIntersection::AddCutSide(
    int sid, const std::vector<int>& nids, DRT::Element::DiscretizationType distype, int mi)
{
  // create side
  return cut_mesh_[mi]->CreateSide(sid, nids, distype);
}

/*----------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic
 * sidehandle for quadratic sides)
 *----------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshIntersection::AddCutSide(int sid, const std::vector<int>& nids,
    const Epetra_SerialDenseMatrix& xyz, DRT::Element::DiscretizationType distype, int mi)
{
  Mesh& cut_mesh = CutMesh(mi);

  int numnode = nids.size();
  if (numnode != xyz.N())
  {
    throw std::runtime_error("node coordiante number mismatch");
  }

  //   PointSet nodalpoints;

  // make sure all nodes are there
  for (int i = 0; i < numnode; ++i)
  {
    cut_mesh.GetNode(nids[i], &xyz(0, i));
    //     nodalpoints.insert( n->point() );
    //     if ( n==NULL )
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
  return cut_mesh_[mi]->CreateSide(sid, nids, distype);
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection in the self cut           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::BuildSelfCutTree()
{
  Mesh& cm = CutMesh();

  cm.BuildSelfCutTree();
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::BuildStaticSearchTree()
{
  Mesh& m = NormalMesh();

  m.BuildStaticSearchTree();
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion via COMBUST-code                        *
 * where dofsets and node positions                                                               *
 * have not to be computed, standard cut for cut_test                                             *
 * !!!IS JUST USED FOR CUT TESTS                                                                  *
 *                                                                                   schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::CutTest_Cut(bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype, INPAR::CUT::BCellGaussPts BCellgausstype,
    bool tetcellsonly, bool screenoutput, bool do_Cut_Positions_Dofsets)
{
  int mypid = 0;
  int mpi_is_running = 0;
  MPI_Initialized(&mpi_is_running);
  if (mpi_is_running)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
  }

  Status();

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

  BuildStaticSearchTree();

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

  Cut_CollisionDetection(include_inner, screenoutput);

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

  Cut_MeshIntersection(screenoutput);

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

    Cut_Positions_Dofsets(include_inner, screenoutput);

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

  // DumpGmshVolumeCells("CUT_vc", true);
  // DumpGmshIntegrationCells("CUT_intcells");

  Status(VCellgausstype);
}

/*------------------------------------------------------------------------------------------------*
 * handles cut sides which cut each other                                                         *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_SelfCut(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 2/6 --- Cut_SelfCut");

  if (myrank_ == 0 and screenoutput) IO::cout << "\t * 2/6 Cut_SelfCut ...";

  SelfCut selfcut(CutMesh());

  selfcut.PerformSelfCut();
}

/*------------------------------------------------------------------------------------------------*
 * detects if a side of the cut mesh possibly collides with an element of the background mesh     *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_CollisionDetection(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 3/6 --- Cut_CollisionDetection");

  if (myrank_ == 0 and screenoutput) IO::cout << "\t * 3/6 Cut_CollisionDetection ...";

  Mesh& m = NormalMesh();

  Mesh& cutmesh = CutMesh();

  m.SearchCollisions(cutmesh);
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI and XFLUIDFLUID where dofsets and node positions        *
 * have to be parallelized                                                           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_MeshIntersection(bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 4/6 --- Cut_MeshIntersection");

  if (myrank_ == 0 and screenoutput) IO::cout << "\t * 4/6 Cut_MeshIntersection ...";

  //----------------------------------------------------------

  Mesh& m = NormalMesh();

  m.FindCutPoints();
  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

}  // GEO::CUT::MeshIntersection::Cut_MeshIntersection

/*------------------------------------------------------------------------------------------------*
 * Routine for deciding the inside-outside position. This creates the dofset data,                *
 * also in parallel                                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_Positions_Dofsets(bool include_inner, bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 5/6 --- Cut_Positions_Dofsets (serial)");

  if (myrank_ == 0 and screenoutput) IO::cout << "\t * 5/6 Cut_Positions_Dofsets ...";

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

}  // GEO::CUT::MeshIntersection::Cut_Positions_Dofsets


/*--------------------------------------------------------------------------------------*
 * get the cut mesh's side based on side id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle* GEO::CUT::MeshIntersection::GetCutSide(int sid, int mi) const
{
  return cut_mesh_[mi]->GetSide(sid);
}


/*--------------------------------------------------------------------------------------*
 * status
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Status(INPAR::CUT::VCellGaussPts gausstype)
{
  // call status of parent intersection
  my::Status(gausstype);

#ifdef DEBUG
  for (std::vector<Teuchos::RCP<MeshHandle>>::iterator i = cut_mesh_.begin(); i != cut_mesh_.end();
       ++i)
  {
    MeshHandle& cut_mesh_handle = **i;
    Mesh& cut_mesh = cut_mesh_handle.LinearMesh();
    cut_mesh.Status();
  }

#ifdef DEBUGCUTLIBRARY
  int count = 0;
  for (std::vector<Teuchos::RCP<MeshHandle>>::iterator i = cut_mesh_.begin(); i != cut_mesh_.end();
       ++i)
  {
    MeshHandle& cut_mesh_handle = **i;
    Mesh& cut_mesh = cut_mesh_handle.LinearMesh();
    std::stringstream str;
    str << "cut_mesh" << count << ".pos";
    cut_mesh.DumpGmsh(str.str().c_str());
    count++;
  }
#endif
#endif
}
