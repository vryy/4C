/*!-----------------------------------------------------------------------------------------------*
\file cut_meshintersection.cpp

\brief provides the specific functionality for cutting a mesh with other meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "cut_integrationcell.H"
#include "cut_selfcut.H"

#include "cut_meshintersection.H"

/*-----------------------------------------------------------------------------------------*
 * add this background element if it falls within the bounding box of cut mesh
 * If it is not within BB, this element is never cut
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle * GEO::CUT::MeshIntersection::AddElement( int eid,
                                                                  const std::vector<int> & nids,
                                                                  const Epetra_SerialDenseMatrix & xyz,
                                                                  DRT::Element::DiscretizationType distype )
{
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();

    // element is to be added only when it falls within bounding box
    // generated over cut mesh. otherwise this element is never cut
    if ( cut_mesh.WithinBB( xyz ) )
    {
      int numnode = nids.size();
      if ( numnode != xyz.N() )
      {
        throw std::runtime_error( "node coordiante number mismatch" );
      }

      // make sure all nodes are there
      for ( int i=0; i<numnode; ++i )
      {
        NormalMesh().GetNode( nids[i], &xyz( 0, i ) );
//         if ( n==NULL )
//         {
//           // if there is no node with that id but a node at the given
//           // location, the element is illegal and cannot be created
//           return;
//         }
      }

      // create element
      return mesh_.CreateElement( eid, nids, distype );
    }
  }
  return NULL;
}

/*-----------------------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for quadratic sides)
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                                               const std::vector<int> & nids,
                                                               DRT::Element::DiscretizationType distype,
                                                               int mi )
{
  // create side
  return cut_mesh_[mi]->CreateSide( sid, nids, distype );
}

/*-----------------------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for quadratic sides)
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                                               const std::vector<int> & nids,
                                                               const Epetra_SerialDenseMatrix & xyz,
                                                               DRT::Element::DiscretizationType distype,
                                                               int mi )
{
  Mesh & cut_mesh = CutMesh( mi );

  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

//   PointSet nodalpoints;

  // make sure all nodes are there
  for ( int i=0; i<numnode; ++i )
  {
    cut_mesh.GetNode( nids[i], &xyz( 0, i ) );
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
  return cut_mesh_[mi]->CreateSide( sid, nids, distype );
}

/*------------------------------------------------------------------------------------------------*
 * build the bounding volume tree for the collision detection in the context of the selfcut       *
 *                                                                                    wirtz 09/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::BuildBVTree()
{

  Mesh & cm = CutMesh();

  cm.BuildBVTree();

}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::BuildStaticSearchTree()
{

  Mesh & m = NormalMesh();

  m.BuildStaticSearchTree();

}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion where dofsets and node positions        *
 * have not to be computed, standard cut for cut_est                                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,
    INPAR::CUT::BCellGaussPts BCellgausstype,
    bool tetcellsonly,
    bool screenoutput)
{
  Status();

  // build the static search tree for the collision detection
  BuildStaticSearchTree();

  // handles cut sides which cut each other
  Cut_SelfCut( include_inner, screenoutput);

  // detects if a side of the cut mesh possibly collides with an element of the background mesh
  Cut_CollisionDetection( include_inner, screenoutput);

  // cut the mesh and create cutlines, facets, volumecells
  Cut_MeshIntersection( include_inner, screenoutput);

  // determine inside-outside position and dofset-data, parallel communication if required
  Cut_Positions_Dofsets( include_inner, screenoutput );

  // create integration points and/or subtetrahedralization
  Cut_Finalize( include_inner, VCellgausstype, BCellgausstype, false, tetcellsonly, screenoutput);

  // DumpGmshVolumeCells("CUT_vc", true);
  // DumpGmshIntegrationCells("CUT_intcells");

  Status(VCellgausstype);
}

/*------------------------------------------------------------------------------------------------*
 * handles cut sides which cut each other                                                         *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_SelfCut( bool include_inner,bool screenoutput)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 2/6 --- Cut_SelfCut" );

  if(myrank_==0 and screenoutput) IO::cout << "\t * 2/6 Cut_SelfCut ...";

  SelfCut selfcut( CutMesh() );

}

/*------------------------------------------------------------------------------------------------*
 * detects if a side of the cut mesh possibly collides with an element of the background mesh     *
 *                                                                                    wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_CollisionDetection( bool include_inner,bool screenoutput)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 3/6 --- Cut_CollisionDetection" );

  if(myrank_==0 and screenoutput) IO::cout << "\t * 3/6 Cut_CollisionDetection ...";

  Mesh & m = NormalMesh();

  Mesh & cutmesh = CutMesh();

  m.SearchCollisions(cutmesh);

}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI and XFLUIDFLUID where dofsets and node positions        *
 * have to be parallelized                                                           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_MeshIntersection( bool include_inner,bool screenoutput)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 4/6 --- Cut_MeshIntersection" );

  if(myrank_==0 and screenoutput) IO::cout << "\t * 4/6 Cut_MeshIntersection ...";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  Mesh & m = NormalMesh();

  m.FindCutPoints();
  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

  //----------------------------------------------------------

//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0  and screenoutput)
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }

} // GEO::CUT::MeshIntersection::Cut_MeshIntersection

/*------------------------------------------------------------------------------------------------*
 * Routine for deciding the inside-outside position. This creates the dofset data,                *
 * also in parallel                                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_Positions_Dofsets( bool include_inner , bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets (serial)" );

  if(myrank_==0 and screenoutput) IO::cout << "\t * 5/6 Cut_Positions_Dofsets ...";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------


  Mesh & m = NormalMesh();

  if ( options_.FindPositions() )
  {

    // find inside and outside positions of nodes, corresponding facets and volumecells using a DFS-algorithm
    m.FindNodePositions();

    // find the positions for all remaining facets ( and points, volumecells)
    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes, also in parallel
    m.FindNodalDOFSets( include_inner );

  }


//  //----------------------------------------------------------
//
//   const double t_diff = Teuchos::Time::wallTime()-t_start;
//   if ( myrank_ == 0 and screenoutput)
//   {
//     IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//   }
} //GEO::CUT::MeshIntersection::Cut_Positions_Dofsets


/*--------------------------------------------------------------------------------------*
 * get the cut mesh's side based on side id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::GetCutSide( int sid, int mi ) const
{
  return cut_mesh_[mi]->GetSide( sid );
}

void GEO::CUT::MeshIntersection::Status(INPAR::CUT::VCellGaussPts gausstype)
{
#ifdef DEBUG
  NormalMesh().Status();
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();
    cut_mesh.Status();
  }

#ifdef DEBUGCUTLIBRARY
  NormalMesh().DumpGmsh( "mesh.pos" );
  int count = 0;
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();
    std::stringstream str;
    str << "cut_mesh" << count << ".pos";
    cut_mesh.DumpGmsh( str.str().c_str() );
    count++;
  }

  //NormalMesh().DumpGmshVolumeCells( "volumecells" );
  if(gausstype==INPAR::CUT::VCellGaussPts_Tessellation)
  {
    DumpGmshIntegrationCells( "integrationcells.pos" );
    DumpGmshVolumeCells("volumecells.pos");
  }
  else
    DumpGmshVolumeCells("volumecells.pos");
#endif
#endif
}
