/*!-----------------------------------------------------------------------------------------------*
\file cut_levelsetintersection.cpp

\brief provides the basic functionality for cutting a mesh with a level set function

<pre>
Maintainer: Benedikt Schott and Magnus Winter
            schott@lnm.mw.tum.de, winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "cut_levelsetside.H"

#include "cut_levelsetintersection.H"

/*-----------------------------------------------------------------------------------------*
 * constructur for LevelSetIntersecton class
 *-----------------------------------------------------------------------------------------*/
  GEO::CUT::LevelSetIntersection::LevelSetIntersection(int myrank)
  :   ParentIntersection( myrank )
  {
    side_ = Teuchos::rcp( new LevelSetSide( 1 ) );
  }

/*-----------------------------------------------------------------------------------------*
 * add this background element if it is cut. Which implies that the level set function of
 * the element has values which are positive and negative.
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::AddElement( int eid,
                                                 const std::vector<int> & nids,
                                                 const Epetra_SerialDenseMatrix & xyz,
                                                 const double * lsv,
                                                 DRT::Element::DiscretizationType distype,
                                                 const bool include_inner_ele)
{
  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

  bool ltz = false;
  bool gtz = false;

  // make sure element has LSV +ve and -ve in one of its nodes
  // ensures this is a cut element
  for ( int i=0; i<numnode; ++i )
  {
    if ( lsv[i] <=   TOLERANCE )
      ltz = true;
    if ( lsv[i] >= - TOLERANCE )
      gtz = true;
  }

  if ( ltz and gtz )
  {
    // add all nodes to mesh
    for ( int i=0; i<numnode; ++i )
    {
      NormalMesh().GetNode( nids[i], &xyz( 0, i ), lsv[i] );
    }

    // create element
    mesh_.CreateElement( eid, nids, distype );
  }
  else if(include_inner_ele and ltz)
  {
    // add all nodes to mesh
    for ( int i=0; i<numnode; ++i )
    {
      NormalMesh().GetNode( nids[i], &xyz( 0, i ), lsv[i] );
    }

    // create element
    mesh_.CreateElement( eid, nids, distype );
  }

}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion where dofsets and node positions        *
 * have not to be computed, standard cut for cut_test (Only used for cut test)       winter 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut( bool include_inner, bool screenoutput )
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut" );

//  if(myrank_==0 and screenoutput) IO::cout << "\n\t * 1/3 Cut ...";

  Mesh & m = NormalMesh();

  //STEP 1/3 CUT THE MESH
//######################################################################################
  //m.Status();
  m.Cut( *side_ );

  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

  //----------------------------------------------------------

//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0  and screenoutput)
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }
//######################################################################################
  //STEP 2/3 ASSIGN DOFS
//######################################################################################
  if ( options_.FindPositions() )
  {
    m.FindLSNodePositions();

    //=====================================================================

    m.FindNodalDOFSets( include_inner );

    //=====================================================================

    }
//######################################################################################
  //STEP 3/3 FINALIZE, ASSIGN INTEGRATIONRULES
//######################################################################################
  m.CreateIntegrationCells( 0, true );
  //m.RemoveEmptyVolumeCells();

#ifdef DEBUGCUTLIBRARY
  m.DumpGmsh( "mesh.pos" );
//  m.DumpGmshVolumeCells( "volumecells" );
  m.DumpGmshIntegrationCells( "integrationcells.pos" );
#endif

  m.TestElementVolume( true );
//######################################################################################

}//GEO::CUT::LevelSetIntersection::Cut

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel Level Set Cut where dofsets and node positions               *
 * have to be parallelized                                                           winter 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut_Mesh( bool include_inner, bool screenoutput )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut" );

  if(myrank_==0 and screenoutput) IO::cout << "\n\t * 1/3 Cut ...";

  Mesh & m = NormalMesh();

  //m.Status();
  m.Cut( *side_ );

  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

  //----------------------------------------------------------

//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0  and screenoutput)
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }

}//GEO::CUT::LevelSetIntersection::Cut_Mesh


void GEO::CUT::LevelSetIntersection::Status(INPAR::CUT::VCellGaussPts gausstype)
{
#ifdef DEBUG
  NormalMesh().Status();

#ifdef DEBUGCUTLIBRARY
  NormalMesh().DumpGmsh( "mesh.pos" );

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
