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

#include "cut_side.H"
#include "cut_levelsetside.H"
#include "cut_levelsetintersection.H"

/*-----------------------------------------------------------------------------------------*
 * constructur for LevelSetIntersection class
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::LevelSetIntersection::LevelSetIntersection(int myrank, bool create_side)
:   ParentIntersection( myrank )
{
  if(create_side)
  {
    side_ = Teuchos::rcp( new LevelSetSide( 1 ) );
  }
  else side_ = Teuchos::null;
}


/*-----------------------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for quadratic sides)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::AddCutSide( int levelset_sid )
{
  if(side_ != Teuchos::null) dserror("currently only one levelset-side supported");

  // create the levelset-side
  side_ = Teuchos::rcp( new LevelSetSide( levelset_sid ) );
}


/*-----------------------------------------------------------------------------------------*
 * add this background element if it is cut. Which implies that the level set function of
 * the element has values which are positive and negative.
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle * GEO::CUT::LevelSetIntersection::AddElement(
    int eid,
    const std::vector<int> & nids,
    const Epetra_SerialDenseMatrix & xyz,
    DRT::Element::DiscretizationType distype,
    const double * lsv,
    const bool lsv_only_plus_domain
)
{
  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

//  bool ltz = false;
//  bool gtz = false;
//
//  // make sure element has LSV +ve and -ve in one of its nodes
//  // ensures this is a cut element
//  for ( int i=0; i<numnode; ++i )
//  {
//    if ( lsv[i] <=   TOLERANCE )
//      ltz = true;
//    if ( lsv[i] >= - TOLERANCE )
//      gtz = true;
//  }

  // add all cut elements (different signs of levelset values) OR
  // if only plus domain is a physical field we have to add also elements with only negative values (as they are not allowed to carry DOFS at the end)
//  if ( (ltz and gtz) or (lsv_only_plus_domain and ltz) )
  {
    // add all nodes to mesh
    for ( int i=0; i<numnode; ++i )
    {
      NormalMesh().GetNode( nids[i], &xyz( 0, i ), lsv[i] );
    }

    // create element
    return mesh_.CreateElement( eid, nids, distype );
  }

  return NULL;
}


/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel Level Set Cut where dofsets and node positions               *
 * have to be parallelized                                                           winter 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut_Mesh( bool screenoutput )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut" );

  if(myrank_==0 and screenoutput) IO::cout << "\n\t * 1/3 Cut ...";

  Mesh & m = NormalMesh();

  //m.Status();
  m.Cut( *side_ );

  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

}//GEO::CUT::LevelSetIntersection::Cut_Mesh



/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion where dofsets and node positions        *
 * have not to be computed, standard cut for cut_test (Only used for cut test)       winter 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut( bool include_inner, bool screenoutput )
{
  //######################################################################################
  //STEP 1/3 CUT THE MESH
  //######################################################################################

  //m.Status();
  Cut_Mesh(screenoutput);


  //######################################################################################
  //STEP 2/3 ASSIGN DOFS
  //######################################################################################

  Mesh & m = NormalMesh();

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
  m.CreateIntegrationCells( 0 );
  //m.RemoveEmptyVolumeCells();

#ifdef DEBUGCUTLIBRARY
  m.DumpGmsh( "mesh.pos" );
//  m.DumpGmshVolumeCells( "volumecells" );
  m.DumpGmshIntegrationCells( "integrationcells.pos" );
#endif

  m.TestElementVolume( true );

#ifdef DEBUGCUTLIBRARY
 // m.TestVolumeSurface(); //Broken test, needs to be fixed for proper usage.
  m.TestFacetArea();
#endif

//######################################################################################

}//GEO::CUT::LevelSetIntersection::Cut


