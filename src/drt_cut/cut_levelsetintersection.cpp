
#include "cut_levelsetintersection.H"

void GEO::CUT::LevelSetIntersection::AddElement( int eid,
                                                 const std::vector<int> & nids,
                                                 const Epetra_SerialDenseMatrix & xyz,
                                                 const double * lsv,
                                                 DRT::Element::DiscretizationType distype )
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
}

void GEO::CUT::LevelSetIntersection::Cut( bool include_inner )
{
  Mesh & m = NormalMesh();

  m.Cut( side_ );

  m.Status();

  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

  if ( options_.FindPositions() )
  {
    m.FindLSNodePositions();
    m.FindNodalDOFSets( include_inner );
  }

  m.CreateIntegrationCells( 0, true );
  //m.RemoveEmptyVolumeCells();

#ifdef DEBUGCUTLIBRARY
  m.DumpGmsh( "mesh.pos" );
  //m.DumpGmshVolumeCells( "volumecells" );
  m.DumpGmshIntegrationCells( "integrationcells.pos" );
#endif

  m.TestElementVolume( true );
}

GEO::CUT::ElementHandle * GEO::CUT::LevelSetIntersection::GetElement( int eid )
{
  return mesh_.GetElement( eid );
}
