
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

  for ( int i=0; i<numnode; ++i )
  {
    if ( lsv[i] <= 0 )
      ltz = true;
    if ( lsv[i] >= 0 )
      gtz = true;
  }

  if ( ltz and gtz )
  {
    // make sure all nodes are there
    for ( int i=0; i<numnode; ++i )
    {
      NormalMesh().GetNode( nids[i], &xyz( 0, i ), lsv[i] );
    }

    // create element
    mesh_.CreateElement( eid, nids, distype );
  }
}

void GEO::CUT::LevelSetIntersection::Cut()
{
  Mesh & m = NormalMesh();

  m.Cut( side_ );
  //m.Status();
  m.MakeFacets();
  m.MakeVolumeCells();

  m.FindLSNodePositions();

  //m.DumpGmsh( "mesh" );

  m.CreateIntegrationCells();
}

GEO::CUT::ElementHandle * GEO::CUT::LevelSetIntersection::GetElement( int eid )
{
  return mesh_.GetElement( eid );
}
