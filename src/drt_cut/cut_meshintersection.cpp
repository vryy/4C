
#include "cut_meshintersection.H"

void GEO::CUT::MeshIntersection::AddElement( int eid,
                                             const std::vector<int> & nids,
                                             const Epetra_SerialDenseMatrix & xyz,
                                             DRT::Element::DiscretizationType distype )
{
  if ( cut_mesh_.WithinBB( xyz ) )
  {
    int numnode = nids.size();
    if ( numnode != xyz.N() )
    {
      throw std::runtime_error( "node coordiante number mismatch" );
    }

    // make sure all nodes are there
    for ( int i=0; i<numnode; ++i )
    {
      mesh_.GetNode( nids[i], &xyz( 0, i ) );
    }

    // create element
    mesh_.CreateElement( eid, nids, distype );
  }
}

void GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                             const std::vector<int> & nids,
                                             const Epetra_SerialDenseMatrix & xyz,
                                             DRT::Element::DiscretizationType distype )
{
  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

  // make sure all nodes are there
  for ( int i=0; i<numnode; ++i )
  {
    cut_mesh_.GetNode( nids[i], &xyz( 0, i ) );
  }

  // create side
  cut_mesh_.CreateSide( sid, nids, distype );
}

void GEO::CUT::MeshIntersection::Cut( CellGenerator * generator )
{
  cut_mesh_.FillComplete();
  mesh_.FillComplete();

  // loop cut sides and cut against elements at the same position in space
  cut_mesh_.Cut( mesh_ );

  mesh_.Status();
  cut_mesh_.Status();

  mesh_.MakeFacets();
  mesh_.FindNodePositions();

  mesh_.DumpGmsh( "mesh" );
  cut_mesh_.DumpGmsh( "cut_mesh" );

  mesh_.GenerateTetgen( generator );
}

GEO::CUT::Side * GEO::CUT::MeshIntersection::GetCutSides( int sid )
{
  const std::vector<GEO::CUT::Side*> & cut_sides = cut_mesh_.GetSides( sid );
  if ( cut_sides.size()==1 )
  {
    return cut_sides[0];
  }
  if ( cut_sides.size()>1 )
  {
    // the first entry has to be the quadratic element, the following entries
    // its linear parts
    return cut_sides[0];
  }
  throw std::runtime_error( "no such side" );
}
