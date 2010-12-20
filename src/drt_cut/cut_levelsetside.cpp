
#include "cut_levelsetside.H"

void GEO::CUT::LevelSetSide::Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
{
  edge.LevelSetCut( mesh, *this, cut_points );
}

void GEO::CUT::LevelSetSide::EdgeAt( double r, double s, std::vector<Edge*> & edges )
{
  throw std::runtime_error( "no edges on level set cut surface" );
}

bool GEO::CUT::LevelSetSide::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  throw std::runtime_error( "no local coordinates on level set cut surface" );
}

void GEO::CUT::LevelSetSide::MakeOwnedSideFacets( Mesh & mesh, const PointLineFilter & filter, std::set<Facet*> & facets )
{
  Side::MakeOwnedSideFacets( mesh, filter, facets );
}

void GEO::CUT::LevelSetSide::MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  Side::MakeSideCutFacets( mesh, element, facets );
}

void GEO::CUT::LevelSetSide::MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  Side::MakeInternalFacets( mesh, element, facets );
}

bool GEO::CUT::LevelSetSide::IsCut()
{
  return Side::IsCut();
}

void GEO::CUT::LevelSetSide::AddLine( Line* cut_line )
{
  Side::AddLine( cut_line );
}

GEO::CUT::Facet * GEO::CUT::LevelSetSide::FindFacet( const std::vector<Point*> & facet_points )
{
  return Side::FindFacet( facet_points );
}
