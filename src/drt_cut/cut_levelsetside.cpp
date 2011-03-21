
#include "cut_levelsetside.H"
#include "cut_mesh.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

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

bool GEO::CUT::LevelSetSide::FindLevelSetCutLines( Mesh & mesh, Element * element, Side & side, const std::set<Point*> & cut )
{
  switch ( side.Shape() )
  {
  case DRT::Element::tri3:
    return false;
  case DRT::Element::quad4:
  {
    if ( cut.size()!=4 )
      return false;

    std::vector<Point*> edge_points;
    edge_points.reserve( 4 );
    const std::vector<Edge*> & edges = side.Edges();
    std::set<Point*> cut_points( cut );
    for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      for ( std::set<Point*>::iterator i=cut_points.begin(); i!=cut_points.end(); )
      {
        Point * p = *i;
        if ( p->IsCut( e ) )
        {
          edge_points.push_back( p );
          cut_points.erase( i++ );
          break;
        }
        else
        {
          ++i;
        }
      }
    }
    if ( edge_points.size()!=4 or cut_points.size()!=0 )
    {
      throw std::runtime_error( "failed to associate cut points with edges" );
    }

    // find levelset value at side center
    LINALG::Matrix<4,1> lsv;
    LINALG::Matrix<4,1> funct;
    DRT::UTILS::shape_function_2D( funct, 0., 0., DRT::Element::quad4 );
    const std::vector<Node*> & nodes = side.Nodes();
    std::vector<int> zero_positions;
    for ( unsigned i=0; i<4; ++i )
    {
      lsv( i ) = nodes[i]->LSV();
      if ( lsv( i )==0 )
        zero_positions.push_back( i );
    }

    LINALG::Matrix<1,1> midlsv;
    midlsv.MultiplyTN( lsv, funct );

    for ( unsigned i=1; i<zero_positions.size(); ++i )
    {
      int diff = zero_positions[i] - zero_positions[i-1];
      if ( diff == 1 or diff == 3 )
      {
        throw std::runtime_error( "cannot have adjacent zeros" );
      }
    }

    bool negativemiddle;

    if ( midlsv( 0 ) < 0 )
    {
      negativemiddle = true;
    }
    else if ( midlsv( 0 ) > 0 )
    {
      negativemiddle = false;
    }
    else
    {
      //throw std::runtime_error( "side center at interface on multiple cuts: undefined" );
      return false;
    }

    if ( lsv( 0 ) <= 0 and
         lsv( 1 ) >= 0 and
         lsv( 2 ) <= 0 and
         lsv( 3 ) >= 0 )
    {
      if ( negativemiddle )
      {
        mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
        mesh.NewLine( edge_points[2], edge_points[3], &side, this, element );
      }
      else
      {
        mesh.NewLine( edge_points[0], edge_points[3], &side, this, element );
        mesh.NewLine( edge_points[2], edge_points[1], &side, this, element );
      }
      return true;
    }
    else if ( lsv( 0 ) >= 0 and
              lsv( 1 ) <= 0 and
              lsv( 2 ) >= 0 and
              lsv( 3 ) <= 0 )
    {
      if ( negativemiddle )
      {
        mesh.NewLine( edge_points[0], edge_points[3], &side, this, element );
        mesh.NewLine( edge_points[2], edge_points[1], &side, this, element );
      }
      else
      {
        mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
        mesh.NewLine( edge_points[2], edge_points[3], &side, this, element );
      }
      return true;
    }
    else
    {
      throw std::runtime_error( "illegal levelset pattern" );
    }

    return false;
  }
  default:
    throw std::runtime_error( "unsupported side shape" );
  }
}
