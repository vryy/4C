
//#include "../drt_geometry/intersection_templates.H"

#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"

#include <string>
#include <stack>

#include "cut_edge.H"
#include "cut_levelsetside.H"


bool GEO::CUT::Edge::FindCutPoints( Mesh & mesh,
                                    Element * element,
                                    Side & side,
                                    Side & other )
{
#if 1
  bool cut = false;
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( &other ) )
    {
      cut = true;
      p->AddElement( element );
    }
  }
  if ( cut )
  {
    return true;
  }
#endif

  // test for the cut of edge and side

  std::set<Point*, PointPidLess> cut_points;
  other.Cut( mesh, *this, cut_points );
  for ( std::set<Point*, PointPidLess>::iterator i=cut_points.begin();
        i!=cut_points.end();
        ++i )
  {
    Point * p = *i;

    p->AddEdge( this );
    p->AddElement( element );

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->AddSide( &side );
    p->AddSide( &other );
    AddPoint( p );
  }

  return cut_points.size()>0;
}

void GEO::CUT::Edge::GetCutPoints( Element * element,
                                   Side & side,
                                   Side & other,
                                   std::set<Point*> & cuts )
{
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( &other ) and p->IsCut( element ) )
    {
      cuts.insert( p );
    }
  }
}


void GEO::CUT::Edge::AddPoint( Point* cut_point )
{
  cut_points_.insert( cut_point );
}

void GEO::CUT::Edge::CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points )
{
  if ( BeginNode()==edge_start and EndNode()==edge_end )
  {
    std::set<Point*, PointPositionLess>::iterator bi = cut_points_.begin();
    std::set<Point*, PointPositionLess>::iterator ei = cut_points_.end();

    edge_points.clear();
    edge_points.assign( bi, ei );
  }
  else if ( EndNode()==edge_start and BeginNode()==edge_end )
  {
    std::set<Point*, PointPositionLess>::reverse_iterator bi = cut_points_.rbegin();
    std::set<Point*, PointPositionLess>::reverse_iterator ei = cut_points_.rend();

    edge_points.clear();
    edge_points.assign( bi, ei );
  }
  else
  {
    throw std::runtime_error( "not a valid edge" );
  }
}

void GEO::CUT::Edge::CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points )
{
  SideCutFilter filter( side );
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    std::set<Line*> cut_lines;
    p->CutLines( filter, cut_lines );
    if ( cut_lines.size()>0 )
    {
      cut_points.insert( p );
    }
  }
}

void GEO::CUT::Edge::CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line )
{
  std::set<Point*, PointPositionLess>::iterator bi = cut_points_.find( begin );
  std::set<Point*, PointPositionLess>::iterator ei = cut_points_.find( end );

  double bt = begin->t( this );
  double et = end->t( this );

  if ( bt < et )
  {
    if ( bi!=cut_points_.end() )
    {
      ++bi;
    }
    std::copy( bi, ei, std::back_inserter( line ) );
  }
  else if ( bt > et )
  {
    if ( ei!=cut_points_.end() )
    {
      ++ei;
    }
    std::copy( ei, bi, std::back_inserter( line ) );
  }
  else
  {
    if ( begin!=end )
    {
      throw std::runtime_error( "different points at the same place" );
    }
  }
}

void GEO::CUT::Edge::CutPointsInside( Element * element, std::vector<Point*> & line )
{
  for ( std::set<Point*, PointPositionLess>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( element ) )
    {
      line.push_back( p );
    }
  }
}

bool GEO::CUT::Edge::IsCut( Side * side )
{
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( side ) )
    {
      return true;
    }
  }
  return false;
}

GEO::CUT::Point* GEO::CUT::Edge::NodeInElement( Element * element, Point * other )
{
  Point * p = BeginNode()->point();
  if ( p!=other and p->IsCut( element ) )
  {
    return p;
  }
  p = EndNode()->point();
  if ( p!=other and p->IsCut( element ) )
  {
    return p;
  }
  return NULL;
}


void GEO::CUT::Edge::Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts )
{
  Intersection<DRT::Element::line2, DRT::Element::tri3> inter( mesh, *this, side );
  inter.Intersect( cuts );
}

void GEO::CUT::Edge::Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts )
{
  Intersection<DRT::Element::line2, DRT::Element::quad4> inter( mesh, *this, side );
  inter.Intersect( cuts );
}

void GEO::CUT::Edge::LevelSetCut( Mesh & mesh, LevelSetSide & side, std::set<Point*, PointPidLess> & cuts )
{
  double blsv = BeginNode()->LSV();
  double elsv = EndNode()  ->LSV();

  if ( ( blsv < 0.0 and elsv > 0.0 ) or
       ( blsv > 0.0 and elsv < 0.0 ) )
  {
    double t = blsv / ( blsv-elsv );

    LINALG::Matrix<3,1> x1;
    LINALG::Matrix<3,1> x2;
    BeginNode()->Coordinates( x1.A() );
    EndNode()  ->Coordinates( x2.A() );

    LINALG::Matrix<3,1> x;
    x.Update( -1., x1, 1., x2, 0. );
    x.Update( 1., x1, t );
    Point * p = Point::NewPoint( mesh, x.A(), 2.*t-1., this, &side );
    cuts.insert( p );
  }
  else
  {
    if ( fabs( blsv ) < std::numeric_limits<double>::min() )
    {
      cuts.insert( Point::InsertCut( this, &side, BeginNode() ) );
    }
    if ( fabs( elsv ) < std::numeric_limits<double>::min() )
    {
      cuts.insert( Point::InsertCut( this, &side, EndNode() ) );
    }
  }
}
