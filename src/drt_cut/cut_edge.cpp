
#include "../drt_geometry/intersection_templates.H"

#ifdef QHULL
#undef PI
#ifdef TETGENINCLUDED
#include "../tetgen/tetgen.h"
#else
#include <tetgen.h>
#endif
#undef DOT
#endif

#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"

#include <string>
#include <stack>

#include "cut_edge.H"

void GEO::CUT::ConcreteEdge<DRT::Element::line2>::FillComplete( Mesh & mesh )
{
}

bool GEO::CUT::ConcreteEdge<DRT::Element::line2>::Cut( Mesh & mesh, Side & side, std::set<Point*, PointPidLess> & cuts )
{
#if 0
  // see if the cut is already there
  // this shadows double cuts since the cuts are shared between sides
  bool found_cut = false;
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point* cut_point = *i;
    if ( cut_point->IsCut( &side ) )
    {
      cuts.insert( cut_point );
      found_cut = true;
    }
  }
  if ( found_cut )
  {
    return true;
  }
#endif

  // test for the cut of edge and side

  std::set<Point*, PointPidLess> cut_points;
  side.Cut( mesh, *this, cut_points );
  for ( std::set<Point*, PointPidLess>::iterator i=cut_points.begin();
        i!=cut_points.end();
        ++i )
  {
    Point * cut_point = *i;
    cut_point->AddEdge( this );
    cuts.insert( cut_point );
  }

  return cut_points.size()>0;
}


void GEO::CUT::ConcreteEdge<DRT::Element::line2>::AddPoint( Point* cut_point )
{
  cut_points_.insert( cut_point );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line2>::CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points )
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

void GEO::CUT::ConcreteEdge<DRT::Element::line2>::CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points )
{
  for ( std::set<Point*>::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    std::set<Line*> cut_lines;
    p->CutLines( side, cut_lines );
    if ( cut_lines.size()>0 )
    {
      cut_points.insert( p );
    }
  }
}

void GEO::CUT::ConcreteEdge<DRT::Element::line2>::CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line )
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

GEO::CUT::Point* GEO::CUT::Edge::NodeInElement( Element * element, Point * other )
{
  Point * p = BeginNode()->point();
  if ( p==other or not element->PointInside( p ) )
  {
    p = EndNode()->point();
  }
  if ( p==other or not element->PointInside( p ) )
  {
    //throw std::runtime_error( "no node inside element" );
    return NULL;
  }
  return p;
}


void GEO::CUT::ConcreteEdge<DRT::Element::line2>::Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts )
{
  Intersection<DRT::Element::line2, DRT::Element::tri3> inter( mesh, *this, side );
  inter.Intersect( cuts );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line2>::Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts )
{
  Intersection<DRT::Element::line2, DRT::Element::quad4> inter( mesh, *this, side );
  inter.Intersect( cuts );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::FillComplete( Mesh & mesh )
{
  subedge1_ = mesh.GetEdge( BeginNode(), MiddleNode() );
  subedge2_ = mesh.GetEdge( MiddleNode(), EndNode() );
}

bool GEO::CUT::ConcreteEdge<DRT::Element::line3>::Cut( Mesh & mesh, Side & side, std::set<Point*, PointPidLess> & cuts )
{
  bool cut1 = subedge1_->Cut( mesh, side, cuts );
  bool cut2 = subedge1_->Cut( mesh, side, cuts );
  return cut1 or cut2;
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts )
{
  //subedge1_->Cut( mesh, side, cuts );
  //subedge2_->Cut( mesh, side, cuts );
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts )
{
  //subedge1_->Cut( mesh, side, cuts );
  //subedge2_->Cut( mesh, side, cuts );
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::AddPoint( Point* cut_point )
{
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points )
{
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points )
{
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteEdge<DRT::Element::line3>::CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line )
{
  throw std::runtime_error( "not supposed to end up here" );
}
