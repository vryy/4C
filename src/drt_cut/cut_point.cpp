
#include <string>

#include "cut_point.H"
#include "cut_point_impl.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_line.H"
#include "cut_facet.H"
#include "cut_mesh.H"

GEO::CUT::Point * GEO::CUT::Point::NewPoint( Mesh & mesh, const double * x, double t, Edge * cut_edge, Side * cut_side )
{
  Point * p = mesh.NewPoint( x, cut_edge, cut_side );
  p->Position( Point::oncutsurface );
  p->t( cut_edge, t );
  return p;
}

GEO::CUT::Point * GEO::CUT::Point::InsertCut( Edge * cut_edge, Side * cut_side, Node * n )
{
  Point * p = n->point();
  const std::set<Edge*> & edges = n->Edges();
  for ( std::set<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    p->AddEdge( e );
  }
  p->AddEdge( cut_edge );
  p->AddSide( cut_side );
  p->Position( Point::oncutsurface );
  return p;
}

GEO::CUT::Point::Point( unsigned pid, const double * x, Edge * cut_edge, Side * cut_side )
  : pid_( pid ),
    position_( undecided )
{
  std::copy( x, x+3, x_ );

  if ( cut_edge!=NULL )
  {
    cut_edges_.insert( cut_edge );

    // copy all sides at the edge to the set of cutted sides
    std::copy( cut_edge->Sides().begin(), cut_edge->Sides().end(),
               std::inserter( cut_sides_, cut_sides_.begin() ) );
  }
  if ( cut_side!=NULL )
  {
    cut_sides_.insert( cut_side );
  }
}

void GEO::CUT::Point::AddEdge( Edge* cut_edge )
{
  cut_edges_.insert( cut_edge );

  // revers add
  cut_edge->AddPoint( this );

  // copy all sides at the edge to the set of cutted sides
  std::copy( cut_edge->Sides().begin(), cut_edge->Sides().end(),
             std::inserter( cut_sides_, cut_sides_.begin() ) );
}

GEO::CUT::Edge * GEO::CUT::Point::CommonEdge( Point * other )
{
  Edge * found = NULL;
  for ( std::set<Edge*>::iterator i=cut_edges_.begin(); i!=cut_edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( other->IsCut( e ) )
    {
      if ( found==NULL )
      {
        found = e;
      }
      else
      {
        throw std::runtime_error( "not unique" );
      }
    }
  }
  return found;
}

// std::vector<GEO::CUT::Edge*> GEO::CUT::Point::CutEdges( Point * other )
// {
//   std::vector<Edge*> matches;
//   for ( std::set<Edge*>::iterator i=cut_edges_.begin(); i!=cut_edges_.end(); ++i )
//   {
//     Edge * e = *i;
//     if ( other->IsCut( e ) )
//     {
//       matches.push_back( e );
//     }
//   }
//   return matches;
// }

void GEO::CUT::Point::CutEdge( Side * side, Line * other_line, std::vector<Edge*> & matches )
{
  for ( std::set<Edge*>::iterator i=cut_edges_.begin(); i!=cut_edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( e->AtSide( side ) )
    {
      if ( e->BeginNode()->point()==this or e->EndNode()->point()==this )
      {
        if ( not other_line->OnEdge( e ) )
        {
          matches.push_back( e );
        }
      }
      else
      {
        matches.push_back( e );
      }
    }
  }
}

GEO::CUT::Line * GEO::CUT::Point::CutLine( const PointLineFilter & filter, bool unique )
{
  Line * line_found = NULL;
  for ( std::set<Line*>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    Line * line = *i;
    if ( filter( line ) )
    {
      if ( line_found==NULL )
      {
        line_found = line;
        if ( not unique )
        {
          break;
        }
      }
      else
      {
        throw std::runtime_error( "not unique" );
      }
    }
  }
  return line_found;
}

GEO::CUT::Line * GEO::CUT::Point::CutLine( Line * line, const PointLineFilter & filter, bool unique )
{
  ExcludeLineFilter f( line, filter );
  return CutLine( f, unique );
}

void GEO::CUT::Point::CutLines( const PointLineFilter & filter, std::set<Line*> & cut_lines )
{
  for ( std::set<Line*>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    Line * line = *i;
    if ( filter( line ) )
    {
      cut_lines.insert( line );
    }
  }
}

void GEO::CUT::Point::CutLines( Side * side, std::set<Line*> & cut_lines )
{
  SideCutFilter filter( side );
  CutLines( filter, cut_lines );
}

double GEO::CUT::Point::t( Edge* edge )
{
  std::map<Edge*, double>::iterator i=t_.find( edge );
  if ( i==t_.end() )
  {
    Point * p1 = edge->BeginNode()->point();
    Point * p2 = edge->EndNode()->point();

    LINALG::Matrix<3, 1> x;
    LINALG::Matrix<3, 1> x1;
    LINALG::Matrix<3, 1> x2;

    Coordinates( x.A() );
    p1->Coordinates( x1.A() );
    p2->Coordinates( x2.A() );

    x.Update( -1, x1, 1 );
    x2.Update( -1, x1, 1 );

    double l1 = x.Norm2();
    double l2 = x2.Norm2();

    if ( fabs( l2 )<TOLERANCE )
    {
      throw std::runtime_error( "edge with no length" );
    }

    double z = l1/l2;

    x.Update( -z, x2, 1 );
    if ( x.Norm2() > MINIMALTOL )
    {
      throw std::runtime_error( "point not on edge, no edge position" );
    }

    double t = 2*z - 1;
    t_[edge] = t;
    return t;
  }
  return i->second;
}

void GEO::CUT::Point::Intersection( std::set<Side*> & sides )
{
  std::set<Side*> intersection;
  std::set_intersection( cut_sides_.begin(), cut_sides_.end(),
                         sides.begin(), sides.end(),
                         std::inserter( intersection, intersection.begin() ) );
  std::swap( sides, intersection );
}

void GEO::CUT::Point::Intersection( std::set<Facet*> & facets )
{
  std::set<Facet*> intersection;
  std::set_intersection( facets_.begin(), facets_.end(),
                         facets.begin(), facets.end(),
                         std::inserter( intersection, intersection.begin() ) );
  std::swap( facets, intersection );
}

bool GEO::CUT::Point::NodalPoint( const std::vector<Node*> & nodes ) const
{
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    if ( n->point()==this )
    {
      return true;
    }
  }
  return false;
}

GEO::CUT::Node * GEO::CUT::Point::CutNode()
{
  for ( std::set<Edge*>::iterator i=cut_edges_.begin(); i!=cut_edges_.end(); ++i )
  {
    Edge * e = *i;
    const std::vector<Node*> & nodes = e->Nodes();
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      if ( n->point()==this )
      {
        return n;
      }
    }
  }
  return NULL;
}

void GEO::CUT::Point::Position( Point::PointPosition pos )
{
  if ( position_ != pos )
  {
    position_ = pos;
    if ( pos==Point::outside or pos==Point::inside )
    {
      for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
      {
        Facet * f = *i;
        f->Position( pos );
      }
    }
  }
}

GEO::CUT::Side * GEO::CUT::Point::CutSide( Side * side, Point * other )
{
  Side * found_side = NULL;
  for ( std::set<Side*>::iterator i=cut_sides_.begin(); i!=cut_sides_.end(); ++i )
  {
    Side * s = *i;
    if ( s != side and other->IsCut( s ) )
    {
      if ( found_side==NULL )
      {
        found_side = s;
      }
      else
      {
        throw std::runtime_error( "side not unique" );
      }
    }
  }
  return found_side;
}

bool GEO::CUT::SideCutFilter::operator()( Line * line ) const
{
  return line->IsInternalCut( side_ );
}

bool GEO::CUT::SideElementCutFilter::operator()( Line * line ) const
{
  //return line->IsCut( side_ ) and line->IsCut( element_ );
  return line->IsInternalCut( side_ ) and line->IsCut( element_ );
}

bool GEO::CUT::SideSideCutFilter::operator()( Line * line ) const
{
  return line->IsCut( side1_ ) and line->IsCut( side2_ );
}

GEO::CUT::Edge * GEO::CUT::Point::CommonCutEdge( Side * side )
{
  for ( std::set<Edge*>::iterator i=cut_edges_.begin(); i!=cut_edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( e->IsCut( side ) )
    {
      return e;
    }
  }
  return NULL;
}

