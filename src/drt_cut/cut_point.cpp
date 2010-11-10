
#include <string>

#include "cut_point.H"
#include "cut_point_impl.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_line.H"

GEO::CUT::Point::Point( unsigned pid, const double * x, Edge * cut_edge, Side * cut_side, bool nodalpoint )
  : pid_( pid ),
    nodalpoint_( nodalpoint )
{
  std::copy( x, x+3, x_ );

  if ( not nodalpoint )
  {
    if ( cut_edge==NULL or cut_side==NULL )
    {
      //throw std::runtime_error( "not a valid cut point" );
      return;
    }
    cut_edges_.insert( cut_edge );
    cut_sides_.insert( cut_side );

    // copy all sides at the edge to the set of cutted sides
    std::copy( cut_edge->Sides().begin(), cut_edge->Sides().end(),
               std::inserter( cut_sides_, cut_sides_.begin() ) );
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

void GEO::CUT::Point::CutLines( Side * side, std::set<Line*> & cut_lines )
{
  for ( std::set<Line*>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    Line * line = *i;
    if ( line->IsInternalCut( side ) )
    {
      cut_lines.insert( line );
    }
  }
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
