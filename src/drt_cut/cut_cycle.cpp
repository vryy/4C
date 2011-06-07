
#include "cut_cycle.H"
#include "cut_point.H"
#include "cut_edge.H"


bool GEO::CUT::Cycle::IsValid() const
{
  if ( points_.size() < 3 )
    return false;

  // ignore cycles with all points on one and the same edge
  {
    std::set<Edge*> edges;
    CommonEdges( edges );
    if ( edges.size() > 0 )
    {
      return false;
    }
  }

  return true;
}

bool GEO::CUT::Cycle::IsCut( Element * element ) const
{
  for ( std::vector<Point*>::const_iterator i=points_.begin(); i!=points_.end(); ++i )
  {
    Point * p = *i;
    if ( not p->IsCut( element ) )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::Cycle::CommonEdges( std::set<Edge*> & edges ) const
{
  std::vector<Point*>::const_iterator i = points_.begin();
  if ( i!=points_.end() )
  {
    edges = ( *i )->CutEdges();
    for ( ++i; i!=points_.end(); ++i )
    {
      Point * p = *i;
      p->Intersection( edges );
      if ( edges.size()==0 )
      {
        break;
      }
    }
  }
  else
  {
    edges.clear();
  }
}

void GEO::CUT::Cycle::CommonSides( std::set<Side*> & sides ) const
{
  std::vector<Point*>::const_iterator i = points_.begin();
  if ( i!=points_.end() )
  {
    sides = ( *i )->CutSides();
    for ( ++i; i!=points_.end(); ++i )
    {
      Point * p = *i;
      p->Intersection( sides );
      if ( sides.size()==0 )
      {
        break;
      }
    }
  }
  else
  {
    sides.clear();
  }
}

void GEO::CUT::Cycle::Intersection( std::set<Side*> & sides ) const
{
  for ( std::vector<Point*>::const_iterator i=points_.begin(); i!=points_.end(); ++i )
  {
    Point * p = *i;
    p->Intersection( sides );
    if ( sides.size()==0 )
      break;
  }
}

bool GEO::CUT::Cycle::Equals( const Cycle & other )
{
  if ( size()!=other.size() )
  {
    return false;
  }

  for ( std::vector<Point*>::const_iterator i=other.points_.begin(); i!=other.points_.end(); ++i )
  {
    Point * p = *i;
    //if ( not std::binary_search( sorted.begin(), sorted.end(), p ) )

    if ( std::count( points_.begin(), points_.end(), p )!=1 )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::Cycle::DropPoint( Point * p )
{
  std::vector<Point*>::iterator j = std::find( points_.begin(), points_.end(), p );
  if ( j!=points_.end() )
  {
    std::vector<Point*>::iterator prev = j==points_.begin() ? points_.end() : j;
    std::advance( prev, -1 );
    std::vector<Point*>::iterator next = j;
    std::advance( next, 1 );
    if ( next==points_.end() )
      next = points_.begin();
    if ( *prev == *next )
    {
      if ( next > j )
      {
        points_.erase( next );
        points_.erase( j );
      }
      else
      {
        points_.erase( j );
        points_.erase( next );
      }
    }
    else
    {
      points_.erase( j );
    }
  }
}

void GEO::CUT::Cycle::TestUnique()
{
  PointSet c_copy;
  c_copy.insert( points_.begin(), points_.end() );
  if ( points_.size()!=c_copy.size() )
    throw std::runtime_error( "double point in cycle" );
}

namespace GEO
{
  namespace CUT
  {
    std::ostream & operator<<( std::ostream & stream, const Cycle & cycle )
    {
      std::copy( cycle.points_.begin(), cycle.points_.end(), std::ostream_iterator<Point*>( stream, " " ) );
      return stream << "\n";
    }
  }
}

void GEO::CUT::Cycle::GnuplotDump( std::ostream & stream ) const
{
  for ( unsigned i=0; i!=points_.size(); ++i )
  {
    Point * p1 = points_[i];
    Point * p2 = points_[( i+1 ) % points_.size()];

    p1->Plot( stream );
    p2->Plot( stream );
    stream << "\n\n";
  }
}
