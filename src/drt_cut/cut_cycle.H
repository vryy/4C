#ifndef CUT_CYCLE_H
#define CUT_CYCLE_H

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <set>

#include "cut_utils.H"

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Edge;
    class Side;
    class Element;

#ifdef CUT_USE_SORTED_VECTOR
    typedef sorted_vector<std::pair<Point*, Point*> > point_line_set;
#else
    typedef std::set<std::pair<Point*, Point*> > point_line_set;
#endif

    class Cycle
    {
    public:

      Cycle()
      {
      }

      Cycle( const std::vector<Point*> & points )
        : points_( points )
      {
      }

      static bool MakeCycle( const point_line_set & lines, Cycle & cycle );

      bool IsValid() const;

      bool IsCut( Element * element ) const;

      const std::vector<Point*> & operator()() const
      {
        return points_;
      }

      void Add( point_line_set & lines ) const;

      void CommonEdges( plain_edge_set & edges ) const;

      void CommonSides( plain_side_set & sides ) const;

      void Intersection( plain_side_set & sides ) const;

      bool Equals( const Cycle & other );

      void DropPoint( Point * p );

      void TestUnique();

      void push_back( Point * p ) { points_.push_back( p ); }

      void reserve( unsigned s ) { points_.reserve( s ); }

      void clear() { points_.clear(); }

      unsigned size() const { return points_.size(); }

      void swap( Cycle & other ) { std::swap( points_, other.points_ ); }

      friend std::ostream & operator<<( std::ostream & stream, const Cycle & cycle );

      void GnuplotDump( std::ostream & stream ) const;

    private:
      std::vector<Point*> points_;
    };

  }
}

#endif
