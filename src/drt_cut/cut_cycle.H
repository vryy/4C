#ifndef CUT_CYCLE_H
#define CUT_CYCLE_H

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <set>

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Edge;
    class Side;
    class Element;

    class Cycle
    {
    public:

      bool IsValid() const;

      bool IsCut( Element * element ) const;

      const std::vector<Point*> & operator()() const
      {
        return points_;
      }

      void CommonEdges( std::set<Edge*> & edges ) const;

      void CommonSides( std::set<Side*> & sides ) const;

      void Intersection( std::set<Side*> & sides ) const;

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
