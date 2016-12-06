/*---------------------------------------------------------------------*/
/*!
\file cut_cycle.H

\brief a cylcle of points (basic to create facets)

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_CYCLE_H
#define CUT_CYCLE_H

#include "cut_utils.H"

// forward declaration for the operator
namespace GEO
{
  namespace CUT
  {
    class Cycle;
  }
}

// must be here, as it is supposed to be a friend of the Cycle-class
std::ostream & operator<<( std::ostream & stream, const GEO::CUT::Cycle & cycle );

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Edge;
    class Side;
    class Element;

    /*!
    \brief Contains closed cycle of points. A utility class for facet creation from this cycle of points
     */
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

      /*!
      \brief Returns true if the set of lines given a input contains the cycle of points.
      This means that a facet can be created form these lines
       */
      static bool MakeCycle( const point_line_set & lines, Cycle & cycle );

      /*!
      \brief Returns true if the cycle of points are suitable for creating a facet
       */
      bool IsValid() const;

      /*!
      \brief Returns true if all the points in the cylce falls within the element
       */
      bool IsCut( Element * element ) const;

      /*!
      \brief Returns all the points in the cycle
       */
      const std::vector<Point*> & operator()() const
      {
        return points_;
      }
      /*!
      \brief add all lines from this cycle to the set of lines
      */
      void Add( point_line_set & lines ) const;

      /*!
      \brief Find the common edges among the cycle of points
       */
      void CommonEdges( plain_edge_set & edges ) const;

      /*!
      \brief Find the common edges among the cycle of points. NOT USED NOW
       */
      void CommonSides( plain_side_set & sides ) const;

      /*!
      \brief Get only the sides that are cut by any one of the points in the cycle
       */
      void Intersection( plain_side_set & sides ) const;

      /*!
      \brief Check whether these two cycles are one and the same
       */
      bool Equals( const Cycle & other );

      /*!
      \brief Delete the specified point from the cycle, and make the resulting cycle valid
       */
      void DropPoint( Point * p );

      void TestUnique();

      /*!
      \brief Add the specified point to the cycle
        */
      void push_back( Point * p ) { points_.push_back( p ); }

      /*!
      \brief Reserve "s" number of points in the cycle
        */
      void reserve( unsigned s ) { points_.reserve( s ); }

      /*!
      \brief Delete all the points in the cycle
        */
      void clear() { points_.clear(); }

      /*!
      \brief Returns the number of points
        */
      unsigned size() const { return points_.size(); }

      /*!
      \brief Reverse the order of storing the points
       */
      void reverse();

      void swap( Cycle & other ) { std::swap( points_, other.points_ ); }

      friend std::ostream & ::operator<<( std::ostream & stream, const Cycle & cycle );

      void GnuplotDump( std::ostream & stream ) const;

    private:
      std::vector<Point*> points_;
    };

  }
}

#endif
