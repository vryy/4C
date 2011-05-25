#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

#include "cut_graph.H"
#include "cut_point.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Line;
  class Side;
  class Element;
  class Mesh;

class PointGraph
{
  class Graph : public GRAPH::Graph
  {
  public:

    void AddAll( Point * p1, Point * p2 );

    Point * GetPoint( int i ) { return all_points_[i]; }

    void PlotAllPoints( std::ostream & stream=std::cout );

    void PlotPoints( Element * element );

  private:

    std::map<int, Point*> all_points_;
  };

public:

  class PointGraphIterator
  {
  public:
    PointGraphIterator( Graph & graph, GRAPH::CycleList & facet_cycles, const GRAPH::CycleListIterator & i, Element * element )
      : graph_( graph ),
        facet_cycles_( facet_cycles ),
        i_( i ),
        element_( element )
    {
      NextActive();
    }

    PointGraphIterator & operator++()
    {
      ++i_;
      NextActive();
      return *this;
    }

    const std::vector<Point*> & operator*()
    {
      if ( points_.size() == 0 )
        throw std::runtime_error( "empty point list" );
      return points_;
    }

    bool operator!=( const PointGraphIterator & other )
    {
      return i_ != other.i_;
    }

  private:

    void NextActive()
    {
      while ( i_ != facet_cycles_.end() )
      {
        if ( FillPoints() )
        {
          return;
        }
        ++i_;
      }
    }

    bool FillPoints()
    {
      points_.clear();
      const std::vector<int> & facet_points = *i_;
      if ( facet_points.size() > 2 )
      {
        points_.reserve( facet_points.size() );
        for ( std::vector<int>::const_iterator i=facet_points.begin(); i!=facet_points.end(); ++i )
        {
          points_.push_back( graph_.GetPoint( *i ) );
          if ( not points_.back()->IsCut( element_ ) )
          {
            points_.clear();
            return false;
          }
        }
        return true;
      }
      return false;
    }

    Graph & graph_;
    GRAPH::CycleList & facet_cycles_;
    GRAPH::CycleListIterator i_;
    Element * element_;
    std::vector<Point*> points_;
  };

  typedef PointGraphIterator iterator;

  PointGraph( Mesh & mesh, Element * element, Side * side, bool inner );

  iterator begin() { return PointGraphIterator( graph_, facet_cycles_, facet_cycles_.begin(), element_ ); }

  iterator end()   { return PointGraphIterator( graph_, facet_cycles_, facet_cycles_.end(), element_ ); }

  unsigned size() { return facet_cycles_.size(); }

  unsigned ActiveCount() { return facet_cycles_.ActiveCount(); }

  Point * GetPoint( int i ) { return graph_.GetPoint( i ); }

  void Print()
  {
    graph_.Print();
    facet_cycles_.Print();
  }

private:

  void FillGraph( Element * element, Side * side, std::vector<int> & cycle );

  void AddFacetPoints( std::vector<int> & cycle, std::set<int> & free, bool inner );

  Element * element_;
  Side * side_;
  Graph graph_;
  GRAPH::CycleList facet_cycles_;
};

}
}


#endif
