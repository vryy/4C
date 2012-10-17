#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

#include "cut_point.H"
#include "cut_cycle.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Line;
  class Side;
  class Element;
  class Mesh;

/// implementation details
namespace IMPL
{

/// a planar graph that is used to create facets out of lines and points
class PointGraph
{
public:

  enum Location
  {
    element_side,
    cut_side
  };

  enum Strategy
  {
    all_lines,
    own_lines     // used in levelset cut to create internal facets
  };

private:

  struct Graph
  {
    void AddEdge( int row, int col );

    void AddEdge( Point * p1, Point * p2 );

    Point * GetPoint( int i );

    void Print( std::ostream & stream=std::cout );

    void PlotAllPoints( std::ostream & stream=std::cout );

    void PlotPoints( Element * element );

    void FindCycles( Element * element, Side * side, Cycle & cycle, Location location, Strategy strategy );

    /*!
    \brief Any edge with single point in the graph is deleted
     */
    void FixSinglePoints( Cycle & cycle );

    bool HasSinglePoints();

    void GnuplotDumpCycles( const std::string & filename, const std::vector<Cycle> & cycles );

    std::map<int, plain_int_set> graph_;
    std::map<int, Point*> all_points_;
    std::vector<Cycle> main_cycles_;
    std::vector<std::vector<Cycle> > hole_cycles_;
  };

public:

  typedef std::vector<Cycle>::iterator facet_iterator;
  typedef std::vector<std::vector<Cycle> >::iterator hole_iterator;

  PointGraph( Mesh & mesh, Element * element, Side * side, Location location, Strategy strategy );

  facet_iterator fbegin() { return graph_.main_cycles_.begin(); }

  facet_iterator fend()   { return graph_.main_cycles_.end(); }

  hole_iterator hbegin() { return graph_.hole_cycles_.begin(); }

  hole_iterator hend()   { return graph_.hole_cycles_.end(); }

  void Print()
  {
    graph_.Print();
  }

private:

  /*!
  \brief Graph is filled with all edges that are created due to additional cut points and cut lines
   */
  void FillGraph( Element * element, Side * side, Cycle & cycle, Strategy strategy );

  Element * element_;
  Side * side_;
  Graph graph_;
};

}
}
}


#endif
