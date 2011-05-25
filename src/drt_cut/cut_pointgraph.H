#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

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
public:

  enum Location
  {
    element_side,
    cut_side
  };

  enum Strategy
  {
    all_lines,
    own_lines
  };

private:

  struct Graph
  {
    void AddEdge( int row, int col );

    void AddEdge( Point * p1, Point * p2 );

    Point * GetPoint( int i ) { return all_points_[i]; }

    void Print( std::ostream & stream=std::cout );

    void PlotAllPoints( std::ostream & stream=std::cout );

    void PlotPoints( Element * element );

    void FindCycles( Side * side, std::vector<Point*> & cycle, Location location, Strategy strategy );

    void FixSinglePoints();

    bool HasSinglePoints();

    void GnuplotDumpCycles( const std::string & filename, const std::vector<std::vector<Point*> > & cycles );

    std::map<int, std::set<int> > graph_;
    std::map<int, Point*> all_points_;
    std::vector<std::vector<Point*> > main_cycles_;
    std::vector<std::vector<std::vector<Point*> > > hole_cycles_;
  };

public:

  typedef std::vector<std::vector<Point*> >::iterator facet_iterator;
  typedef std::vector<std::vector<std::vector<Point*> > >::iterator hole_iterator;

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

  void FillGraph( Element * element, Side * side, std::vector<Point*> & cycle, Strategy strategy );

  Element * element_;
  Side * side_;
  Graph graph_;
};

}
}


#endif
