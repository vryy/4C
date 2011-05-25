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
  class Graph
  {
  public:

    void AddEdge( int row, int col );

    void AddEdge( Point * p1, Point * p2 );

    Point * GetPoint( int i ) { return all_points_[i]; }

    void Print( std::ostream & stream=std::cout );

    void PlotAllPoints( std::ostream & stream=std::cout );

    void PlotPoints( Element * element );

    void FindCycles( const std::vector<int> & cycle, std::vector<std::vector<Point*> > & cycles, bool inner );

    void FixSinglePoints();

    bool HasSinglePoints();

  private:

    std::map<int, std::set<int> > graph_;
    std::map<int, Point*> all_points_;
  };

public:

  typedef std::vector<std::vector<Point*> >::iterator iterator;

  PointGraph( Mesh & mesh, Element * element, Side * side, bool inner );

  iterator begin() { return cycles_.begin(); }

  iterator end()   { return cycles_.end(); }

  unsigned size() { return cycles_.size(); }

  void Print()
  {
    graph_.Print();
  }

private:

  void FillGraph( Element * element, Side * side, std::vector<int> & cycle );

  void AddFacetPoints( std::vector<int> & cycle, std::set<int> & free, bool inner );

  Element * element_;
  Side * side_;
  Graph graph_;
  std::vector<std::vector<Point*> > cycles_;
};

}
}


#endif
