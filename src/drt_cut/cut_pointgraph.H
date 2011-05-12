#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

#include "cut_graph.H"

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

//     void FindSideLines( Mesh & mesh, Element * element, Side * side );

  private:

    std::map<int, Point*> all_points_;
  };

public:

  typedef GRAPH::CycleListIterator iterator;

  PointGraph( Mesh & mesh, Element * element, Side * side, bool inner );

  iterator begin() { return facet_cycles_.begin(); }

  iterator end()   { return facet_cycles_.end(); }

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

  Side * side_;
  Graph graph_;
  GRAPH::CycleList facet_cycles_;
};

}
}


#endif
