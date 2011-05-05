#ifndef CUT_FACETGRAPH_H
#define CUT_FACETGRAPH_H

#include <algorithm>

#include "cut_coloredgraph.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Facet;
  class Side;
  class Mesh;
  class Element;
  class VolumeCell;

class FacetGraph
{
public:

  FacetGraph( const std::vector<Side*> & sides, const std::set<Facet*> & facets );

  void Print()
  {
    graph_.Print();
    cycle_list_.Print();
  }

  void CreateVolumeCells( Mesh & mesh, Element * element, std::set<VolumeCell*> & cells );

#ifdef DEBUGCUTLIBRARY
  bool InCollectedFacets( const std::set<Facet*> & collected_facets );

  void PrintAllCollected();

#endif

private:

//   void Add( GRAPH::Graph & graph,
//             std::map<std::pair<Point*, Point*>, int> & all_facets,
//             std::map<std::pair<Point*, Point*>, std::set<Facet*> > & facets
//             );

  int FacetId( Facet * f )
  {
    return std::lower_bound( all_facets_.begin(), all_facets_.end(), f ) - all_facets_.begin();
  }

  COLOREDGRAPH::Graph graph_;
  COLOREDGRAPH::CycleList cycle_list_;
  std::vector<Facet*> all_facets_;
  std::vector<std::pair<Point*, Point*> > all_lines_;
#ifdef DEBUGCUTLIBRARY
  std::vector<std::set<Facet*> > all_collected_facets_;
#endif
};

}
}

#endif
