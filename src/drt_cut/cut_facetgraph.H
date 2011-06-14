#ifndef CUT_FACETGRAPH_H
#define CUT_FACETGRAPH_H

#include <algorithm>

#include "cut_coloredgraph.H"
#include "cut_utils.H"

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

  FacetGraph( const std::vector<Side*> & sides, const plain_facet_set & facets );

  void Print()
  {
    graph_.Print();
    cycle_list_.Print();
  }

  void CreateVolumeCells( Mesh & mesh, Element * element, plain_volumecell_set & cells );

#ifdef DEBUGCUTLIBRARY
  bool InCollectedFacets( const plain_facet_set & collected_facets );

  void PrintAllCollected();

#endif

private:

//   void Add( GRAPH::Graph & graph,
//             std::map<std::pair<Point*, Point*>, int> & all_facets,
//             std::map<std::pair<Point*, Point*>, plain_facet_set> & facets
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
  std::vector<plain_facet_set> all_collected_facets_;
#endif
};

}
}

#endif
