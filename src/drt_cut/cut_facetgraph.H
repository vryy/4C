/*---------------------------------------------------------------------*/
/*!
\file cut_facetgraph.H

\brief graph to create volumecells from facets and lines

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_FACETGRAPH_H
#define CUT_FACETGRAPH_H

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

  /*!
  \brief Interface to the algorithm that creates volume cells.
  The connection between lines and facets forms a (particular) graph. Volume
  cells are minimal but complete loops in this graph.

  The current algorithm is purely topological and does not consider
  geometrical information. Therefore, cuts that cannot be decided on pure
  topological ground will not work properly. This is one of the many
  shortcomings of this library.
 */
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
