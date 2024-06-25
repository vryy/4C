/*---------------------------------------------------------------------*/
/*! \file

\brief graph to create volumecells from facets and lines

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_FACETGRAPH_HPP
#define FOUR_C_CUT_FACETGRAPH_HPP

#include "4C_config.hpp"

#include "4C_cut_coloredgraph.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Point;
    class Facet;
    class Side;
    class Mesh;
    class Element;
    class VolumeCell;

    /*! \brief Interface to the algorithm that creates volume cells.
     *  The connection between lines and facets forms a (particular) graph. Volume
     *  cells are minimal but complete loops in this graph.
     *
     *  The current algorithm is purely topological and does not consider
     *  geometrical information. Therefore, cuts that cannot be decided on pure
     *  topological ground will not work properly. This is one of the many
     *  shortcomings of this library. */
    class FacetGraph
    {
     public:
      /// create the FacetGraph object for the given element dimension
      static Teuchos::RCP<Core::Geo::Cut::FacetGraph> Create(
          const std::vector<Side*>& sides, const plain_facet_set& facets);

     public:
      FacetGraph(const std::vector<Side*>& sides, const plain_facet_set& facets);

      /// destructor
      virtual ~FacetGraph() = default;

      void print() const
      {
        std::cout << "\n=== FacetGraph ===\n";
        std::cout << "--- Graph ---\n";
        graph_.print();
        std::cout << "\n--- CycleList ---\n";
        cycle_list_.print();
      }

      virtual void CreateVolumeCells(Mesh& mesh, Element* element, plain_volumecell_set& cells);

     protected:
      /// empty constructor for derived classes only
      FacetGraph(const unsigned& graph_size = 0)
          : graph_(graph_size){/* intentionally left blank */};

      void collect_volume_lines(plain_facet_set& collected_facets,
          std::map<std::pair<Point*, Point*>, plain_facet_set>& volume_lines) const;

      void add_to_volume_cells(Mesh& mesh, Element* element, std::vector<plain_facet_set>& volumes,
          plain_volumecell_set& cells) const;

      int facet_id(Facet* f)
      {
        return std::lower_bound(all_facets_.begin(), all_facets_.end(), f) - all_facets_.begin();
      }

      std::vector<Facet*> all_facets_;
      ColoredGraph::Graph graph_;
      ColoredGraph::CycleList cycle_list_;
      std::vector<std::pair<Point*, Point*>> all_lines_;
    };  // class FacetGraph

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
