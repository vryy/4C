/*----------------------------------------------------------------------------*/
/*! \file

\brief Create a simple facet graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )

\date Nov 14, 2016

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_FACETGRAPH_SIMPLE_HPP
#define FOUR_C_CUT_FACETGRAPH_SIMPLE_HPP

#include "4C_config.hpp"

#include "4C_cut_facetgraph.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    /*--------------------------------------------------------------------------*/
    class SimpleFacetGraph1D : public FacetGraph
    {
     public:
      /// constructor
      SimpleFacetGraph1D(const std::vector<Side*>& sides, const plain_facet_set& facets);

      void create_volume_cells(Mesh& mesh, Element* element, plain_volumecell_set& cells) override;

     private:
      void sort_facets(const Element* element, std::map<double, Facet*>& sorted_facets) const;

      void combine_facets_to_line_volumes(const std::map<double, Facet*>& sorted_facets,
          std::vector<plain_facet_set>& volumes) const;

    };  // class  SimpleFacetGraph_1D

    /*--------------------------------------------------------------------------*/
    class SimpleFacetGraph2D : public FacetGraph
    {
     public:
      /// constructor
      SimpleFacetGraph2D(const std::vector<Side*>& sides, const plain_facet_set& facets);

      void create_volume_cells(Mesh& mesh, Element* element, plain_volumecell_set& cells) override;

    };  // class  SimpleFacetGraph_2D
  }     // namespace Cut
}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
