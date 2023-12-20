/*----------------------------------------------------------------------------*/
/*! \file

\brief Create a simple facet graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )

\date Nov 14, 2016

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef BACI_CUT_FACETGRAPH_SIMPLE_H
#define BACI_CUT_FACETGRAPH_SIMPLE_H

#include "baci_config.H"

#include "baci_cut_facetgraph.H"

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    /*--------------------------------------------------------------------------*/
    class SimpleFacetGraph_1D : public FacetGraph
    {
     public:
      /// constructor
      SimpleFacetGraph_1D(const std::vector<Side*>& sides, const plain_facet_set& facets);

      void CreateVolumeCells(Mesh& mesh, Element* element, plain_volumecell_set& cells) override;

     private:
      void SortFacets(const Element* element, std::map<double, Facet*>& sorted_facets) const;

      void CombineFacetsToLineVolumes(const std::map<double, Facet*>& sorted_facets,
          std::vector<plain_facet_set>& volumes) const;

    };  // class  SimpleFacetGraph_1D

    /*--------------------------------------------------------------------------*/
    class SimpleFacetGraph_2D : public FacetGraph
    {
     public:
      /// constructor
      SimpleFacetGraph_2D(const std::vector<Side*>& sides, const plain_facet_set& facets);

      void CreateVolumeCells(Mesh& mesh, Element* element, plain_volumecell_set& cells) override;

    };  // class  SimpleFacetGraph_2D
  }     // namespace CUT
}  // namespace CORE::GEO


BACI_NAMESPACE_CLOSE

#endif  // CUT_FACETGRAPH_SIMPLE_H
