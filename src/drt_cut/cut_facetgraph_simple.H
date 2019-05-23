/*----------------------------------------------------------------------------*/
/*!

\brief Create a simple facet graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )

\date Nov 14, 2016

\level 2

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#ifndef SRC_DRT_CUT_CUT_FACETGRAPH_SIMPLE_H_
#define SRC_DRT_CUT_CUT_FACETGRAPH_SIMPLE_H_

#include "cut_facetgraph.H"

namespace GEO
{
  namespace CUT
  {
    /*--------------------------------------------------------------------------*/
    class SimpleFacetGraph_1D : public FacetGraph
    {
     public:
      /// constructor
      SimpleFacetGraph_1D(const std::vector<Side*>& sides, const plain_facet_set& facets);

      /// destructor
      virtual ~SimpleFacetGraph_1D(){};

      virtual void CreateVolumeCells(Mesh& mesh, Element* element, plain_volumecell_set& cells);

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

      /// destructor
      virtual ~SimpleFacetGraph_2D(){};

      virtual void CreateVolumeCells(Mesh& mesh, Element* element, plain_volumecell_set& cells);

    };  // class  SimpleFacetGraph_2D
  }     // namespace CUT
}  // namespace GEO


#endif /* SRC_DRT_CUT_CUT_FACETGRAPH_SIMPLE_H_ */
