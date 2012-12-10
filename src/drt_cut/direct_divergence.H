/*!-----------------------------------------------------------------------------------------------*
\file direct_divergence.H

\brief Generate main Gauss points when using "DirectDivergence" approach.
equations
 *------------------------------------------------------------------------------------------------*/
#ifndef DIRECT_DIVERGENCE_H
#define DIRECT_DIVERGENCE_H

#include "cut_element.H"
#include "cut_volumecell.H"

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;
    class VolumeCell;
    class Mesh;

    /*!
    \brief A class to construct Gaussian rule for volumecell by direct application of divergence theorem.
    This generate only the integration points on the facets.
    */
    class DirectDivergence
    {
      public:
        DirectDivergence(VolumeCell* volcell, Element* elem, const GEO::CUT::Point::PointPosition posi, Mesh &mesh)
        :volcell_(volcell),
         elem1_(elem),
         position_(posi),
         mesh_(mesh)
        {
        }

        ~DirectDivergence(){};

        /*!
        \brief Generate integration points on the facets of the volumecell
        */
        Teuchos::RCP<DRT::UTILS::GaussPoints> VCIntegrationRule(std::vector<double> &RefPlaneEqn);

        /*!
        \brief Compute and set correspondingly the volume of the considered volumecell from the generated integration
        rule and compare it with full application of divergence theorem
         */
        void DebugVolume( const DRT::UTILS::GaussIntegration & gpv,
                          const std::vector<double> &RefPlaneEqn,
                          const std::vector<DRT::UTILS::GaussIntegration> intGRule );

      private:

        /*!
        \brief Identify the list of facets which need to be triangulated, and also get the reference facet that will be used in xfluid part
         */
        void ListFacets( std::vector<plain_facet_set::const_iterator>& facetIterator,
                         std::vector<double>& RefPlaneEqn,
                         plain_facet_set::const_iterator& IteratorRefFacet,
                         bool & IsRefFacet );

        /*!
        \brief Geometry of volumecell, reference facet and main gauss points for gmsh output.
         */
        void DivengenceCellsGMSH( plain_facet_set::const_iterator& IteratorRefFacet,
                                  std::vector<plain_facet_set::const_iterator>& facetIterator,
                                  const DRT::UTILS::GaussIntegration & gpv,
                                  bool & IsRefFacet );

        /*!
        \brief Integration of specific predefined functions using the Gaussian rule generated from DirectDivergence
        */
        void IntegrateSpecificFuntions( const DRT::UTILS::GaussIntegration & gpv,
                                        const std::vector<double> &RefPlaneEqn  );



        VolumeCell* volcell_;
        Element* elem1_;
        const GEO::CUT::Point::PointPosition position_;
        Mesh& mesh_;
    };
  }
}

#endif
