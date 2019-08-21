/*---------------------------------------------------------------------*/
/*! \file

\brief Generate main Gauss points when using "DirectDivergence" approach.

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/
#ifndef DIRECT_DIVERGENCE_H
#define DIRECT_DIVERGENCE_H

#include "cut_element.H"
#include "cut_volumecell.H"

// Choose whether to output divirgence cell information in global coordinates (for easier comparison
// to Tesselation).
#define OUTPUT_GLOBAL_DIVERGENCE_CELLS

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;
    class VolumeCell;
    class Mesh;

    /*!
    \brief A class to construct Gaussian rule for volumecell by direct application of divergence
    theorem. This generate only the integration points on the facets.
    */
    class DirectDivergence
    {
     public:
      DirectDivergence(
          VolumeCell* volcell, Element* elem, const GEO::CUT::Point::PointPosition posi, Mesh& mesh)
          : volcell_(volcell), elem1_(elem), position_(posi), mesh_(mesh)
      {
      }

      ~DirectDivergence(){};

      /*!
      \brief Generate integration points on the facets of the volumecell
      */
      Teuchos::RCP<DRT::UTILS::GaussPoints> VCIntegrationRule(std::vector<double>& RefPlaneEqn);

      /*!
      \brief Compute and set correspondingly the volume of the considered volumecell from the
      generated integration rule and compare it with full application of divergence theorem
       */
      void DebugVolume(const DRT::UTILS::GaussIntegration& gpv, bool& isNeg);

      /*!
      \brief Integration of specific predefined functions using the Gaussian rule generated from
      DirectDivergence
      */
      void IntegrateSpecificFuntions(const DRT::UTILS::GaussIntegration& gpv);

      /*!
      \brief Geometry of volumecell, reference facet, main and internal gauss points for gmsh
      output.
       */
      void DivengenceCellsGMSH(
          const DRT::UTILS::GaussIntegration& gpv, Teuchos::RCP<DRT::UTILS::GaussPoints>& gpmain);

     private:
      /*!
      \brief Identify the list of facets which need to be triangulated, and also get the reference
      facet that will be used in xfluid part
       */
      void ListFacets(std::vector<plain_facet_set::const_iterator>& facetIterator,
          std::vector<double>& RefPlaneEqn, plain_facet_set::const_iterator& IteratorRefFacet,
          bool& IsRefFacet);



      //! volumecell over which we construct integration scheme
      VolumeCell* volcell_;

      //! background element that contains this volumecell
      Element* elem1_;

      //! position of this volumecell
      const GEO::CUT::Point::PointPosition position_;

      //! mesh that contains the background element
      Mesh& mesh_;

      //! reference facet identified for this volumecell
      Facet* refFacet_;

      //! true if the reference plane is on a facet of volumecell
      bool isRef_;

      //! Points that define the reference plane used for this volumecell
      std::vector<Point*> refPtsGmsh_;
    };
  }  // namespace CUT
}  // namespace GEO

#endif
