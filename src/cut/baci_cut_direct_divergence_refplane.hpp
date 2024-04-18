/*----------------------------------------------------------------------*/
/*! \file

\brief Construct reference plane for direct divergence method when used in global
coordinate system

\level 2

 *------------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_CUT_DIRECT_DIVERGENCE_REFPLANE_HPP
#define FOUR_C_CUT_DIRECT_DIVERGENCE_REFPLANE_HPP

#include "baci_config.hpp"

#include "baci_cut_element.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Options;
    class Point;

    // When direct divergence is used in local coordinate system, a reference plane can be easily
    // chosen to ensure that all the integration points are completely within the background
    // element. However, ALE methods requires the integration to be constructed in global coordinate
    // system in which obtaining a proper reference plane is not trivial. This class handles all
    // possible choices of choosing a reference plane in such a way that all integration points are
    // ensured to be within the background element

    class DirectDivergenceGlobalRefplane
    {
     public:
      DirectDivergenceGlobalRefplane(Element* elem, VolumeCell* vc, Options& options)
          : elem1_(elem), volcell_(vc), options_(options)
      {
      }

      /*!
      \Compute the reference plane for this element
       */
      std::vector<double> GetReferencePlane();

      /*!
      \brief Get the reference points that are used to define the reference plane. This is used in
      gmsh output of volume cells
       */
      std::vector<Point*> GetReferencePointGmsh() { return ref_pts_gmsh_; }

     private:
      /*!
      \brief Compute reference plane based on the diagonals of the element
       */
      bool DiagonalBasedRef(
          std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol);

      /*!
      \brief Compute reference plane based on the facets of the volumecell
       */
      bool FacetBasedRef(std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol);
      /*!
      \brief Compute reference plane based on the sides of the element
       */
      bool SideBasedRef(std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol);

      /*!
      \brief Returns true if all the projected points are within the element
       */
      bool isAllProjectedCornersInsideEle(
          std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol = 1e-8);

      /*!
      \brief Scale the given equation of plane
       */
      void scaleEquationOfPlane(std::vector<double>& RefPlaneEqn);

      /*!
      \brief Class used to sort maps in descending order
       */
      struct CompareClass
      {
        bool operator()(const double& left, const double& right) const { return left > right; }
      };

      //! background element that contains the volumecell
      Element* elem1_;

      //! volumecell over which we construct integration scheme
      VolumeCell* volcell_;

      //! options container
      Options& options_;

      //! Points that define the reference plane used for this volumecell
      std::vector<Point*> ref_pts_gmsh_;

      //! Triangular diagonals in an hex8 element
      static const unsigned tri_diags_[24][3];

      //! Split of a side
      static const unsigned side_split_[4][3];
    };

  } /* namespace CUT */
}  // namespace CORE::GEO
FOUR_C_NAMESPACE_CLOSE

#endif
