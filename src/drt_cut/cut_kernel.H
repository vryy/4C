/*---------------------------------------------------------------------*/
/*!
\file cut_kernel.H

\brief The cut kernel computes basic geometric operation, implemented are
    - Intersection of Surface and line or line and line
    - Calculate local coordinates inside an element
    - Compute Distance from a point to an embedded geometrical object
      ( surface or line )

\level 1

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_KERNEL_H
#define CUT_KERNEL_H

#include "cut_tolerance.H"
#include "cut_utils.H"
#include "cut_enum.H"
#include "cut_clnwrapper.H"
#include "cut_memory_manager.H"

#include "../linalg/linalg_gauss.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

// activate debug intersection output in cut_kernel
// #define DEBUG_CUTKERNEL_OUTPUT


#ifdef DEBUG_CUTKERNEL_OUTPUT
#include "cut_output.H"
#endif

namespace GEO
{
  namespace CUT
  {
    /*! \brief Core geometrical operations like finding the cut points,
     *   distance calculations takes place here */
    namespace KERNEL
    {
      // functions to compare determinant ot zero, e.g when computing the ConditionNumber
      bool closeToZero(const double a);

      bool closeToZero(const GEO::CUT::ClnWrapper& a);

      // Class to collects statistics about runs on double and cln in the cut intersection
      class CutKernelStatistics
      {
       public:
        static CutKernelStatistics& getCutKernelStatistics()
        {
          static CutKernelStatistics intersection_counter_;
          return intersection_counter_;
        }
        void DoubleIntersectionCounter() { double_int_++; };
        void DoubleDistanceCounter() { double_dist_++; };
        void ClnIntersectionCounter() { cln_int_++; };
        void ClnDistanceCounter() { cln_dist_++; };
        ~CutKernelStatistics()
        {
          std::cout << "\n\n =====================INTERSECTION "
                       "STATISTICS====================================\n\n";
          std::cout << "During compute intersection " << double_int_ << "/"
                    << double_int_ + cln_int_ << " was done on double only " << std::endl;
          std::cout << "During compute distance     " << double_dist_ << "/"
                    << double_dist_ + cln_dist_ << " was done on double only " << std::endl;
        }

       private:
        CutKernelStatistics(){/* blank */};
        unsigned long long int double_int_;
        unsigned long long int double_dist_;
        unsigned long long int cln_int_;
        unsigned long long int cln_dist_;
      };

      /// Information about the location of the point on the surface
      class PointOnSurfaceLoc
      {
       public:
        PointOnSurfaceLoc() : within_side_(false), on_side_(false) {}

        PointOnSurfaceLoc(bool within_side, bool on_side)
            : within_side_(within_side), on_side_(on_side)
        {
        }

        bool WithinSide() { return within_side_; }

        bool OnSide() { return on_side_; }

       private:
        bool within_side_;
        bool on_side_;
      };

      void FindCornerPoints(const std::vector<Point*>& points, std::vector<Point*>& corner_points);

      unsigned FindNextCornerPoint(const std::vector<Point*>& points, LINALG::Matrix<3, 1>& x1,
          LINALG::Matrix<3, 1>& x2, LINALG::Matrix<3, 1>& x3, LINALG::Matrix<3, 1>& b1,
          LINALG::Matrix<3, 1>& b2, LINALG::Matrix<3, 1>& b3, unsigned i);

      /// check if the given corner points can belong to a point1 "element"
      bool IsValidPoint1(const std::vector<Point*>& corner_points);

      /// check if the given corner points can belong to a line2 element
      bool IsValidLine2(const std::vector<Point*>& corner_points);

      /// check if the given corner points can belong to a tri3 element
      inline bool IsValidTri3(const std::vector<Point*>& points) { return points.size() == 3; }

      /// check if the given corner points can belong to a quad4 element
      bool IsValidQuad4(const std::vector<Point*>& points);

      DRT::Element::DiscretizationType CalculateShape(
          const std::vector<Point*>& points, std::vector<Point*>& line_points);

      template <int points>
      double FindL2Scaling(LINALG::Matrix<3, points>& xyze)
      {
        double scale = 0;
        LINALG::Matrix<3, 1> d;
        for (unsigned int i = 0; i < points; ++i)
        {
          LINALG::Matrix<3, 1> x1(&xyze(0, i), true);
          LINALG::Matrix<3, 1> x2(&xyze(0, (i + 1) % points), true);
          d.Update(1, x2, -1, x1, 0);
          scale += d.Norm2();
        }
        scale /= points;
        return scale;
      }

      /// get all edges adjacent to given local coordinates
      template <DRT::Element::DiscretizationType sideType, class T, class floatType>
      void GetEdgesAt(const T& xsi, std::vector<int>& edges_id, const floatType& tol)
      {
        switch (sideType)
        {
          case DRT::Element::tri3:
          {
            if (std::abs(xsi(1)) <= std::abs(tol(1))) edges_id.push_back(0);
            if (std::abs(xsi(0) + xsi(1) - 1.0) <=
                std::sqrt(tol(0) * tol(0) * 0.25 + tol(1) * tol(1) * 0.25))
              edges_id.push_back(1);
            if (std::abs(xsi(0)) <= std::abs(tol(0))) edges_id.push_back(2);
            break;
          }
          case DRT::Element::quad4:
          {
            if (std::abs(xsi(1) + 1.0) <= std::abs(tol(1))) edges_id.push_back(0);
            if (std::abs(xsi(0) - 1.0) <= std::abs(tol(0))) edges_id.push_back(1);
            if (std::abs(xsi(1) - 1.0) <= std::abs(tol(1))) edges_id.push_back(2);
            if (std::abs(xsi(0) + 1.0) <= std::abs(tol(0))) edges_id.push_back(3);
            break;
          }
          case DRT::Element::line2:
          {
            if (std::abs(xsi(0) + 1.0) <= std::abs(tol(0))) edges_id.push_back(0);
            if (std::abs(xsi(0) - 1.0) <= std::abs(tol(0))) edges_id.push_back(1);
            break;
          }
          default:
          {
            throw "Unknown/unsupported side type! \n";
            break;
          }
        }  // switch (sidetype)
      }

      /// get all nodes adjacent to given local coordinates
      template <DRT::Element::DiscretizationType sideType, class T, class floatType>
      void GetNodesAt(const T& xsi, std::vector<int>& nodes_id, const floatType& tol)
      {
        switch (sideType)
        {
          case DRT::Element::tri3:
          {
            if ((std::abs(xsi(1)) <= std::abs(tol(1))) and (std::abs(xsi(0)) <= std::abs(tol(0))))
              nodes_id.push_back(0);
            if ((std::abs(xsi(1) <= std::abs(tol(1)))) and
                (std::abs(1.0 - xsi(0)) <= std::abs(tol(0))))
              nodes_id.push_back(1);
            if (std::abs(1.0 - xsi(1)) <= std::abs(tol(1))) nodes_id.push_back(2);
            break;
          }
          case DRT::Element::quad4:
          {
            if ((std::abs(xsi(1) + 1.0) <= std::abs(tol(1))) and
                (std::abs(xsi(0) + 1.0) <= std::abs(tol(0))))
              nodes_id.push_back(0);
            if ((std::abs(xsi(1) + 1.0) <= std::abs(tol(1))) and
                (std::abs(xsi(0) - 1.0) <= std::abs(tol(0))))
              nodes_id.push_back(1);
            if ((std::abs(xsi(1) - 1.0) <= std::abs(tol(1))) and
                (std::abs(xsi(0) - 1.0) <= std::abs(tol(0))))
              nodes_id.push_back(2);
            if ((std::abs(xsi(0) - 1.0) <= std::abs(tol(1))) and
                (std::abs(xsi(0) + 1.0) <= std::abs(tol(0))))
              nodes_id.push_back(3);
            break;
          }
          case DRT::Element::line2:
          {
            if (std::abs(xsi(0) + 1.0) <= std::abs(tol(0))) nodes_id.push_back(0);
            if (std::abs(xsi(0) - 1.0) <= std::abs(tol(0))) nodes_id.push_back(1);
            break;
          }
          default:
          {
            throw "Unknown/unsupported side type! \n";
            break;
          }
        }  // switch (sidetype)
      }

      /** \brief check if the given local point \c xsi of the side is on an edge
       *
       *  The function returns \TRUE if the given point lies on one of the edges
       *  of the side.
       *
       *  \param xsi (in) : local parameter space coordinates
       *  \param tol (in) : tolerance for the check (default: \c tol=REFERENCETOL) */
      template <DRT::Element::DiscretizationType elementType, class T>
      bool AtEdge(const T& xsi, const double& tol)
      {
        // sanity check
        if (xsi.Rows() < DRT::UTILS::DisTypeToDim<elementType>::dim)
          dserror(
              "The given local coordinate has the wrong dimension!\n"
              "xsi.Rows() < Dim <===> %d<%d",
              xsi.Rows(), DRT::UTILS::DisTypeToDim<elementType>::dim);

        switch (elementType)
        {
          case DRT::Element::line2:
          {
            return (std::abs(xsi(0) + 1) < tol or std::abs(xsi(0) - 1) < tol);
            break;
          }
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::quad9:
          {
            return (std::abs(xsi(0) + 1) < tol or std::abs(xsi(1) + 1) < tol or
                    std::abs(xsi(0) - 1) < tol or std::abs(xsi(1) - 1) < tol);
            break;
          }
          case DRT::Element::tri3:
          case DRT::Element::tri6:
          {
            return (std::abs(xsi(0) + 0) < tol or std::abs(xsi(1) + 0) < tol or
                    std::abs(xsi(1) + xsi(0) - 1) < tol);
            break;
          }
          default:
            dserror("unsupported element type: %i | %s", elementType,
                DRT::DistypeToString(elementType).c_str());
            break;
        }
      }
      template <DRT::Element::DiscretizationType elementType, class T>
      bool AtEdge(const T& xsi)
      {
        return AtEdge<elementType>(xsi, REFERENCETOL);
      }


      /// withinlimits with 3 manually specified (independent from each other) tolerances
      template <DRT::Element::DiscretizationType elementType, class T, class floatType,
          unsigned int dim>
      bool WithinLimitsSplittedQuad(const T& xsi, const LINALG::Matrix<dim, 1, floatType>& tol)
      {
        if (elementType == DRT::Element::tri3)
        {
          return xsi(0) >= 0.0 - std::abs(tol(0)) and xsi(1) >= 0.0 - std::abs(tol(1)) and
                 xsi(0) <= 1.0 + std::abs(tol(0)) and xsi(1) <= 1.0 + std::abs(tol(1)) and
                 xsi(1) <= (1.0 - xsi(0)) + tol(2);
        }
        else
          dserror("This function is only called for triangulated elements");
      }


      template <DRT::Element::DiscretizationType elementType, class T, class floatType,
          unsigned int dim>
      bool WithinLimits(const T& xsi, const LINALG::Matrix<dim, 1, floatType>& tol)
      {
        if (xsi.M() < DRT::UTILS::DisTypeToDim<elementType>::dim)
          dserror(
              "The given local coordinate has the wrong dimension!\n"
              "xsi.Rows() < eleDim <===> %d<%d",
              xsi.M(), DRT::UTILS::DisTypeToDim<elementType>::dim);

        switch (elementType)
        {
          // ---------------------------------------------------------------
          // 1-D elements
          // ---------------------------------------------------------------
          case DRT::Element::line2:
          {
            return (xsi(0) >= -1.0 - std::abs(tol(0)) and xsi(0) <= 1.0 + std::abs(tol(0)));
            break;
          }
          // ---------------------------------------------------------------
          // 2-D elements
          // ---------------------------------------------------------------
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::quad9:
          {
            return xsi(0) >= -1.0 - std::abs(tol(0)) and xsi(1) >= -1.0 - std::abs(tol(1)) and
                   xsi(0) <= 1.0 + std::abs(tol(0)) and xsi(1) <= 1.0 + std::abs(tol(1));
            break;
          }
          case DRT::Element::tri3:
          case DRT::Element::tri6:
          {
            return xsi(0) >= 0.0 - std::abs(tol(0)) and xsi(1) >= 0.0 - std::abs(tol(1)) and
                   xsi(0) <= 1.0 + std::abs(tol(0)) and xsi(1) <= 1.0 + std::abs(tol(1)) and
                   xsi(1) <= (1.0 - xsi(0)) + std::sqrt(std::abs(tol(0)) * std::abs(tol(0)) * 0.25 +
                                                        std::abs(tol(1)) * std::abs(tol(1)) * 0.25);
            break;
          }
          // NOTE: Not sure if all the tolerances are set up correctly, especially points close to
          // the "\" diagonal line
          // ---------------------------------------------------------------
          // 3-D elements
          // ---------------------------------------------------------------
          case DRT::Element::hex8:
          case DRT::Element::hex16:
          case DRT::Element::hex20:
          case DRT::Element::hex27:
          {
            return (xsi(0) >= -1.0 - std::abs(tol(0)) and xsi(1) >= -1.0 - std::abs(tol(1)) and
                    xsi(2) >= -1.0 - std::abs(tol(2)) and xsi(0) <= 1.0 + std::abs(tol(0)) and
                    xsi(1) <= 1.0 + std::abs(tol(1)) and xsi(2) <= 1.0 + std::abs(tol(2)));
            break;
          }
          case DRT::Element::tet4:
          case DRT::Element::tet10:
          {
            return (xsi(0) >= -std::abs(tol(0)) and xsi(1) >= -std::abs(tol(1)) and
                    xsi(2) >= -std::abs(tol(2)) and
                    xsi(0) + xsi(1) + xsi(2) <= 1.0 + std::abs(tol(0)));
            break;
          }
          case DRT::Element::pyramid5:
          {
            return (xsi(0) >= -1.0 - std::abs(tol(0)) and xsi(1) >= -1.0 - std::abs(tol(1)) and
                    xsi(2) >= -std::abs(tol(2)) and xsi(0) <= 1.0 + std::abs(tol(0)) and
                    xsi(1) <= 1.0 + std::abs(tol(1)) and
                    (std::abs(xsi(0)) > std::abs(xsi(1))
                            ? (xsi(2) + std::abs(xsi(0)) <= 1.0 + std::abs(tol(0)))
                            : (xsi(2) + std::abs(xsi(1)) <= 1.0 + std::abs(tol(0)))));
            break;
          }
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          {
            return (xsi(0) >= -std::abs(tol(0)) and xsi(1) >= -std::abs(tol(1)) and
                    xsi(2) >= -1.0 - std::abs(tol(0)) and xsi(2) <= 1.0 + std::abs(tol(0)) and
                    xsi(0) + xsi(1) <= 1.0 + std::abs(tol(0)));
            break;
          }
          default:
            dserror("unsupported element type: %i | %s", elementType,
                DRT::DistypeToString(elementType).c_str());
            break;
        }  // switch ( elementType )
        return false;
      }


      // template< DRT::Element::DiscretizationType elementType, class T>

      /** \brief check if the local coordinate is within the specified limits
       *
       *  Returns \TRUE if the given point lies inside the element.
       *
       *  \param xsi (in) : local parameter space coordinates
       *  \param tol (in) : tolerance for the check (default: \c tol=REFERENCETOL) */
      template <DRT::Element::DiscretizationType elementType, class T, class floatType>
      bool WithinLimits(const T& xsi, const floatType& tol)
      {
        if (xsi.M() < DRT::UTILS::DisTypeToDim<elementType>::dim)
          dserror(
              "The given local coordinate has the wrong dimension!\n"
              "xsi.Rows() < eleDim <===> %d<%d",
              xsi.M(), DRT::UTILS::DisTypeToDim<elementType>::dim);

        switch (elementType)
        {
          // ---------------------------------------------------------------
          // 1-D elements
          // ---------------------------------------------------------------
          case DRT::Element::line2:
          {
            return (xsi(0) >= -1.0 - tol and xsi(0) <= 1.0 + tol);
            break;
          }
          // ---------------------------------------------------------------
          // 2-D elements
          // ---------------------------------------------------------------
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::quad9:
          {
            return xsi(0) >= -1.0 - tol and xsi(1) >= -1.0 - tol and xsi(0) <= 1.0 + tol and
                   xsi(1) <= 1.0 + tol;
            break;
          }
          case DRT::Element::tri3:
          case DRT::Element::tri6:
          {
            return xsi(0) >= 0.0 - tol and xsi(1) >= 0.0 - tol and xsi(0) <= 1.0 + tol and
                   xsi(1) <= 1.0 + tol and xsi(1) <= (1.0 - xsi(0)) + tol;
            break;
          }
          // ---------------------------------------------------------------
          // 3-D elements
          // ---------------------------------------------------------------
          case DRT::Element::hex8:
          case DRT::Element::hex16:
          case DRT::Element::hex20:
          case DRT::Element::hex27:
          {
            return (xsi(0) >= -1.0 - tol and xsi(1) >= -1.0 - tol and xsi(2) >= -1.0 - tol and
                    xsi(0) <= 1.0 + tol and xsi(1) <= 1.0 + tol and xsi(2) <= 1.0 + tol);
            break;
          }
          case DRT::Element::tet4:
          case DRT::Element::tet10:
          {
            return (xsi(0) >= -tol and xsi(1) >= -tol and xsi(2) >= -tol and
                    xsi(0) + xsi(1) + xsi(2) <= 1.0 + tol);
            break;
          }
          case DRT::Element::pyramid5:
          {
            return (
                xsi(0) >= -1.0 - tol and xsi(1) >= -1.0 - tol and xsi(2) >= -tol and
                xsi(0) <= 1.0 + tol and xsi(1) <= 1.0 + tol and
                (std::abs(xsi(0)) > std::abs(xsi(1)) ? (xsi(2) + std::abs(xsi(0)) <= 1.0 + tol)
                                                     : (xsi(2) + std::abs(xsi(1)) <= 1.0 + tol)));
            break;
          }
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          {
            return (xsi(0) >= -tol and xsi(1) >= -tol and xsi(2) >= -1.0 - tol and
                    xsi(2) <= 1.0 + tol and xsi(0) + xsi(1) <= 1.0 + tol);
            break;
          }
          default:
            dserror("unsupported element type: %i | %s", elementType,
                DRT::DistypeToString(elementType).c_str());
            break;
        }  // switch ( elementType )
        return false;
      }


      template <DRT::Element::DiscretizationType elementType, class T>
      bool WithinLimits(const T& xsi)
      {
        return WithinLimits<elementType>(xsi, REFERENCETOL);
      }

      /** \brief check if the local augmented coordinates are within the specified limits
       *
       *  \author hiermeier */
#if __cplusplus < 201103L

      template <unsigned probDim, DRT::Element::DiscretizationType sideType>
      bool WithinLimitsEmbeddedManifold(
          const Epetra_SerialDenseMatrix& xsi_aug, double tol, bool allow_dist, double tol2)
      {
        const LINALG::Matrix<probDim, 1> xsi_aug_linalg(xsi_aug.A(), true);
        return WithinLimitsEmbeddedManifold<probDim, sideType>(
            xsi_aug_linalg, tol, allow_dist, tol2);
      }

      template <unsigned probDim, DRT::Element::DiscretizationType sideType>
      bool WithinLimitsEmbeddedManifold(
          const LINALG::Matrix<probDim, 1>& xsi_aug, double tol, bool allow_dist, double tol2)
      {
        const unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim;
#else
      template <unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim, class T>
      bool WithinLimitsEmbeddedManifold(const T& xsi_aug, double tol, bool allow_dist, double tol2)
      {
#endif
        if (probDim == 1) dserror("probDim is not allowed to be equal 1!");
        if (xsi_aug.M() != probDim)
          dserror(
              "Wrong dimension of the augmented xsi input variable. \n"
              "The dimension is supposed to be equal to the problem dimension. \n"
              "The last coordinates are holding the distances from the global \n"
              "point to the side element. ");

        bool check = true;
        // check the distances if desired
        if (not allow_dist)
          for (unsigned r = dimSide; r < probDim; ++r)
            check = (check and (xsi_aug(r, 0) >= -tol2 and xsi_aug(r, 0) <= tol2) ? true : false);

        /* Since the first dimSide vector entries are the important ones,
         * we do not have to perform any extraction. */
        return (check and WithinLimits<sideType>(xsi_aug, tol));
      }
      template <unsigned probDim, DRT::Element::DiscretizationType sideType, class T>
      bool WithinLimitsEmbeddedManifold(
          const T& xsi_aug, const bool& allow_dist, const double& tol2)
      {
        return KERNEL::WithinLimitsEmbeddedManifold<probDim, sideType>(
            xsi_aug, REFERENCETOL, allow_dist, tol2);
      }

      /*! \brief Check whether three points are on the same line
       *  by checking whether the cross product is zero
       *
       *  \author Sudhakar
       *  \date 04/12 */
      bool IsOnLine(Point*& pt1, Point*& pt2, Point*& pt3, bool DeleteInlinePts = false);

      /*! \brief Check whether the polygon defined by the set of points is convex */
      std::vector<int> CheckConvexity(const std::vector<Point*>& ptlist,
          GEO::CUT::FacetShape& geomType, bool InSplit = true, bool DeleteInlinePts = false);

      // std::vector<double> EqnPlanePolygon( const std::vector<Point*>& ptlist, bool
      // DeleteInlinePts = false );


      /*! \brief Compute the equation of plane of this polygon using Newell's method */
      std::vector<double> EqnPlaneOfPolygon(const std::vector<Point*>& ptlist);

      /*! \brief Compute the equation of plane of this polygon using Newell's method */
      std::vector<double> EqnPlaneOfPolygon(const std::vector<std::vector<double>>& vertices);

      /*! \brief Find the equation of plane that contains these non-collinear points
       *
       *  It must be noted while using this function to find equation of facets,
       *  none of these 3 points must be a reflex (concave) point */
      std::vector<double> EqnPlane(Point*& pt1, Point*& pt2, Point*& pt3);

      /*! \brief Check whether the point named \c check is inside the triangle
       *  formed by tri */
      bool PtInsideTriangle(std::vector<Point*> tri, Point* check, bool DeleteInlinePts = false);

      /*! \brief Check whether the point \c check is inside the Quad */
      bool PtInsideQuad(std::vector<Point*> quad, Point* check);

      /*! \brief Return \TRUE if the points of the polygon are ordered clockwise
       *
       *  Polygon in 3D space is first projected into 2D plane, and the plane
       *  of projection is returned in projType.
       *
       *  \author sudhakar
       *  \date 05/12 */
      bool IsClockwiseOrderedPolygon(std::vector<Point*> polyPoints, std::string& projType);

      /*! \brief Delete unnecessary in-line points.
       *
       *  If more than two points are on a line, all points except the end points
       *  are deleted. This is checked for all the lines for a facet. So once this is
       *  called the facet is free of more than 2 inline points
       *
       *  \author sudhakar
       *  \date 06/12 */
      void DeleteInlinePts(std::vector<Point*>& poly);

      /*! \brief Returns true if at least 3 points are collinear
       *
       * \author wirtz
       * \date 05/13 */
      bool HaveInlinePts(std::vector<Point*>& poly);

      /*! \brief Finds tree points of the polygon which are not collinear
       *
       *  \author wirtz
       *  \date 05/13 */
      std::vector<Point*> Get3NoncollinearPts(std::vector<Point*>& polyPoints);

      /*! \brief Find appropriate projection plane
       *
       *  In several cases, it is appropriate to project the surface in 3D space into appropriate
       *  coordinate plane. It is better to project over the plane which is max normal component
       *  because this will reduce the round-off error in further calculations.
       *
       *  \author sudhakar
       *  \date 06/12 */
      void FindProjectionPlane(std::string& projPlane, const std::vector<double>& eqn);

      /*! \brief Split a QUAD4 element into two TRI3 elements
       *
       *  Here it is important that the triangle is created in the same rotation
       *  as the QUAD4 is, to get normal in the same direction and therefore the
       *  same signed distance:
       *
       *  tri3_id = 0 ---> QUAD4 nodes = {0 1 2}
       *  tri3_id = 1 ---> QUAD4 nodes = {2 3 0}
       *
       *  \param xyze_quad4 (in) : nodal coordinates of a QUAD4 element
       *  \param tri3_id    (in) : id of the desired TRI3 element \f& id \in \{0,\; 1\}\f$
       *  \param xyze_tri3  (out): created/filled TRI3 nodal coordinate matrix
       *
       *  \author hiermeier \date 11/16 */
      template <class T1, class T2>
      void SplitQuad4IntoTri3(const T1& xyze_quad4, const unsigned& tri3_id, T2& xyze_tri3)
      {
        if (tri3_id > 1)
          dserror(
              "A QUAD4 is supposed to be split into 2 TRI3 elements, \n"
              "therefore you have the choice between the tri3_id's 0 and 1. \n"
              "Nevertheless, you tried to access the id %d.",
              tri3_id);

        unsigned n0 = 2 * tri3_id;
        for (unsigned r = 0; r < static_cast<unsigned>(xyze_quad4.M()); ++r)
          for (unsigned c = 0; c < 3; ++c) xyze_tri3(r, c) = xyze_quad4(r, ((c + n0) % 4));
      };

      /*! \brief Get area of triangle in 3D space
       *
       *  \author sudhakar
       *  \date 11/14 */
      double getAreaTri(const std::vector<Point*>& poly, LINALG::Matrix<3, 1>* normalvec = NULL);
      double getAreaTri(const double* a, const double* b, const double* c,
          LINALG::Matrix<3, 1>* normalvec = NULL);
      template <class T>
      double getAreaTri(const T& xyze, LINALG::Matrix<3, 1>* normalvec = NULL)
      {
        if (xyze.M() != 3) run_time_error("Currently unsupported element dimension!");
        if (xyze.N() != 3) run_time_error("Wrong node number!");

        return getAreaTri(&xyze(0, 0), &xyze(0, 1), &xyze(0, 2), normalvec);
      }

      /*! \brief Get area of convex Quad in 3D space
       *
       *  Quad is split into two triangles and area of tri are summed up
       *
       *  \param poly (in): vector containting four pointers to the quad points
       *
       *  \author sudhakar
       *  \date 11/14 */
      double getAreaConvexQuad(std::vector<Point*>& poly);

      enum NewtonStatus
      {
        converged = 1,
        unconverged = 0,
        failed = -1
      };

      /** \brief Build an adaptive combined Newton tolerance
       *
       *  First a scaled absolute tolerance is calculated, based on the
       *  inf-norm of the absolute coordinate values of the considered side
       *  element and point.
       *
       *  In a second step a relative tolerance is added, which accounts for the
       *  fact, that the point px_ is maybe far outside of the considered side.
       *  If the point is inside or near the side, the additional term will stay
       *  very low and will not distort the used tolerance too much.
       *
       *  \author hiermeier \date 02/17 */
      template <class T1, class T2, class T3>
      double AdaptiveCombinedNewtonTolerance(

          const T1& xyze, const T2& px, const T3& initial_rhs)
      {
        /* --- Build the absolute tolerance */
        double tol = xyze.NormInf();
        double linescale = px.NormInf();
        if (linescale > tol) tol = linescale;
        tol *= LINSOLVETOL;

        /* --- Add the relative tolerance */
        tol += LINSOLVETOL * initial_rhs.NormInf();

        return tol;
      }

      template <class T1, class T2, class T3>
      double AdaptiveCombinedNewtonTolerance(

          const T1& xyze, const T2& px, const T3& initial_rhs, double)
      {
        return AdaptiveCombinedNewtonTolerance(
            xyze, px, initial_rhs);  // forwarding to normal function
      }

#ifdef CLN_CALC

      /// computes adaptive precision for AdaptivePrecision strategies
      template <class T1, class T2, class T3>
      ClnWrapper AdaptiveCombinedNewtonTolerance(
          const T1& xyze, const T2& px, const T3& initial_rhs, ClnWrapper&)
      {
        /* --- Build the absolute tolerance */
        // cln does not have built it pow function and using double will cause floating point
        // overflow, so we using manual cln conversion
        std::stringstream string_buffer;
        int nsize = ClnWrapper::GetPrecision();
        string_buffer << nsize;
        // construct string that represent 1e^-(current_cln_precision - 2) explicitely
        // otherwise cln might thinkg it too close to zero, and convert to short float zero and mess
        // the whole computation up
        std::string clnumstr =
            "0." + std::string(nsize - 2, '0') + "1" + "e+1_" + string_buffer.str();
        cln::cl_F cln_tol = clnumstr.c_str();
        // in order to make sure that we maintain proper precision
        ClnWrapper cln_tol_cln =
            cln::cl_float(cln_tol, cln::float_format(ClnWrapper::GetPrecision()));
        ClnWrapper cln_linsolvetol = cln_tol_cln / cln::sqrt(ClnWrapper(3.0)) *
                                     ClnWrapper(10.0);  //  similarly as for double

        ClnWrapper tol = xyze.NormInf();
        ClnWrapper linescale = px.NormInf();

        if (linescale > tol) tol = linescale;

        tol *= cln_linsolvetol;
        /* --- Add the relative tolerance */
        tol += cln_linsolvetol * initial_rhs.NormInf();

        if (tol == 0.0)
        {
          dserror(
              "Newton tolerance is cut_kernel for CLN is equal to zero! This should not happen!");
        }

        return tol;
      }

#endif

      /// Convert Newton Status enumerator to string
      static inline std::string NewtonStatus2String(const enum NewtonStatus& status)
      {
        switch (status)
        {
          case converged:
            return "CONVERGED";
          case unconverged:
            return "UNCONVERGED";
          case failed:
            return "FAILED";
          default:
            dserror("Unknown Newton status!");
            exit(EXIT_FAILURE);
        }
        exit(EXIT_FAILURE);
      }

      /*--------------------------------------------------------------------------*/
      /// generic Newton algorithm
      template <class Strategy, unsigned dim = 3, unsigned maxiter = 50>
      class NewtonSolve : public Strategy
      {
       public:
        /// constructor
        NewtonSolve(LINALG::Matrix<dim, 1>& xsi, bool checklimits) : Strategy(xsi, checklimits) {}
#ifdef CLN_CALC
        /// required constructor for CLN
        NewtonSolve(LINALG::Matrix<dim, 1, ClnWrapper>& xsi, bool checklimits)
            : Strategy(xsi, checklimits)
        {
        }
#endif

        bool Solve()
        {
          this->SetupSolve();

          for (unsigned iter = 0; iter < maxiter; ++iter)
          {
            this->SetupStep(iter);

            switch (this->TestConverged(iter))
            {
              // Newton iteration was successful
              case converged:
              {
                return true;
              }
              // Newton did not yet converge
              case unconverged:
              {
                break;
              }
              // Newton failed, thus we can stop here
              case failed:
              {
                return this->NewtonFailed();
              }
            }

            if (not this->LinearSolve(iter))
            {
              return false;
            }

            if (not this->Update(iter))
            {
              return false;
            }
          }

          return this->NewtonFailed();
        }

      };  // class class NewtonSolve

      /*--------------------------------------------------------------------------*/
      /// empty strategy for generic Newton algorithm
      class EmptyNewtonStrategy
      {
       public:
        void SetupSolve() {}

        void SetupStep(int iter) {}

        enum NewtonStatus TestConverged(int iter) { return unconverged; }

        bool LinearSolve(int iter) { return true; }

        bool Update(int iter) { return true; }

        bool NewtonFailed() { return false; }

        double GetTolerance()
        {
          dserror("Try to get tolerance from EmptyNewtonStrategy!");
          return 0.0;  // just to make compiler happy!
        }

        bool ZeroArea()
        {
          dserror("Try to get ZeroArea Information from EmptyNewtonStrategy!");
          return false;  // just to make compiler happy!
        }

        void WritetoGmsh(std::ofstream& file)
        {
          dserror("Try WritetoGmsh() from EmptyNewtonStrategy!");
        }
      };  // class EmptyNewtonStrategy

      /*--------------------------------------------------------------------------*/
      /// Debug helper strategy for generic Newton algorithm
      template <class Strategy, unsigned dim = 3>
      class DebugNewtonStrategy : public Strategy
      {
       public:
        /// constructor
        DebugNewtonStrategy(LINALG::Matrix<dim, 1>& xsi, bool checklimits)
            : Strategy(xsi, checklimits)
        {
        }

        void SetupSolve()
        {
          std::cout << "SetupSolve()\n";
          Strategy::SetupSolve();
        }

        void SetupStep(int iter)
        {
          std::cout << "SetupStep( iter = " << std::setw(2) << iter << " )\n";
          Strategy::SetupStep(iter);
        }

        enum NewtonStatus TestConverged(int iter)
        {
          enum NewtonStatus status = Strategy::TestConverged(iter);
          std::cout << "TestConverged( iter = " << std::setw(2) << iter
                    << " ) = " << NewtonStatus2String(status) << "\n"
                    << std::flush;
          return status;
        }

        bool LinearSolve(int iter)
        {
          bool res = Strategy::LinearSolve(iter);
          std::cout << "LinearSolve( iter = " << std::setw(2) << iter
                    << " )   = " << (res ? "SUCCESS" : "FAILED") << "\n"
                    << std::flush;
          return res;
        }

        bool Update(int iter)
        {
          bool res = Strategy::Update(iter);
          std::cout << "Update( iter = " << std::setw(2) << iter
                    << " )        = " << (res ? "SUCCESS" : "FAILED") << "\n"
                    << std::flush;
          return res;
        }

        bool NewtonFailed()
        {
          bool res = Strategy::NewtonFailed();
          std::cout << "NewtonFailed() = " << (res ? "SUCCESS" : "FAILED") << "\n" << std::flush;
          return res;
        }
      };  // class DebugNewtonStrategy

      // Static data storage of ComputePosition data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim, typename floatType = double>
      struct ComputePositionStaticMembers
      {
        /// nodal shape function values at the position xsi_
        static LINALG::Matrix<numNodesElement, 1, floatType> funct_;
        /** nodal first derivative shape function values at the position xsi_ */
        static LINALG::Matrix<probDim, numNodesElement, floatType> deriv1_;

        /** \brief (extended) jacobian matrix
         *
         *  Keep in mind, that the jacobian has to be extended if dim < probDim! */
        static LINALG::Matrix<probDim, probDim, floatType> A_;
        /// right hand side vector b_ = -(px_ - x_(xsi_))
        static LINALG::Matrix<probDim, 1, floatType> b_;
        /// newton increment (extended) parameter space coordinates
        static LINALG::Matrix<probDim, 1, floatType> dx_;
      };

      // Data storage of ComputePosition data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim, typename floatType = double>
      struct ComputePositionNoStaticMembers
      {
        LINALG::Matrix<numNodesElement, 1, floatType> funct_;

        LINALG::Matrix<probDim, numNodesElement, floatType> deriv1_;

        LINALG::Matrix<probDim, probDim, floatType> A_;

        LINALG::Matrix<probDim, 1, floatType> b_;

        LINALG::Matrix<probDim, 1, floatType> dx_;
      };



      /*--------------------------------------------------------------------------*/
      /** \brief strategy for position of point within element
       *
       *  inheritance diagram:
       *
       *  ComputePosition --> GenericComputePosition --> NewtonSolve
       *  --> ComputePositionStrategy --> EmptyNewtonStrategy
       *      ^^^^^^^^^^^^^^^^^^^^^^^
       */
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim, typename floatType = double,
          template <bool, unsigned, DRT::Element::DiscretizationType, unsigned, unsigned, typename>
          class MemberStoragePolicy = ComputePositionStaticMembers>
      class ComputePositionStrategy : public EmptyNewtonStrategy,
                                      private MemberStoragePolicy<debug, probDim, elementType,
                                          numNodesElement, dim, floatType>
      {
       public:
        /// constructor
        ComputePositionStrategy(LINALG::Matrix<dim, 1, floatType>& xsi, bool checklimits)
            : xsi_(xsi), xyze_(NULL), px_(NULL), tol_(0.0)
        {
        }

        /** setup the computation
         *
         *  \param xyze (in) : global nodal coordinates
         *  \param px   (in) : global position vector of the searched parameter space point */
        void Setup(const LINALG::Matrix<probDim, numNodesElement, floatType>& xyze,
            const LINALG::Matrix<probDim, 1, floatType>& px)
        {
          // std::cout << "Entering setup for the position strategy" << std::endl;
          xyze_ = &xyze;
          px_ = &px;

          /* use tol_ first just to store the scaling for the tolerance ..
           * set it to the maximal absolute coordinate value of the involved element
           * or point!
           *
           * --> this defines accuracy of the intersection (for good conditioned
           *     linear systems) */

          tol_ = xyze_->NormInf();
          floatType linescale = px_->NormInf();

          if (linescale > tol_) tol_ = linescale;

          tol_ *= LINSOLVETOL;

          if (debug)
          {
            std::cout << "\n=== ComputePosition ===\n";
            std::cout << "Element    = " << DRT::DistypeToString(elementType) << "\n";
            std::cout << "xyze_      = " << *xyze_;
            std::cout << "px_        = " << *px_;
            std::cout << "\n";
          }
        }

        /// get the local solution coordinates
        const LINALG::Matrix<dim, 1, floatType>& LocalCoordinates() { return xsi_; }
        /// get the solution coordaintes in the CLN format
        /// initialize the solution variable to zero
        void SetupSolve()
        {
          xsi_ = 0.0;
          b_ = 0.0;
          PositionRHS(*xyze_, *px_, b_);
          tol_ = AdaptiveCombinedNewtonTolerance(*xyze_, *px_, b_, xsi_(0, 0));
        }

        /** \brief evaluate the current right hand side
         *
         *  \f[
         *      b_ = - (px_ - x(xsi_))
         *  \f] */
        void SetupStep(int iter)
        {
          if (iter > 0)
          {
            PositionRHS(*xyze_, *px_, b_);
          }
        }

        /// calculate the position rhs
        void PositionRHS(const LINALG::Matrix<probDim, numNodesElement, floatType>& xyze,
            const LINALG::Matrix<probDim, 1, floatType>& px,
            LINALG::Matrix<probDim, 1, floatType>& b)
        {
          DRT::UTILS::shape_function<elementType>(xsi_, funct_);
          b = *px_;
          b.Multiply(-1.0, xyze, funct_, 1.0);
        }

        /// check for convergence
        enum NewtonStatus TestConverged(int iter)
        {
          floatType residual = b_.Norm2();

          if (debug)
          {
            std::cout << "WithinLimits             = "
                      << (WithinLimits<elementType>(xsi_) ? "TRUE" : "FALSE") << "\n"
                      << "rsd_.Norm2()             = " << xsi_.Norm2() << "\n"
                      << "dx_.Norm2() / rsd.Norm2()= "
                      << (iter != 0 ? (dx_.Norm2() / xsi_.Norm2()) : -1.0) << "\n"
                      << "b_.Norm2() ( RHS )       = " << b_.Norm2() << "\n"
                      << "b_.Norm2() / rsd.Norm2() = "
                      << (iter != 0 ? (b_.Norm2() / xsi_.Norm2()) : -1.0) << "\n"
                      << "tol_                     = " << tol_ << "\n"
                      << "rsd_ ( aka xsi_ )        = " << xsi_
                      << "b_ (RHS)                 = " << b_ << "\n"
                      << std::flush;
          }

          return (std::abs(residual) < tol_ ? converged : unconverged);
        }

        /** \brief Build and solve the linear system
         *
         *  \param iter (in): Currently unused. Maximal number of linear iterations,
         *                    but we do a direct solve here.
         *
         * \f[
         * \underline{\underline{J}} \Delta \underline{\xi} = - (\underline{x}_{p} -
         * \underline{x(\underline{\xi})},\\ \intertext{where} \underline{\underline{J}} =
         * \begin{pmatrix}
         * \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial
         * x(\underline{\xi})}{\partial \eta} & \frac{\partial x(\underline{\xi})}{\partial \zeta}
         * \\
         * \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial
         * y(\underline{\xi})}{\partial \eta} & \frac{\partial y(\underline{\xi})}{\partial \zeta}
         * \\ \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial
         * z(\underline{\xi})}{\partial \eta} & \frac{\partial z(\underline{\xi})}{\partial \zeta}
         * \end{pmatrix}.
         * \f]
         *
         * For the special case that the element dimension is smaller than the
         * problem dimension (i.e. manifold), see the ComputeDistance classes. */
        bool LinearSolve(int iter)
        {
          floatType det = 0.0;
          // ToDo check if it is possible to switch completely to the manifold case,
          // without too much loss in performance.
          DRT::UTILS::shape_function_deriv1<elementType>(xsi_, deriv1_);
          A_.MultiplyNT(*xyze_, deriv1_);
          dx_ = 0.0;
          det = LINALG::gaussElimination<true, probDim, floatType>(A_, b_, dx_);

          if (debug)
          {
            std::cout << "det                        = " << det << "\n";
            std::cout << "dx_ ( solution increment ) = " << dx_;
            std::cout << "dx_.Norm2()                = " << dx_.Norm2() << "\n";
          }

          return std::abs(det) >= LINSOLVETOL;
        }

        /// Update the solution parameter space coordinates
        bool Update(int iter)
        {
          // update the solution variables (element dimension)
          for (unsigned i = 0; i < dim; ++i) xsi_(i, 0) += dx_(i, 0);

          return true;
        }

        bool NewtonFailed() { return false; }

        std::pair<bool, floatType> ConditionNumber()
        {
          LINALG::Matrix<probDim, probDim, floatType> A_inv;
          floatType det = A_.Determinant();
          if (KERNEL::closeToZero(det))
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            std::cout << " WARNING: Condition number is equal to infinity in the problem. Stopping "
                         "the increase in precision"
                      << std::endl;
#endif
            return std::make_pair(false, -1.0);
          }
          A_inv.Invert(A_);
          floatType cond = A_inv.Norm2() * A_.Norm2();
          return std::make_pair(true, cond);
        }
        /// get the Newton tolerance
        floatType GetTolerance() const { return tol_; }

        floatType GetLocalTolerance(const floatType& real_tolerance)
        {
          dserror(
              "Not impelemented! Have a look at  GetLocalTolerance for ComputeDistance strategy to "
              "find out the idea of implementation!");
          return real_tolerance;
        }

       private:
        /** parameter space coordinates (final result of the calculation)
         * corresponding to the global position px_. */
        LINALG::Matrix<dim, 1, floatType>& xsi_;

        /// pointer to the global nodal positions
        const LINALG::Matrix<probDim, numNodesElement, floatType>* xyze_;
        /// global position vector (we are looking for the corresponding parameter coordinates)
        const LINALG::Matrix<probDim, 1, floatType>* px_;

        using MemberStoragePolicy<debug, probDim, elementType, numNodesElement, dim,
            floatType>::funct_;

        using MemberStoragePolicy<debug, probDim, elementType, numNodesElement, dim,
            floatType>::deriv1_;

        using MemberStoragePolicy<debug, probDim, elementType, numNodesElement, dim, floatType>::A_;

        using MemberStoragePolicy<debug, probDim, elementType, numNodesElement, dim, floatType>::b_;

        using MemberStoragePolicy<debug, probDim, elementType, numNodesElement, dim,
            floatType>::dx_;

        /// adapted Newton tolerance
        floatType tol_;

        // precision for cln
      };  // class ComputePositionStrategy


#ifdef DEBUG_MEMORY_ALLOCATION
      // shared data between all the template instantiations, for debugging memory allocations
      class AdaptiveKernelSharedData
      {
       public:
        static bool custom_allocator_run_;

       protected:
        static size_t cln_byte_size_[CLN_LIMIT_ITER];

        static boost::unordered_map<size_t, int> memory_allocations_;

        static std::map<size_t, int> memory_allocations_intersection_;

        static std::map<size_t, int> memory_allocations_distance_;

        static std::map<size_t, int> memory_allocations_position_;

        static bool all_intersections_done_once_;

        static bool all_distance_done_once_;

        static bool all_position_done_once_;
        // specifies whether to run on the custor allocator now

        void UpdateMemoryAllocations(const boost::unordered_map<size_t, int>& allocations)
        {
          for (boost::unordered_map<size_t, int>::const_iterator it = allocations.begin();
               it != allocations.end(); ++it)
          {
            if (memory_allocations_[it->first] < (it->second))
            {
              std::cout << "Size " << (*it).first << " was allocated " << (*it).second << " times "
                        << std::endl;
              memory_allocations_[it->first] = it->second;
            }
          }
        }

        static int& GetPositionCount()
        {
          static int int_count = 0;
          return int_count;
        }
        static int& GetDistanceCount()
        {
          static int int_count = 0;
          return int_count;
        }
        static int& GetIntersectionCount()
        {
          static int int_count = 0;
          return int_count;
        }

        static int& PositionCounter()
        {
          static int counter = 0;
          return counter;
        }
        static int IncreasePositionCounter() { return (++PositionCounter()); }

        static int& DistanceCounter()
        {
          static int counter = 0;
          return counter;
        }
        static int IncreaseDistanceCounter() { return (++DistanceCounter()); }

        static int& IntersectionCounter()
        {
          static int counter = 0;
          return counter;
        }
        static int IncreaseIntersectionCounter() { return (++IntersectionCounter()); }

        static void ReportIntersectionAllocated()
        {
          size_t total = 0;
          std::cout << "REPORTING MAXIMUM MEMORY ALLOCATION IN INTERSECTION" << std::endl;
          for (std::map<size_t, int>::iterator it = memory_allocations_intersection_.begin();
               it != memory_allocations_intersection_.end(); ++it)
          {
            std::cout << "Size " << (*it).first << " was allocated " << (*it).second << " times "
                      << std::endl;
            if (memory_allocations_[it->first] < it->second)
            {
              memory_allocations_[it->first] = it->second;
            }
            total += (*it).second * (*it).first;
          }
          std::cout << "Totally maximum required allocationd " << total << " bytes" << std::endl;
        }

        static void ReportPositionAllocated()
        {
          size_t total = 0;
          std::cout << "REPORTING MAXIMUM MEMORY ALLOCATION IN THE POSITION" << std::endl;
          for (std::map<size_t, int>::iterator it = memory_allocations_position_.begin();
               it != memory_allocations_position_.end(); ++it)
          {
            std::cout << "Size " << (*it).first << " was allocated " << (*it).second << " times "
                      << std::endl;
            if (memory_allocations_[it->first] < it->second)
            {
              memory_allocations_[it->first] = it->second;
            }
            total += (*it).second * (*it).first;
          }
          std::cout << "Totally maximum required allocationd " << total << " bytes" << std::endl;
        }

        static void ReportDistanceAllocated()
        {
          size_t total = 0;
          std::cout << "REPORTING MAXIMUM MEMORY ALLOCATION IN THE DISTANCE" << std::endl;
          for (std::map<size_t, int>::iterator it = memory_allocations_distance_.begin();
               it != memory_allocations_distance_.end(); ++it)
          {
            std::cout << "Size " << (*it).first << " was allocated " << (*it).second << " times "
                      << std::endl;
            if (memory_allocations_[it->first] < it->second)
            {
              memory_allocations_[it->first] = it->second;
            }
            total += (*it).second * (*it).first;
          }
          std::cout << "Totally maximum required allocationd " << total << " bytes" << std::endl;
        }

        static void ReportTotalAllocated()
        {
          size_t total = 0;
          std::cout << "REPORTING MAXIMUM MEMORY ALLOCATION" << std::endl;
          for (boost::unordered_map<size_t, int>::iterator it = memory_allocations_.begin();
               it != memory_allocations_.end(); ++it)
          {
            std::cout << "Size " << (*it).first << " was allocated " << (*it).second << " times "
                      << std::endl;
            total += (*it).second * (*it).first;
          }
          std::cout << "Totally maximum required allocationd " << total << " bytes" << std::endl;
        }
      };

#endif

#ifdef CLN_CALC

      /*--------------------------------------------------------------------------*/
      /** \brief Strategy for position of point within element for arbitrary precision
       *
       *  inheritance diagram:
       *
       *  ComputePosition --> GenericComputePosition --> NewtonSolve
       *  --> ComputePositionAdaptivePrecision --> EmptyNewtonStrategy
       *      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType elementType,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
      class ComputePositionAdaptivePrecision : Strategy
#ifdef DEBUG_MEMORY_ALLOCATION
          ,
                                               AdaptiveKernelSharedData
#endif
      {
       public:
        ComputePositionAdaptivePrecision(LINALG::Matrix<dim, 1>& xsi)
            : Strategy(clnxsi_, true), xsi_(xsi), clnxsi_(true), cond_infinity_(false)
        {
          dserror(
              "ComputePositionAdaptivePrecision: This function is not tested! Most probably it "
              "will work, but be careful");
        }

        bool operator()(const LINALG::Matrix<probDim, numNodesElement>& xyze,
            const LINALG::Matrix<probDim, 1>& px)
        {
          if (!((elementType == DRT::Element::hex8) || (elementType == DRT::Element::tet4) ||
                  (elementType == DRT::Element::pyramid5) || (elementType == DRT::Element::quad4) ||
                  (elementType == DRT::Element::hex20) || (elementType == DRT::Element::tri3) ||
                  (elementType == DRT::Element::line2) || (elementType == DRT::Element::wedge6)))
          {
            dserror(
                "This type of element %s  is not yet implemented for the CLN calculation. You are "
                "welcome to edit fem_shapefunctions.H file to fix it",
                DRT::DistypeToString(elementType).c_str());
          }

          bool conv;
          int iter = 0;
          int prec = CLN_START_PRECISION;  // standart initial precision
          ClnWrapper err;
          ClnWrapper cond_number;
          cond_infinity_ = false;  // test check if condition number is equal to infinity

#if DEBUG_MEMORY_ALLOCATION
          // Report total memory allocation done, for all possible template instantiations
          if (all_distance_done_once_ and all_intersections_done_once_ and
              all_position_done_once_ and (!custom_allocator_run_))
          {
            GEO::CUT::MemorySingleton::getInstance().ReportAllocated();
            ReportIntersectionAllocated();
            ReportPositionAllocated();
            ReportDistanceAllocated();
            ReportTotalAllocated();
            GEO::CUT::MemorySingleton::getInstance().SetState(1, memory_allocations_);
            custom_allocator_run_ = true;
          }
          if (first_run_) std::cout << "In cut position statistics is" << std::endl;
#endif
          do
          {
            ClnWrapper::SetPrecision(prec);
#ifdef CUSTOM_MEMORY_ALLOCATOR
#if DEBUG_MEMORY_ALLOCATION
            // set current constainer for most frequenty byte size allocation for this iteration
            if (custom_allocator_run_)
#endif
              GEO::CUT::MemorySingleton::getInstance().GetMemoryPoolAllocator().SetCurrent(
                  cln_byte_size_[iter]);
#endif

#if DEBUG_MEMORY_ALLOCATION
            // Report memory allocation of the class if doing this for the first time
            if (first_run_)
            {
              if (GEO::CUT::MemorySingleton::getInstance().IsRecording())
              {
                UpdateMemoryAllocations(
                    GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern());
                GEO::CUT::MemorySingleton::getInstance().StopRecord();
              }
              GEO::CUT::MemorySingleton::getInstance().StartRecord();
            }
#endif
            // convertion to correspondent precison
            ConvDoubleCLN(xyze, clnxyze_, prec);
            ConvDoubleCLN(px, clnpx_, prec);

            this->Setup(clnxyze_, clnpx_);
            conv = this->Solve();

            // safety check of precision loss
            cln::float_format_t prec_beg = cln::float_format(clnxyze_(0, 0));
            for (unsigned int i = 0; i < dim; ++i)
            {
              cln::float_format_t prec_end = cln::float_format(clnxsi_(i, 0));
              if (prec_beg != prec_end)
              {
                dserror(
                    "There is a loss of cln-precision during computation of Intersection. "
                    "Something is wrong with the conversion");
              }
            }

            std::pair<bool, ClnWrapper> cond_pair = this->ConditionNumber();
            cond_number = cond_pair.second;
            if (not cond_pair.first)
            {
#ifdef DEBUG_CUTKERNEL_OUTPUT
              std::cout << "Condition number of this problem is equal to infinity" << std::endl;
#endif
              cond_infinity_ = true;
            }

            err =
                ComputeError(xyze, px, this->LocalCoordinates(), CLN_REFERENCE_PREC) * cond_number;
            // get location taking error into account
            //  NOTE: Here one might need to get tolerance in the local coordinates, similarly as
            //  done
            // in the ComputeDistance and ComputeIntersection
            //  location_ = GetPreciseLocation<elementType>( clnxsi_, err );

#if DEBUG_MEMORY_ALLOCATION
            UpdateMemoryUsage(prec, iter);
#endif

            prec += CLN_INCREMENT_STEP;
            iter++;

          } while (((!conv) || (err > CLN_LIMIT_ERROR) || (cond_infinity_)) &&
                   (iter < CLN_LIMIT_ITER) && (prec < CLN_LIMIT_PREC));


          // safely convert the values back
          ConvClnDouble(clnxsi_, xsi_);
          if ((not conv) && (err < CLN_LIMIT_ERROR) && (err > 0)) conv = true;

#ifdef CUSTOM_MEMORY_ALLOCATOR
          if (first_run_)
          {
            const int num_inter = 3;
            if ((++GetPositionCount()) == num_inter)
            {
              all_position_done_once_ = true;
              std::cout << "All possible template position are done!" << std::endl;
            }
            first_run_ = false;
          }
#endif
          return conv;
        }

        const LINALG::Matrix<dim, 1, ClnWrapper>& LocalCoordinates()
        {
          return Strategy::LocalCoordinates();
        }

        std::pair<bool, ClnWrapper> ConditionNumber() { return Strategy::ConditionNumber(); }

        ClnWrapper GetTolerance() const { return Strategy::GetTolerance(); }

        PointOnSurfaceLoc GetSideLocation() { return location_; }

       private:
        // Evaluate difference between inital global coordinates passed to ComputePosition and
        // and global coordinates based on loc coordinates calculated in the ComputePosition,
        // using conversion to very high reference precision
        ClnWrapper ComputeError(
            const LINALG::Matrix<probDim,
                DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement>& ref_shape_xyz,
            const LINALG::Matrix<probDim, 1>& glob_init,
            const LINALG::Matrix<dim, 1, ClnWrapper>& loc_calc, int prec)
        {
          unsigned int prev_prec = ClnWrapper::GetPrecision();
          ClnWrapper::SetPrecision(prec);
          // Converting input data into higher precision floating point format
          LINALG::Matrix<probDim,
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement, ClnWrapper>
              cln_ref_shape_xyz;
          LINALG::Matrix<probDim, 1, ClnWrapper> cln_glob_init;
          LINALG::Matrix<probDim, 1, ClnWrapper> cln_loc_calc;

          ConvDoubleCLN(loc_calc, cln_loc_calc, prec);
          ConvDoubleCLN(ref_shape_xyz, cln_ref_shape_xyz, prec);
          ConvDoubleCLN(glob_init, cln_glob_init, prec);

          LINALG::Matrix<probDim, 1, ClnWrapper> cln_glob_calc;
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement, 1,
              ClnWrapper>
              shapeFunct;
          DRT::UTILS::shape_function<elementType>(cln_loc_calc, shapeFunct);
          for (unsigned int inode = 0;
               inode < DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement; ++inode)
          {
            for (unsigned int isd = 0; isd < probDim; ++isd)
            {
              cln_glob_calc(isd) += cln_ref_shape_xyz(isd, inode) * shapeFunct(inode);
            }
          }
          // evaluating error by comparison with global position, that was set up initially and the
          // one that we obtained from the ComputePosition evaluation and the interpolation
          LINALG::Matrix<probDim, 1, ClnWrapper> diff_vec;
          for (unsigned int i = 0; i < probDim; ++i)
            diff_vec(i) = cln_glob_calc(i) - cln_glob_init(i);
          // resetting precision to previous
          ClnWrapper::SetPrecision(prev_prec);
          return diff_vec.Norm2();
        }


        // Update global maximum number of allocations ( number of containers of particular byte
        // size) if allocations on this precisions were more. Also set what was the most frequent
        // size of allocation for the current iteration
        void UpdateMemoryUsage(int prec, int iter)
        {
#if DEBUG_MEMORY_ALLOCATION
          if (first_run_)
          {
            GEO::CUT::MemorySingleton::getInstance().StopRecord();

            boost::unordered_map<size_t, int>& allocation_map =
                GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern();

            // update the global number of allocations to the maximum
            int max_num = 0;
            int size_max_num = 0;
            for (boost::unordered_map<size_t, int>::iterator it = allocation_map.begin();
                 it != allocation_map.end(); ++it)
            {
              size_t size = it->first;
              if (memory_allocations_position_[size] < allocation_map[size])
              {
                memory_allocations_position_[size] = allocation_map[size];
              }
              if (allocation_map[size] > max_num)
              {
                size_max_num = size;
                max_num = allocation_map[size];
              }
            }

            if (size_max_num != 0)
            {
              cln_sizes_[iter] = size_max_num;
              std::cout << "Size for the precision " << prec << " is " << size_max_num << std::endl;
              std::cout << "It is allocated " << max_num << " times " << std::endl;
            }
            else
              dserror("This should not be possible!");

            GEO::CUT::MemorySingleton::getInstance().ResetAllocated();
          }
#endif
        }

       private:
        // variable for cln computation
        LINALG::Matrix<probDim, numNodesElement, ClnWrapper> clnxyze_;

        LINALG::Matrix<probDim, 1, ClnWrapper> clnpx_;

        // reference to the double TMATRIX
        LINALG::Matrix<dim, 1>& xsi_;

        LINALG::Matrix<dim, 1, ClnWrapper> clnxsi_;

        // if condition number if infinity or not
        bool cond_infinity_;

        // determines location of this point with respect to the edge and side that cut it
        PointOnSurfaceLoc location_;

        static bool first_run_;

        // sizes of cln variables with different preciion
        static size_t cln_sizes_[CLN_LIMIT_ITER];
      };

#endif

      /*--------------------------------------------------------------------------*/
      /** \brief generic strategy for position of point within element
       *
       *  inheritance diagram:
       *
       *  ComputePosition --> GenericComputePosition --> NewtonSolve
       *                      ^^^^^^^^^^^^^^^^^^^^^^
       *  --> ComputePositionStrategy --> EmptyNewtonStrategy
       */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType elementType,
          bool computeCln = false,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
      class GenericComputePosition : Strategy
      {
       public:
        GenericComputePosition(LINALG::Matrix<dim, 1>& xsi)
            : Strategy(xsi, true)
#ifdef CLN_CALC
              ,
              xsi_(xsi)
#endif
        {
        }

        bool operator()(const LINALG::Matrix<probDim, numNodesElement>& xyze,
            const LINALG::Matrix<probDim, 1>& px)
        {
          bool conv;
          this->Setup(xyze, px);
          conv = this->Solve();
#ifdef CLN_CALC
          if (computeCln)  // this is not used up to now!
          {
            {
              dserror("Position Strategy with CLN not tested!");
              ComputePositionAdaptivePrecision<
                  NewtonSolve<ComputePositionStrategy<false, probDim, elementType, numNodesElement,
                                  dim, ClnWrapper, ComputePositionNoStaticMembers>,
                      dim>,
                  probDim, elementType>
                  clncalc(xsi_);

              bool clnsolver = clncalc(xyze, px);
              conv = clnsolver;
              location_ = clncalc.GetSideLocation();
              ClnWrapper::ResetPrecision();
            }
#ifdef CUSTOM_MEMORY_ALLOCATOR
            GEO::CUT::MemorySingleton::getInstance().Finalize();
#endif
          }
#endif
          return conv;
        }

        const LINALG::Matrix<dim, 1>& LocalCoordinates() { return Strategy::LocalCoordinates(); }

        std::pair<bool, double> ConditionNumber() { return Strategy::ConditionNumber(); }

        double GetTolerance() const { return Strategy::GetTolerance(); }
#ifdef CLN_CALC

        PointOnSurfaceLoc GetSideLocation() { return location_; }

       private:
        PointOnSurfaceLoc location_;

        // touched edges ids
        static std::vector<int> touched_edges_ids_;

        // to hold the reference to cln calculation object
        LINALG::Matrix<dim, 1>& xsi_;
#endif
      };  // class GenericComputePosition

      /*--------------------------------------------------------------------------*/
      /** \brief most derived strategy class for position of point within element
       *
       *  inheritance diagram:
       *
       *  ComputePosition --> GenericComputePosition --> NewtonSolve
       *  ^^^^^^^^^^^^^^^
       *  --> ComputePositionStrategy --> EmptyNewtonStrategy
       */
      template <unsigned probDim, DRT::Element::DiscretizationType elementType,
          bool computeCln = false,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
      class ComputePosition
          : public GenericComputePosition<
                NewtonSolve<ComputePositionStrategy<false, probDim, elementType>, dim>, probDim,
                elementType, computeCln>
      {
       public:
        ComputePosition(LINALG::Matrix<dim, 1>& xsi)
            : GenericComputePosition<
                  NewtonSolve<ComputePositionStrategy<false, probDim, elementType, numNodesElement>,
                      dim>,
                  probDim, elementType, computeCln>(xsi)
        {
          if (dim != probDim)
            dserror(
                "You called the wrong object! Use the ComputeDistance class "
                "for the embedded case.");
        }
      };  // class ComputePosition

      /*--------------------------------------------------------------------------*/
      template <unsigned probDim, DRT::Element::DiscretizationType elementType,
          unsigned numNodesElement =
              DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
          unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
      class DebugComputePosition
          : public GenericComputePosition<
                NewtonSolve<
                    DebugNewtonStrategy<
                        ComputePositionStrategy<true, probDim, elementType, numNodesElement>, dim>,
                    dim>,
                probDim, elementType>
      {
       public:
        DebugComputePosition(LINALG::Matrix<dim, 1>& xsi)
            : GenericComputePosition<
                  NewtonSolve<DebugNewtonStrategy<ComputePositionStrategy<true, probDim,
                                                      elementType, numNodesElement>,
                                  dim>,
                      dim>,
                  probDim, elementType>(xsi)
        {
        }
      };  // class DebugComputePosition

      // Static storage class of ComputeDistance data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double>
      struct ComputeDistanceStaticMembers
      {
        /// nodal shape function values at \c xsi
        static LINALG::Matrix<numNodesSide, 1, floatType> sideFunct_;

        /** nodal 1-st derivative values at \c xsi
         *
         *  \remark This is an augmented derivatives matrix. The actual
         *  dimension is \f$ dim \times numNodesSide \f$. The remaining
         *  entries are filled with zeros. This is due to the use of an
         *  UTILS function, which expects this input. */
        static LINALG::Matrix<probDim, numNodesSide, floatType> sideDeriv1_;

        /** \brief nodal 2-nd derivatives at \c xsi
         *
         *  1-D case: \f$ x_{,\xi \xi}
         *                \;\in \mathbb{R}^{1 \times N}\f$,
         *
         *  2-D case: \f$ \left( x_{,\xi \xi}, \;x_{,\eta \eta}, \;
         *                       x_{,\xi\eta} \right)
         *                \;\in \mathbb{R}^{3 \times N} \f$ */
        static LINALG::Matrix<2 * dimSide - 1, numNodesSide, floatType> sideDeriv2_;

        /// complete linearization matrix
        static LINALG::Matrix<probDim, probDim, floatType> A_;
        /// auxiliary matrix
        static LINALG::Matrix<probDim, probDim, floatType> B_;
        /** auxiliary matrix holding 2-nd derivatives at \c xsi_
         *  w.r.t. the parameter space coordinates */
        static LINALG::Matrix<probDim, 2 * dimSide - 1, floatType> C_;
        /// right hand side vector
        static LINALG::Matrix<probDim, 1, floatType> b_;
        /** \brief solution increment
         *
         * (parameter space coordinates of the side + distance increment) */
        static LINALG::Matrix<probDim, 1, floatType> dx_;
        /// (unscaled) normal vectors
        static LINALG::Matrix<probDim, 2, floatType> N_;

        static LINALG::Matrix<probDim, 2, floatType> nvec_;  // proper normal vector
      };

      // Static storage class of ComputeDistance data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double>
      struct ComputeDistanceNoStaticMembers
      {
        LINALG::Matrix<numNodesSide, 1, floatType> sideFunct_;

        LINALG::Matrix<probDim, numNodesSide, floatType> sideDeriv1_;

        LINALG::Matrix<2 * dimSide - 1, numNodesSide, floatType> sideDeriv2_;

        LINALG::Matrix<probDim, probDim, floatType> A_;

        LINALG::Matrix<probDim, probDim, floatType> B_;

        LINALG::Matrix<probDim, 2 * dimSide - 1, floatType> C_;

        LINALG::Matrix<probDim, 1, floatType> b_;

        LINALG::Matrix<probDim, 1, floatType> dx_;

        LINALG::Matrix<probDim, 2, floatType> N_;

        LINALG::Matrix<probDim, 2, floatType> nvec_;
      };

      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double,
          template <bool, unsigned, DRT::Element::DiscretizationType, unsigned, unsigned, typename>
          class MemberStoragePolicy = ComputeDistanceStaticMembers>
      class ComputeDistanceStrategy
          : public EmptyNewtonStrategy,
            private MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>
      {
       public:
        ComputeDistanceStrategy(LINALG::Matrix<probDim, 1, floatType>& xsi, bool checklimits)
            : xsi_(xsi),
              // we do initializing of distance in SetupSolve, because xsi could be from derived
              // class and not initialize at this point yet
              distance_(NULL),
              xyze_side_(NULL),
              px_(NULL),
              tol_(0.0),
              zeroarea_(false)
        {
        }

        /** \brief Setup the distance calculation problem
         *
         *  \param xyze_side (in) : global coordinates of the considered side
         *  \param px        (in) : global coordinates of the considered point */
        void Setup(const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            const LINALG::Matrix<probDim, 1, floatType>& px, bool zeroarea)
        {
          zeroarea_ = zeroarea;
          xyze_side_ = &xyze_side;
          px_ = &px;

          if (debug)
          {
            std::cout << "\n=== ComputeDistance ===\n";
            std::cout << "Side       = " << DRT::DistypeToString(sideType) << "\n";
            std::cout << "xyze_side_ = " << *xyze_side_;
            std::cout << "px_        = " << *px_;
            std::cout << "\n";
          }
        }
        /// return the absolute value of the distance of the point \c px_ to the given side
        floatType Distance() const
        {
          switch (probDim - dimSide)
          {
            case 1:
              return std::abs(distance_[0]);
            case 2:
              return std::sqrt(distance_[0] * distance_[0] + distance_[1] * distance_[1]);
            default:
              dserror("Unsupported probDim and dimSide combination!");
              exit(EXIT_FAILURE);
          }
          exit(EXIT_FAILURE);
        }

        /// return the signed distance of the point px_ to the given side
        const floatType* SignedDistance() const { return distance_; }

        /// return the local solution vector (parameter space coordinates + distance)
        const LINALG::Matrix<probDim, 1, floatType>& LocalCoordinates() const { return xsi_; }

        /// transform tolerance from global to local coordinates
        void GetLocalTolerance(const floatType& global_tolerance,
            LINALG::Matrix<dimSide, 1, floatType>& scaled_tolerance)

        {
          LINALG::Matrix<probDim, 1, floatType> real_tolerance;
          // set all values of the vector to be equal to the global tolerance
          real_tolerance = global_tolerance;
          // extract part of the side from the Jacobian
          LINALG::Matrix<probDim, dimSide, floatType> A;
          for (unsigned int row = 0; row < probDim; ++row)
          {
            for (unsigned int col = 0; col < dimSide; ++col) A(row, col) = A_(row, col);
          }

          LINALG::Matrix<dimSide, dimSide, floatType> TN_inv;
          LINALG::Matrix<dimSide, probDim, floatType> aux;
          TN_inv.MultiplyTN(A, A);

          TN_inv.Invert();
          aux.MultiplyNT(TN_inv, A);
          scaled_tolerance.Multiply(aux, real_tolerance);
        }

        /// return isInfinite and possibly floatType
        std::pair<bool, floatType> ConditionNumber()
        {
          LINALG::Matrix<probDim, probDim, floatType> A_inv;
          floatType det = A_.Determinant();
          if (KERNEL::closeToZero(det))
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            std::cout << " WARNING: Condition number is equal to infinity in the problem. Stopping "
                         "the increase in precision"
                      << std::endl;
#endif
            return std::make_pair(false, -1.0);
          }
          A_inv.Invert(A_);
          floatType cond = A_inv.Norm2() * A_.Norm2();
          return std::make_pair(true, cond);
        }

        /// set the initial solution vector to zero
        void SetupSolve()
        {
          xsi_ = 0.0;
          distance_ = xsi_.A() + dimSide;
          // evaluate initial rhs value (w/o distance contribution)
          DistanceRHS(*xyze_side_, *px_, b_);
          tol_ = AdaptiveCombinedNewtonTolerance(*xyze_side_, *px_, b_, xsi_(0, 0));
        }

        /** \brief Setup routine for a new Newton step
         *
         *  Setup the system of equations. */
        void SetupStep(int iter)
        {
          if (iter > 0)
          {
            DistanceRHS(*xyze_side_, *px_, b_);
          }

          // build the linear system of equations

          DistanceSystem(*xyze_side_, *px_, distance_, A_, B_, C_, N_, b_);
        }
        /// compute the rhs value ( w/o distance contributions )
        void DistanceRHS(const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            const LINALG::Matrix<probDim, 1, floatType>& px,
            LINALG::Matrix<probDim, 1, floatType>& b)
        {
          LINALG::Matrix<dimSide, 1, floatType> xsi_side(xsi_.A(), true);
          DRT::UTILS::shape_function<sideType>(xsi_side, sideFunct_);
          b = px;
          b.Multiply(-1.0, xyze_side, sideFunct_, 1.0);
        }

        /// Test the stop / convergence criterion of the Newton scheme
        enum NewtonStatus TestConverged(int iter)
        {
          if (debug)
          {
            std::cout << "WithinLimits             = "
                      << (WithinLimits<sideType>(xsi_) ? "TRUE" : "FALSE") << "\n"
                      << "rsd_.Norm2()             = " << xsi_.Norm2() << "\n"
                      << "dx_.Norm2() / rsd.Norm2()= "
                      << (iter != 0 ? (dx_.Norm2() / xsi_.Norm2()) : -1.0) << "\n"
                      << "b_.Norm2() ( RHS )       = " << b_.Norm2() << "\n"
                      << "b_.Norm2() / rsd.Norm2() = "
                      << (iter != 0 ? (b_.Norm2() / xsi_.Norm2()) : -1.0) << "\n"
                      << "tol_                     = " << tol_ << "\n"
                      << "rsd_ ( aka xsi_ )        = " << xsi_
                      << "b_ (RHS)                 = " << b_ << "\n"
                      << std::flush;
          }

          if (zeroarea_)
          {
            return failed;
          }

          floatType residual = b_.Norm2();

          if (residual < tol_)
          {
            // actually I'm not sure why we are doing this here... ( hiermeier 08/16 )
            tol_ *= std::sqrt(static_cast<double>(probDim));
            return converged;
          }
          else
            return unconverged;
        }

        /** solve the linear system:
         *  \f[
         *      A \; \Delta \xi = b
         *  \f] */
        bool LinearSolve(int iter)
        {
          dx_ = 0.0;
          if (zeroarea_) return false;
          if (debug)
          {
            std::cout << "Matrix A  is " << A_ << std::endl;
            std::cout << "Matrix dx is  " << dx_ << std::endl;
            std::cout << "Matrix b  is  " << b_ << std::endl;
          }

          floatType det = LINALG::gaussElimination<true, probDim, floatType>(A_, b_, dx_);

          if (debug)
          {
            std::cout << "det                        = " << det << "\n";
            std::cout << "dx_ ( solution increment ) = " << dx_;
            std::cout << "dx_.Norm2()                = " << dx_.Norm2() << "\n";
          }

          return (det != 0.0);  // here det = 0 will just happen if the surface element is
                                // distorted!
        }

        /// update the local coordinates solution vector
        bool Update(int iter)
        {
          xsi_ += dx_;
          return true;
        }

        /// fall back routine, if the Newton failed
        bool NewtonFailed()
        {
          tol_ = b_.Norm2();
#ifdef DEBUG_CUTKERNEL_OUTPUT
          std::stringstream str;
          str << "ComputeDistance: Newton scheme did not converge:\n"
              << std::setprecision(16) << (*xyze_side_) << (*px_) << xsi_;

          std::string filename = OUTPUT::GenerateGmshOutputFilename(".NewtonFailed_distance.pos");
          std::ofstream file(filename.c_str());
          WritetoGmsh(file);
          file.close();
#endif
          return false;
        }

        /// Gmsh debug output
        void WritetoGmsh(std::ofstream& file)
        {
          file.precision(32);  // higher precision!
          char elementType;
          switch (numNodesSide)
          {
            case 3:
              elementType = 'T';
              break;
            case 4:
              elementType = 'Q';
              break;
            default:
              dserror(
                  "unsupported element type in GmshSideDump."
                  " Please feel free to extend the functionality if necessary.");
          }

          file << "View \""
               << "Side"
               << "\" {\n";
          {
            file << "S" << elementType << "(";
            for (uint i = 0; i < numNodesSide; ++i)
            {
              if (i != 0) file << ",";
              file << (*xyze_side_)(0, i) << "," << (*xyze_side_)(1, i) << ","
                   << (*xyze_side_)(2, i);
            }
            file << "){";
            for (uint i = 0; i < numNodesSide; ++i)
            {
              if (i != 0) file << ",";
              file << "1";
            }
            file << "};\n";
          }
          file << "};\n";

          // calculate the points
          LINALG::Matrix<numNodesSide, 1, floatType> surface;

          DRT::UTILS::shape_function_2D(surface, xsi_(0), xsi_(1), sideType);

          LINALG::Matrix<probDim, 1, floatType> x1;
          ;

          x1.Multiply(*xyze_side_, surface);

          file << "View \""
               << "Point"
               << "\" {\n";
          {
            file << "SP (";
            file << (*px_)(0, 0) << "," << (*px_)(1, 0) << "," << (*px_)(2, 0);
            file << "){";
            file << 1;
            file << "};\n";

            file << "SP (";
            file << x1(0, 0) << "," << x1(1, 0) << "," << x1(2, 0);
            file << "){";
            file << 2;
            file << "};\n";
          }
          file << "};";
        }

        /// get normal vectors
        const LINALG::Matrix<probDim, 2, floatType>& GetNormalVector() { return nvec_; }

        /// get the adapted Newton tolerance
        floatType GetTolerance() const { return tol_; }

        /// get actual reached residual
        floatType GetResidualL2Norm() const { return b_.Norm2(); }

        /// return \TRUE if the side metric vectors is zero
        bool ZeroArea() const { return zeroarea_; }

       private:
        enum NormalPlane
        {
          normal_in_xy_plane,
          normal_in_yz_plane,
          normal_plane_undefined
        };

        enum NormalPlane DetectNormalPlane(const LINALG::Matrix<probDim, 1, floatType>& n) const
        {
          if (probDim < 3) dserror("This function makes only sense for the 3-D case!");

          if (n(2) == 0.0)
            return normal_in_xy_plane;
          else if (n(0) == 0.0)
            return normal_in_yz_plane;

          dserror("Couldn't detect a feasible plane for the given normal vector!");
          exit(EXIT_FAILURE);
        }

        /** \build Build the linear system of equations for the Newton scheme
         *
         *  Right-hand-side:
         *  \f[
         *    b = \hat{x} - x(\xi) - \frac{\tilde{n}}{\| \tilde{n} \|} d,
         *  \f]
         *  where \f$ \hat{x} \in \mathbb{R}^{d}\f$ is the given point,
         *        \f$ \xi \in \mathbb{R}^{d-1}\f$ is the projected in-plane
         *        parameter space coordinate vector,
         *        \f$ \tilde{n} \f$ is the non-unit normal vector on the side
         *        in the point \f$ \xi \f$,
         *        \f$ d \f$ is the distance in normal direction between the point
         *        \f$\hat{x}\f$ and the side.
         *
         *  The matrix \f$ A \f$ holds the consistent linearization.
         *
         *  \author hiermeier
         *  \date 08/16    */
        bool DistanceSystem(const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            const LINALG::Matrix<probDim, 1, floatType>& px, floatType* distance,
            LINALG::Matrix<probDim, probDim, floatType>& A,
            LINALG::Matrix<probDim, probDim, floatType>& B,
            LINALG::Matrix<probDim, 2 * dimSide - 1, floatType>& C,
            LINALG::Matrix<probDim, 2, floatType>& N, LINALG::Matrix<probDim, 1, floatType>& b)
        {
          LINALG::Matrix<probDim, 1, floatType> n1(&N(0, 0), true);
          LINALG::Matrix<probDim, 1, floatType> n2(&N(0, 1), true);
          /* get only the local side coordinates and use them to evaluate
           * the 1-st and 2-nd derivatives */

          LINALG::Matrix<dimSide, 1, floatType> xsi_side(xsi_.A(), true);

          // pre-evaluate some important variables
          B = 0.0;

          floatType det = EvalDerivsInParameterSpace<probDim, sideType, floatType>(
              xyze_side, xsi_side, sideDeriv1_, B, NULL, &n1, &n2, false);

          A.UpdateT(B);


          if (debug)
          {
            std::cout << "det-metric               = " << det << std::endl;
            if (det < 1.0e-16)
              std::cout << "!!! Determinant of jacobian is smaller than 1.0e-16 !!!\n";
          }

          if (det < 1.0e-16 && det > 1.0e-16)
#if EXTENDED_CUT_DEBUG_OUTPUT
            std::cout << "Determinant in  compute position is very close to zero" << std::endl;
#endif
          if (std::abs(det) == 0.0)  // here might lie problem for the cln
          {
            /* then calculation of a normal to a line like this makes no sense at
             * all! */
#if EXTENDED_CUT_DEBUG_OUTPUT
            std::cout << "Absolute value of the determinant is zero " << std::endl;
#endif
            zeroarea_ = true;
            return false;
          }

          // inverse of the normal direction norm
          floatType n1norm_inv = 1.0 / n1.Norm2();
          // --- add the distance part to the right-hand-side
          b.Update(-distance[0] * n1norm_inv, n1, 1.0);

          for (unsigned r = 0; r < probDim; ++r)  // store normal vector
            nvec_(r, 0) = n1(r);
          // --- evaluate the linearization
          // scale the 1-st normal in A to unit length
          for (unsigned r = 0; r < probDim; ++r) A(r, dimSide) *= n1norm_inv;
          // --- evaluate 2-nd derivatives at xsi w.r.t. the parameter space coordinates
          C = 0.0;
          DRT::UTILS::shape_function_deriv2<sideType>(xsi_side, sideDeriv2_);
          C.MultiplyNT(xyze_side, sideDeriv2_);

          // reset B-matrix
          B = 0.0;
          // 1-dimensional side ( a.k.a. edge or line ) embedded in 3-dimensional space
          if (probDim == 3 and dimSide == 1)
          {
            floatType n2norm_inv = 1.0 / n2.Norm2();

            // scale the 2-nd normal in the last column of A to unit length
            for (unsigned r = 0; r < probDim; ++r) A(r, 2) *= n2norm_inv;

            for (unsigned r = 0; r < probDim; ++r) nvec_(r, 1) = n2(r);

            // the normal linearization can be skipped for this case
            if (Distance() == 0.0) return true;

            // add extra term to the right hand side
            b.Update(-distance[1] * n2norm_inv, n2, 1.0);
            // update stored normal vector

            // linearization of the 1-st unscaled normal vector
            switch (DetectNormalPlane(n1))
            {
              case normal_in_xy_plane:
              {
                // d(n1_0) / dr =   d^2(y) / dr^2
                B(0, 0) = C(1, 0);
                // d(n1_1) / dr = - d^2(x) / dr^2
                B(1, 0) = -C(0, 0);
                // d(n1_2) / dr =   0.0
                break;
              }
              case normal_in_yz_plane:
              {
                // d(n1_0) / dr =    0.0
                // d(n1_1) / dr =  - d^2(z) / dr^2
                B(1, 0) = -C(2, 0);
                // d(n1_2) / dr =    d^2(y) / dr^2
                // u
                B(2, 0) = C(1, 0);
                break;
              }
              default:
              {
                dserror("Shouldn't happen!");
                exit(EXIT_FAILURE);
              }
            }

            // linearization of the 2-nd unscaled normal vector
            // d(n2_0) / dr
            B(0, 1) = (C(1, 0) * n1(2) + A(1, 0) * B(2, 0) - (C(2, 0) * n1(1) + A(2, 0) * B(1, 0)));
            // d(n2_1) / dr
            B(1, 1) = (C(2, 0) * n1(0) + A(2, 0) * B(0, 0) - (C(0, 0) * n1(2) + A(0, 0) * B(2, 0)));
            // d(n2_2) / dr
            B(2, 1) = (C(0, 0) * n1(1) + A(0, 0) * B(1, 0) - (C(1, 0) * n1(0) + A(1, 0) * B(0, 0)));

            // linearization of the scaling factor ( 1-st normal )
            floatType lin_nnorm = 0.0;

            lin_nnorm =
                (n1(0) * B(0, 0) + n1(1) * B(1, 0) + n1(2) * B(2, 0)) * n1norm_inv * n1norm_inv;
            for (unsigned r = 0; r < probDim; ++r) B(r, 0) -= lin_nnorm * A(r, 1);

            // linearization of the scaling factor ( 2-nd normal )
            lin_nnorm =
                (n2(0) * B(0, 1) + n2(1) * B(1, 1) + n2(2) * B(2, 1)) * n2norm_inv * n2norm_inv;
            for (unsigned r = 0; r < probDim; ++r) B(r, 1) -= lin_nnorm * A(r, 2);

            // final scaling and add
            for (unsigned r = 0; r < probDim; ++r)
            {
              A(r, 0) += distance[0] * B(r, 0);
              A(r, 0) += distance[1] * B(r, 1);
            }
          }
          // 2-dimensional side embedded in 3-dimensional space
          else if (probDim == 3 and dimSide == 2)
          {
            // the normal linearization can be skipped for this case
            if (Distance() == 0.0) return true;

            // linearization of the unscaled normal vector
            // d(n_0) / dr
            B(0, 0) =
                (C(1, 0) * A(2, 1) + A(1, 0) * C(2, 2) - (C(2, 0) * A(1, 1) + A(2, 0) * C(1, 2)));
            // d(n_0) / ds
            B(0, 1) =
                (C(1, 2) * A(2, 1) + A(1, 0) * C(2, 1) - (C(2, 2) * A(1, 1) + A(2, 0) * C(1, 1)));

            // d(n_1) / dr
            B(1, 0) =
                (C(2, 0) * A(0, 1) + A(2, 0) * C(0, 2) - (C(0, 0) * A(2, 1) + A(0, 0) * C(2, 2)));
            // d(n_1) / ds
            B(1, 1) =
                (C(2, 2) * A(0, 1) + A(2, 0) * C(0, 1) - (C(0, 2) * A(2, 1) + A(0, 0) * C(2, 1)));

            // d(n_2) / dr
            B(2, 0) =
                (C(0, 0) * A(1, 1) + A(0, 0) * C(1, 2) - (C(1, 0) * A(0, 1) + A(1, 0) * C(0, 2)));
            // d(n_2) / ds
            B(2, 1) =
                (C(0, 2) * A(1, 1) + A(0, 0) * C(1, 1) - (C(1, 2) * A(0, 1) + A(1, 0) * C(0, 1)));

            // linearization of the scaling factor fact
            floatType lin_nnorm = 0.0;

            for (unsigned c = 0; c < dimSide; ++c)
            {
              lin_nnorm =
                  (n1(0) * B(0, c) + n1(1) * B(1, c) + n1(2) * B(2, c)) * n1norm_inv * n1norm_inv;
              for (unsigned r = 0; r < probDim; ++r) B(r, c) -= lin_nnorm * A(r, 2);
            }

            // final scaling and add
            A.Update(distance[0], B, 1.0);
          }
          // 1-dimensional side/line element embedded in 2-dimensional space
          else if (probDim == 2 and dimSide == 1)
          {
            // the normal linearization can be skipped for this case
            if (Distance() == 0.0) return true;

            // linearization of the unscaled normal vector
            B(0, 0) = C(1, 0);
            B(1, 0) = -C(0, 0);

            // linearization of the scaling factor fact
            floatType lin_nnorm = (C(1, 0) + C(0, 0)) * n1norm_inv * n1norm_inv;
            for (unsigned r = 0; r < probDim; ++r) B(r, 0) -= lin_nnorm * A(r, 1);

            // final scaling and add
            A.Update(distance[0], B, 1.0);
          }
          else
            dserror("Unsupported dim <--> probDim relation!");
          return true;
        }

        /// solution vector (parameter space coordinates + distance values)
        LINALG::Matrix<probDim, 1, floatType>& xsi_;

        /// distance values
        floatType* distance_;

        const LINALG::Matrix<probDim, numNodesSide, floatType>* xyze_side_;
        /// given global point coordinates
        const LINALG::Matrix<probDim, 1, floatType>* px_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide,
            floatType>::sideFunct_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide,
            floatType>::sideDeriv1_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide,
            floatType>::sideDeriv2_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::A_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::B_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::C_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::b_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::dx_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide, floatType>::N_;

        /// adapted Newton tolerance
        floatType tol_;

        /** \brief  indicates that the surface element has a area value of zero
         *
         *  --> cross product for normal is zero!!! */
        bool zeroarea_;

        using MemberStoragePolicy<debug, probDim, sideType, dimSide, numNodesSide,
            floatType>::nvec_;
      };  // class ComputeDistanceStrategy


#ifdef CLN_CALC
      /*--------------------------------------------------------------------------*/
      /** \brief Strategy class for  distance between side and point for arbitrary precision
       *
       *  inheritance diagram:
       *
       *  ComputeDistance --> ComputeDistanceAdaptivePrecision --> NewtonSolve
       *                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy  */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class ComputeDistanceAdaptivePrecision : Strategy
#ifdef DEBUG_MEMORY_ALLOCATION
          ,
                                               AdaptiveKernelSharedData
#endif
      {
       public:
        /// constructor
        ComputeDistanceAdaptivePrecision(LINALG::Matrix<probDim, 1>& xsi, bool checklimits)
            : Strategy(clnxsi_, checklimits), xsi_(xsi), clnxsi_(true), cond_infinity_(false)
        {
        }

        /// start and solve the distance calculation
        bool operator()(const LINALG::Matrix<probDim, numNodesSide>& xyze_side,
            const LINALG::Matrix<probDim, 1>& px, double& distance, bool signeddistance = false)
        {
          if (!((sideType == DRT::Element::hex8) || (sideType == DRT::Element::quad4) ||
                  (sideType == DRT::Element::hex8) || (sideType == DRT::Element::tri3) ||
                  (sideType == DRT::Element::line2)) ||
              (sideType == DRT::Element::hex20) || (sideType == DRT::Element::tet4) ||
              (sideType == DRT::Element::pyramid5))
          {
            dserror(
                "This type of element (%s)  is not tested for the CLN calculation. You are welcome "
                "to edit ../drt_fem_general/drt_utils_fem_shapefunctions.H file to fix it. Just "
                "change all the integers occuring there double for CLN to work.",
                DRT::DistypeToString(sideType).c_str());
          }

          int prec = CLN_START_PRECISION;  // standart initial precision
          bool conv;
          ClnWrapper err;
          ClnWrapper cond_number;
          int iter = 0;
          cond_infinity_ = false;

          // converting everythign to the cln precision
#ifdef CUSTOM_MEMORY_ALLOCATOR
          if (all_distance_done_once_ and all_intersections_done_once_ and
              all_position_done_once_ and (!custom_allocator_run_))
          {
#if DEBUG_MEMORY_ALLOCATION
            GEO::CUT::MemorySingleton::getInstance().ReportAllocated();
            ReportIntersectionAllocated();
            ReportPositionAllocated();
            ReportDistanceAllocated();
            ReportTotalAllocated();
            // start setting up the memory container
            GEO::CUT::MemorySingleton::getInstance().SetState(1, memory_allocations_);
#else
            GEO::CUT::MemorySingleton::getInstance().SwitchState();
#endif
            custom_allocator_run_ = true;
          }
#if DEBUG_MEMORY_ALLOCATION
          if (first_run_) std::cout << "In cut distance statistics is" << std::endl;
#endif
#endif
          do
          {
            ClnWrapper::SetPrecision(prec);
#ifdef CUSTOM_MEMORY_ALLOCATOR
#if DEBUG_MEMORY_ALLOCATION
            if (custom_allocator_run_)
#endif
              GEO::CUT::MemorySingleton::getInstance().GetMemoryPoolAllocator().SetCurrent(
                  cln_byte_size_[iter]);
#endif
#if DEBUG_MEMORY_ALLOCATION

            if (first_run_)
            {
              if (GEO::CUT::MemorySingleton::getInstance().IsRecording())
              {
                UpdateMemoryAllocations(
                    GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern());
                GEO::CUT::MemorySingleton::getInstance().StopRecord();
              }
              GEO::CUT::MemorySingleton::getInstance().StartRecord();
            }

#endif
            ConvDoubleCLN(xyze_side, clnxyze_side_, prec);
            ConvDoubleCLN(px, clnpx_, prec);
            //  no need to convert, we don't actually compute anything
            // clndistance_ = cln::cl_float(distance, cln::float_format(prec));
            this->Setup(clnxyze_side_, clnpx_, false);

            conv = this->Solve();

#if CUT_DEVELOP
            // safety check for precision loss across the array
            cln::float_format_t prec_beg = cln::float_format(clnxyze_side_(0, 0));
            for (unsigned int i = 0; i < probDim; ++i)
            {
              cln::float_format_t prec_end = cln::float_format(clnxsi_(i, 0));
              if (prec_beg != prec_end)
              {
                dserror(
                    "There is a loss of cln-precision during computation of Intersection. "
                    "Something is wrong with the conversion");
              }
            }
#endif
            if (not signeddistance)
              clndistance_ = this->Distance();
            else
            {
              switch (probDim - dimSide)
              {
                case 1:
                  clndistance_ = this->SignedDistance()[0];
                  break;
                default:
                  dserror("A scalar signed distance value is not available!");
                  exit(EXIT_FAILURE);
              }
            }
            std::pair<bool, ClnWrapper> cond_pair = this->ConditionNumber();
            cond_number = cond_pair.second;

            if (not cond_pair.first)
            {
              cond_infinity_ = true;
            }
            err = ComputeError(xyze_side, px, this->LocalCoordinates(), this->GetNormalVector(),
                      this->SignedDistance(), CLN_REFERENCE_PREC) *
                  cond_number;
#if DEBUG_MEMORY_ALLOCATION
            UpdateMemoryUsage(prec, iter);
#endif
            prec += CLN_INCREMENT_STEP;
            iter++;

          } while (((!conv) || (err > CLN_LIMIT_ERROR) || (cond_infinity_)) &&
                   (prec < CLN_LIMIT_PREC) && (iter < CLN_LIMIT_ITER));

          // now since we converged or went over maximum allowed number of loops  we convert the
          // values back if we did not converge, by newton tolerance, but the error is small enough
          // treat as converged
          if ((not conv) && (err < CLN_LIMIT_ERROR) && (err > 0)) conv = true;

          if (cond_infinity_)
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            std::cout << "Condition number of this problem is equal to infinity" << std::endl;
#endif
          }

          GetTopologyInformation();
          // NOTE: Might be not needed later
          FixCornerCase();
          // converting all the values back
          ConvClnDouble(clnxsi_, xsi_);
          distance = cln::double_approx(clndistance_);

#ifdef CUSTOM_MEMORY_ALLOCATOR
          if (first_run_)
          {
            const int num_inter = 3;
            if ((++GetDistanceCount()) == num_inter)
            {
              all_distance_done_once_ = true;
              std::cout << "All possible template distance are done!" << std::endl;
            }
            first_run_ = false;
          }
#endif
          return conv;
        }
        // get the local coordinates
        const LINALG::Matrix<probDim, 1, ClnWrapper>& LocalCoordinates() const
        {
          return Strategy::LocalCoordinates();
        }

        const ClnWrapper* SignedDistance() const { return Strategy::SignedDistance(); }

        ClnWrapper Distance() const { return Strategy::Distance(); }

        std::pair<bool, ClnWrapper> ConditionNumber() { return Strategy::ConditionNumber(); }

        const LINALG::Matrix<probDim, 2, ClnWrapper> GetNormalVector()
        {
          return Strategy::GetNormalVector();
        }

        /// access the Newton tolerance
        ClnWrapper GetTolerance() const { return Strategy::GetTolerance(); }
        PointOnSurfaceLoc GetSideLocation() { return location_; }
        const std::vector<int>& GetTouchedSideEdges() { return touched_edges_ids_; }

        const std::vector<int>& GetTouchedNodes()
        {
          return touched_nodes_ids_;
          ;
        }

        void WritetoGmsh(std::ofstream& file) { Strategy::WritetoGmsh(file); }

        bool IsConditionInfinity() { return cond_infinity_; }

       private:
        // Evaluate difference between inital global coordinates passed to ComputeDistance and
        // and global coordinates based on loc coordinates and distance  calculated in the
        // ComputeDistance, using conversion to very high reference precision
        ClnWrapper ComputeError(
            const LINALG::Matrix<probDim,
                DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>&
                refshape_xyze,                    // referenape
            const LINALG::Matrix<probDim, 1>& p,  // global position of the point
            const LINALG::Matrix<probDim, 1, ClnWrapper>&
                loc_calc,  // calculated wrt reference shape local position
            const LINALG::Matrix<probDim, 2, ClnWrapper>& nvec,  // normal vector
            const ClnWrapper*
                distance,  // distance from the projection on the reference shape  to the point
            int prec       // precision for calculation
        )
        {
          unsigned int prev_prec = ClnWrapper::GetPrecision();
          ClnWrapper::SetPrecision(prec);
          // Converting input arrays  to higher precision floating points
          ClnWrapper clndistance;
          LINALG::Matrix<probDim, DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
              ClnWrapper>
              xyze_side;
          LINALG::Matrix<probDim, 1, ClnWrapper> clnxi;
          LINALG::Matrix<probDim, 1, ClnWrapper> clnpx;
          LINALG::Matrix<probDim, 1, ClnWrapper> n1;

          clndistance = cln::cl_float(distance[0], cln::float_format(prec));
          ConvDoubleCLN(refshape_xyze, xyze_side, prec);
          ConvDoubleCLN(p, clnpx, prec);
          ConvDoubleCLN(loc_calc, clnxi, prec);

          for (unsigned int i = 0; i < probDim; ++i)
          {
            n1(i) = cln::cl_float(nvec(i, 0), cln::float_format(prec));
          }

          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement, 1,
              ClnWrapper>
              surfaceFunct;
          LINALG::Matrix<probDim, 1, ClnWrapper> b;

          LINALG::Matrix<dimSide, 1, ClnWrapper> clnxiside(clnxi.A(), true);
          DRT::UTILS::shape_function<sideType>(clnxiside, surfaceFunct);
          b = clnpx;
          b.Multiply(-1.0, xyze_side, surfaceFunct, 1.0);

          // inverse of the normal direction norm
          ClnWrapper n1norm_inv = 1.0 / n1.Norm2();
          // --- add the distance part to the right-hand-side
          b.Update(-clndistance * n1norm_inv, n1, 1.0);
          // in a special case add more of the distance
          if (probDim == 3 and dimSide == 1)
          {
            LINALG::Matrix<probDim, 1, ClnWrapper> n2;
            for (unsigned int i = 0; i < probDim; ++i)
              n2(i) = cln::cl_float(nvec(i, 1), cln::float_format(prec));

            ClnWrapper n2norm_inv = 1.0 / n2.Norm2();
            b.Update(-distance[1] * n2norm_inv, n2, 1.0);
          }
          ClnWrapper res = b.Norm2();
          // resetting precision to previous
          ClnWrapper::SetPrecision(prev_prec);
          return res;
        }

        enum PointOnSurfacePlane
        {
          above = 0,
          below = 1,
          on = 2
        };

        // Transform tolerances into local coordinate system and get and
        // get location of the point on the surface  as well as touched edges
        void GetTopologyInformation()
        {
          LINALG::Matrix<dimSide, 1, ClnWrapper> scaled_tolerance_side_touched_edges;
          this->GetLocalTolerance(SIDE_DETECTION_TOLERANCE, scaled_tolerance_side_touched_edges);
          LINALG::Matrix<dimSide, 1, ClnWrapper> zero_tolerance_side;
          bool is_inside = WithinLimits<sideType>(clnxsi_, zero_tolerance_side);
          PointOnSurfacePlane point_on_surface = GetLocation(SIDE_DETECTION_TOLERANCE);

          location_ = PointOnSurfaceLoc(is_inside, point_on_surface == on);
          touched_edges_ids_.clear();
          touched_nodes_ids_.clear();
          // remove edges, that are too far away
          if (location_.OnSide())
          {
            GetEdgesAt<sideType>(clnxsi_, touched_edges_ids_, scaled_tolerance_side_touched_edges);
            GetNodesAt<sideType>(clnxsi_, touched_nodes_ids_, scaled_tolerance_side_touched_edges);
          }
        }

        /// try to fix the case that is close to the corner and outside - in that case point
        /// is close to corner and outside might be touching two edges but difference between it
        /// and corner point > tolerance - hence this must be corrected
        void FixCornerCase()
        {
          if (touched_edges_ids_.size() >= 2)
          {
            // if we are inside the "tolerance wrapper" around the side
            if (location_.WithinSide())
            {
              LINALG::Matrix<dimSide, 1, ClnWrapper> zero_tolerance_side;
              // this means we are on the outer boundary of the tolerance
              if (not WithinLimits<sideType>(clnxsi_, zero_tolerance_side))
              {
                // we transform our point to global tolerance and check if it is close enough to the
                // nodal point of the side, that is common point of the touching edges
                const LINALG::Matrix<probDim, numNodesSide, ClnWrapper>& clnxyze_side_latest =
                    clnxyze_side_;
                const LINALG::Matrix<probDim, 1, ClnWrapper>& xsi_global = clnpx_;

                // find minimum distance and check if it is closer than 1e-14
                LINALG::Matrix<probDim, 1, ClnWrapper> dist;
                ClnWrapper min_dist;
                for (unsigned int inode = 0; inode < numNodesSide; ++inode)
                {
                  // acess raw responsible for that node
                  LINALG::Matrix<probDim, 1, ClnWrapper> coord(
                      clnxyze_side_latest.A() + inode * probDim, true);
                  dist.Update(1.0, coord, -1.0, xsi_global);
                  ClnWrapper tmp_dist = dist.Norm2();
                  if ((inode == 0) or (tmp_dist < min_dist))
                  {
                    min_dist = tmp_dist;
                  }
                }

                if (min_dist > NON_TOPOLOGICAL_TOLERANCE)
                {
                  // then it will not be merged to this point in the pointpool. Hence we remove the
                  // topological connection and change location of the point to outisde
                  touched_edges_ids_.clear();
                  // toggle the bit off,
                  location_ = PointOnSurfaceLoc(false, location_.OnSide());
                }
              }
            }
          }
        }

        // Update global maximum number of allocations ( number of containers of particular byte
        // size) if allocations on this precisions were more. Also set what was the most frequent
        // size of allocation for the current iteration
        void UpdateMemoryUsage(int prec, int iter)
        {
#if DEBUG_MEMORY_ALLOCATION

          if (first_run_)
          {
            GEO::CUT::MemorySingleton::getInstance().StopRecord();

            boost::unordered_map<size_t, int>& allocation_map =
                GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern();
            int max_num = 0;
            int size_max_num = 0;
            // update the global number of allocations to the maximum
            for (boost::unordered_map<size_t, int>::iterator it = allocation_map.begin();
                 it != allocation_map.end(); ++it)
            {
              size_t size = it->first;
              if (memory_allocations_distance_[size] < allocation_map[size])
              {
                memory_allocations_distance_[size] = allocation_map[size];
              }
              if (allocation_map[size] > max_num)
              {
                size_max_num = size;
                max_num = allocation_map[size];
              }
            }

            if (size_max_num != 0)
            {
              std::cout << "Size for the precision " << prec << " is " << size_max_num << std::endl;
              std::cout << "It is allocated " << max_num << " times " << std::endl;
              cln_sizes_[iter] = size_max_num;
            }
            else
              dserror("This should not be possible");

            GEO::CUT::MemorySingleton::getInstance().ResetAllocated();
          }
#endif
        }

        // Get location of the point with respect to the plain ( above/on/below )
        PointOnSurfacePlane GetLocation(const ClnWrapper& err)
        {
          switch (probDim - dimSide)
          {
            case 1:
            {
              PointOnSurfacePlane result;
              const ClnWrapper* signed_distance = Strategy::SignedDistance();

              if ((signed_distance[0] + err < 0.0) && (signed_distance[0] - err < 0.0))
              {
                // it is below the side (distance is negative and bigger than tolerance)
                result = below;
              }
              else if ((signed_distance[0] + err > 0.0) && (signed_distance[0] - err > 0.0))
              {
                // it is above the side (distance is negative and higher than tolerance)
                result = above;
              }
              else
              {
                // else it is on surface
                result = on;
              }
              return result;
              break;
            }
            default:
            {
              // we cat only get signed distance for edge edge intersection
              const ClnWrapper distance = Strategy::Distance();
              // some safety checks
              if (distance < 0) dserror("Not possible");
              if (err <= 0) dserror("Error should be equal to basic tolerance!");
              if ((distance - err <= 0.0) and (distance + err > 0.0))
                return on;
              else
                return below;
            }
          }
        }

        // touched edges ids
        static std::vector<int> touched_edges_ids_;

        // touched edges ids
        static std::vector<int> touched_nodes_ids_;

        // variable for cln computation
        LINALG::Matrix<probDim, numNodesSide, ClnWrapper> clnxyze_side_;

        // global position of the point converted to cln precision
        LINALG::Matrix<probDim, 1, ClnWrapper> clnpx_;

        ClnWrapper clndistance_;

        // reference to the double TMatrix
        LINALG::Matrix<probDim, 1>& xsi_;

        LINALG::Matrix<probDim, 1, ClnWrapper> clnxsi_;

        // determines position of this point with respect to the side limits (values:
        // inside/outsaide)
        PointOnSurfaceLoc location_;

        // is condition number infinity or not
        bool cond_infinity_;

        static bool first_run_;

        // sizes of cln variables with different preciion
        static size_t cln_sizes_[CLN_LIMIT_ITER];
      };  // class ComputeDistanceAdaptivePrecision

#endif

      /*--------------------------------------------------------------------------*/
      /** \brief generic distance between side and point
       *
       *  inheritance diagram:
       *
       *  ComputeDistance --> GenericComputeDistance --> NewtonSolve
       *                      ^^^^^^^^^^^^^^^^^^^^^^
       *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy  */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
          bool computeCln = false, unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class GenericComputeDistance : Strategy
      {
       public:
        /// constructor
        GenericComputeDistance(LINALG::Matrix<probDim, 1>& xsi, bool checklimits)
            : Strategy(xsi, checklimits),
              xsi_ref_(xsi),
              cond_infinity_(false)
#ifdef CLN_CALC
              ,
              checklimits_ref_(checklimits)
#endif
        {
        }

        /// start and solve the distance calculation
        bool operator()(const LINALG::Matrix<probDim, numNodesSide>& xyze_side,
            const LINALG::Matrix<probDim, 1>& px, double& distance, bool signeddistance = false)
        {
          // now we are doing basic iteration on double precision
          bool conv;
          this->Setup(xyze_side, px, false);
          conv = this->Solve();

          if (not signeddistance)
            distance = this->Distance();
          else
          {
            switch (probDim - dimSide)
            {
              case 1:
                distance = this->SignedDistance()[0];
                break;
              default:
                dserror("A scalar signed distance value is not available!");
                exit(EXIT_FAILURE);
            }
          }
          GetTopologyInformation();
#ifdef CLN_CALC
          if (computeCln)
          {
#if DOUBLE_PLUS_CLN_COMPUTE
            if ((ComputeError(xyze_side, px) > DOUBLE_LIMIT_ERROR) or (not conv) or
                (cond_infinity_))
            {
              // CutKernelStatistics::getCutKernelStatistics().ClnDistanceCounter();

#endif
              {
                ComputeDistanceAdaptivePrecision<
                    NewtonSolve<ComputeDistanceStrategy<false, probDim, sideType, dimSide,
                                    numNodesSide, ClnWrapper, ComputeDistanceNoStaticMembers>,
                        probDim>,
                    probDim, sideType>
                    cln_calc(xsi_ref_, checklimits_ref_);

                bool distanceworked = cln_calc(xyze_side, px, distance, signeddistance);
                conv = distanceworked;
                location_ = cln_calc.GetSideLocation();
                // assign other location
                touched_edges_ids_.clear();
                touched_nodes_ids_.clear();
                const std::vector<int>& touched_edges_ids_cln = cln_calc.GetTouchedSideEdges();
                const std::vector<int>& touched_nodes_ids_cln = cln_calc.GetTouchedNodes();
                touched_edges_ids_ = touched_edges_ids_cln;
                touched_nodes_ids_ = touched_nodes_ids_cln;
                cond_infinity_ = cln_calc.IsConditionInfinity();
                ClnWrapper::ResetPrecision();
              }
              // finalize memory allocator
#ifdef CUSTOM_MEMORY_ALLOCATOR
              GEO::CUT::MemorySingleton::getInstance().Finalize();
#endif
#if DOUBLE_PLUS_CLN_COMPUTE
            }
            else
            {
              // CutKernelStatistics::getCutKernelStatistics().DoubleDistanceCounter();
            }
#endif
          }
#endif
          return conv;
        }
        // get the local coordinates
        const LINALG::Matrix<probDim, 1>& LocalCoordinates() const
        {
          return Strategy::LocalCoordinates();
        }

        const double* SignedDistance() const { return Strategy::SignedDistance(); }

        double Distance() const { return Strategy::Distance(); }

        std::pair<bool, double> ConditionNumber() { return Strategy::ConditionNumber(); }

        const LINALG::Matrix<probDim, 2> GetNormalVector() { return Strategy::GetNormalVector(); }

        bool SurfaceWithinLimits(double tolerance = REFERENCETOL) const
        {
          return GEO::CUT::KERNEL::WithinLimits<sideType>(LocalCoordinates(), tolerance);
        }

        /// access the Newton tolerance
        double GetTolerance() const { return Strategy::GetTolerance(); }

        /// access the l2 norm of the reached residual
        double GetResidualL2Norm() const { return Strategy::GetResidualL2Norm(); }

        bool ZeroArea() const { return Strategy::ZeroArea(); }
        const std::vector<int>& GetTouchedSideEdges() { return touched_edges_ids_; }

        const std::vector<int>& GetTouchedNodes() { return touched_nodes_ids_; }

        void GetTouchedSideEdges(std::vector<int>& inout)
        {
          inout.clear();
          inout = touched_edges_ids_;
        }

        bool IsConditionInfinity() { return cond_infinity_; }


        // Get location of the point, based on the extended tolerance of 1e-14 for the diagonal
        // edge of the splitted triangle, and tolerance 0.0 for the normal edges, used
        // to avoid possibility of middle point of the quad4 split ot be missed because of
        // numerical errors
        PointOnSurfaceLoc GetSideLocationTriangleSplit()
        {
          if (sideType != DRT::Element::tri3)
            dserror("This method only works for tri3 side. Current side is %s",
                DRT::DistypeToString(sideType).c_str());

          LINALG::Matrix<dimSide, 1> scaled_tolerance;
          double distance_tolerance = TOPOLOGICAL_TOLERANCE;
          // get tolerance with 1e-14 around the triangle
          this->GetLocalTolerance(distance_tolerance, scaled_tolerance);
          // Diagonal one, corresponds to the middle one, tolerance there should be original
          LINALG::Matrix<3, 1> real_tolerance;
          real_tolerance(0) = 0.0;
          real_tolerance(1) = scaled_tolerance(1);
          real_tolerance(2) = 0.0;
          if (WithinLimitsSplittedQuad<sideType>(xsi_ref_, real_tolerance))
            sideLocation_triangle_split_ = PointOnSurfaceLoc(true, true);
          else
            sideLocation_triangle_split_ = PointOnSurfaceLoc(false, true);
          return sideLocation_triangle_split_;
        }

        PointOnSurfaceLoc GetSideLocation() { return location_; }

        //  Converts tolerance in the local coordinates and computes topological information, such
        //  as touched_edges and  location on the surface
        void GetTopologyInformation()
        {
          LINALG::Matrix<probDim, 1> xsi = this->LocalCoordinates();
          // set up required tolerances
          LINALG::Matrix<dimSide, 1> scaled_tolerance;
          double distance_tolerance = SIDE_DETECTION_TOLERANCE;
          this->GetLocalTolerance(distance_tolerance, scaled_tolerance);
          LINALG::Matrix<dimSide, 1> zero_tolerance_side;
          zero_tolerance_side = 0.0;
          // find location
          bool is_inside = WithinLimits<sideType>(xsi, zero_tolerance_side);
          PointOnSurfacePlane point_on_surface = GetLocation(distance_tolerance);
          location_ = PointOnSurfaceLoc(is_inside, point_on_surface == on);
          // clear previous and get edges nearby
          touched_edges_ids_.clear();
          touched_nodes_ids_.clear();
          // if on side, get neigboring edges
          if (location_.OnSide())
          {
            GetEdgesAt<sideType>(xsi, touched_edges_ids_, scaled_tolerance);
            GetNodesAt<sideType>(xsi, touched_nodes_ids_, scaled_tolerance);
          }
        }

        void WritetoGmsh(std::ofstream& file) { Strategy::WritetoGmsh(file); }

        // Evaluate difference between inital global coordinates passed to ComputeDistance and
        // and global coordinates based on loc coordinates and distance  calculated in the
        // ComputeDistance
        double ComputeError(
            const LINALG::Matrix<probDim,
                DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>& xyze_side,
            const LINALG::Matrix<probDim, 1>& px)
        {
          const LINALG::Matrix<probDim, 1>& xsi = this->LocalCoordinates();
          const LINALG::Matrix<probDim, 2>& n_vec = this->GetNormalVector();
          const double* distance = this->SignedDistance();
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement, 1>
              surfaceFunct;
          LINALG::Matrix<probDim, 1> b;
          LINALG::Matrix<probDim, 1> n1(n_vec.A(), true);

          LINALG::Matrix<dimSide, 1> xsi_side(xsi.A(), true);
          DRT::UTILS::shape_function<sideType>(xsi_side, surfaceFunct);
          b = px;

          b.Multiply(-1.0, xyze_side, surfaceFunct, 1.0);

          double n1norm_inv = 1.0 / n1.Norm2();
          b.Update(-distance[0] * n1norm_inv, n1, 1.0);
          if (probDim == 3 and dimSide == 1)
          {
            LINALG::Matrix<probDim, 1> n2(n_vec.A() + probDim, true);
            double n2norm_inv = 1.0 / n2.Norm2();
            b.Update(-distance[1] * n2norm_inv, n2, 1.0);
          }
          double res = b.Norm2();
          if (res == 0) res = MINIMUM_DOUBLE_ERROR;
          std::pair<bool, double> cond_pair = this->ConditionNumber();
          if (not cond_pair.first)
          {
            cond_infinity_ = true;
          }
          return (res * cond_pair.second);
        }

       private:
        enum PointOnSurfacePlane
        {
          above = 0,
          below = 1,
          on = 2
        };

        // Get location of the point with respect to the plain ( above/on/below )
        PointOnSurfacePlane GetLocation(const double& err)
        {
          switch (probDim - dimSide)
          {
            case 1:
            {
              PointOnSurfacePlane result;
              // use it anywhere for now.
              const double* signed_distance = Strategy::SignedDistance();

              if ((signed_distance[0] + err < 0.0) && (signed_distance[0] - err < 0.0))
              {
                /// it is below the side (distance is negative and bigger than tolerance)
                result = below;
              }
              else if ((signed_distance[0] + err > 0.0) && (signed_distance[0] - err > 0.0))
              {
                /// it is above the side (distance is negative and higher than tolerance)
                result = above;
              }
              else
              {
                // else it is on surface
                result = on;
              }
              return result;
              break;
            }
            default:
            {
              // we cat only get signed distance for edge edge intersection
              const double distance = Strategy::Distance();
              // some safety checks
              if (distance < 0) dserror("Not possible");
              if (err <= 0) dserror("Error should be equal to basic tolerance!");
              if ((distance - err <= 0.0) and (distance + err > 0.0))
                return on;
              else
                return below;
            }
          }
        }
        // That encodes the location, last two bits represent whether point is on the surface and
        // within its limits
        PointOnSurfaceLoc location_;

        // touched edges ids
        static std::vector<int> touched_edges_ids_;

        // touched nodes id
        static std::vector<int> touched_nodes_ids_;

        // determines location of this point, with extended tolerance over side, that is shared for
        // another triangle (only for tri3)
        PointOnSurfaceLoc sideLocation_triangle_split_;
        // we need to do detection of location and touching edge
        LINALG::Matrix<probDim, 1>& xsi_ref_;

        // if condition number is equal to infinity
        bool cond_infinity_;
#ifdef CLN_CALC
        // Need to remember it in order to call to ComputeDistanceAdaptivePrecision
        bool checklimits_ref_;
#endif
      };  // class GenericComputeDistance

      /*--------------------------------------------------------------------------*/
      /** \brief most derived strategy class for distance between side and point
       *
       *  inheritance diagram:
       *
       *  ComputeDistance --> GenericComputeDistance --> NewtonSolve
       *  ^^^^^^^^^^^^^^^
       *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy  */
      template <unsigned probDim, DRT::Element::DiscretizationType sideType,
          bool computeCln = false, unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class ComputeDistance
          : public GenericComputeDistance<
                NewtonSolve<ComputeDistanceStrategy<false, probDim, sideType>, probDim>, probDim,
                sideType, computeCln>
      {
       public:
        /// constructor
        ComputeDistance(LINALG::Matrix<probDim, 1>& xsi, bool checklimits = true)
            : GenericComputeDistance<
                  NewtonSolve<ComputeDistanceStrategy<false, probDim, sideType>, probDim>, probDim,
                  sideType, computeCln>(xsi, checklimits)
        {
        }
      };  // class ComputeDistance

      /*--------------------------------------------------------------------------*/
      template <unsigned probDim, DRT::Element::DiscretizationType sideType,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class DebugComputeDistance
          : public GenericComputeDistance<
                NewtonSolve<
                    DebugNewtonStrategy<ComputeDistanceStrategy<true, probDim, sideType>, probDim>,
                    probDim>,
                probDim, sideType>
      {
       public:
        /// constructor
        DebugComputeDistance(LINALG::Matrix<probDim, 1>& xsi, bool checklimits = true)
            : GenericComputeDistance<
                  NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true, probDim, sideType>,
                                  probDim>,
                      probDim>,
                  probDim, sideType>(xsi, checklimits)
        {
        }
      };  // class DebugComputeDistance

      // Static storage class of ComputeIntersection data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double>
      struct ComputeIntersectionStaticMembers
      {
        static LINALG::Matrix<numNodesSide, 1, floatType> sideFunct_;
        static LINALG::Matrix<numNodesEdge, 1, floatType> edgeFunct_;
        static LINALG::Matrix<dimSide, numNodesSide, floatType> sideDeriv1_;
        static LINALG::Matrix<dimEdge, numNodesEdge, floatType> edgeDeriv1_;

        static LINALG::Matrix<dimEdge + dimSide, dimEdge + dimSide, floatType> A_;
        static LINALG::Matrix<probDim, dimEdge + dimSide, floatType> B_;
        static LINALG::Matrix<dimEdge + dimSide, 1, floatType> b_;
        static LINALG::Matrix<probDim, 1, floatType> c_;

        /// increment in local coordinates d(xi1_side, xi2_side, xi_line)
        static LINALG::Matrix<dimEdge + dimSide, 1, floatType> dx_;
      };


      // Storage class of ComputeIntersection data members to inherit from
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double>
      struct ComputeIntersectionNoStaticMembers
      {
        LINALG::Matrix<numNodesSide, 1, floatType> sideFunct_;
        LINALG::Matrix<numNodesEdge, 1, floatType> edgeFunct_;
        LINALG::Matrix<dimSide, numNodesSide, floatType> sideDeriv1_;
        LINALG::Matrix<dimEdge, numNodesEdge, floatType> edgeDeriv1_;

        LINALG::Matrix<dimEdge + dimSide, dimEdge + dimSide, floatType> A_;
        LINALG::Matrix<probDim, dimEdge + dimSide, floatType> B_;
        LINALG::Matrix<dimEdge + dimSide, 1, floatType> b_;
        LINALG::Matrix<probDim, 1, floatType> c_;

        LINALG::Matrix<dimEdge + dimSide, 1, floatType> dx_;
      };

      /*--------------------------------------------------------------------------*/
      /** \brief strategy class for intersection of side and edge/line
       *
       *  inheritance diagram:
       *
       *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
       *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy
       *      ^^^^^^^^^^^^^^^^^^^^^^^^^^^  */
      template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
          typename floatType = double,
          template <bool, unsigned, DRT::Element::DiscretizationType,
              DRT::Element::DiscretizationType, unsigned, unsigned, unsigned, unsigned, typename>
          class MemberStoragePolicy = ComputeIntersectionStaticMembers>
      class ComputeIntersectionStrategy
          : public EmptyNewtonStrategy,
            MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge,
                numNodesSide, floatType>
      {
       public:
        /// constructor
        ComputeIntersectionStrategy(
            LINALG::Matrix<dimEdge + dimSide, 1, floatType>& xsi, bool checklimits)
            : xsi_(xsi),
              xyze_side_(NULL),
              xyze_edge_(NULL),
              residual_(0.0),
              off_count_(0),
              check_limits_(checklimits),
              tol_(0.0)
        {
        }

        void Setup(const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            const LINALG::Matrix<probDim, numNodesEdge, floatType>& xyze_edge)
        {
          xyze_side_ = &xyze_side;
          xyze_edge_ = &xyze_edge;

          if (debug)
          {
            std::cout << "\n\n === ComputeIntersection ===\n";
            std::cout << "--- Setup()\n";
            std::cout << "  Edge = " << DRT::DistypeToString(edgeType) << "\n";
            std::cout << "  Side = " << DRT::DistypeToString(sideType) << "\n";
            std::cout << "  xyze_side = " << std::setprecision(15) << xyze_side;
            std::cout << "  xyze_edge = " << std::setprecision(15) << xyze_edge;
          }
        }

        std::pair<bool, floatType> ConditionNumber()
        {
          LINALG::Matrix<dimEdge + dimSide, dimEdge + dimSide, floatType> A_inv;
          floatType det = A_.Determinant();
          if (KERNEL::closeToZero(det))
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            std::cout << " WARNING: Condition number is equal to infinity in the problem. "
                         "Continuing increase in precision "
                      << std::endl;
#endif
            return std::make_pair(false, -1.0);
          }
          A_inv.Invert(A_);
          floatType cond = A_inv.Norm2() * A_.Norm2();
          return std::make_pair(true, cond);
        }

        void GetLocalToleranceEdge(const floatType& global_tolerance,
            LINALG::Matrix<dimEdge, 1, floatType>& scaled_tolerance)
        {
          LINALG::Matrix<probDim, 1, floatType> real_tolerance;
          real_tolerance =
              global_tolerance;  // set all values of the vector to be equal to the global tolerance
          // extract part of the side from the Jacobian
          LINALG::Matrix<probDim, dimEdge, floatType> A;
          if (dimEdge + dimSide == probDim)
          {
            for (unsigned int row = 0; row < probDim; ++row)
            {
              for (unsigned int col = 0; col < dimEdge; ++col) A(row, col) = A_(row, dimSide + col);
            }
          }
          else
          {
            for (unsigned int row = 0; row < probDim; ++row)
            {
              for (unsigned int col = 0; col < dimEdge; ++col) A(row, col) = B_(row, dimSide + col);
            }
          }

          LINALG::Matrix<dimEdge, dimEdge, floatType> TN_inv;
          LINALG::Matrix<dimEdge, probDim, floatType> aux;
          TN_inv.MultiplyTN(A, A);

          TN_inv.Invert();
          aux.MultiplyNT(TN_inv, A);
          scaled_tolerance.Multiply(aux, real_tolerance);
        }

        void GetLocalToleranceSide(const floatType& global_tolerance,
            LINALG::Matrix<dimSide, 1, floatType>& scaled_tolerance)
        {
          LINALG::Matrix<probDim, 1, floatType> real_tolerance;
          real_tolerance =
              global_tolerance;  // set all values of the vector to be equal to the global tolerance
          // extract part of the side from the Jacobian
          LINALG::Matrix<probDim, dimSide, floatType> A;
          if (dimEdge + dimSide == probDim)
          {
            for (unsigned int row = 0; row < probDim; ++row)
            {
              for (unsigned int col = 0; col < dimSide; ++col) A(row, col) = A_(row, col);
            }
          }
          else
          {
            for (unsigned int row = 0; row < probDim; ++row)
            {
              for (unsigned int col = 0; col < dimSide; ++col) A(row, col) = B_(row, col);
            }
          }

          LINALG::Matrix<dimSide, dimSide, floatType> TN_inv;
          LINALG::Matrix<dimSide, probDim, floatType> aux;
          TN_inv.MultiplyTN(A, A);

          TN_inv.Invert();
          aux.MultiplyNT(TN_inv, A);
          scaled_tolerance.Multiply(aux, real_tolerance);
        }

        const LINALG::Matrix<dimEdge + dimSide, 1, floatType>& LocalCoordinates() const
        {
          return xsi_;
        }

        floatType DistanceBetween() const { return c_.Norm2(); }

        void SetupSolve()
        {
          xsi_ = 0.0;
          dx_ = 0.0;
          c_ = 0.0;


          off_count_ = 0;
          residual_ = 0.0;

          // evaluate initial rhs value
          const LINALG::Matrix<dimSide, 1, floatType> xsi_side(xsi_.A(), true);
          const LINALG::Matrix<dimEdge, 1, floatType> xsi_edge(xsi_.A() + dimSide, true);
          IntersectionRHS(xsi_edge, xsi_side, *xyze_edge_, *xyze_side_, c_);
          tol_ = AdaptiveCombinedNewtonTolerance(*xyze_side_, *xyze_edge_, c_, xsi_(0, 0));
        }

        void SetupStep(int iter)
        {
          // build linear system of equations
          const LINALG::Matrix<dimSide, 1, floatType> xsi_side(xsi_.A(), true);
          const LINALG::Matrix<dimEdge, 1, floatType> xsi_edge(xsi_.A() + dimSide, true);
          if (iter > 0)
          {
            IntersectionRHS(xsi_edge, xsi_side, *xyze_edge_, *xyze_side_, c_);
          }
          IntersectionSystem(xsi_edge, xsi_side, *xyze_edge_, *xyze_side_, c_, A_, B_, b_);
        }

        /// compute the right-hand-side
        void IntersectionRHS(const LINALG::Matrix<dimEdge, 1, floatType>& xsi_edge,
            const LINALG::Matrix<dimSide, 1, floatType>& xsi_side,
            const LINALG::Matrix<probDim, numNodesEdge, floatType>& xyze_edge,
            const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            LINALG::Matrix<probDim, 1, floatType>& c)
        {
          DRT::UTILS::shape_function<sideType>(xsi_side, sideFunct_);
          DRT::UTILS::shape_function<edgeType>(xsi_edge, edgeFunct_);

          c.Multiply(xyze_edge, edgeFunct_);
          c.Multiply(-1.0, xyze_side, sideFunct_, 1.0);
        }


        enum NewtonStatus TestConverged(int iter)
        {
          residual_ = b_.Norm2();


          if (debug)
          {
            std::cout << "--- TestConverged\n";
            floatType incr = dx_.Norm2();
            std::cout << "  residual = " << std::setprecision(15) << residual_
                      << "  incr = " << std::setprecision(15) << incr << "\n";
          }

          if (c_.Norm2() < tol_ && iter > 0)
          {
            tol_ *= std::sqrt(static_cast<double>(probDim));
            return converged;
          }

          /* If the norm of c_ does not but the norm of b_ does fulfill the criterion,
           * the least square scheme failed! */
          if (residual_ < tol_ && iter > 0)
          {
            return failed;
          }

          return unconverged;
        }

        bool LinearSolve(int iter)
        {
          dx_ = 0.0;
          floatType det = LINALG::gaussElimination<true, dimSide + dimEdge, floatType>(A_, b_, dx_);

          if (debug)
          {
            std::cout << "--- LinearSolve\n";
            std::cout << "  A_ = " << A_;
            std::cout << "  b_ = " << b_;
            std::cout << "  det = " << std::setprecision(15) << det << "\n";
            std::cout << "  dxi_ = " << dx_;
          }

          return (det != 0.0);
        }

        bool Update(int iter)
        {
          xsi_.Update(1.0, dx_, 1.0);


#ifndef CLN_CALC
          if (debug)
          {
            std::cout << "--- Update\n";
            const LINALG::Matrix<dimSide, 1> xsi_side(xsi_.A(), true);
            const LINALG::Matrix<dimEdge, 1> xsi_edge(xsi_.A() + dimSide, true);
            std::cout << "  off_count=" << off_count_ << "\n  Side within limits = "
                      << (WithinLimits<sideType>(xsi_side) ? "TRUE" : "FALSE")
                      << "\n  Edge within limits = "
                      << (WithinLimits<edgeType>(xsi_edge) ? "TRUE" : "FALSE") << "\n";
            std::cout << "  xsi = " << xsi_;
            std::cout << "  xsi_side = " << xsi_side;
            std::cout << "  xsi_edge = " << xsi_edge;

            LINALG::Matrix<numNodesSide, 1> sideFunct;
            LINALG::Matrix<numNodesEdge, 1> edgeFunct;

            DRT::UTILS::shape_function<sideType>(xsi_side, sideFunct);
            DRT::UTILS::shape_function<edgeType>(xsi_edge, edgeFunct);

            LINALG::Matrix<probDim, 1> x1;
            LINALG::Matrix<probDim, 1> x2;

            x1.Multiply(*xyze_side_, sideFunct);
            x2.Multiply(*xyze_edge_, edgeFunct);

            std::cout << "  x_side =" << x1;
            std::cout << "  x_edge =" << x2;

            x2.Update(1, x1, -1);
            std::cout << "  (x_side - x_edge) =" << x2;
          }
#endif

          return true;
        }

        bool NewtonFailed()
        {
#ifdef DEBUG_CUTKERNEL_OUTPUT
          // If it is not edge-edge interseciton we explicitely notify, that Newton's method failed
          if (probDim == dimEdge + dimSide)
          {
            std::cout << "Current type is " << typeid(floatType).name() << std::endl;
            std::cout << "Newton method failed" << std::endl;
            std::cout << "Residual is" << residual_ << std::endl;
            std::cout << "Tolerance is" << tol_ << std::endl;
            std::cout << "Precision is" << ClnWrapper::precision_ << "decimal points" << std::endl;
            std::cout << "Local coordinates are " << xsi_ << std::endl;
          }
#endif

          tol_ = c_.Norm2() * std::sqrt(3.0);

#ifdef DEBUG_CUTKERNEL_OUTPUT
          if (probDim == dimEdge + dimSide)
          {
            std::stringstream str;
            str << "Newton scheme did not converge:\n  "
#ifndef CLN_CALC
                << (*xyze_side_) << (*xyze_edge_)
#endif
                ;
#ifndef CLN_CALC
            DumpDoubles(str, xyze_side_->A(), xyze_side_->M() * xyze_side_->N());
            str << "\n";
            DumpDoubles(str, xyze_edge_->A(), xyze_edge_->M() * xyze_edge_->N());
#endif


            std::stringstream f_str;
            f_str << ".NewtonFailed_intersection.pos";
            std::string filename(GEO::CUT::OUTPUT::GenerateGmshOutputFilename(f_str.str()));
            std::ofstream file(filename.c_str());
            WritetoGmsh(file);
          }
#endif
          return false;
        }

        floatType GetTolerance() const { return tol_; }

        void WritetoGmsh(std::ofstream& file)
        {
          file.precision(32);  // higher precision!
          char elementType = '\0';
          switch (numNodesSide)
          {
            case 3:
              elementType = 'T';
              break;
            case 4:
              elementType = 'Q';
              break;
            default:
              dserror(
                  "unsupported element type ( % s ) in WritetoGmsh."
                  " Please feel free to extend the functionality if necessary.",
                  DRT::DistypeToString(sideType).c_str());
          }

          file << "View \""
               << "Side"
               << "\" {\n";
          {
            file << "S" << elementType << "(";
            for (uint i = 0; i < numNodesSide; ++i)
            {
              if (i != 0) file << ",";
              file << (*xyze_side_)(0, i) << "," << (*xyze_side_)(1, i) << ","
                   << (*xyze_side_)(2, i);
            }
            file << "){";
            for (uint i = 0; i < numNodesSide; ++i)
            {
              if (i != 0) file << ",";
              file << "1";
            }
            file << "};\n";
          }
          file << "};\n";

          file << "View \""
               << "Line"
               << "\" {\n";
          for (unsigned iline = 0; iline < numNodesEdge - 1; ++iline)
          {
            file << "SL (";
            file << (*xyze_edge_)(0, iline) << "," << (*xyze_edge_)(1, iline) << ","
                 << (*xyze_edge_)(2, iline) << ",";
            file << (*xyze_edge_)(0, iline + 1) << "," << (*xyze_edge_)(1, iline + 1) << ","
                 << (*xyze_edge_)(2, iline + 1);
            file << "){";
            file << 1 << ",";
            file << 1;
            file << "};\n";
          }
          file << "};\n";


          // calculate the points
          LINALG::Matrix<numNodesSide, 1, floatType> surface;
          LINALG::Matrix<numNodesEdge, 1, floatType> line;

          DRT::UTILS::shape_function_2D(surface, xsi_(0), xsi_(1), sideType);
          DRT::UTILS::shape_function_1D(line, xsi_(2), edgeType);

          LINALG::Matrix<probDim, 1, floatType> x1;
          LINALG::Matrix<probDim, 1, floatType> x2;

          x1.Multiply(*xyze_side_, surface);
          x2.Multiply(*xyze_edge_, line);


          file << "View \""
               << "SidePoint"
               << "\" {\n";
          {
            file << "SP (";
            file << x1(0, 0) << "," << x1(1, 0) << "," << x1(2, 0);
            file << "){";
            file << 1;
            file << "};\n";
            file << "};\n";
          }
          file << "View \""
               << "EdgePoint"
               << "\" {\n";
          {
            file << "SP (";
            file << x2(0, 0) << "," << x2(1, 0) << "," << x2(2, 0);
            file << "){";
            file << 2;
            file << "};\n";
            file << "};\n";
          }
        }

       private:
        /** \brief Build the linear system of equations
         *
         *  Right-hand-side:
         *  \f[
         *    b = x(\xi) - x(\tilde{\xi}),
         *  \f]
         *  where \f$\xi \in \mathbb{R}^{dimEgde}\f$ is the parameter space
         *  coordinate of the edge-element and \f$ \tilde{\xi} \in \mathbb{R}^{dimSide} \f$
         *  are the parameter space coordinates of the side element.
         *
         *  The matrix \f$ A \f$ holds the consistent linearization. If the problem
         *  dimension is larger than the accumulated dimensions of the side and edge
         *  elements, a least square approach will be used.
         *
         *  \author hiermeier
         *  \date 08/16 */
        void IntersectionSystem(const LINALG::Matrix<dimEdge, 1, floatType>& xsi_edge,
            const LINALG::Matrix<dimSide, 1, floatType>& xsi_side,
            const LINALG::Matrix<probDim, numNodesEdge, floatType>& xyze_edge,
            const LINALG::Matrix<probDim, numNodesSide, floatType>& xyze_side,
            LINALG::Matrix<probDim, 1, floatType>& c,
            LINALG::Matrix<dimEdge + dimSide, dimEdge + dimSide, floatType>& A,
            LINALG::Matrix<probDim, dimEdge + dimSide, floatType>& B,
            LINALG::Matrix<dimEdge + dimSide, 1, floatType>& b)
        {
          // --- compute the linearization
          DRT::UTILS::shape_function_deriv1<sideType>(xsi_side, sideDeriv1_);
          DRT::UTILS::shape_function_deriv1<edgeType>(xsi_edge, edgeDeriv1_);

          A = 0.0;
          B = 0.0;

          // linearization w.r.t. the parameter space coordinate(s) of the side
          for (unsigned inode = 0; inode < numNodesSide; ++inode)
            for (unsigned row = 0; row < probDim; ++row)
              for (unsigned col = 0; col < dimSide; ++col)
                B(row, col) += xyze_side(row, inode) * sideDeriv1_(col, inode);

          // linearization w.r.t. the parameter space coordinate(s) of the edge
          for (unsigned inode = 0; inode < numNodesEdge; ++inode)
            for (unsigned row = 0; row < probDim; ++row)
              for (unsigned col = 0; col < dimEdge; ++col)
                B(row, dimSide + col) -= xyze_edge(row, inode) * edgeDeriv1_(col, inode);

          // default case
          if (probDim == dimSide + dimEdge)
          {
            A.SetView(B.A());
            b.SetView(c.A());
          }
          /* Actually there is only one relevant case, which comes into my mind
           * and that is the intersection of two lines in 3 dimensions. Anyway
           * we use a more general if-statement here.
           *
           * If we have to handle such a configuration we use a Least-Squares
           * approach. */
          else if (probDim > dimSide + dimEdge)
          {
            A.MultiplyTN(B, B);
            b.MultiplyTN(B, c);
          }
          else
            dserror(
                "The problem dimension is smaller than the combination of"
                "the side and edge element dimensions. This case is currently "
                "unsupported!");
        }

       private:
        /// local coordinates (xi1_side, xi2_side, xi_line)
        LINALG::Matrix<dimEdge + dimSide, 1, floatType>& xsi_;

        const LINALG::Matrix<probDim, numNodesSide, floatType>* xyze_side_;

        const LINALG::Matrix<probDim, numNodesEdge, floatType>* xyze_edge_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::sideFunct_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::edgeFunct_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::sideDeriv1_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::edgeDeriv1_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::A_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::B_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::b_;

        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::c_;
        /// increment in local coordinates d(xi1_side, xi2_side, xi_line)
        using MemberStoragePolicy<debug, probDim, edgeType, sideType, dimEdge, dimSide,
            numNodesEdge, numNodesSide, floatType>::dx_;

        /// residual norm
        floatType residual_;

        /** count how many iterations during the Newton the local coordinates
         *  are not within limits! Stop the newton scheme if the iteration is
         *  too often outside the limits */
        int off_count_;

        /** shall we check the limits? ( seems out-dated, since we donot check
         *  limits by default ) */
        bool check_limits_;
        // basic tolerance!!!
        floatType tol_;
      };  // class ComputeIntersectionStrategy


#ifdef CLN_CALC

      /*--------------------------------------------------------------------------*/
      /** \brief Strategy class for intersection of side and line/edge for arbitrary precision
       *
       *  inheritance diagram:
       *
       *  ComputeIntersection --> ComputeIntersectionAdaptivePrecision --> NewtonSolve
       *                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy  */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class ComputeIntersectionAdaptivePrecision : Strategy
#ifdef DEBUG_MEMORY_ALLOCATION
          ,
                                                   AdaptiveKernelSharedData
#endif
      {
       public:
        //! constructor
        ComputeIntersectionAdaptivePrecision(
            LINALG::Matrix<dimEdge + dimSide, 1>& xsi, bool checklimits)
            : Strategy(clnxsi_, checklimits), xsi_(xsi), clnxsi_(true), cond_infinity_(false)
        {
        }

        bool operator()(const LINALG::Matrix<probDim, numNodesSide>& xyze_side,
            const LINALG::Matrix<probDim, numNodesEdge>& xyze_edge)
        {
          if (!(((edgeType == DRT::Element::tet4) || (edgeType == DRT::Element::pyramid5) ||
                    (edgeType == DRT::Element::hex8) || (edgeType == DRT::Element::quad4) ||
                    (edgeType == DRT::Element::hex20) || (edgeType == DRT::Element::tri3) ||
                    (edgeType == DRT::Element::hex8) || (edgeType == DRT::Element::line2)) &&
                  ((sideType == DRT::Element::hex8) || (sideType == DRT::Element::quad4) ||
                      (sideType == DRT::Element::tri3) || (sideType == DRT::Element::line2) ||
                      (sideType == DRT::Element::hex20) || (sideType == DRT::Element::tet4) ||
                      (sideType == DRT::Element::pyramid5) || (sideType == DRT::Element::hex8))))
          {
            dserror(
                "This type of element (%s)  is not tested for the CLN calculation. You are welcome "
                "to edit ../drt_fem_general/drt_utils_fem_shapefunctions.H file to fix it. Just "
                "change all the integers occuring there double for CLN to work.",
                DRT::DistypeToString(sideType).c_str());
          }

          // Converting values for the calculation on the CLN
          int prec = CLN_START_PRECISION;
          int iter = 0;
          bool conv;
          ClnWrapper err;
          ClnWrapper cond_number;
          cond_infinity_ = false;

#ifdef CUSTOM_MEMORY_ALLOCATOR
          if (all_distance_done_once_ and all_intersections_done_once_ and
              all_position_done_once_ and (!custom_allocator_run_))
          {
#if DEBUG_MEMORY_ALLOCATION
            GEO::CUT::MemorySingleton::getInstance().ReportAllocated();
            ReportIntersectionAllocated();
            ReportPositionAllocated();
            ReportDistanceAllocated();
            ReportTotalAllocated();
            GEO::CUT::MemorySingleton::getInstance().SetState(1, memory_allocations_);
#else
            GEO::CUT::MemorySingleton::getInstance().SwitchState();
#endif
            custom_allocator_run_ = true;
          }
#if DEBUG_MEMORY_ALLOCATION
          if (first_run_) std::cout << "In cut intersection statistics is" << std::endl;
#endif
#endif

          do
          {
            ClnWrapper::SetPrecision(prec);
#ifdef CUSTOM_MEMORY_ALLOCATOR
#if DEBUG_MEMORY_ALLOCATION
            if (custom_allocator_run_)
#endif
              GEO::CUT::MemorySingleton::getInstance().GetMemoryPoolAllocator().SetCurrent(
                  cln_byte_size_[iter]);
#endif
#if DEBUG_MEMORY_ALLOCATION

            if (first_run_)
            {
              if (GEO::CUT::MemorySingleton::getInstance().IsRecording())
              {
                UpdateMemoryAllocations(
                    GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern());
                GEO::CUT::MemorySingleton::getInstance().StopRecord();
              }
              GEO::CUT::MemorySingleton::getInstance().StartRecord();
            }

#endif
            ConvDoubleCLN(xyze_side, clnxyze_side_, prec);
            ConvDoubleCLN(xyze_edge, clnxyze_edge_, prec);
            this->Setup(clnxyze_side_, clnxyze_edge_);

            conv = this->Solve();

#if CUT_DEVELOP
            // safety check for precision loss
            cln::float_format_t prec_beg = cln::float_format(clnxyze_side_(0, 0));
            for (unsigned int i = 0; i < dimEdge + dimSide; ++i)
            {
              cln::float_format_t prec_end = cln::float_format(clnxsi_(i, 0));
              if (prec_beg != prec_end)
              {
                dserror(
                    "There is a loss of cln-precision during computation of Intersection. "
                    "Something is wrong with the conversion");
              }
            }
#endif

            std::pair<bool, ClnWrapper> cond_pair = this->ConditionNumber();
            cond_number = cond_pair.second;
            if (not cond_pair.first)
            {
              cond_infinity_ = true;
            }

            err = cond_number *
                  ComputeError(xyze_side, xyze_edge, this->LocalCoordinates(), CLN_REFERENCE_PREC);

            // compute scaled value of the error with respect to global coordinate
#if DEBUG_MEMORY_ALLOCATION
            UpdateMemoryUsage(prec, iter);
#endif
            prec += CLN_INCREMENT_STEP;
            iter++;
            // increase the precision
          } while (((!conv) || (err > CLN_LIMIT_ERROR) || (cond_infinity_)) &&
                   (prec < CLN_LIMIT_PREC) && (iter < CLN_LIMIT_ITER));


          if ((not conv) && (err < CLN_LIMIT_ERROR) && (err > 0) &&
              ((dimEdge + dimSide) == probDim))
            conv = true;
          else if (not conv)
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            if (probDim == dimSide + dimEdge)
              std::cout << "Newton method failed and the residual is " << err << " / "
                        << CLN_LIMIT_ERROR << std::endl;
#endif
          }

          if (cond_infinity_)
          {
#ifdef DEBUG_CUTKERNEL_OUTPUT
            std::cout << "Condition number of this problem is equal to infinity" << std::endl;
#endif
          }

          GetTopologyInformation();
          // NOTE: Might be not needed later
          FixCornerCase();
#ifdef CUSTOM_MEMORY_ALLOCATOR
          // We need it because on the first run, we try to do the try run with allocation of the
          // required non-reusable doubles, e.g. double->cln compile time constant and other global
          // variables. For compile time constants it does not really matter, but for the global
          // variables we need somehow to allocate them in the not reusable container
          if (first_run_)
          {
            const int num_inter = 3;
            if ((++GetIntersectionCount()) == num_inter)
            {
              all_intersections_done_once_ = true;
              std::cout << "All possible template intersection are done!" << std::endl;
            }
            first_run_ = false;
          }
#endif
          ConvClnDouble(clnxsi_, xsi_);
          clnxsi_ = 0.0;  // resetting clnxsi_
          return conv;
        }

        const LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper>& LocalCoordinates() const
        {
          return Strategy::LocalCoordinates();
        }

        std::pair<bool, ClnWrapper> ConditionNumber() { return Strategy::ConditionNumber(); }

        bool SurfaceWithinLimits(double tolerance = REFERENCETOL) const
        {
          return WithinLimits<sideType>(LocalCoordinates(), tolerance);
        }

        bool LineWithinLimits(double tolerance = REFERENCETOL) const
        {
          const LINALG::Matrix<dimEdge, 1> xsi_line(LocalCoordinates().A() + dimSide, true);
          return WithinLimits<edgeType>(xsi_line, tolerance);
        }

        ClnWrapper GetTolerance() const { return Strategy::GetTolerance(); }

        void WritetoGmsh(std::ofstream& file) { Strategy::WritetoGmsh(file); }

        bool IsConditionInfinity() { return cond_infinity_; }

        PointOnSurfaceLoc GetSideLocation() { return sideLocation_; }

        PointOnSurfaceLoc GetEdgeLocation() { return edgeLocation_; }

        const std::vector<int>& GetTouchedSideEdges() { return touched_edges_ids_; }

        ClnWrapper DistanceBetween() const { return Strategy::DistanceBetween(); }

       private:
        //  Converts tolerance in the local coordinates and computes topological information, such
        //  as touched_edges and  location on the surface
        void GetTopologyInformation()
        {
          // tolerance of getting nearby edges
          LINALG::Matrix<dimSide, 1, ClnWrapper> scaled_tolerance_side_touched_edges;
          // tolerances of determining inside/outside surface limits
          LINALG::Matrix<dimSide, 1, ClnWrapper> scaled_tolerance_side;
          LINALG::Matrix<dimEdge, 1, ClnWrapper> scaled_tolerance_edge;
          // compute location
          this->GetLocalToleranceSide(INSIDE_OUTSIDE_TOLERANCE, scaled_tolerance_side);
          this->GetLocalToleranceSide(
              SIDE_DETECTION_TOLERANCE, scaled_tolerance_side_touched_edges);
          this->GetLocalToleranceEdge(INSIDE_OUTSIDE_TOLERANCE, scaled_tolerance_edge);
          // compute inside-outside relation
          if (WithinLimits<sideType>(clnxsi_, scaled_tolerance_side))
            sideLocation_ = PointOnSurfaceLoc(true, true);
          else
            sideLocation_ = PointOnSurfaceLoc(false, true);
          const LINALG::Matrix<dimEdge, 1, ClnWrapper> clnxsi_line(clnxsi_.A() + dimSide, true);
          if (WithinLimits<edgeType>(clnxsi_line, scaled_tolerance_edge))
            edgeLocation_ = PointOnSurfaceLoc(true, true);
          else
            edgeLocation_ = PointOnSurfaceLoc(false, true);
          // NOTE: In compute intersection we don't actually use touched_edges now. left for the
          // future
          touched_edges_ids_.clear();
          if (sideLocation_.WithinSide())
          {
            GetEdgesAt<sideType>(clnxsi_, touched_edges_ids_, scaled_tolerance_side_touched_edges);
          }
        }

        // Evaluate difference between global coordinates in the intersection when computed based on
        // side and based on edge using conversion to very high reference precision
        ClnWrapper ComputeError(
            const LINALG::Matrix<probDim,
                DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>& refside_xyz,
            const LINALG::Matrix<probDim,
                DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement>& refedge_xyz,
            const LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper>& loc_calc, int prec)
        {
          unsigned int prev_prec = ClnWrapper::GetPrecision();
          ClnWrapper::SetPrecision(prec);

          // Converting input arrays  to higher precision floating points
          LINALG::Matrix<probDim, DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement,
              ClnWrapper>
              cln_refside_xyz;
          LINALG::Matrix<probDim, DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
              ClnWrapper>
              cln_refedge_xyz;
          LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper> cln_loc_calc;

          ConvDoubleCLN(refside_xyz, cln_refside_xyz, prec);
          ConvDoubleCLN(refedge_xyz, cln_refedge_xyz, prec);
          ConvDoubleCLN(loc_calc, cln_loc_calc, prec);

          // Calculating interpolation from the shapefunction of the side
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement, 1,
              ClnWrapper>
              sideFunct;
          LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper> cln_glob_calc_side;
          LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper> cln_loc_calc_side(
              cln_loc_calc.A(), true);
          DRT::UTILS::shape_function<sideType>(cln_loc_calc_side, sideFunct);

          for (unsigned int inode = 0;
               inode < DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement; ++inode)
          {
            for (unsigned int isd = 0; isd < (dimEdge + dimSide); ++isd)
            {
              cln_glob_calc_side(isd) += cln_refside_xyz(isd, inode) * sideFunct(inode);
            }
          }

          // Calculating interpolation from the shapefunction of the edge
          LINALG::Matrix<probDim, 1, ClnWrapper> cln_glob_calc_edge;
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement, 1,
              ClnWrapper>
              edgeFunct;
          LINALG::Matrix<probDim, 1, ClnWrapper> cln_loc_calc_edge(
              cln_loc_calc.A() + DRT::UTILS::DisTypeToDim<sideType>::dim, true);

          DRT::UTILS::shape_function<edgeType>(cln_loc_calc_edge, edgeFunct);
          for (unsigned int inode = 0;
               inode < DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement; ++inode)
          {
            for (unsigned int isd = 0; isd < (dimEdge + dimSide); ++isd)
            {
              cln_glob_calc_edge(isd) += cln_refedge_xyz(isd, inode) * edgeFunct(inode);
            }
          }

          // Evaluating error
          LINALG::Matrix<probDim, 1, ClnWrapper>
              diff_vec;  // vector of difference between interpolation
          for (unsigned int i = 0; i < (dimEdge + dimSide); ++i)
            diff_vec(i) = cln_glob_calc_edge(i) - cln_glob_calc_side(i);

          // resetting precision to previous
          ClnWrapper::SetPrecision(prev_prec);
          return diff_vec.Norm2();
        }

        /// Try to fix the case that is close to the corner and outside - in that case point
        /// is close to corner and outside might be touching two edges but difference between it
        /// and corner point > tolerance - hence this must be corrected
        void FixCornerCase()
        {
          // if touches more than one edge
          if (touched_edges_ids_.size() >= 2)
          {
            if (sideLocation_.WithinSide())
            {
              LINALG::Matrix<dimSide, 1, ClnWrapper> zero_tolerance_side;
              // this means we are on the outer boundary of the tolerance
              if (not WithinLimits<sideType>(clnxsi_, zero_tolerance_side))
              {
                // we transform our point to global tolerance and check if it is close enough to the
                // nodal point of the side, that is common point of the touching edges
                LINALG::Matrix<numNodesSide, 1, ClnWrapper> sideFunct;
                LINALG::Matrix<probDim, 1, ClnWrapper> xsi_global;
                xsi_global = 0.0;
                const LINALG::Matrix<dimSide, 1> xsi_side(xsi_.A(), true);
                DRT::UTILS::shape_function<sideType>(xsi_side, sideFunct);
                // taking latest cln based coordinate
                const LINALG::Matrix<probDim, numNodesSide, ClnWrapper>& clnxyze_side_latest =
                    clnxyze_side_;
                for (unsigned int inode = 0; inode < numNodesSide; ++inode)
                  for (unsigned int isd = 0; isd < probDim; ++isd)
                    xsi_global(isd) += clnxyze_side_latest(isd, inode) * sideFunct(inode);

                // lame method for now. find minimum distance and check if it is closer than 1e-14
                LINALG::Matrix<probDim, 1, ClnWrapper> dist;
                ClnWrapper min_dist;
                for (unsigned int inode = 0; inode < numNodesSide; ++inode)
                {
                  // acess raw responsible for that node
                  LINALG::Matrix<probDim, 1, ClnWrapper> coord(
                      clnxyze_side_latest.A() + inode * probDim, true);
                  dist.Update(1.0, coord, -1.0, xsi_global);
                  ClnWrapper tmp_dist = dist.Norm2();
                  if ((inode == 0) or (tmp_dist < min_dist)) min_dist = tmp_dist;
                }

                // this should not be possible now
                if (min_dist > NON_TOPOLOGICAL_TOLERANCE)
                {
                  touched_edges_ids_.clear();
                  sideLocation_ = PointOnSurfaceLoc(false, true);
                  std::stringstream msg;
                  msg << "NOTICE: Intersection point is close to 2 side edges, but too far from "
                         "the corner point"
                      << "\n Distance is " << min_dist;
                  dserror(msg.str());
                }
              }
            }
            else
              dserror("This should not be possible!");
          }
        }

        // Update global maximum number of allocations ( number of containers of particular byte
        // size) if allocations on this precisions were more. Also set what was the most frequent
        // size of allocation for the current iteration
        void UpdateMemoryUsage(int prec, int iter)
        {
#if DEBUG_MEMORY_ALLOCATION
          // if first run need also to run until the end
          if (first_run_)
          {
            GEO::CUT::MemorySingleton::getInstance().StopRecord();

            boost::unordered_map<size_t, int>& allocation_map =
                GEO::CUT::MemorySingleton::getInstance().GetMemoryPattern();

            int max_num = 0;
            int size_max_num = 0;
            // update the global number of allocations to the maximum
            for (boost::unordered_map<size_t, int>::iterator it = allocation_map.begin();
                 it != allocation_map.end(); ++it)
            {
              size_t size = it->first;
              if (memory_allocations_intersection_[size] < allocation_map[size])
              {
                memory_allocations_intersection_[size] = allocation_map[size];
              }
              if (allocation_map[size] > max_num)
              {
                size_max_num = size;
                max_num = allocation_map[size];
              }
            }

            if (size_max_num != 0)
            {
              std::cout << "Size for the precision " << prec << " is " << size_max_num << std::endl;
              std::cout << "It is allocated " << max_num << " times " << std::endl;
              cln_sizes_[iter] = size_max_num;
            }
            else
              dserror("This should not be possible!");

            GEO::CUT::MemorySingleton::getInstance().ResetAllocated();
          }
#endif
        }

        static std::vector<int> touched_edges_ids_;

        LINALG::Matrix<probDim, numNodesSide, ClnWrapper> clnxyze_side_;

        LINALG::Matrix<probDim, numNodesEdge, ClnWrapper> clnxyze_edge_;

        // holds reference to original value
        LINALG::Matrix<dimEdge + dimSide, 1>& xsi_;

        LINALG::Matrix<dimEdge + dimSide, 1, ClnWrapper> clnxsi_;

        static bool first_run_;

        // if condition number is infinity or not
        bool cond_infinity_;

        // determines location of this point with respect to the edge that cuts it
        PointOnSurfaceLoc sideLocation_;

        // determines location of this point with respect to the side that cuts it
        PointOnSurfaceLoc edgeLocation_;

        // sizes of cln variables with different preciion
        static size_t cln_sizes_[CLN_LIMIT_ITER];
      };  // class ComputeIntersectionAdaptivePrecision

#endif

      /*--------------------------------------------------------------------------*/
      /** \brief generic strategy class for intersection of side and line/edge
       *
       *  inheritance diagram:
       *
       *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
       *                          ^^^^^^^^^^^^^^^^^^^^^^^^^^
       *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy  */
      template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType, bool computeCln = false,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class GenericComputeIntersection : Strategy
      {
       public:
        //! constructor
        GenericComputeIntersection(LINALG::Matrix<dimEdge + dimSide, 1>& xsi, bool checklimits)
            : Strategy(xsi, checklimits),
              cond_infinity_(false),
              distance_between_(0.0),
              xsi_(xsi)
#ifdef CLN_CALC
              ,
              checklimits_(checklimits)
#endif
        {
        }

        bool operator()(const LINALG::Matrix<probDim, numNodesSide>& xyze_side,
            const LINALG::Matrix<probDim, numNodesEdge>& xyze_edge)
        {
          bool conv;
          this->Setup(xyze_side, xyze_edge);
          conv = this->Solve();
          GetTopologyInformation();

#ifdef CLN_CALC
          if (computeCln)
          {
#if DOUBLE_PLUS_CLN_COMPUTE
            if ((ComputeError(xyze_side, xyze_edge, xsi_) > DOUBLE_LIMIT_ERROR) or (not conv) or
                (cond_infinity_))
            {
//             CutKernelStatistics::getCutKernelStatistics().ClnIntersectionCounter();
#endif
              {
                ComputeIntersectionAdaptivePrecision<
                    NewtonSolve<ComputeIntersectionStrategy<false, probDim, edgeType, sideType,
                                    dimEdge, dimSide, numNodesEdge, numNodesSide, ClnWrapper,
                                    ComputeIntersectionNoStaticMembers>,
                        dimEdge + dimSide>,
                    probDim, edgeType, sideType>
                    clncalc(xsi_, checklimits_);
                bool clnsolved = clncalc(xyze_side, xyze_edge);
                conv = clnsolved;
                edgeLocation_ = clncalc.GetEdgeLocation();
                sideLocation_ = clncalc.GetSideLocation();
                cond_infinity_ = clncalc.IsConditionInfinity();
                touched_edges_ids_.clear();
                const std::vector<int>& touched_edges_ids_cln = clncalc.GetTouchedSideEdges();
                touched_edges_ids_ = touched_edges_ids_cln;

                if (probDim > dimEdge + dimSide)
                  distance_between_ = cln::double_approx(clncalc.DistanceBetween());

                ClnWrapper::ResetPrecision();
              }
#ifdef CUSTOM_MEMORY_ALLOCATOR
              GEO::CUT::MemorySingleton::getInstance().Finalize();
#endif
#if DOUBLE_PLUS_CLN_COMPUTE
            }
            else
            {
              //             CutKernelStatistics::getCutKernelStatistics().DoubleIntersectionCounter();
            }
#endif
          }
#endif
          return conv;
        }

        const LINALG::Matrix<dimEdge + dimSide, 1>& LocalCoordinates() const
        {
          return Strategy::LocalCoordinates();
        }

        std::pair<bool, double> ConditionNumber() { return Strategy::ConditionNumber(); }

        bool SurfaceWithinLimits(double tolerance = REFERENCETOL) const
        {
          return WithinLimits<sideType>(LocalCoordinates(), tolerance);
        }

        bool LineWithinLimits(double tolerance = REFERENCETOL) const
        {
          const LINALG::Matrix<dimEdge, 1> xsi_line(LocalCoordinates().A() + dimSide, true);
          return WithinLimits<edgeType>(xsi_line, tolerance);
        }

        double GetTolerance() const { return Strategy::GetTolerance(); }

        void WritetoGmsh(std::ofstream& file) { Strategy::WritetoGmsh(file); }

        PointOnSurfaceLoc GetSideLocation() { return sideLocation_; }

        PointOnSurfaceLoc GetEdgeLocation() { return edgeLocation_; }

        bool IsConditionInfinity() { return cond_infinity_; }

        const std::vector<int>& GetTouchedSideEdges() { return touched_edges_ids_; }

        void GetTouchedSideEdges(std::vector<int>& inout)
        {
          inout.clear();
          inout = touched_edges_ids_;
        }

        // Get location of the point, based on the extended tolerance of 1e-14 for the diagonal
        // edge of the splitted triangle, and tolerance 0.0 for the normal edges.
        // Used to avoid possibility of middle point of the quad4 split ot be missed because of
        // numerical errors
        PointOnSurfaceLoc GetSideLocationTriangleSplit()
        {
          if (sideType != DRT::Element::tri3)
            dserror("This method only works for tri3 side. Current side is %s",
                DRT::DistypeToString(sideType).c_str());

          double distance_tolerance = TOPOLOGICAL_TOLERANCE;
          // get tolerance with 1e-14 around the triangle, but in local coordinates
          LINALG::Matrix<dimSide, 1> scaled_tolerance;
          this->GetLocalToleranceSide(distance_tolerance, scaled_tolerance);
          // Diagonal one, corresponds to the middle one for both triangles, tolerance there should
          // be 1e-14
          LINALG::Matrix<3, 1> real_tolerance;
          real_tolerance(0) = 0.0;
          real_tolerance(1) = scaled_tolerance(1);
          real_tolerance(2) = 0.0;
          if (WithinLimitsSplittedQuad<sideType>(xsi_, real_tolerance))
            sideLocation_triangle_split_ = PointOnSurfaceLoc(true, true);
          else
            sideLocation_triangle_split_ = PointOnSurfaceLoc(false, true);
          return sideLocation_triangle_split_;
        }

        // Distance between edges in edge-edge intersection
        double DistanceBetween() const
        {
          if (probDim > dimEdge + dimSide)
            return distance_between_;
          else
            dserror("This method should only be used for edge-edge intersection");
        }

       private:
        // Get local tolerance and based on it get location of the point on the surface
        // as well as touched edges
        void GetTopologyInformation()
        {
          LINALG::Matrix<dimEdge + dimSide, 1> xsi = Strategy::LocalCoordinates();
          // tolerance of getting nearby edges
          LINALG::Matrix<dimSide, 1> scaled_tolerance_side_touched_edges;
          // tolerances of determining inside/outside surface limits
          LINALG::Matrix<dimEdge, 1> scaled_tolerance_edge;
          LINALG::Matrix<dimSide, 1> scaled_tolerance_side;
          // compute location finally
          this->GetLocalToleranceEdge(INSIDE_OUTSIDE_TOLERANCE, scaled_tolerance_edge);
          this->GetLocalToleranceSide(INSIDE_OUTSIDE_TOLERANCE, scaled_tolerance_side);
          this->GetLocalToleranceSide(
              SIDE_DETECTION_TOLERANCE, scaled_tolerance_side_touched_edges);
          // compute inside-outside relation
          if (WithinLimits<sideType>(xsi, scaled_tolerance_side))
            sideLocation_ = PointOnSurfaceLoc(true, true);
          else
            sideLocation_ = PointOnSurfaceLoc(false, true);
          const LINALG::Matrix<dimEdge, 1> xsi_line(xsi.A() + dimSide, true);
          if (WithinLimits<edgeType>(xsi_line, scaled_tolerance_edge))
            edgeLocation_ = PointOnSurfaceLoc(true, true);
          else
            edgeLocation_ = PointOnSurfaceLoc(false, true);
          if (probDim > dimEdge + dimSide) distance_between_ = Strategy::DistanceBetween();

          touched_edges_ids_.clear();
          if (sideLocation_.WithinSide())
          {
            GetEdgesAt<sideType>(xsi, touched_edges_ids_, scaled_tolerance_side_touched_edges);
          }
        }


        // Evaluate difference between global coordinates in the intersection when computed based on
        // side and based on edge
        double ComputeError(const LINALG::Matrix<probDim, numNodesSide>& side_xyz,
            const LINALG::Matrix<probDim, numNodesEdge>& edge_xyz,
            const LINALG::Matrix<dimEdge + dimSide, 1>& locxyz)
        {
          // Calculating interpolation from the shapefunction of the side
          LINALG::Matrix<numNodesSide, 1> sideFunct;
          LINALG::Matrix<dimEdge + dimSide, 1> globxyz_side;
          LINALG::Matrix<dimEdge + dimSide, 1> locxyz_side(locxyz.A(), true);
          DRT::UTILS::shape_function<sideType>(locxyz_side, sideFunct);
          for (unsigned int inode = 0; inode < numNodesSide; ++inode)
          {
            for (unsigned int isd = 0; isd < dimEdge + dimSide; ++isd)
            {
              globxyz_side(isd) += side_xyz(isd, inode) * sideFunct(inode);
            }
          }

          // Calculating interpolation from the shapefunction of the edge

          LINALG::Matrix<probDim, 1> globxyz_edge;
          ;
          LINALG::Matrix<numNodesEdge, 1> edgeFunct;
          LINALG::Matrix<probDim, 1> locxyz_edge(locxyz.A() + dimSide, true);

          DRT::UTILS::shape_function<edgeType>(locxyz_edge, edgeFunct);
          for (unsigned int inode = 0; inode < numNodesEdge; ++inode)
          {
            for (unsigned int isd = 0; isd < dimEdge + dimSide; ++isd)
            {
              globxyz_edge(isd) += edge_xyz(isd, inode) * edgeFunct(inode);
            }
          }

          // Evaluating error
          LINALG::Matrix<probDim, 1> diffVec;
          for (unsigned int i = 0; i < (dimEdge + dimSide); ++i)
            diffVec(i) = globxyz_edge(i) - globxyz_side(i);
          double error = diffVec.Norm2();
          if (error == 0) error = MINIMUM_DOUBLE_ERROR;
          std::pair<bool, double> cond_pair = this->ConditionNumber();
          if (not cond_pair.first)
          {
            cond_infinity_ = true;
          }
          return (error * cond_pair.second);
        }

        // touched edges ids
        static std::vector<int> touched_edges_ids_;

        // determines location of this point with respect to the edge that cuts it
        PointOnSurfaceLoc sideLocation_;

        // determines location of this point with respect to the side that cuts it
        PointOnSurfaceLoc edgeLocation_;
        // determines location of this point, with extended tolerance over side, that is shared for
        // another triangle (only for tri3)
        PointOnSurfaceLoc sideLocation_triangle_split_;

        // if condition number is infinity or not
        bool cond_infinity_;

        // holds distance between lines in case of edge-edge interesection
        double distance_between_;

        // holds reference to original value passed
        LINALG::Matrix<dimEdge + dimSide, 1>& xsi_;

#ifdef CLN_CALC
        bool checklimits_;
#endif
      };  // class GenericComputeIntersection


      /*--------------------------------------------------------------------------*/
      /** \brief most derived strategy class for intersection of side and line
       *
       *  inheritance diagram:
       *
       *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
       *  ^^^^^^^^^^^^^^^^^^^
       *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy  */
      template <unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType, bool computeCln = false,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class ComputeIntersection
          : public GenericComputeIntersection<
                NewtonSolve<ComputeIntersectionStrategy<false, probDim, edgeType, sideType>,
                    dimEdge + dimSide>,
                probDim, edgeType, sideType, computeCln>
      {
       public:
        //! constructor
        ComputeIntersection(LINALG::Matrix<dimEdge + dimSide, 1>& xsi, bool checklimits = true)
            : GenericComputeIntersection<
                  NewtonSolve<ComputeIntersectionStrategy<false, probDim, edgeType, sideType>,
                      dimEdge + dimSide>,
                  probDim, edgeType, sideType, computeCln>(xsi, checklimits)
        {
        }
      };  // class ComputeIntersection

      /*--------------------------------------------------------------------------*/
      /** \brief most derived strategy class for debugging intersection of side
       *  and line */
      template <unsigned probDim, DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType,
          unsigned dimEdge = DRT::UTILS::DisTypeToDim<edgeType>::dim,
          unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
          unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
          unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
      class DebugComputeIntersection
          : public GenericComputeIntersection<
                NewtonSolve<DebugNewtonStrategy<
                                ComputeIntersectionStrategy<true, probDim, edgeType, sideType>,
                                dimEdge + dimSide>,
                    dimEdge + dimSide>,
                probDim, edgeType, sideType>
      {
       public:
        //! constructor
        DebugComputeIntersection(LINALG::Matrix<dimEdge + dimSide, 1>& xsi, bool checklimits = true)
            : GenericComputeIntersection<
                  NewtonSolve<DebugNewtonStrategy<
                                  ComputeIntersectionStrategy<true, probDim, edgeType, sideType>,
                                  dimEdge + dimSide>,
                      dimEdge + dimSide>,
                  probDim, edgeType, sideType>(xsi, checklimits)
        {
        }

      };  // class DebugComputeIntersection

    }  // namespace KERNEL
  }    // namespace CUT
}  // namespace GEO


// static members in KERNEL
// in compute position strategy
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim, typename floatType>
LINALG::Matrix<numNodesElement, 1, floatType> GEO::CUT::KERNEL::ComputePositionStaticMembers<debug,
    probDim, elementType, numNodesElement, dim, floatType>::funct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim, typename floatType>
LINALG::Matrix<probDim, numNodesElement, floatType> GEO::CUT::KERNEL::ComputePositionStaticMembers<
    debug, probDim, elementType, numNodesElement, dim, floatType>::deriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim, typename floatType>
LINALG::Matrix<probDim, probDim, floatType> GEO::CUT::KERNEL::ComputePositionStaticMembers<debug,
    probDim, elementType, numNodesElement, dim, floatType>::A_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim, typename floatType>
LINALG::Matrix<probDim, 1, floatType> GEO::CUT::KERNEL::ComputePositionStaticMembers<debug, probDim,
    elementType, numNodesElement, dim, floatType>::b_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim, typename floatType>
LINALG::Matrix<probDim, 1, floatType> GEO::CUT::KERNEL::ComputePositionStaticMembers<debug, probDim,
    elementType, numNodesElement, dim, floatType>::dx_;

// in compute distance strategy
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<numNodesSide, 1, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug,
    probDim, sideType, dimSide, numNodesSide, floatType>::sideFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, numNodesSide, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<
    debug, probDim, sideType, dimSide, numNodesSide, floatType>::sideDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<2 * dimSide - 1, numNodesSide, floatType>
    GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug, probDim, sideType, dimSide, numNodesSide,
        floatType>::sideDeriv2_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, probDim, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug,
    probDim, sideType, dimSide, numNodesSide, floatType>::A_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, probDim, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug,
    probDim, sideType, dimSide, numNodesSide, floatType>::B_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 2 * dimSide - 1, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<
    debug, probDim, sideType, dimSide, numNodesSide, floatType>::C_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 1, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug, probDim,
    sideType, dimSide, numNodesSide, floatType>::b_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 1, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug, probDim,
    sideType, dimSide, numNodesSide, floatType>::dx_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 2, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug, probDim,
    sideType, dimSide, numNodesSide, floatType>::N_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide,
    unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 2, floatType> GEO::CUT::KERNEL::ComputeDistanceStaticMembers<debug, probDim,
    sideType, dimSide, numNodesSide, floatType>::nvec_;

// in compute intersection strategy
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<numNodesSide, 1, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::sideFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<numNodesEdge, 1, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::edgeFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<dimSide, numNodesSide, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::sideDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<dimEdge, numNodesEdge, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::edgeDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<dimEdge + dimSide, dimEdge + dimSide, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::A_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, dimEdge + dimSide, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::B_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<dimEdge + dimSide, 1, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::b_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<probDim, 1, floatType> GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug,
    probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide, floatType>::c_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide, typename floatType>
LINALG::Matrix<dimEdge + dimSide, 1, floatType>
    GEO::CUT::KERNEL::ComputeIntersectionStaticMembers<debug, probDim, edgeType, sideType, dimEdge,
        dimSide, numNodesEdge, numNodesSide, floatType>::dx_;


// for generic compute intersection
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, bool computeCln, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::GenericComputeIntersection<Strategy, probDim, edgeType, sideType,
    computeCln, dimEdge, dimSide, numNodesEdge, numNodesSide>::touched_edges_ids_;

// for generic compute distance
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    bool computeCln, unsigned dimSide, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::GenericComputeDistance<Strategy, probDim, sideType, computeCln,
    dimSide, numNodesSide>::touched_edges_ids_;

template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    bool computeCln, unsigned dimSide, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::GenericComputeDistance<Strategy, probDim, sideType, computeCln,
    dimSide, numNodesSide>::touched_nodes_ids_;

#ifdef CLN_CALC
// for generic compute position
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType elementType,
    bool computeCln, unsigned numNodesElement, unsigned dim>
std::vector<int> GEO::CUT::KERNEL::GenericComputePosition<Strategy, probDim, elementType,
    computeCln, numNodesElement, dim>::touched_edges_ids_;
#endif



#ifdef CLN_CALC

// Initializing for static arrays for AdaptivePrecision strategies


// ComputePosition
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim>
bool GEO::CUT::KERNEL::ComputePositionAdaptivePrecision<Strategy, probDim, elementType,
    numNodesElement, dim>::first_run_ = true;
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType elementType,
    unsigned numNodesElement, unsigned dim>
size_t GEO::CUT::KERNEL::ComputePositionAdaptivePrecision<Strategy, probDim, elementType,
    numNodesElement, dim>::cln_sizes_[];

// GEO::CUT::KERNEL::ComputeDistance

template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    unsigned dimSide, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::ComputeDistanceAdaptivePrecision<Strategy, probDim, sideType,
    dimSide, numNodesSide>::touched_edges_ids_;
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    unsigned dimSide, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::ComputeDistanceAdaptivePrecision<Strategy, probDim, sideType,
    dimSide, numNodesSide>::touched_nodes_ids_;
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    unsigned dimSide, unsigned numNodesSide>
bool GEO::CUT::KERNEL::ComputeDistanceAdaptivePrecision<Strategy, probDim, sideType, dimSide,
    numNodesSide>::first_run_ = true;
template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType sideType,
    unsigned dimSide, unsigned numNodesSide>
size_t GEO::CUT::KERNEL::ComputeDistanceAdaptivePrecision<Strategy, probDim, sideType, dimSide,
    numNodesSide>::cln_sizes_[];

// GEO::CUT::KERNEL::ComputeIntersection

template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide>
std::vector<int> GEO::CUT::KERNEL::ComputeIntersectionAdaptivePrecision<Strategy, probDim, edgeType,
    sideType, dimEdge, dimSide, numNodesEdge, numNodesSide>::touched_edges_ids_;

template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::KERNEL::ComputeIntersectionAdaptivePrecision<Strategy, probDim, edgeType, sideType,
    dimEdge, dimSide, numNodesEdge, numNodesSide>::first_run_ = true;

template <class Strategy, unsigned probDim, DRT::Element::DiscretizationType edgeType,
    DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide,
    unsigned numNodesEdge, unsigned numNodesSide>
size_t GEO::CUT::KERNEL::ComputeIntersectionAdaptivePrecision<Strategy, probDim, edgeType, sideType,
    dimEdge, dimSide, numNodesEdge, numNodesSide>::cln_sizes_[];

#endif

#endif
