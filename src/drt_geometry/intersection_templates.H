/*!----------------------------------------------------------------------

\brief collection of templated methods in intersection computations
      !WARNING: not used at the moment
      (remove this comment and change level as soon as this functionality is tested again!)

\level 3

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#ifndef INTERSECTION_TEMPLATES_H_
#define INTERSECTION_TEMPLATES_H_


#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/intersection_math.H"


namespace GEO
{
  /*!
      \brief tests if one entry of a fixed size vector is NaN using std::isnan function
   */
  template <unsigned length>
  bool isNaN(const LINALG::Matrix<length, 1>& xsi)
  {
    for (unsigned i = 0; i < length; i++)
    {
      if (std::isnan(xsi(i))) return true;
    }
    return false;
  }

  /*!
      \brief checks if intersection point lies within a surface and line element
   */
  template <DRT::Element::DiscretizationType surftype, class V1, class V2, class V3>
  static inline bool intersectionPointInElements(
      const V1& xsi, const V2& upLimit, const V3& loLimit)
  {
    switch (surftype)
    {
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        if ((xsi(0) > (upLimit(0) + GEO::TOL7)) || (xsi(1) > (upLimit(1) + GEO::TOL7)) ||
            (xsi(2) > (upLimit(2) + GEO::TOL7)) || (xsi(0) < (loLimit(0) - GEO::TOL7)) ||
            (xsi(1) < (loLimit(1) - GEO::TOL7)) || (xsi(2) < (loLimit(2) - GEO::TOL7)))
          return false;
        else
          return true;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        if ((xsi(0) > (upLimit(0) + GEO::TOL7)) || (xsi(1) > (upLimit(1) + GEO::TOL7)) ||
            (xsi(2) > (upLimit(2) + GEO::TOL7)) || (xsi(0) < (loLimit(0) - GEO::TOL7)) ||
            (xsi(1) < (loLimit(1) - GEO::TOL7)) || (xsi(2) < (loLimit(2) - GEO::TOL7)) ||
            (xsi(1) > (1.0 - xsi(0) + GEO::TOL7)))  // checks if in triangle
          return false;
        else
          return true;
      }
      default:
        dserror("surface type not implemented");
        return false;
    };
  }



  /*!
      \brief updates the systemmatrix at the corresponding element coordinates
             for the computation of curve surface intersections
   */
  template <DRT::Element::DiscretizationType surftype, DRT::Element::DiscretizationType linetype,
      class M1, class M2, class V>
  static inline void updateAForCSI(LINALG::Matrix<3, 3>& A,  ///< system matrix
      const V& xsi,                   ///< vector of element coordinates (r,s,t)
      const M1& xyze_surfaceElement,  ///< nodal positions of surface element
      const M2& xyze_lineElement      ///< nodal positions of line element
  )
  {
    const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
    const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;

    A.Clear();
    static LINALG::Matrix<2, numNodesSurface> surfaceDeriv1;
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);

    for (int inode = 0; inode < numNodesSurface; inode++)
      for (int isd = 0; isd < 3; isd++)
      {
        A(isd, 0) += xyze_surfaceElement(isd, inode) * surfaceDeriv1(0, inode);
        A(isd, 1) += xyze_surfaceElement(isd, inode) * surfaceDeriv1(1, inode);
      }


    static LINALG::Matrix<2, numNodesLine> lineDeriv1;
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), linetype);

    for (int inode = 0; inode < numNodesLine; inode++)
      for (int isd = 0; isd < 3; isd++)
        A(isd, 2) -= xyze_lineElement(isd, inode) * lineDeriv1(0, inode);
  }



  /*!
      \brief updates the rhs at the corresponding element coordinates
             for the computation of curve surface intersections
   */
  template <DRT::Element::DiscretizationType surftype, DRT::Element::DiscretizationType linetype,
      class M1, class M2, class V>
  static inline void updateRHSForCSI(LINALG::Matrix<3, 1>& b,  ///< right-hand-side
      const V& xsi,                   ///< vector of element coordinates (r,s,t)
      const M1& xyze_surfaceElement,  ///< nodal positions of surface element
      const M2& xyze_lineElement      ///< nodal positions of line element
  )
  {
    const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
    const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;

    b.Clear();

    static LINALG::Matrix<numNodesSurface, 1> surfaceFunct;
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);
    for (int i = 0; i < numNodesSurface; i++)
    {
      for (int dim = 0; dim < 3; dim++) b(dim) -= xyze_surfaceElement(dim, i) * surfaceFunct(i);
    }

    static LINALG::Matrix<numNodesLine, 1> lineFunct;
    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);
    for (int i = 0; i < numNodesLine; i++)
    {
      for (int dim = 0; dim < 3; dim++) b(dim) += xyze_lineElement(dim, i) * lineFunct(i);
    }
  }



  /*!
      \brief checks if intersection point lies within a surface element
   */
  template <DRT::Element::DiscretizationType surftype, class V>
  static inline bool intersectionPointInElementsRCI(const V& xsi)
  {
    switch (surftype)
    {
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        if ((xsi(0) > (1.0 + GEO::TOL7)) || (xsi(1) > (1.0 + GEO::TOL7)) ||
            (xsi(0) < (-1.0 - GEO::TOL7)) || (xsi(1) < (-1.0 - GEO::TOL7)))
          return false;
        else
          return true;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        // checks if in triangle
        if ((xsi(0) > (1.0 + GEO::TOL7)) || (xsi(1) > (1.0 - xsi(0) + GEO::TOL7)) ||
            (xsi(0) < ((-1.0) * GEO::TOL7)) || (xsi(1) < ((-1.0) * GEO::TOL7)))
          return false;
        else
          return true;
      }
      default:
        dserror("surface type not implemented");
        break;
    };
    return false;
  }



  /*!
      \brief  updates the system matrix for the intersection computation
              of an average normal and the corresponding cutter element
              to lift a point of the linearized interface onto the curved interface

      \param A                    (out)       :   system matrix
      \param xsi                  (in)        :   intersection point
      \param normal               (in)        :   nodes of the intersection normal
      \param surfaceElement       (in)        :   surface element
      \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
      \param onBoundary           (in)        :   true, if intersection normal lies on the
                                                  element boundary
   */
  template <DRT::Element::DiscretizationType surftype, class M, class V>
  static inline void updateAForRCINormal(LINALG::Matrix<3, 3>& A, const V& xsi,
      const std::vector<LINALG::Matrix<3, 1>>& normal, const M& xyze_surfaceElement,
      const bool onBoundary)
  {
    const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;

    A.Clear();
    static LINALG::Matrix<2, numNodesSurface> surfaceDeriv1;
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);

    if (!onBoundary)
    {
      for (int i = 0; i < numNodesSurface; i++)
        for (int dim = 0; dim < 3; dim++)
        {
          A(dim, 0) += xyze_surfaceElement(dim, i) * surfaceDeriv1(0, i);
          A(dim, 1) += xyze_surfaceElement(dim, i) * surfaceDeriv1(1, i);
        }

      for (int dim = 0; dim < 3; dim++) A(dim, 2) -= 0.5 * (-normal[0](dim) + normal[1](dim));
    }
    else
    {
      const int numNodesLine = 3;
      static LINALG::Matrix<1, 3> lineDeriv1;
      // LINALG::SerialDenseMatrix lineDeriv1(1,numNodesLine);
      DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), DRT::Element::line3);

      for (int i = 0; i < numNodesSurface; i++)
        for (int dim = 0; dim < 3; dim++)
        {
          A(dim, 0) += xyze_surfaceElement(dim, i) * surfaceDeriv1(0, i);
          A(dim, 1) += xyze_surfaceElement(dim, i) * surfaceDeriv1(1, i);
        }

      for (int i = 0; i < numNodesLine; i++)
        for (int dim = 0; dim < 3; dim++)
        {
          int index = i;
          if (i > 1) index = 4;
          A(dim, 2) -= normal[index](dim) * lineDeriv1(0, i);
        }
    }
  }



  /*!
      \brief  updates the right-hand-side for the intersection computation
              of an average normal and the corresponding cutter element
              to lift a point of the linearized interface onto the curved interface

      \param b                    (out)       :   right-hand-side
      \param xsi                  (in)        :   intersection point
      \param normal               (in)        :   nodes of the intersection normal
      \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
      \param onBoundary           (in)        :   true, if intersection normal lies on the
                                                  element boundary
   */
  template <DRT::Element::DiscretizationType surftype, class M, class V>
  static inline void updateRHSForRCINormal(LINALG::Matrix<3, 1>& b, const V& xsi,
      const std::vector<LINALG::Matrix<3, 1>>& normal, const M& xyze_surfaceElement,
      const bool onBoundary)
  {
    const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
    static LINALG::Matrix<numNodesSurface, 1> surfaceFunct;
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);

    b.Clear();
    if (!onBoundary)
    {
      for (int i = 0; i < numNodesSurface; i++)
        for (int dim = 0; dim < 3; dim++) b(dim) -= xyze_surfaceElement(dim, i) * surfaceFunct(i);


      for (int dim = 0; dim < 3; dim++)
        b(dim) += 0.5 * (normal[0](dim) * (1.0 - xsi(2)) + normal[1](dim) * (1.0 + xsi(2)));
    }
    else
    {
      const int numNodesLine = 3;
      static LINALG::Matrix<3, 1> lineFunct;
      DRT::UTILS::shape_function_1D(lineFunct, xsi(2), DRT::Element::line3);

      for (int i = 0; i < numNodesSurface; i++)
        for (int dim = 0; dim < 3; dim++) b(dim) -= xyze_surfaceElement(dim, i) * surfaceFunct(i);



      for (int i = 0; i < numNodesLine; i++)
        for (int dim = 0; dim < 3; dim++)
        {
          int index = i;
          if (i > 1) index = 4;

          b(dim) += normal[index](dim) * lineFunct(i);
        }
    }
  }



  /*!
      \brief  updates the system matrix for the intersection computation
              of an average normal and the corresponding cutter element
              to lift a point of the linearized interface onto the curved interface

      \param A                        (out)       :   system matrix
      \param xsi                      (in)        :   intersection point
      \param plane                    (in)        :   nodes of the intersection normal
      \param xyze_lineElement         (in)        :   nodal coordinates of line element
   */
  template <DRT::Element::DiscretizationType linetype, class V, class M>
  static inline void updateAForRCIPlane(LINALG::Matrix<3, 3>& A, const V& xsi,
      const std::vector<LINALG::Matrix<3, 1>>& plane, const M& xyze_lineElement)
  {
    const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
    const int numNodesSurface = 4;

    static LINALG::Matrix<2, numNodesSurface> surfaceDeriv;
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv, xsi(0), xsi(1), DRT::Element::quad4);
    static LINALG::Matrix<1, numNodesLine> lineDeriv;
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv, xsi(2), linetype);

    dsassert((int)plane.size() >= numNodesSurface,
        "plane array has to have size numNodesSurface ( = 4)!");

    A.Clear();
    for (int dim = 0; dim < 3; dim++)
      for (int i = 0; i < numNodesSurface; i++)
      {
        A(dim, 0) += plane[i](dim) * surfaceDeriv(0, i);
        A(dim, 1) += plane[i](dim) * surfaceDeriv(1, i);
      }

    for (int i = 0; i < numNodesLine; i++)
      for (int dim = 0; dim < 3; dim++) A(dim, 2) -= xyze_lineElement(dim, i) * lineDeriv(0, i);
  }



  /*!
      \brief  updates the right-hand-side for the intersection computation
              of the cutter element edge with a plane through the midpoint
              of the face edge to lift a point of the linearized interface
              onto the curved interface

      \param b                        (out)       :   right-hand-side
      \param xsi                      (in)        :   intersection point
      \param plane                    (in)        :   nodes of the intersection normal
      \param xyze_lineElement         (in)        :   nodal coordinates of line element
   */
  template <DRT::Element::DiscretizationType linetype, class V, class M>
  static inline void updateRHSForRCIPlane(LINALG::Matrix<3, 1>& b, const V& xsi,
      const std::vector<LINALG::Matrix<3, 1>>& plane, const M& xyze_lineElement)
  {
    const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
    const int numNodesSurface = 4;

    static LINALG::Matrix<numNodesSurface, 1> surfaceFunct;
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), DRT::Element::quad4);
    static LINALG::Matrix<numNodesLine, 1> lineFunct;
    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);

    dsassert((int)plane.size() >= numNodesSurface,
        "plane array has to have size numNodesSurface ( = 4)!");

    b.Clear();
    for (int dim = 0; dim < 3; dim++)
      for (int i = 0; i < numNodesSurface; i++) b(dim) -= plane[i](dim) * surfaceFunct(i);

    for (int i = 0; i < numNodesLine; i++)
      for (int dim = 0; dim < 3; dim++) b(dim) += xyze_lineElement(dim, i) * lineFunct(i);
  }



  /*!
      \brief  computes the recovery of the curved interface
              for one xfem element after a Contrained Delaunay
              Tetrahedralization.

      \param lineIndex                (in)    :   line index
      \param currentcutterpositions   (in)    :   current cutter node positions
      \param xsi                      (out)   :   intersection point
      \param plane                    (in)    :   nodes of the intersection normal
      \param cutterElement            (in)    :   cutter element
      \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element

      \return true if plane intersects with the line element or false otherwise
   */
  template <DRT::Element::DiscretizationType linetype, class V, class M>
  bool computeRecoveryPlaneT(
      const M& xyze_lineElement, const std::vector<LINALG::Matrix<3, 1>>& plane, V& xsi)
  {
    int iter = 0;
    const int maxiter = 20;
    double residual = 1.0;

    static LINALG::Matrix<3, 3> A;
    static LINALG::Matrix<3, 1> b;
    static LINALG::Matrix<3, 1> dx;

    updateRHSForRCIPlane<linetype>(b, xsi, plane, xyze_lineElement);

    while (residual > GEO::TOL13)
    {
      updateAForRCIPlane<linetype>(A, xsi, plane, xyze_lineElement);

      // check if singular computatio necessary
      double det = LINALG::gaussElimination<true, 3>(A, b, dx);
      if (fabs(det) > TOL14) return false;

      if (iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi))
        return false;

      xsi += dx;
      updateRHSForRCIPlane<linetype>(b, xsi, plane, xyze_lineElement);
      residual = b.Norm2();
      iter++;
      // printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol =
      // %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
    }

    if ((xsi(2) > (1.0 + GEO::TOL7)) ||
        (xsi(2) < (-1.0 - GEO::TOL7)))  // planes coordinate may be bigger than 1
      return false;

    // printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol =
    // %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
    return true;
  }



  /*!
      \brief  computes the recovery of the curved interface
              for one xfem element after a Contrained Delaunay
              Tetrahedralization.

      \param lineIndex                (in)    :   line index
      \param currentcutterpositions   (in)    :   current cutter node positions
      \param xsi                      (out)   :   intersection point
      \param plane                    (in)    :   nodes of the intersection normal
      \param cutterElement            (in)    :   cutter element
      \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element

      \return true if plane intersects with the line element or false otherwise
   */
  template <class M, class V>
  bool computeRecoveryPlane(const DRT::Element* lineElement, const M& xyze_lineElement,
      const std::vector<LINALG::Matrix<3, 1>>& plane, V& xsi)
  {
    switch (lineElement->Shape())
    {
      case DRT::Element::line2:
        return computeRecoveryPlaneT<DRT::Element::line2>(xyze_lineElement, plane, xsi);
      case DRT::Element::line3:
        return computeRecoveryPlaneT<DRT::Element::line3>(xyze_lineElement, plane, xsi);
      default:
        dserror("linetype not yet implemented");
        return false;
    };
    return true;
  }


} /*namespace GEO*/



#endif /*INTERSECTION_TEMPLATES_H_*/
