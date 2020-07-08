/*----------------------------------------------------------------------*/
/*! \file

\brief calculate second global derivatives w.r.t x/y/z at point r,s,t

\level 1
*/
#ifndef UTILS_GDER2_H
#define UTILS_GDER2_H

#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
  namespace UTILS
  {
    /*---------------------------------------------------------------------*/
    /*!
     * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
     */
    /*                                                          gammi 07/07
    |
    | To compute the second derivative, the first derivative is derived once more.
    | This is an application of the product rule and the chain rule.
    | The trick to compute it is to realize that the transformation dr/dx must be
    | derived w.r.t. dx (but dr/dx = (dx/dr)^-1), which gives the term -(dx/dr)^{-2}
    | plus some inner derivatives.
    |
    | This formula is valid in all dimension if one interprets x and r as tensors
    | and interprets multiplication with (dx/dr)^-1 from the right is as a matrix
    | multiplication with (dx/dr)^-T.
    |
    |             +-            -+
    |  d^2N    d  |  /dx\-1   dN |   d   /dx\-1   dN     /dx\    d  dN
    |  ----  = -- | | -- |  * -- | = -- | -- |  * --  + | -- | * -- --
    |  dx^2    dx |  \dr/     dr |   dx  \dr/     dr     \dr/    dx dr
    |             +-            -+
    |
    |             /dx\-2 dN d  dx    /dx\-1   d^2N    /dx\-1
    |        = - | -- | *--*-- -- + | -- |  * ---- * | -- |
    |             \dr/   dr dx dr    \dr/     dr^2    \dr/

    |             /dx\-1 dN d^2x  /dx\-1     /dx\-1 d^2N  /dx\-1
    |        = - | -- | *--*----*| -- |   + | -- | *----*| -- |
    |             \dr/   dr dr^2  \dr/       \dr/   dr^2  \dr/

    |                    +-                -+
    |           /dx\-1   |   dN d^2x   d^2N |    /dx\-1
    |        = | -- |  * | - --*---- + ---- | * | -- |
    |           \dr/     |   dx dr^2   dr^2 |    \dr/
    |                    +-                -+
    *----------------------------------------------------------------------*/
    // for enriched elements (e.g. xwall), NumNode may be larger than the number of element nodes
    // for all other elements, NumNode==numnode
    template <DRT::Element::DiscretizationType DISTYPE, int NumNode, int ProbDim>
    void gder2(const LINALG::Matrix<ProbDim, ProbDim>& xjm,
        const LINALG::Matrix<ProbDim, NumNode>& derxy,
        const LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2, NumNode>& deriv2,
        const LINALG::Matrix<ProbDim, NumNode>& xyze,
        LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2, NumNode>& derxy2)
    {
      // some numbers already known during compilation
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      dsassert(numnode <= NumNode, "Expect at least numNodePerElement matrix columns");
      const int nsd = ProbDim;
      const int numderiv2 = DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2;

      // compute d^2x/dr^2
      double xder2[numderiv2 * nsd];
      LINALG::DENSEFUNCTIONS::multiplyNT<double, numderiv2, numnode, nsd>(
          xder2, deriv2.A(), xyze.A());

      // compute -(dN/dx)*(d^2x/dr^2)
      LINALG::DENSEFUNCTIONS::multiply<double, numderiv2, nsd, numnode>(
          derxy2.A(), -1.0, xder2, derxy.A());

      // compute -(dN/dx)*(d^2x/dr^2) + (d^2N/dr^2)
      derxy2.Update(1.0, deriv2, 1.0);

      // finally multiply by (dx/dr)^-1 from the left and (dx/dr)^-T from the right
      // write out the products by hand because derxy2 only stores the symmetric part
      // of the second derivative tensor (which means we cannot use matrix-matrix products)
      double xjiData[9];
      double* xji[3];
      for (int i = 0; i < nsd; ++i) xji[i] = &xjiData[nsd * i];
      LINALG::DENSEFUNCTIONS::invert<double, nsd, nsd>(xjiData, xjm.A());
      for (int node = 0; node < numnode; ++node)
      {
        double tmp[3][3];
        switch (nsd)
        {
          case 3:
            for (int j = 0; j < nsd; ++j)
            {
              tmp[0][j] = derxy2(0, node) * xji[0][j] + derxy2(3, node) * xji[1][j] +
                          derxy2(4, node) * xji[2][j];
              tmp[1][j] = derxy2(3, node) * xji[0][j] + derxy2(1, node) * xji[1][j] +
                          derxy2(5, node) * xji[2][j];
              tmp[2][j] = derxy2(4, node) * xji[0][j] + derxy2(5, node) * xji[1][j] +
                          derxy2(2, node) * xji[2][j];
            }
            derxy2(0, node) = xji[0][0] * tmp[0][0] + xji[1][0] * tmp[1][0] + xji[2][0] * tmp[2][0];
            derxy2(1, node) = xji[0][1] * tmp[0][1] + xji[1][1] * tmp[1][1] + xji[2][1] * tmp[2][1];
            derxy2(2, node) = xji[0][2] * tmp[0][2] + xji[1][2] * tmp[1][2] + xji[2][2] * tmp[2][2];
            derxy2(3, node) = xji[0][0] * tmp[0][1] + xji[1][0] * tmp[1][1] + xji[2][0] * tmp[2][1];
            derxy2(4, node) = xji[0][0] * tmp[0][2] + xji[1][0] * tmp[1][2] + xji[2][0] * tmp[2][2];
            derxy2(5, node) = xji[0][1] * tmp[0][2] + xji[1][1] * tmp[1][2] + xji[2][1] * tmp[2][2];
            break;
          case 2:
            for (int j = 0; j < nsd; ++j)
            {
              tmp[0][j] = derxy2(0, node) * xji[0][j] + derxy2(2, node) * xji[1][j];
              tmp[1][j] = derxy2(2, node) * xji[0][j] + derxy2(1, node) * xji[1][j];
            }
            derxy2(0, node) = xji[0][0] * tmp[0][0] + xji[1][0] * tmp[1][0];
            derxy2(1, node) = xji[0][1] * tmp[0][1] + xji[1][1] * tmp[1][1];
            derxy2(2, node) = xji[0][0] * tmp[0][1] + xji[1][0] * tmp[1][1];
            break;
          case 1:
            derxy2(0, node) *= xji[0][0] * xji[0][0];
            break;
          default:
            dserror("Illegal number of space dimensions: %d", nsd);
            break;
        }
      }

      return;

    }  // DRT::UTILS::gder2

    template <DRT::Element::DiscretizationType DISTYPE, int NumNode>
    void gder2(const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
                   DRT::UTILS::DisTypeToDim<DISTYPE>::dim>& xjm,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim, NumNode>& derxy,
        const LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2, NumNode>& deriv2,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim, NumNode>& xyze,
        LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2, NumNode>& derxy2)
    {
      gder2<DISTYPE, NumNode, DRT::UTILS::DisTypeToDim<DISTYPE>::dim>(
          xjm, derxy, deriv2, xyze, derxy2);
    }


  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifndef UTILS_GDER2_H
