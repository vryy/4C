/*!

\brief methods for the integration over boundary elements

1) computation of kovariant metric tensor for surface element
2) mapping of gausspoints on surface element to 3d space of parent element
   (required for integrations of parent-element shape functions
    over boundary elements, required for example in weak
    dirichlet boundary conditions).

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089-289-15235
*/

#ifndef UTILS_BOUNDARY_INTEGRATION_H
#define UTILS_BOUNDARY_INTEGRATION_H


#include "drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_utils_integration.H"
#include "drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_lib/drt_discret.H"

namespace DRT
{
  namespace UTILS
  {
    struct IntegrationPoints1D;
    struct IntegrationPoints2D;

    /*-----------------------------------------------------------------

      \brief Transform Gausspoints on surface element to 3d space of
             parent element (required for integrations of parent-element
             shape functions over boundary elements, required for example
             in weak dirichlet boundary conditions).

    \param pqxg       (o) transformed integration points of surface
                          elements in parent elements reference
                          coordinates
    \param derivtrafo (o) transformation matrix between parent and boundary local coordinates
                          for derivatives w.r.t. local coordinates
    \param intpoints  (i) integration points of surface element in its
                          reference coordinates
    \param pdistype   (i) discretisation type of parent element
    \param distype    (i) discretisation type of boundary element
    \param lineid     (i) local id of boundary element

    \author gammi (gamnitzer@lnm.mw.tum.de), 05/09

      -----------------------------------------------------------------*/
    template <class V, class W, typename IntegrationPoints>
    void SurfaceGPToParentGP(V& pqxg, W& derivtrafo, const IntegrationPoints& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid)
    {
      if ((distype == DRT::Element::quad4 && pdistype == DRT::Element::hex8) or
          (distype == DRT::Element::quad9 && pdistype == DRT::Element::hex27))
      {
        switch (surfaceid)
        {
          case 0:
          {
            // t=-1
            /*
                      parent               surface

                       r|                    s|
                        |                     |
                   1         2           3         2
                    +-------+             +-------+
                    |   |   |      s      |   |   |      r
                    |   +-- |  -----      |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   0         3           0         1
            */

            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[1];
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = -1.0;
            }
            derivtrafo(0, 1) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(2, 2) = -1.0;
            break;
          }
          case 1:
          {
            // s=-1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                   4         5           3         2
                    +-------+             +-------+
                    |   |   |      r      |   |   |      r
                    |   +-- |  -----      |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   0         1           0         1
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = -1.0;
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 2) = -1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 2:
          {
            // r= 1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                   5         6           3         2
                    +-------+             +-------+
                    |   |   |      s      |   |   |      r
                    |   +-- |  -----      |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   1         2           0         1
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 2) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 3:
          {
            // s= 1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                   6         7           3         2
                    +-------+             +-------+
              r     |   |   |             |   |   |      r
              ----  |   +-- |             |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   2         3           0         1
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = -intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = 1.0;
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = -1.0;
            derivtrafo(1, 2) = 1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 4:
          {
            // r=-1
            /*
                      parent               surface

                       s|                    s|
                        |                     |
                   3         7           3         2
                    +-------+             +-------+
                    |   |   |      t      |   |   |      r
                    |   +-- |  -----      |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   0         4           0         1
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = intpoints.Point(iquad)[1];
              pqxg(iquad, 2) = intpoints.Point(iquad)[0];
            }
            derivtrafo(0, 2) = -1.0;
            derivtrafo(1, 1) = 1.0;
            derivtrafo(2, 0) = 1.0;
            break;
          }
          case 5:
          {
            // t=1
            /*
                      parent               surface

                       s|                    s|
                        |                     |
                   7         6           3         2
                    +-------+             +-------+
                    |   |   |      r      |   |   |      r
                    |   +-- |  -----      |   +-- |  -----
                    |       |             |       |
                    +-------+             +-------+
                   4         5           0         1
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = intpoints.Point(iquad)[1];
              pqxg(iquad, 2) = 1.0;
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 1) = 1.0;
            derivtrafo(2, 2) = 1.0;
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else if (distype == DRT::Element::nurbs9 && pdistype == DRT::Element::nurbs27)
      {
        switch (surfaceid)
        {
          case 0:
          {
            // t=-1
            /*
                      parent               surface

                       s|                    s|
                        |                     |
                    +---+---+             +---+---+
                   6|  7|  8|      r     6|  7|  8|      r
                    +   +-- +  -----      +   +-- +  -----
                   3|  4   5|            3|  4   5|
                    +---+---+             +---+---+
                   0   1   2             0   1   2
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = intpoints.Point(iquad)[1];
              pqxg(iquad, 2) = -1.0;
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 1) = 1.0;
            derivtrafo(2, 2) = -1.0;
            break;
          }
          case 1:
          {
            // t=+1
            /*
                      parent               surface

                       s|                    s|
                        |                     |
                    +---+---+             +---+---+
                  24| 25| 26|      r     6|  7|  8|      r
                    +   +-- +  -----      +   +-- +  -----
                  21| 22  23|            3|  4   5|
                    +---+---+             +---+---+
                  18  19  20             0   1   2
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = intpoints.Point(iquad)[1];
              pqxg(iquad, 2) = 1.0;
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 1) = 1.0;
            derivtrafo(2, 2) = 1.0;
            break;
          }
          case 2:
          {
            // s=-1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                    +---+---+             +---+---+
                  18| 19| 20|      r     6|  7|  8|      r
                    +   +-- +  -----      +   +-- +  -----
                   9| 10  11|            3|  4   5|
                    +---+---+             +---+---+
                   0   1   2             0   1   2
            */

            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = -1.0;
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 2) = -1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 3:
          {
            // s=+1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                    +---+---+             +---+---+
                  24| 25| 26|    r       6|  7|  8|      r
                    +   +-- + ----        +   +-- +  -----
                  15| 16  17|            3|  4   5|
                    +---+---+             +---+---+
                   6   7   8             0   1   2
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = 1.0;
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 2) = 1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 4:
          {
            // r=+1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                    +---+---+             +---+---+
                  20| 23| 26|      s     6|  7|  8|      r
                    +   +-- +  -----      +   +-- +  -----
                  11| 14  17|            3|  4   5|
                    +---+---+             +---+---+
                   2   5   8             0   1   2
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 2) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 5:
          {
            // r=-1
            /*
                      parent               surface

                       t|                    s|
                        |                     |
                    +---+---+             +---+---+
                  18| 21| 24|      s     6|  7|  8|      r
                    +   +-- +  -----      +   +-- +  -----
                   9| 12  15|            3|  4   5|
                    +---+---+             +---+---+
                   0   3   6             0   1   2
            */
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 2) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else if ((distype == DRT::Element::tri3 && pdistype == DRT::Element::tet4) or
               (distype == DRT::Element::tri6 && pdistype == DRT::Element::tet10))
      {
        switch (surfaceid)
        {
          case 0:
          {
            // s=0
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[0];
              pqxg(iquad, 1) = 0.0;
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = 1.0;
            derivtrafo(1, 2) = -1.0;
            derivtrafo(2, 1) = 1.0;
            break;
          }
          case 1:
          {
            // r+s+t=1
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0 - intpoints.Point(iquad)[0] - intpoints.Point(iquad)[1];
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = intpoints.Point(iquad)[1];
            }
            derivtrafo(0, 0) = -1.0;
            derivtrafo(0, 1) = -1.0;
            derivtrafo(0, 2) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(1, 2) = 1.0;
            derivtrafo(2, 1) = 1.0;
            derivtrafo(2, 2) = 1.0;
            break;
          }
          case 2:
          {
            // r=0
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.0;
              pqxg(iquad, 1) = intpoints.Point(iquad)[1];
              pqxg(iquad, 2) = intpoints.Point(iquad)[0];
            }
            derivtrafo(0, 2) = -1.0;
            derivtrafo(1, 1) = 1.0;
            derivtrafo(2, 0) = 1.0;
            break;
          }
          case 3:
          {
            // t=0
            for (int iquad = 0; iquad < pqxg.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints.Point(iquad)[1];
              pqxg(iquad, 1) = intpoints.Point(iquad)[0];
              pqxg(iquad, 2) = 0.0;
            }
            derivtrafo(0, 1) = 1.0;
            derivtrafo(1, 0) = 1.0;
            derivtrafo(2, 2) = -1.0;
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else
      {
        dserror(
            "only quad4/hex8, quad9/hex27, tri3/tet4 and nurbs9/nurbs27 mappings of surface "
            "gausspoint to parent element implemented up to now\n");
      }

      return;
    }

    /*-----------------------------------------------------------------

    \brief Transform Gausspoints on line element to 2d space of
           parent element (required for integrations of parent-element
           shape functions over boundary elements, for example
           in weak dirichlet boundary conditions).


    \param pqxg       (o) transformed integration points of surface
                          elements in parent elements reference
                          coordinates
    \param derivtrafo (o) transformation matrix between parent and boundary local coordinates
                          for derivatives w.r.t. local coordinates
    \param intpoints  (i) integration points of surface element in its
                          reference coordinates
    \param pdistype   (i) discretisation type of parent element
    \param distype    (i) discretisation type of boundary element
    \param lineid     (i) local id of boundary element

    \author gammi (gamnitzer@lnm.mw.tum.de), 05/09

      -----------------------------------------------------------------*/
    template <class V, class W, typename IntegrationPoints>
    void LineGPToParentGP(V& pqxg, W& derivtrafo, const IntegrationPoints& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int lineid);

    /*-----------------------------------------------------------------

    \brief Template version of transformation of Gausspoints on boundary element to space of
           parent element

    \param pqxg       (o) transformed integration points of surface
                          elements in parent elements reference
                          coordinates
    \param derivtrafo (o) transformation matrix between parent and boundary local coordinates
                          for derivatives w.r.t. local coordinates
    \param intpoints  (i) integration points of surface element in its
                          reference coordinates
    \param pdistype   (i) discretisation type of parent element
    \param distype    (i) discretisation type of boundary element
    \param lineid     (i) local id of boundary element

      ----------------------------------------------------------------------------------*/
    template <const int NSD>
    void BoundaryGPToParentGP(LINALG::SerialDenseMatrix& pqxg, Epetra_SerialDenseMatrix& derivtrafo,
        const DRT::UTILS::IntPointsAndWeights<NSD - 1>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 3D
    template <>
    void BoundaryGPToParentGP<3>(LINALG::SerialDenseMatrix& pqxg,
        Epetra_SerialDenseMatrix& derivtrafo, const DRT::UTILS::IntPointsAndWeights<2>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 2D
    template <>
    void BoundaryGPToParentGP<2>(LINALG::SerialDenseMatrix& pqxg,
        Epetra_SerialDenseMatrix& derivtrafo, const DRT::UTILS::IntPointsAndWeights<1>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);



    template <const int NSD>
    void BoundaryGPToParentGP(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<NSD, NSD>& derivtrafo,
        const DRT::UTILS::IntPointsAndWeights<NSD - 1>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 3D
    template <>
    void BoundaryGPToParentGP<3>(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<3, 3>& derivtrafo,
        const DRT::UTILS::IntPointsAndWeights<2>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 2D
    template <>
    void BoundaryGPToParentGP<2>(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<2, 2>& derivtrafo,
        const DRT::UTILS::IntPointsAndWeights<1>& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);



    template <const int NSD>
    void BoundaryGPToParentGP(LINALG::SerialDenseMatrix& pqxg, Epetra_SerialDenseMatrix& derivtrafo,
        const DRT::UTILS::GaussPoints& intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 3D
    template <>
    void BoundaryGPToParentGP<3>(LINALG::SerialDenseMatrix& pqxg,
        Epetra_SerialDenseMatrix& derivtrafo, const DRT::UTILS::GaussPoints& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 2D
    template <>
    void BoundaryGPToParentGP<2>(LINALG::SerialDenseMatrix& pqxg,
        Epetra_SerialDenseMatrix& derivtrafo, const DRT::UTILS::GaussPoints& intpoints,
        const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    template <const int NSD>
    void BoundaryGPToParentGP(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<NSD, NSD>& derivtrafo,
        const DRT::UTILS::GaussPoints& intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 3D
    template <>
    void BoundaryGPToParentGP<3>(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<3, 3>& derivtrafo,
        const DRT::UTILS::GaussPoints& intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    //! specialization for 2D
    template <>
    void BoundaryGPToParentGP<2>(LINALG::SerialDenseMatrix& pqxg, LINALG::Matrix<2, 2>& derivtrafo,
        const DRT::UTILS::GaussPoints& intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int surfaceid);

    /*----------------------------------------------------------------------------------------*/

    //! compute kovariant metric tensor for surface element
    void ComputeMetricTensorForSurface(const Epetra_SerialDenseMatrix& xyze,
        const Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& metrictensor,
        double* sqrtdetg);

    //! compute covariant metric tensor for surface/line element and optionally, the normalized
    //! normal vector at the Gausspoint (template)
    template <DRT::Element::DiscretizationType DISTYPE, int probdim, typename valueType = double>
    void ComputeMetricTensorForBoundaryEle(
        const LINALG::Matrix<probdim,
            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement, valueType>& xyze,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement, valueType>& deriv,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
            DRT::UTILS::DisTypeToDim<DISTYPE>::dim, valueType>& metrictensor,
        valueType& sqrtdetg, LINALG::Matrix<probdim, 1, valueType>* normalvec = NULL,
        const bool& unit_normal = true, const bool& throw_error = true)
    {
      /* 2D boundary Element
      |                                              0 1 2
      |                                             +-+-+-+
      |       0 1 2              0...iel-1          | | | | 0
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |      | | | | 1           | | | | | 0        | | | | .
      |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
      |      | | | | 2           | | | | | 1        | | | | .
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |                                             | | | | iel-1
      |                                             +-+-+-+
      |
      |       dxyzdrs             deriv              xyze^T
      |
      |
      |                                 +-            -+
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | dr   dr   dr |
      |     yields           dxyzdrs =  |              |
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | ds   ds   ds |
      |                                 +-            -+
      |
      */

      /* 1D boundary Element
      |
      |   dxyzdrs(1,2) = deriv(1,iel) * xyze(2, iel)^T
       */
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim, probdim, valueType> dxyzdrs;
      dxyzdrs.MultiplyNT(deriv, xyze);

      /* 2D boundary Element
      |
      |      +-           -+    +-            -+   +-            -+ T
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
      |      |             |    | dr   dr   dr |   | dr   dr   dr |
      |      |             |  = |              | * |              |
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
      |      |             |    | ds   ds   ds |   | ds   ds   ds |
      |      +-           -+    +-            -+   +-            -+
      |
      | the calculation of g21 is redundant since g21=g12
      */

      /* 1D boundary Element
      |
      |   metrictensor(1,1) = dxyzdrs(1,2) * dxyzdrs(1,2)^T
       */

      metrictensor.Clear();
      metrictensor.MultiplyNT(dxyzdrs, dxyzdrs);

      /*
                                +--------------+
                               /               |
                 sqrtdetg =   /  g11*g22-g12^2
                            \/
      */
      sqrtdetg = metrictensor.Determinant();
      if (sqrtdetg > 0.0)
        sqrtdetg = std::sqrt(sqrtdetg);
      else if (throw_error)
      {
        std::ostringstream err_msg;
        err_msg << "--- ERROR DETECTED ---\n";
        err_msg << "xyze = \n" << xyze;
        err_msg << "The determinant of the metric-tensor is equal zero or negative! ( det( g ) = "
                << sqrtdetg << " )";
        dserror(err_msg.str().c_str());
      }

      // Calculate outward pointing unit normal vector
      if (normalvec)
      {
        if (probdim == 3 and DRT::UTILS::DisTypeToDim<DISTYPE>::dim == 2)
        {
          (*normalvec)(0) = dxyzdrs(0, 1) * dxyzdrs(1, 2) - dxyzdrs(1, 1) * dxyzdrs(0, 2);
          (*normalvec)(1) = dxyzdrs(0, 2) * dxyzdrs(1, 0) - dxyzdrs(1, 2) * dxyzdrs(0, 0);
          (*normalvec)(2) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(1, 0) * dxyzdrs(0, 1);
          // compute unit normal (outward pointing)
          if (unit_normal)
          {
            const valueType norm2 = normalvec->Norm2();
            if (norm2 < 0.0)
              dserror(
                  "The l2-norm of the normal vector is smaller than 0.0!"
                  " ( norm2 = %e )",
                  (const double)norm2);
            normalvec->Scale(1.0 / norm2);
          }
        }
        else if (probdim == 2 and DRT::UTILS::DisTypeToDim<DISTYPE>::dim == 1)
        {
          (*normalvec)(0) = dxyzdrs(0, 1);
          (*normalvec)(1) = -dxyzdrs(0, 0);
          // compute unit normal (outward pointing)
          if (unit_normal)
          {
            const valueType norm2 = normalvec->Norm2();
            if (norm2 < 0.0)
              dserror(
                  "The l2-norm of the normal vector is smaller than 0.0!"
                  " ( norm2 = %e )",
                  (const double)norm2);
            normalvec->Scale(1.0 / norm2);
          }
        }
        else if (probdim == 3 and DRT::UTILS::DisTypeToDim<DISTYPE>::dim == 1)
        {
          if (std::abs(dxyzdrs(0, 0)) > 1.0e-6 or std::abs(dxyzdrs(0, 1)) > 1.0e-6)
          {
            (*normalvec)(0) = dxyzdrs(0, 1);
            (*normalvec)(1) = -dxyzdrs(0, 0);
            (*normalvec)(2) = 0.0;
          }
          // to handle small normal vectors with dxyzdrs(0,0) < 1e-6
          // and other componets  = 0.0
          else if (dxyzdrs(0, 2) != 0.0 or dxyzdrs(0, 1) != 0.0)
          {
            (*normalvec)(0) = 0.0;
            (*normalvec)(1) = -dxyzdrs(0, 2);
            (*normalvec)(2) = dxyzdrs(0, 1);
          }
          else
          {
            (*normalvec)(0) = dxyzdrs(0, 1);
            (*normalvec)(1) = -dxyzdrs(0, 0);
            (*normalvec)(2) = 0.0;
          }

          // compute unit normal (outward pointing)
          if (unit_normal)
          {
            const valueType norm2 = normalvec->Norm2();
            if (norm2 < 0.0)
              dserror(
                  "The l2-norm of the normal vector is smaller than 0.0!"
                  " ( norm2 = %e )",
                  (const double)norm2);
            normalvec->Scale(1.0 / norm2);
          }
        }
        else
          dserror("There are only 2D and 1D boundary elements");
      }

      return;
    }

    //! compute kovariant metric tensor for surface/line element and optionally, the normalized
    //! normal vector at the Gausspoint (template)
    template <DRT::Element::DiscretizationType DISTYPE>
    void ComputeMetricTensorForBoundaryEle(
        const LINALG::Matrix<(1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim),
            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& xyze,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& deriv,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
            DRT::UTILS::DisTypeToDim<DISTYPE>::dim>& metrictensor,
        double& sqrtdetg,
        LINALG::Matrix<(1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim), 1>* normalvec = NULL)
    {
      ComputeMetricTensorForBoundaryEle<DISTYPE, (1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim)>(
          xyze, deriv, metrictensor, sqrtdetg, normalvec);
    }

    /*-----------------------------------------------------------------

    \brief Transform Gausspoints on boundary element to space of
           parent element (required for integrations of parent-element
           shape functions over boundary elements, for example
           in weak dirichlet boundary conditions).
           2D-version
      -----------------------------------------------------------------*/
    inline void BoundaryGPToParentGP2(Epetra_SerialDenseMatrix& pqxg,
        const Epetra_SerialDenseMatrix intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int beleid)
    {
      // resize output array
      pqxg.Shape(intpoints.M(), 2);

      if ((distype == DRT::Element::line2 && pdistype == DRT::Element::quad4) or
          (distype == DRT::Element::line3 && pdistype == DRT::Element::quad8) or
          (distype == DRT::Element::line3 && pdistype == DRT::Element::quad9))
      {
        switch (beleid)
        {
          case 0:
          {
            /*                s|
                               |

                               3                   2
                               +-----------------+
                               |                 |
                               |                 |
                               |                 |
                               |        |        |             r
                               |        +--      |         -----
                               |                 |
                               |                 |
                               |                 |
                               |                 |
                               +-----------*-----+
                               0                   1
                                     -->|gp|<--               */


            // s=-1
            /*

            parent                line

            r                     r
            +---+---+  -----      +---+---+ ------
            0   1   2             0   1   2

            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = -1.0;
            }
            break;
          }
          case 1:
          {
            /*                s|
                               |

                               3                   2
                               +-----------------+
                               |                 | |
                               |                 | v
                               |                 *---
                               |        |        | gp          r
                               |        +--      |---      -----
                               |                 | ^
                               |                 | |
                               |                 |
                               |                 |
                               +-----------------+
                               0                   1
             */

            // r=+1
            /*
              parent               surface

             s|                        r|
              |                         |
              +                     +
             8|                    2|
              +                     +
             5|                    1|
              +                     +
             2                     0
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
            }
            break;
          }
          case 2:
          {
            /*                s|
                               |

                              3   -->|gp|<--
                               +-----*-----------+
                               |                 |
                               |                 |
                               |                 |
                               |        |        |             r
                               |        +--      |         -----
                               |                 |
                               |                 |
                               |                 |
                               |                 |
                               +-----------------+
                              0                   1
            */

            // s=+1
            /*

            parent                line

            r                           r
            +---+---+  -----             +---+---+ -----
            6   7   8                    0   1   2

            */

            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -intpoints(iquad, 0);
              pqxg(iquad, 1) = 1.0;
            }
            break;
          }
          case 3:
          {
            /*                s|
                               |

                               3
                               +-----*-----------+
                               |                 |
                               |                 |
                             | |                 |
                             v |        |        |             r
                            ---|        +--      |         -----
                             gp|                 |
                            ---*                 |
                             ^ |                 |
                             | |                 |
                               +-----------------+
                              0                   1
             */

            // r=-1
            /*
              parent               surface

             s|                           r|
              |                            |
              +                            +
             6|                           2|
              +                            +
             3|                           1|
              +                            +
             0                            0
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = -intpoints(iquad, 0);
            }
            break;
          }
          default:
            dserror("invalid number of lines, unable to determine intpoint in parent");
        }
      }
      else if ((distype == DRT::Element::line2 && pdistype == DRT::Element::tri3) or
               (distype == DRT::Element::line3 && pdistype == DRT::Element::tri6))
      {
        switch (beleid)
        {
          case 0:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.5 + 0.5 * intpoints(iquad, 0);
              pqxg(iquad, 1) = 0.0;
            }
            break;
          }
          case 1:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.5 - 0.5 * intpoints(iquad, 0);
              pqxg(iquad, 1) = 0.5 + 0.5 * intpoints(iquad, 0);  // r+s=1.0 -> s= 1.0-r
            }
            break;
          }
          case 2:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.0;
              pqxg(iquad, 1) = 0.5 - 0.5 * intpoints(iquad, 0);
            }
            break;
          }
          default:
            dserror("invalid number of lines, unable to determine intpoint in parent");
        }
      }
      else if (distype == DRT::Element::nurbs3 && pdistype == DRT::Element::nurbs9)
      {
        switch (beleid)
        {
          case 0:
          {
            // s=-1
            /*

            parent                line

            r                     r
            +---+---+  -----      +---+---+ ------
            0   1   2             0   1   2

            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = -1.0;
            }
            break;
          }
          case 1:
          {
            // r=+1
            /*
              parent               surface

              s|                    r|
               |                     |
               +                     +
              8|                    2|
               +                     +
              5|                    1|
               +                     +
              2                     0
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
            }
            break;
          }
          case 2:
          {
            // s=+1
            /*

            parent                line

            r                           r
            +---+---+  -----             +---+---+ -----
            6   7   8                    0   1   2

            */

            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = 1.0;
            }
            break;
          }
          case 3:
          {
            // r=-1
            /*
              parent               surface

             s|                           r|
              |                            |
              +                            +
             6|                           2|
              +                            +
             3|                           1|
              +                            +
             0                            0
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
            }
            break;
          }
          default:
            dserror("invalid number of lines, unable to determine intpoint in parent");
        }
      }
      else
      {
        dserror(
            "only line2/quad4 and nurbs3/nurbs9 mappings of surface gausspoint to parent element "
            "implemented up to now\n");
      }
    };
    /*-----------------------------------------------------------------
    \brief Transform Gausspoints on boundary element to space of
           parent element (required for integrations of parent-element
           shape functions over boundary elements, for example
           in weak dirichlet boundary conditions).
           3D-version
      -----------------------------------------------------------------*/
    inline void BoundaryGPToParentGP3(Epetra_SerialDenseMatrix& pqxg,
        const Epetra_SerialDenseMatrix intpoints, const DRT::Element::DiscretizationType pdistype,
        const DRT::Element::DiscretizationType distype, const int beleid)
    {
      // resize output array
      pqxg.Shape(intpoints.M(), 3);

      if ((distype == DRT::Element::quad4 && pdistype == DRT::Element::hex8) or
          (distype == DRT::Element::quad8 && pdistype == DRT::Element::hex20) or
          (distype == DRT::Element::quad9 && pdistype == DRT::Element::hex27))
      {
        switch (beleid)
        {
          case 0:
          {
            // t=-1
            /*
                    parent               surface

                     r|                    s|
                      |                     |
                 1         2           3         2
                  +-------+             +-------+
                  |   |   |      s      |   |   |      r
                  |   +-- |  -----      |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 0         3           0         1
            */

            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 1);
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = -1.0;
            }
            break;
          }
          case 1:
          {
            // s=-1
            /*
                    parent               surface
                     t|                    s|
                      |                     |
                 4         5           3         2
                  +-------+             +-------+
                  |   |   |      r      |   |   |      r
                  |   +-- |  -----      |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 0         1           0         1
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = -1.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 2:
          {
            // r= 1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                 5         6           3         2
                  +-------+             +-------+
                  |   |   |      s      |   |   |      r
                  |   +-- |  -----      |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 1         2           0         1
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 3:
          {
            // s= 1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                 6         7           3         2
                  +-------+             +-------+
            r     |   |   |             |   |   |      r
            ----  | --+   |             |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 2         3           0         1
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -intpoints(iquad, 0);
              pqxg(iquad, 1) = 1.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 4:
          {
            // r=-1
            /*
                    parent               surface

                     s|                    s|
                      |                     |
                 3         7           3         2
                  +-------+             +-------+
                  |   |   |      t      |   |   |      r
                  |   +-- |  -----      |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 0         4           0         1
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = intpoints(iquad, 0);
            }
            break;
          }
          case 5:
          {
            // t=1
            /*
                    parent               surface

                     s|                    s|
                      |                     |
                 7         6           3         2
                  +-------+             +-------+
                  |   |   |      r      |   |   |      r
                  |   +-- |  -----      |   +-- |  -----
                  |       |             |       |
                  +-------+             +-------+
                 4         5           0         1
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = 1.0;
            }
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else if ((distype == DRT::Element::tri3 && pdistype == DRT::Element::tet4) or
               (distype == DRT::Element::tri6 && pdistype == DRT::Element::tet10))
      {
        switch (beleid)
        {
          case 0:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = 0.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 1:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0 - intpoints(iquad, 0) - intpoints(iquad, 1);
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 2:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.0;
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = intpoints(iquad, 0);
            }
            break;
          }
          case 3:
          {
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 1);
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = 0.0;
            }
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else if ((distype == DRT::Element::quad4 && pdistype == DRT::Element::wedge6) or
               (distype == DRT::Element::quad8 && pdistype == DRT::Element::wedge15))
      {
        switch (beleid)
        {
          case 0:
          {
            // r+s=1
            /*
                          parent               surface

                       t|                        s|
                        |                         |
                       3|        4           3         2
                        +-------+             +-------+
                        |       |      s      |   |   |      r
                        +-------|  -----      |   +-- |  -----
                        |       |             |       |
                        +-------+             +-------+
                       0         1           0         1
             */

            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -0.5 * intpoints(iquad, 0) + 0.5;
              pqxg(iquad, 1) = 0.5 * intpoints(iquad, 0) + 0.5;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 1:
          {
            // r=0
            /*
                          parent               surface
                               t|                s|
                                |                 |
                       4        |5           3         2
                        +-------+             +-------+
                   s    |       |             |   |   |      r
                   -----|-------+             |   +-- |  -----
                        |       |             |       |
                        +-------+             +-------+
                       1         2           0         1
             */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.0;
              pqxg(iquad, 1) = -0.5 * intpoints(iquad, 0) + 0.5;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 2:
          {
            // s=0
            /*
                          parent               surface

                       t|                        s|
                        |                         |
                       5|        3           3    |    2
                        +-------+             +-------+
                        |       |      r      |   |   |      r
                        +       |  -----      |   +-- |  -----
                        |       |             |       |
                        +-------+             +-------+
                       2         0           0         1
             */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 0.5 * intpoints(iquad, 0) + 0.5;
              pqxg(iquad, 1) = 0.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          default:
            dserror(
                "invalid number of quad4 surface for wedge element, unable to determine intpoint "
                "in parent");
        }
      }
      else if ((distype == DRT::Element::tri3 && pdistype == DRT::Element::wedge6) or
               (distype == DRT::Element::tri6 && pdistype == DRT::Element::wedge15))
      {
        switch (beleid)
        {
          case 3:
          {
            // t=-1.0
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 1);
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = -1.0;
            }
            break;
          }
          case 4:
          {
            // t=1.0
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = 1.0;
            }
            break;
          }
          default:
            dserror(
                "invalid number of tri3 surface for wedge element, unable to determine intpoint in "
                "parent");
        }
      }
      else if (distype == DRT::Element::nurbs9 && pdistype == DRT::Element::nurbs27)
      {
        switch (beleid)
        {
          case 0:
          {
            // t=-1
            /*
                    parent               surface

                     s|                    s|
                      |                     |
                  +---+---+             +---+---+
                 6|  7|  8|      r     6|  7|  8|      r
                  +   +-- +  -----      +   +-- +  -----
                 3|  4   5|            3|  4   5|
                  +---+---+             +---+---+
                 0   1   2             0   1   2
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = -1.0;
            }
            break;
          }
          case 1:
          {
            // t=+1
            /*
                    parent               surface

                     s|                    s|
                      |                     |
                  +---+---+             +---+---+
                24| 25| 26|      r     6|  7|  8|      r
                  +   +-- +  -----      +   +-- +  -----
                21| 22  23|            3|  4   5|
                  +---+---+             +---+---+
                18  19  20             0   1   2
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = intpoints(iquad, 1);
              pqxg(iquad, 2) = 1.0;
            }
            break;
          }
          case 2:
          {
            // s=-1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                  +---+---+             +---+---+
                18| 19| 20|      r     6|  7|  8|      r
                  +   +-- +  -----      +   +-- +  -----
                 9| 10  11|            3|  4   5|
                  +---+---+             +---+---+
                 0   1   2             0   1   2
            */

            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = -1.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 3:
          {
            // s=+1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                  +---+---+             +---+---+
                24| 25| 26|    r       6|  7|  8|      r
                  +   +-- + ----        +   +-- +  -----
                15| 16  17|            3|  4   5|
                  +---+---+             +---+---+
                 6   7   8             0   1   2
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = intpoints(iquad, 0);
              pqxg(iquad, 1) = 1.0;
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 4:
          {
            // r=+1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                  +---+---+             +---+---+
                20| 23| 26|      s     6|  7|  8|      r
                  +   +-- +  -----      +   +-- +  -----
                11| 14  17|            3|  4   5|
                  +---+---+             +---+---+
                 2   5   8             0   1   2
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = 1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          case 5:
          {
            // r=-1
            /*
                    parent               surface

                     t|                    s|
                      |                     |
                  +---+---+             +---+---+
                18| 21| 24|      s     6|  7|  8|      r
                  +   +-- +  -----      +   +-- +  -----
                 9| 12  15|            3|  4   5|
                  +---+---+             +---+---+
                 0   3   6             0   1   2
            */
            for (int iquad = 0; iquad < intpoints.M(); ++iquad)
            {
              pqxg(iquad, 0) = -1.0;
              pqxg(iquad, 1) = intpoints(iquad, 0);
              pqxg(iquad, 2) = intpoints(iquad, 1);
            }
            break;
          }
          default:
            dserror("invalid number of surfaces, unable to determine intpoint in parent");
        }
      }
      else
      {
        dserror(
            "only quad4/hex8 and nurbs9/nurbs27 mappings of surface gausspoint to parent element "
            "implemented up to now\n");
      }
    };

    // Calculation of the following things for the boundary element (all outputs):
    // - Shape function 'funct' and its derivative 'deriv' at the Gauss point
    //   as well as the integration factor 'fac'
    // - Unit normal vector 'unitnormal' (pointing out of the domain) and the
    //   infinitesimal area element 'drs' at the Gauss points
    // - Coordinates of the current integration point in reference coordinates
    //   'xsi' and the node coordinates for the boundary element 'xyze'
    // Inputs:
    // - intpoints: integration points
    // - gpid: actual Gauss point
    // - myknots: nurbs specific knotes
    // - isnurbs: if true, some NURBS specific stuff is done
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    template <DRT::Element::DiscretizationType distype>
    void EvalShapeFuncAtBouIntPoint(
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>& funct,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv,
        double& fac, LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim + 1, 1>& unitnormal,
        double& drs, LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>& xsi,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim + 1,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& xyze,
        const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim>& intpoints,
        const int gpid, const std::vector<Epetra_SerialDenseVector>* myknots,
        const Epetra_SerialDenseVector* weights, const bool isnurbs)
    {
      // Obtain number of spatial dimensions of boundary element
      const int bdrynsd = DRT::UTILS::DisTypeToDim<distype>::dim;

      // local coordinates of the current integration point
      const double* gpcoord = (intpoints.IP().qxg)[gpid];
      for (int idim = 0; idim < bdrynsd; ++idim)
      {
        xsi(idim) = gpcoord[idim];
      }

      // get shape functions and derivatives in the plane of the element
      if (not isnurbs)
      {
        // shape functions and their first derivatives of boundary element
        DRT::UTILS::shape_function<distype>(xsi, funct);
        DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
      }
      // only for NURBS!!!
      else
      {
        DRT::NURBS::UTILS::nurbs_get_funct_deriv(funct, deriv, xsi, *myknots, *weights, distype);
      }

      // compute measure tensor for surface element, infinitesimal area element drs
      // and (outward-pointing) unit normal vector
      LINALG::Matrix<bdrynsd, bdrynsd> metrictensor(true);
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(
          xyze, deriv, metrictensor, drs, &unitnormal);

      // compute integration factor
      fac = intpoints.IP().qwgt[gpid] * drs;

      return;
    };  // EvalShapeFuncAtBouIntPoint

    // Calculate mass-consistent node normals for the boundary element
    // according to dissertation of Prof. Wall, equation (7.13)
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    template <DRT::Element::DiscretizationType distype>
    void ElementNodeNormal(
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>& funct,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv,
        double& fac, LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim + 1, 1>& unitnormal,
        double& drs, LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>& xsi,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim + 1,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& xyze,
        DRT::Element* ele, DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1,
        const std::vector<double>& edispnp, const bool isnurbs, const bool isale)
    {
      /*----------------------------------------------------------------------*
       |                          Initialization                              |
       *----------------------------------------------------------------------*/
      // Number of spatial dimensions of boundary element
      const int bdrynsd = DRT::UTILS::DisTypeToDim<distype>::dim;

      // Number of spatial dimensions of parent element
      static const int nsd = bdrynsd + 1;

      // Number of nodes of boundary element
      static const int bdrynen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // Number of degrees of freedom per node
      int numdofpernode = nsd + 1;  // standard case is 'fluid'
      if ((discretization.Name() == "ale") or (discretization.Name() == "structure"))
        numdofpernode = nsd;

      // This functionality has only been tested for types fluid, ale and structure
      if (not((discretization.Name() == "fluid") or (discretization.Name() == "inflow") or
              (discretization.Name() == "ale") or (discretization.Name() == "structure")))
        dserror(
            "ElementNodeNormal: The mass-consistent-node-normal calculation can currently only be "
            "performed for discretization types 'fluid', 'inflow', 'ale' and 'structure'.");

      // Get gaussrule
      const DRT::UTILS::IntPointsAndWeights<bdrynsd> intpoints(
          DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

      // Get node coordinates
      GEO::fillInitialPositionArray<distype, nsd, LINALG::Matrix<nsd, bdrynen>>(ele, xyze);

      // Add displacements to reference coordinates, if an ALE description is used
      if (isale)
      {
        dsassert(edispnp.size() != 0, "paranoid");

        for (int inode = 0; inode < bdrynen; ++inode)
        {
          for (int idim = 0; idim < (nsd); ++idim)
          {
            xyze(idim, inode) += edispnp[numdofpernode * inode + idim];
          }
        }
      }

      /*----------------------------------------------------------------------*
       |               start loop over integration points                     |
       *----------------------------------------------------------------------*/
      for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
      {
        // Computation of the integration factor & shape function at the Gauss point & derivative of
        // the shape function at the Gauss point Computation of the unit normal vector at the Gauss
        // points Computation of nurb specific stuff is not activated here
        DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(
            funct, deriv, fac, unitnormal, drs, xsi, xyze, intpoints, gpid, NULL, NULL, isnurbs);

        for (int inode = 0; inode < bdrynen; ++inode)
        {
          for (int idim = 0; idim < nsd; ++idim)
          {
            elevec1(inode * numdofpernode + idim) += unitnormal(idim) * funct(inode) * fac;
          }
          if (numdofpernode > nsd)
          {  // this is the case for the fluid
            // pressure dof is set to zero
            elevec1(inode * numdofpernode + (nsd)) = 0.0;
          }
        }
      } /* end of loop over integration points gpid */

      return;
    }  // ElementNodeNormal

  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifndef UTILS_BOUNDARY_INTEGRATION_H
