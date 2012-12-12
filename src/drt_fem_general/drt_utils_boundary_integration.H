/*!
\file drt_utils_boundary_integration.H

\brief methods for the integration over boundary elements

1) computation of kovariant metric tensor for surface element
2) mapping of gausspoints on surface element to 3d space of parent element
   (required for integrations of parent-element shape functions
    over boundary elements, required for example in weak
    dirichlet boundary conditions).

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/

#ifndef UTILS_BOUNDARY_INTEGRATION_H
#define UTILS_BOUNDARY_INTEGRATION_H


#include "drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_utils_integration.H"

namespace DRT
{
namespace UTILS
{
  class IntegrationPoints1D;
  class IntegrationPoints2D;

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
void SurfaceGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  Epetra_SerialDenseMatrix              & derivtrafo,
  const DRT::UTILS::IntegrationPoints2D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               surfaceid);

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
void LineGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  Epetra_SerialDenseMatrix              & derivtrafo,
  const DRT::UTILS::IntegrationPoints1D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               lineid   );

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
template<const int NSD>
void BoundaryGPToParentGP(
  Epetra_SerialDenseMatrix                                                    & pqxg     ,
  Epetra_SerialDenseMatrix                                                    & derivtrafo,
  const DRT::UTILS::IntPointsAndWeights<NSD-1> &                                intpoints,
  const DRT::Element::DiscretizationType                                        pdistype ,
  const DRT::Element::DiscretizationType                                        distype  ,
  const int                                                                     surfaceid);

//! specialization for 3D
template<>
void BoundaryGPToParentGP<3>(
Epetra_SerialDenseMatrix                                                    & pqxg     ,
Epetra_SerialDenseMatrix                                                    & derivtrafo,
const DRT::UTILS::IntPointsAndWeights<2> &                                    intpoints,
const DRT::Element::DiscretizationType                                        pdistype ,
const DRT::Element::DiscretizationType                                        distype  ,
const int                                                                     surfaceid);

//! specialization for 2D
template<>
void BoundaryGPToParentGP<2>(
Epetra_SerialDenseMatrix                                                    & pqxg     ,
Epetra_SerialDenseMatrix                                                    & derivtrafo,
const DRT::UTILS::IntPointsAndWeights<1> &                                    intpoints,
const DRT::Element::DiscretizationType                                        pdistype ,
const DRT::Element::DiscretizationType                                        distype  ,
const int                                                                     surfaceid);

/*----------------------------------------------------------------------------------------*/

//! compute kovariant metric tensor for surface element
void ComputeMetricTensorForSurface(
  const Epetra_SerialDenseMatrix& xyze,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix&       metrictensor,
  double*                         sqrtdetg
  );

//! compute kovariant metric tensor for surface/line element and optionally, the normalized normal vector at the Gausspoint (template)
template<DRT::Element::DiscretizationType DISTYPE>
void ComputeMetricTensorForBoundaryEle(
  const LINALG::Matrix< (1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim), DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  xyze,
  const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,  DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  deriv,
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,        DRT::UTILS::DisTypeToDim<DISTYPE>::dim>&                    metrictensor,
  double&  sqrtdetg,
  LINALG::Matrix<(1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim),        1>* normalvec = NULL
  )
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
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim> dxyzdrs(true);
  dxyzdrs.MultiplyNT(1.0,deriv,xyze,0.0);

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
  metrictensor.MultiplyNT(1.0,dxyzdrs,dxyzdrs,0.0);

  /*
                            +--------------+
                           /               |
             sqrtdetg =   /  g11*g22-g12^2
                        \/
  */
  sqrtdetg = sqrt(metrictensor.Determinant());

  // Calculate outward pointing unit normal vector
  if (normalvec != NULL)
  {
    if (DRT::UTILS::DisTypeToDim<DISTYPE>::dim == 2)
    {
      (*normalvec)(0) = dxyzdrs(0,1)*dxyzdrs(1,2)-dxyzdrs(1,1)*dxyzdrs(0,2);
      (*normalvec)(1) = dxyzdrs(0,2)*dxyzdrs(1,0)-dxyzdrs(1,2)*dxyzdrs(0,0);
      (*normalvec)(2) = dxyzdrs(0,0)*dxyzdrs(1,1)-dxyzdrs(1,0)*dxyzdrs(0,1);
      // compute unit normal (outward pointing)
      const double norm2 = normalvec->Norm2();
      normalvec->Scale(1/norm2);

    }
    else if (DRT::UTILS::DisTypeToDim<DISTYPE>::dim == 1)
    {
      (*normalvec)(0) = dxyzdrs(0,1);
      (*normalvec)(1) = -dxyzdrs(0,0);
      // compute unit normal (outward pointing)
      const double norm2 = normalvec->Norm2();
      normalvec->Scale(1/norm2);
    }
    else
      dserror("There are only 2D and 1D boundary elements");
  }

  return;
}

/*-----------------------------------------------------------------

\brief Transform Gausspoints on boundary element to space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, for example
       in weak dirichlet boundary conditions).
       2D-version
  -----------------------------------------------------------------*/
inline void BoundaryGPToParentGP2(
  Epetra_SerialDenseMatrix                   & pqxg     ,
  const Epetra_SerialDenseMatrix               intpoints,
  const DRT::Element::DiscretizationType       pdistype ,
  const DRT::Element::DiscretizationType       distype  ,
  const int                                    beleid   )
{
    // resize output array
    pqxg.Shape(intpoints.M(),2);

    if(distype==DRT::Element::line2 && pdistype==DRT::Element::quad4)
    {
      switch(beleid)
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)=-1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= 1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
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

        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-intpoints(iquad,0);
          pqxg(iquad,1)= 1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-1.0;
          pqxg(iquad,1)=-intpoints(iquad,0);
        }
        break;
      }
      default:
        dserror("invalid number of lines, unable to determine intpoint in parent");
      }

    }
    else if(distype==DRT::Element::nurbs3 && pdistype==DRT::Element::nurbs9)
    {
      switch(beleid)
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)=-1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= 1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
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

        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)= 1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
        }
        break;
      }
      default:
        dserror("invalid number of lines, unable to determine intpoint in parent");
      }
    }
    else
    {
      dserror("only line2/quad4 and nurbs3/nurbs9 mappings of surface gausspoint to parent element implemented up to now\n");
    }
};
/*-----------------------------------------------------------------
\brief Transform Gausspoints on boundary element to space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, for example
       in weak dirichlet boundary conditions).
       3D-version
  -----------------------------------------------------------------*/
inline void BoundaryGPToParentGP3(
  Epetra_SerialDenseMatrix                   & pqxg     ,
  const Epetra_SerialDenseMatrix               intpoints,
  const DRT::Element::DiscretizationType       pdistype ,
  const DRT::Element::DiscretizationType       distype  ,
  const int                                    beleid   )
{

    // resize output array
    pqxg.Shape(intpoints.M(),3);

    if(   (distype==DRT::Element::quad4 && pdistype==DRT::Element::hex8)
       or (distype==DRT::Element::quad8 && pdistype==DRT::Element::hex20)) //schott 05/11
    {
      switch(beleid)
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

        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,1);
          pqxg(iquad,1)= intpoints(iquad,0);
          pqxg(iquad,2)=-1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)=-1.0;
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= 1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-intpoints(iquad,0);
          pqxg(iquad,1)= 1.0;
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-1.0;
          pqxg(iquad,1)= intpoints(iquad,1);
          pqxg(iquad,2)= intpoints(iquad,0);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)= intpoints(iquad,1);
          pqxg(iquad,2)= 1.0;
        }
        break;
      }
      default:
        dserror("invalid number of surfaces, unable to determine intpoint in parent");
      }
    }
    else if(distype==DRT::Element::nurbs9 && pdistype==DRT::Element::nurbs27)
    {
      switch(beleid)
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)= intpoints(iquad,1);
          pqxg(iquad,2)=-1.0;
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)= intpoints(iquad,1);
          pqxg(iquad,2)= 1.0;
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

        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)=-1.0;
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= intpoints(iquad,0);
          pqxg(iquad,1)= 1.0;
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)= 1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
          pqxg(iquad,2)= intpoints(iquad,1);
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
        for (int iquad=0;iquad<intpoints.M();++iquad)
        {
          pqxg(iquad,0)=-1.0;
          pqxg(iquad,1)= intpoints(iquad,0);
          pqxg(iquad,2)= intpoints(iquad,1);
        }
        break;
      }
      default:
        dserror("invalid number of surfaces, unable to determine intpoint in parent");
      }
    }
    else
    {
      dserror("only quad4/hex8 and nurbs9/nurbs27 mappings of surface gausspoint to parent element implemented up to now\n");
    }
};



} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_BOUNDARY_INTEGRATION_H
