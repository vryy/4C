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

#ifdef CCADISCRET
#ifndef UTILS_BOUNDARY_INTEGRATION_H
#define UTILS_BOUNDARY_INTEGRATION_H

#include <Epetra_SerialDenseMatrix.h>
#include "../drt_lib/drt_element.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
namespace UTILS
{
/*-----------------------------------------------------------------

\brief Transform Gausspoints on surface element to 3d space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, required for example
       in weak dirichlet boundary conditions).


\param pqxg      (o) transformed integration points of surface
                     elements in parent elements reference
                     coordinates
\param intpoints (i) integration points of surface element in its
                     reference coordinates
\param pdistype  (i) discretisation type of parent element
\param distype   (i) discretisation type of surface element
\param surfaceid (i) local id of surface element

\author gammi (gamnitzer@lnm.mw.tum.de), 05/09

  -----------------------------------------------------------------*/
void SurfaceGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  const DRT::UTILS::IntegrationPoints2D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               surfaceid);

/*-----------------------------------------------------------------

\brief Transform Gausspoints on line element to 2d space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, for example
       in weak dirichlet boundary conditions).


\param pqxg      (o) transformed integration points of surface
                     elements in parent elements reference
                     coordinates
\param intpoints (i) integration points of surface element in its
                     reference coordinates
\param pdistype  (i) discretisation type of parent element
\param distype   (i) discretisation type of surface element
\param lineid    (i) local id of surface element

\author gammi (gamnitzer@lnm.mw.tum.de), 05/09

  -----------------------------------------------------------------*/
void LineGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  const DRT::UTILS::IntegrationPoints1D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               lineid   );


//! compute kovariant metric tensor for surface element
void ComputeMetricTensorForSurface(
  const Epetra_SerialDenseMatrix& xyze,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix&       metrictensor,
  double*                         sqrtdetg
  );

//! compute kovariant metric tensor for surface/line element (template)
template<DRT::Element::DiscretizationType DISTYPE>
void ComputeMetricTensorForBoundaryEle(
  const LINALG::Matrix< (1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim), DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  xyze,
  const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,  DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  deriv,
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,        DRT::UTILS::DisTypeToDim<DISTYPE>::dim>&                    metrictensor,
  double&  sqrtdetg
  )
{
  /*
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
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,1 + DRT::UTILS::DisTypeToDim<DISTYPE>::dim> dxyzdrs(true);
  dxyzdrs.MultiplyNT(1.0,deriv,xyze,0.0);

  /*
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
  metrictensor.Clear();
  metrictensor.MultiplyNT(1.0,dxyzdrs,dxyzdrs,0.0);

  /*
                            +--------------+
                           /               |
             sqrtdetg =   /  g11*g22-g12^2
                        \/
  */
  sqrtdetg = sqrt(metrictensor.Determinant());

  return;
}

} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_BOUNDARY_INTEGRATION_H
#endif  // #ifdef CCADISCRET
