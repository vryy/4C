/*!
\file drt_utils_boundary_integration.H

\brief methods for the integration over boundary elements

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


namespace DRT
{
namespace UTILS
{

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
