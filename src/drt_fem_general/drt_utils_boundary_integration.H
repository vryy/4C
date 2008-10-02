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


namespace DRT
{
namespace UTILS
{

//! compute kovariant metric tensor for surface element
void ComputeMetricTensorForSurface(
  const Epetra_SerialDenseMatrix& xyze,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix&       metrictensor,
  double                         *sqrtdetg
  );

} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_BOUNDARY_INTEGRATION_H
#endif  // #ifdef CCADISCRET
