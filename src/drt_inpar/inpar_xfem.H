/*----------------------------------------------------------------------*/
/*!
\file inpar_xfem.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_XFEM_H
#define INPAR_XFEM_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for XFEM
  namespace XFEM
  {

    /// which method is used to enforce the boundary condition
    enum BoundaryIntegralType
    {
      BoundaryTypeSigma,
      BoundaryTypeTauPressure,
      BoundaryTypeNitsche,
      BoundaryTypeNeumann
    };

    /// which sided-mortaring is used to enforce the boundary condition
    enum CouplingStrategy
    {
      Xfluid_Sided_Mortaring,
      Embedded_Sided_Mortaring,
      Two_Sided_Mortaring
    };

  }  // namespace XFEM

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
#endif  // #ifdef CCADISCRET
