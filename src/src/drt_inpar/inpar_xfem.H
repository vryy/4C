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

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for XFEM
  namespace XFEM
  {

    /// which algorithm is used to enforce the boundary condition
    enum BoundaryIntegralType
    {
      BoundaryTypeSigma,
      BoundaryTypeTauPressure
    };

  }  // namespace FSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
#endif  // #ifdef CCADISCRET
