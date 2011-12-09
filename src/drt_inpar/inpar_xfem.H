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

    enum MovingBoundary
    {
      XFluidMovingBoundary,
      XFluidStationaryBoundary,
      XFSIMovingBoundary
    };

    /// how to define the interface (solid) velocity
    enum InterfaceInitVel
    {
      interface_vel_init_by_funct,
      interface_vel_init_zero
    };

    /// how to define the interface (solid) velocity
    enum InterfaceVel
    {
      interface_vel_by_disp,
      interface_vel_by_funct,
      interface_vel_zero
    };

    /// how to prescribe the interface (solid) displacement
    enum InterfaceDisplacement
    {
      interface_disp_by_fsi,
      interface_disp_by_funct,
      interface_disp_zero
    };

    /// which sided-mortaring is used to enforce the boundary condition
    enum CouplingStrategy
    {
      Xfluid_Sided_Mortaring,
      Embedded_Sided_Mortaring,
      Two_Sided_Mortaring
    };

    // xfluidfluid-fsi-monolithic approach
    enum Monolithic_xffsi_Approach
    {
      XFFSI_Full_Newton,
      XFFSI_FixedALE_Interpolation,
      XFFSI_FixedALE_Partitioned
    };


  }  // namespace XFEM

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
#endif  // #ifdef CCADISCRET
