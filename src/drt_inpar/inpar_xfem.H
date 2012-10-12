/*----------------------------------------------------------------------*/
/*!
\file inpar_xfem.H

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
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
      BoundaryTypeNitsche,
      BoundaryTypeNeumann
    };

    /// which type of boundary movement
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
      interface_disp_zero,
      interface_disp_by_implementation
    };

    /// which sided-mortaring is used to enforce the boundary condition
    enum CouplingStrategy
    {
      Xfluid_Sided_Mortaring,
      Embedded_Sided_Mortaring,
      Two_Sided_Mortaring
    };

    /// L2 stress projection on whole fluid element or on partial fluid element
    enum MSH_L2_Proj
    {
      MSH_L2_Proj_full,
      MSH_L2_Proj_part
    };

    /// scaling of viscous interface stabilization
    enum ViscStabScaling
    {
      ViscStabScaling_visc_div_by_hk,
      ViscStabScaling_inv_hk,
      ViscStabScaling_const
    };

    /// definition of characteristic element length in cut elements
    enum ViscStab_hk
    {
      ViscStab_hk_vol_equivalent,
      ViscStab_hk_vol_div_by_surf,
      ViscStab_hk_longest_ele_length
    };

    /// scaling for convective/inflow stabilization term
    enum ConvStabScaling
    {
      ConvStabScaling_inflow,
      ConvStabScaling_abs_normal_vel,
      ConvStabScaling_max_abs_normal_vel,
      ConvStabScaling_const,
      ConvStabScaling_averaged,
      ConvStabScaling_none
    };

    /// EdgeBased (EOS) and Ghost Penalty matrix pattern
    enum EOS_GP_Pattern
    {
      EOS_GP_Pattern_uvwp,
      EOS_GP_Pattern_up,
      EOS_GP_Pattern_full
    };

    /// xfluidfluid-fsi-monolithic approach
    enum Monolithic_xffsi_Approach
    {
      XFFSI_Full_Newton,
      XFFSI_FixedALE_Interpolation,
      XFFSI_FixedALE_Partitioned
    };

    /// xfluidfluid time integration approach
    enum XFluidFluidTimeInt
    {
      Xff_TimeInt_FullProj,
      Xff_TimeInt_ProjIfMoved,
      Xff_TimeInt_KeepGhostValues,
      Xff_TimeInt_IncompProj
    };

    /// xfluid time integration approach
    enum XFluidTimeInt
    {
      Xf_TimeInt_SemiLagrange,
      Xf_TimeInt_GhostPenalty,
      Xf_TimeInt_Copy
    };

  }  // namespace XFEM

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
