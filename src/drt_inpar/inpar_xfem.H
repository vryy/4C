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
    enum CouplingMethod
    {
      Hybrid_LM_Cauchy_stress,
      Hybrid_LM_viscous_stress,
      Nitsche
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

    /// enforce weak DBC or coupling conditions with average w.r.t which side?
    enum CouplingStrategy
    {
      Xfluid_Sided_weak_DBC,     /// weak Dirichlet boundary conditions w.r.t fixed background mesh using XFEM
      Xfluid_Sided_Coupling,     /// monolithic coupling between fluid and (fluid/structure), average w.r.t xfluid side
      Embedded_Sided_Coupling,   /// monolithic coupling between fluid and (fluid/structure), average w.r.t embedded fluid/structure side
      Two_Sided_Coupling         /// monolithic coupling between fluid and (fluid/structure), mean average between fluid and fluid/structure
    };

    /// L2 stress projection on whole fluid element or on partial fluid element
    enum Hybrid_LM_L2_Proj
    {
      Hybrid_LM_L2_Proj_full,
      Hybrid_LM_L2_Proj_part
    };

    /// estimate of the scaling of the trace inequality for viscous interface stabilization (Nitsche's method)
    enum ViscStab_TraceEstimate
    {
      ViscStab_TraceEstimate_CT_div_by_hk,  /// estimate of the trace-inequality by trace-constant divided by a characteristic element length
      ViscStab_TraceEstimate_eigenvalue     /// estimate of the trace-inequality by solving an eigenvalue problem
    };

    /// definition of characteristic element length in cut elements
    enum ViscStab_hk
    {
      ViscStab_hk_vol_equivalent,             /// volume equivalent element diameter
      ViscStab_hk_cut_vol_div_by_cut_surf,    /// physical partial/cut volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, not useful for sliver and/or dotted cut situations)
      ViscStab_hk_ele_vol_div_by_cut_surf,    /// full element volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, however, avoids problems with sliver cuts, not useful for dotted cuts)
      ViscStab_hk_ele_vol_div_by_ele_surf,    /// full element volume divided by surface measure ( used for uncut situations, standard weak Dirichlet boundary/coupling conditions)
      ViscStab_hk_ele_vol_div_by_max_ele_surf /// full element volume divided by maximal element surface measure ( used to estimate the trace inequality for stretched elements in combination with ghost-penalties)
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

    /// Xfluidfluid scaling for convective/inflow stabilization term
    enum XFF_ConvStabScaling
    {
      XFF_ConvStabScaling_onesidedinflow,
      XFF_ConvStabScaling_averaged,
      XFF_ConvStabScaling_onesidedinflow_max_penalty,
      XFF_ConvStabScaling_averaged_max_penalty,
      XFF_ConvStabScaling_none
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

    /// xfluid time integration technique
    enum XFluidTimeIntScheme
    {
      //! std-DOFS:   only copying for std-dofs allowed (no semi-lagrangean algo (SL))
      //! ghost-DOFS: copy dofs or ghost-penalty reconstruction
      Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP,
      //! std-DOFS:   copying for std-dofs preferred and semi-lagrangean algo (SL) used for large displacements
      //! ghost-DOFS: copy dofs or ghost-penalty reconstruction
      Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP,
      //! std-DOFS:   semi-lagrangean algo (SL) used for large displacements in the whole cut-zone
      //! ghost-DOFS: only ghost-penalty reconstruction
      Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP
    };

    /// xfluid time integration approach for single dofs
    enum XFluidTimeInt
    {
      Xf_TimeInt_STD_by_SL,                 //! Semi-lagrangean algorithm for standard dof
      Xf_TimeInt_STD_by_COPY_from_STD,      //! copy value for std-dof at t^(n+1) from std-dof at t^n
      Xf_TimeInt_STD_by_COPY_from_GHOST,    //! copy value for std-dof at t^(n+1) from ghost-dof at t^n
      Xf_TimeInt_GHOST_by_GP,               //! Ghost-penalty reconstruction for ghost dof
      Xf_TimeInt_GHOST_by_COPY_from_STD,    //! copy value for ghost-dof at t^(n+1) from std-dof at t^n
      Xf_TimeInt_GHOST_by_COPY_from_GHOST,  //! copy value for ghost-dof at t^(n+1) from ghost-dof at t^n
      Xf_TimeInt_undefined
    };

    /// type of face
    enum FaceType
    {
      face_type_std,
      face_type_ghost_penalty,
      face_type_ghost,
      face_type_double_plus, // combustion only
      face_type_double_minus // combustion only
    };

  }  // namespace XFEM

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
