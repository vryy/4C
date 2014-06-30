/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.H

\brief Implementation class for coupling of two different meshes using Stress/Hybrid method of Nitsche's method to enforce
       interface conditions weakly
<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H
#define FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H

#include "../drt_inpar/inpar_xfem.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {

    /// Side implementation class (SideImpl) for XFEM interface coupling
    /*!
      This class provides a templated Side implementation class (SideImpl).
      Dependent on the template pairs between cut background element, the cutting interface Side and the number
      of nodal dofs the SideImpl object is created, which provides the implementation of coupling matrices for different
      coupling strategies like Nitsche's method or Mixed/Stress/Hybrid method.

      <h3>Usability</h3>

      The calculations are done by the EvaluateXfemInterface...() methods that call this SideInterface factory class.
      The SideInterface class creates the implementation template class.

      \author schott
      \date 04/12
     */
    template<DRT::Element::DiscretizationType distype,
    DRT::Element::DiscretizationType side_distype,
    const int side_numdof>

    class SideImpl : public SideInterface<distype>
    {
    public:

      //! number of boundary element nodes
      static const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;
      //! number of background element nodes (P. Hughes: The Finite Element Method)
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //! number of spatial dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
      //! number of fluid dof
      static const int bg_numdof_ = nsd_ + 1;
      //! number of independent stress-dof
      static const int numstressdof_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      //! Mixed/hybrid stress-based Lagrange multiplier (MHCS or MHVS) coupling between background element and side element
      SideImpl(
          DRT::Element *              side,            ///< side element
          Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
          Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
          Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
          Epetra_SerialDenseMatrix &  Gsui,            ///< interface \f$G_{\sigma u^i}\f$ coupling matrix
          Epetra_SerialDenseMatrix &  Guis,            ///< interface \f$G_{u^i \sigma}\f$ coupling matrix
          Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        G_sui_(Gsui.A(),true),
        G_uis_(Guis.A(),true),
        xyze_(side_xyze.A(),true)
      {
      }

      //! Nitsche (NIT) coupling between background element and side element
      SideImpl(
          DRT::Element *              side,            ///< side element
          Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
          Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
          Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
          Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
          Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        C_uiui_(C_uiui.A(),true),
        xyze_(side_xyze.A(),true)
      {
      }

      //! simple side coupling together with embedded Nitsche coupling
      SideImpl(
          DRT::Element *              side,            ///< side element
          Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
      )
      : xyze_(side_xyze.A(),true)
      {
      }

      //! add side's interface displacements and set current side node coordiantes
      virtual void addeidisp(
          const DRT::Discretization &  cutdis,       ///< cut discretization
          const std::string            state,        ///< state
          const std::vector<int>&      lm            ///< local map
          );

      //! extract/set side's interface velocity
      virtual void eivel(
          const DRT::Discretization &  cutdis,       ///< cut discretization
          const std::string            state,        ///< state
          const std::vector<int>&      lm            ///< local map
          );

      virtual void getivelint(
          LINALG::Matrix<nsd_,1>& ivelint  ///< interface velocity at embedded side
          );

      //! evaluate shape function, derivatives, normal and transformation w.r.t side element at gaussian point
      virtual void Evaluate(
          const LINALG::Matrix<nsd_-1,1>  & eta,     ///< local coordinates w.r.t side element
          LINALG::Matrix<nsd_,1>          & x,       ///< global coordinates of gaussian point
          LINALG::Matrix<nsd_,1>          & normal,  ///< normal vector
          double                          & drs      ///< transformation factor
      );

      //! set prescribed WDBC at Gaussian point
      virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & ivelint )
      {
        ivelint.Clear();
        ivelint.Multiply(eivel_,side_funct_);
      }

      //! compute interface force for side nodes
      virtual void InterfaceForce(
          Epetra_SerialDenseVector &   iforce,     ///< interface force vector
          LINALG::Matrix<nsd_,1> &     traction,   ///< traction vector at gaussian point
          const double &               fac         ///< integration factor
      );

      //! project gaussian point from linearized interfac in normal direction onto corresponding side
      virtual void ProjectOnSide(
          LINALG::Matrix<3,1> & x_gp_lin,            ///< global coordinates of gaussian point w.r.t linearized interface
          LINALG::Matrix<3,1> & x_side,              ///< projected gaussian point on side
          LINALG::Matrix<2,1> & xi_side              ///< local coordinates of projected gaussian point w.r.t side
      );

      //! evaluate interface coupling matrices for mixed/hybrid Cauchy stress-based (MHCS) coupling
      virtual void MHCS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
          const double &                                                fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                funct,      ///< shape function
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s       ///< block rhs vector \f$ rhs_{\sigma} \f$
      );

      //! evaluate interface matrices for mixed/hybrid viscous stress-based (MHVS) coupling
      //! \author kruse \date 06/14
      virtual void MHVS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                    normal,     ///< normal vector
          const double &                                                    fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                    funct,      ///< background element shape functions
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> &     rhs_s,      ///< block rhs vector \f$ rhs_{\sigma}\f$
          const double &                                                    press,      ///< background element pressure
          LINALG::Matrix<nen_,1> &                                          rhs_pui     ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
      );

      //! build the final coupling matrices for mixed/hybrid Cauchy or viscous stress-based coupling (MHCS or MHVS)
      //! \author kruse \date 06/14
      virtual void HybridLM_buildFinalCouplingMatrices(
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,numstressdof_> &     BinvK_ss,   ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,bg_numdof_,numstressdof_> &        BKusInvKss, ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,bg_numdof_> &        BK_su,      ///< block matrix \f$ K_{\sigma u} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>, numstressdof_,1>      &              rhs_s       ///< block rhs vector \f$ rhs_{\sigma}\f$
      );

      //! build stabilizing inflow terms
      virtual void HybridLM_Stab_InflowCoercivity(
          Epetra_SerialDenseMatrix &          C_uu_,                ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &          rhs_Cu_,              ///< standard bg-rhs
          const bool &                        coupling,             ///< indicates the need to assemble coupling terms for the side (yes/no)
          const bool &                        bg_mortaring,         ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
          const LINALG::Matrix<nsd_,1> &      normal,               ///< normal vector
          const double &                      timefacfac,           ///< theta*dt
          const double &                      visceff_1,            ///< viscosity in background fluid
          const double &                      visceff_2,            ///< viscosity in embedded fluid
          const double &                      kappa1,               ///< mortaring weighting
          const double &                      kappa2,               ///< mortaring weighting
          const double &                      stabfac,              ///< Nitsche penalty
          const double &                      stabfac_avg,          ///< Nitsche convective non-dimensionless stabilization factor
          const LINALG::Matrix<nen_,1> &      funct_,               ///< bg shape functions
          const LINALG::Matrix<nsd_,nen_> &   derxy_,               ///< bg deriv
          const LINALG::Matrix<nsd_,nsd_> &   vderxy_,              ///< bg deriv^n
          const LINALG::Matrix<nsd_,1> &      velint,               ///< bg u^n
          const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,    ///< Dirichlet velocity vector or prescribed jump vector
          INPAR::XFEM::ConvStabScaling        conv_stab_scaling,    ///< Inflow term strategies xfluid
          INPAR::XFEM::XFF_ConvStabScaling    xff_conv_stab_scaling ///< Inflow term strategies xfluidfluid
        );

      //! build coupling matrices and assemble terms for Nitsche's (NIT) method
      virtual void NIT_buildCouplingMatrices(
          Epetra_SerialDenseMatrix &          C_uu_,                  ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &          rhs_Cu_,                ///< standard bg-rhs
          const bool &                        coupling,               ///< assemble coupling terms (yes/no)
          const bool &                        bg_mortaring,           ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring, two-sided impl)
          const LINALG::Matrix<nsd_,1> &      normal,                 ///< normal vector
          const double &                      timefacfac,             ///< theta*dt
          const double &                      visceff_1,              ///< viscosity in background fluid
          const double &                      visceff_2,              ///< viscosity in embedded fluid
          const double &                      kappa1,                 ///< mortaring weighting
          const double &                      kappa2,                 ///< mortaring weighting
          const double &                      NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
          const double &                      avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
          const LINALG::Matrix<nen_,1> &      funct_,                 ///< bg shape functions
          const LINALG::Matrix<nsd_,nen_> &   derxy_,                 ///< bg deriv
          const LINALG::Matrix<nsd_,nsd_> &   vderxy_,                ///< bg deriv^n
          const double &                      press,                  ///< bg p^n
          const LINALG::Matrix<nsd_,1> &      velint,                 ///< bg u^n
          const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,      ///< Dirichlet velocity vector or prescribed jump vector
          INPAR::XFEM::ConvStabScaling        conv_stab_scaling,      ///< Inflow term strategies xfluid
          INPAR::XFEM::XFF_ConvStabScaling    xff_conv_stab_scaling   ///< Inflow term strategies xfluidfluid
      );

      //! evaluate pressure-consistency term for Nitsche's method
      virtual void NIT_p_Consistency(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const double &                           press,                            ///< pressure
          const LINALG::Matrix<nen_,1> &           funct_timefacfac_k1,              ///< funct * timefacfac *kappa1
          const LINALG::Matrix<side_nen_,1> &      side_funct_timefacfac_k1,         ///< sidefunct * timefacfac *kappa1
          const LINALG::Matrix<side_nen_,nen_> &   side_funct_dyad_k1_timefacfac,    ///< (sidefunct^T * funct) * timefacfac *kappa1
          const LINALG::Matrix<nen_,nen_> &        funct_dyad_k1_timefacfac,         ///< (funct^T * funct) * timefacfac *kappa1
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

      //! evaluate pressure-adjoint-consistency term for Nitsche's method
      virtual void NIT_p_AdjointConsistency(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          double &                                 velint_normal,                    ///< velocity in normal direction
          double &                                 ivelint_normal,                   ///< interface velocity in normal direction
          double &                                 ivelint_WDBC_JUMP_normal,         ///< prescribed interface velocity or jump vector in normal direction
          const LINALG::Matrix<nen_,1> &           funct_timefacfac_k1,              ///< funct * timefacfac *kappa1
          const LINALG::Matrix<side_nen_,nen_> &   side_funct_dyad_k1_timefacfac,    ///< (sidefunct^T * funct) * timefacfac *kappa1
          const LINALG::Matrix<nen_,nen_> &        funct_dyad_k1_timefacfac,         ///< (funct^T * funct) * timefacfac *kappa1
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

      //! evaluate viscous-consistency term for Nitsche's method
      virtual void NIT_visc_Consistency(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const LINALG::Matrix<nsd_,nen_> &        derxy,                            ///< bg deriv
          const LINALG::Matrix<nsd_,nsd_> &        vderxy,                           ///< bg deriv^n
          const LINALG::Matrix<nen_,1> &           e_funct_visc1_timefacfac,         ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<side_nen_,1> &      s_funct_visc1_timefacfac,         ///< side element funct *mu*timefacfac
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

      //! evaluate viscous-adjoint-consistency term for Nitsche's method
      virtual void NIT_visc_AdjointConsistency(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity
          const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity
          const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
          const LINALG::Matrix<nsd_,nen_> &        derxy,                            ///< bg deriv
          const double &                           visceff_1,                        ///< viscosity in background fluid element
          const double &                           timefacfac,                       ///< timefacfac
          const LINALG::Matrix<nen_,1> &           e_funct_visc1_timefacfac,         ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<side_nen_,1> &      s_funct_visc1_timefacfac,         ///< side element funct *mu*timefacfac
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

      //!  evaluate stabilizing viscous term for Nitsche's method
      virtual void NIT_Stab_ViscCoercivity(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
          const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
          const LINALG::Matrix<side_nen_,1>&       side_funct_timefacfac,            ///< sidefunct^T * timefacfac
          const LINALG::Matrix<side_nen_,nen_>&    side_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
          const double &                           stabfac,                          ///< stabilization factor
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

      virtual void NIT_Stab_ConvAveraged(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
          const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
          const LINALG::Matrix<side_nen_,1>&       side_funct_timefacfac,            ///< sidefunct^T * timefacfac
          const LINALG::Matrix<side_nen_,nen_>&    side_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
          const double &                           stabfac_avg,                      ///< stabilization factor
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );



    private:

      LINALG::Matrix<side_numdof*side_nen_,bg_numdof_*nen_>                   C_uiu_;       ///< coupling matrix row: sidenode1:u1,u2,u3(,p), sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<bg_numdof_*nen_,side_numdof*side_nen_>                   C_uui_;       ///< coupling matrix includes (ui,pi)
      LINALG::Matrix<side_numdof*side_nen_,1>                                 rhC_ui_;      ///< coupling matrix includes (ui,pi)
      LINALG::Matrix<side_numdof*side_nen_,side_numdof*side_nen_>             C_uiui_;      ///< includes (ui,pi) just for Nitsche coupling
      LINALG::Matrix<side_nen_,1>                                             side_funct_;  ///< side's shape functions
      LINALG::Matrix<nsd_-1,side_nen_>                                        side_deriv_;  ///< side's shape derivatives
      LINALG::Matrix<numstressdof_*nen_,side_numdof*side_nen_>                G_sui_;       ///< G_sui coupling matrix
      LINALG::Matrix<side_numdof*side_nen_,numstressdof_*nen_>                G_uis_;       ///< G_uis coupling matrix
      LINALG::Matrix<nsd_,side_nen_>                                          eivel_;       ///< side interface velocity
      LINALG::Matrix<side_nen_,1>                                             eipres_;      ///< side interface pressure
      LINALG::Matrix<nsd_,side_nen_>                                          eidisp_;      ///< side interface displacements
      LINALG::Matrix<nsd_,side_nen_>                                          xyze_;        ///< side's node coordinates
      LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,numstressdof_,nsd_>  BG_sui_;      ///< block G_uis matrix
      LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,nsd_,numstressdof_>  BG_uis_;      ///< block G_sui matrix
      LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,1,nsd_>              BG_pui_;      ///< block G_pui matrix
      LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,nsd_,1>              BG_uip_;      ///< block G_uip matrix

    };


    /// Embedded Element implementation class (EmbImpl) for XFEM interface coupling
    /*!
      This class provides a templated element implementation class (EmbImpl).
      Dependent on the template pairs between cut background element and the cutting embedded element
      the EmbImpl object is created, which provides the implementation of coupling matrices for different
      coupling strategies like Nitsche's method.

      <h3>Usability</h3>

      The calculations are done by the EvaluateXfemInterface...() methods that call this EmbCoupling factory class.
      The EmbCoupling class creates the implementation template class.

      \author schott
      \date 04/12
     */
    template<DRT::Element::DiscretizationType distype,
    DRT::Element::DiscretizationType emb_distype>
    class EmbImpl : public EmbCoupling<distype>
    {
    public:

      //! number of element nodes (P. Hughes: The Finite Element Method)
      static const int emb_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<emb_distype>::numNodePerElement;
      //! number of background element nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //! number of spatial dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
      //! number of fluid dof per node
      static const int bg_numdof_ = nsd_ + 1;


      //! Nitsche (NIT) coupling between background element and embedded element
      EmbImpl(
          DRT::Element *              emb_ele,         ///< side element
          Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
          Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
          Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
          Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
          Epetra_SerialDenseMatrix &  emb_xyze         ///< global node coordinates
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        C_uiui_(C_uiui.A(),true),
        emb_xyze_(emb_xyze.A(),true)
      {
      }

      //!  evaluate shape function, derivatives and transformation w.r.t embedded element at gaussian point
      virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside );

      //! extract/set embedded element velocity
      virtual void emb_vel(
          const DRT::Discretization &  embdis,       ///< embedded discretization
          const std::string            state,        ///< state
          const std::vector<int>&      lm            ///< local map
          );

      //! add embedded element displacements and set current element node coordinates
      virtual void addembdisp(
          const DRT::Discretization &  embdis,       ///< embedded discretization
          const std::string            state,        ///< state
          const std::vector<int>&      lm            ///< local map
          );

      virtual void getembpress(
        double& press
        );

      virtual void getembvelgradint(
        LINALG::Matrix<nsd_,nsd_>& velgradint    ///< interface velocity gradients at embedded side
        );

      //! compute embedded element's element length
      virtual void element_length( double & hk_emb );

      //! build coupling matrices and assemble terms for two-sided mortaring Nitsche's (NIT2) method
      virtual void NIT2_buildCouplingMatrices(
          Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
          Epetra_SerialDenseVector &    rhs_Cu_,        // standard bg-rhs
          bool &                        coupling,       // assemble coupling terms (yes/no)
          bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
          LINALG::Matrix<nsd_,1> &      normal,         // normal vector
          const double                  timefacfac,     // theta*dt
          const double                  visceff_1,      // viscosity in background fluid
          const double                  visceff_2,      // viscosity in embedded fluid
          double &                      kappa1,         // mortaring weighting
          double &                      kappa2,         // mortaring weighting
          const double &                NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
          const double &                avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
          bool &                        velgrad_interface_stab,// penalty term for velocity gradients at the interface
          double &                      velgrad_interface_fac, //stabilization fac for velocity gradients at the interface
          bool &                        presscoupling_interface_stab,// penalty term for pressure coupling at the interface
          double &                      presscoupling_interface_fac,//stabilization fac for pressure coupling at the interfac
          LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
          LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
          LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
          double &                      press,          // bg p^n
          LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
          LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP, // Dirichlet velocity vector or prescribed jump vector
          INPAR::XFEM::XFF_ConvStabScaling  xff_conv_stab_scaling // Inflow term strategies
        );

       virtual void NIT2_Stab_ConvAveraged(
          Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
          Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
          const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
          const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
          const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
          const LINALG::Matrix<emb_nen_,1>&        emb_funct_timefacfac,            ///< sidefunct^T * timefacfac
          const LINALG::Matrix<emb_nen_,nen_>&     emb_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
          const double &                           stabfac_avg,                          ///< stabilization factor
          const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
          const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
      );

    private:

      LINALG::Matrix<bg_numdof_*emb_nen_ ,bg_numdof_*nen_>      C_uiu_;       ///< coupling matrix row: sidenode1:u1,u2,u3,p, sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<bg_numdof_*nen_,bg_numdof_*emb_nen_>       C_uui_;       ///< coupling matrix includes (ui,pi)
      LINALG::Matrix<bg_numdof_*emb_nen_,1>                     rhC_ui_;      /// includes (ui,pi)
      LINALG::Matrix<bg_numdof_*emb_nen_,bg_numdof_*emb_nen_>   C_uiui_;      ///< coupling matrix includes (ui,pi)
      LINALG::Matrix<emb_nen_,1>                                emb_funct_;   ///< element's shape functions
      LINALG::Matrix<nsd_,emb_nen_>                             emb_deriv_;   ///< element's shape derivatives
      LINALG::Matrix<nsd_,nsd_>                                 emb_vderxy_;  ///< velocity derivatives
      LINALG::Matrix<nsd_,emb_nen_>                             emb_derxy_;   ///< element's global derivatives
      LINALG::Matrix<nsd_,emb_nen_>                             emb_vel_;     ///< element node velocity
      LINALG::Matrix<emb_nen_,1>                                emb_pres_;    ///< element node pressure
      LINALG::Matrix<nsd_,emb_nen_>                             emb_disp_;    ///< element node displacements
      LINALG::Matrix<nsd_,emb_nen_>                             emb_xyze_;    ///< element node coordinates

      LINALG::Matrix<nsd_,1>                                    emb_prederxy_; ///< pressure derivatives at GP

    };





    } // XFLUID
  } //ELEMENTS
} //DRT




#endif

