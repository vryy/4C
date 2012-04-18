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


namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {


    template<DRT::Element::DiscretizationType distype,
    DRT::Element::DiscretizationType side_distype,
    const int numdof>

    class SideImpl : public SideInterface<distype>
    {
    public:
      SideImpl();

      //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
      static const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! Mixed/Stress/Hybrid (MSH) coupling
      SideImpl(
          DRT::Element *              side,            ///< side element
          Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
          Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
          Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
          Epetra_SerialDenseMatrix &  Gsui,            ///< interface sigma-u_interface coupling
          Epetra_SerialDenseMatrix &  Guis,            ///< interface u_interface-sigma coupling
          Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        K_sui_(Gsui.A(),true),
        K_uis_(Guis.A(),true),
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
          const vector<int>&           lm,           ///< local map
          Epetra_SerialDenseMatrix  &  side_xyze     ///< side node coordinates
          );

      //! extract/set side's interface velocity
      virtual void eivel(
          const DRT::Discretization &  cutdis,       ///< cut discretization
          const std::string            state,        ///< state
          const vector<int>&           lm            ///< local map
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

      //! build coupling matrices for Mixed/Stress/Hybrid (MSH) method
      virtual void MSH_buildCouplingMatrices(
          LINALG::Matrix<nsd_,1> &                             normal,     ///< normal vector
          const double                                         fac,        ///< integration factor
          LINALG::Matrix<nen_,1> &                             funct,      ///< shape function
          LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs         ///< rhs block matrix
      );

      //! build final coupling matrices for Mixed/Stress/Hybrid (MSH) method
      virtual void MSH_buildFinalCouplingMatrices(
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  BinvK_ss,     ///< block inverse K_sigma_sigma matrix
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,         ///< block K_iK matrix
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,         ///< block K_su matrix
          LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &   rhs           ///< block rhs vector
      );

      //! build coupling matrices and assemble terms for Nitsche's (NIT) method
      virtual void NIT_buildCouplingMatrices(
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
          double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
          double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
          LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
          LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
          LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
          double &                      press,          // bg p^n
          LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
          LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
      );




      LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,6,3> BK_sui_;  ///< block K_sui matrix
      LINALG::BlockMatrix<LINALG::Matrix<nen_,   1>,6,1>      rhsi_;    ///< block rhs ui matrix
      LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,3,6> BK_uis_;  ///< block K_sui_ matrix

      LINALG::Matrix<numdof*side_nen_,4*nen_>            C_uiu_;       ///< coupling matrix row: sidenode1:u1,u2,u3(,p), sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<4*nen_,numdof*side_nen_>            C_uui_;       ///< coupling matrix includes ui(,pi)
      LINALG::Matrix<numdof*side_nen_,1>                 rhC_ui_;      ///< coupling matrix includes (ui,pi)
      LINALG::Matrix<numdof*side_nen_,numdof*side_nen_>  C_uiui_;      ///< includes (ui,pi) just for Nitsche coupling
      LINALG::Matrix<side_nen_,1>                        side_funct_;  ///< side's shape functions
      LINALG::Matrix<nsd_-1,side_nen_>                   side_deriv_;  ///< side's shape derivatives
      LINALG::Matrix<6*nen_,numdof*side_nen_>            K_sui_;       ///< G_sui coupling matrix
      LINALG::Matrix<numdof*side_nen_,6*nen_>            K_uis_;       ///< G_uis coupling matrix
      LINALG::Matrix<3,side_nen_>                        eivel_;       ///< side interface velocity
      LINALG::Matrix<side_nen_,1>                        eipres_;      ///< side interface pressure
      LINALG::Matrix<3,side_nen_>                        eidisp_;      ///< side interface displacements
      LINALG::Matrix<3,side_nen_>                        xyze_;        ///< side's node coordinates
    };



    template<DRT::Element::DiscretizationType distype,
    DRT::Element::DiscretizationType emb_distype>
    class EmbImpl : public EmbCoupling<distype>
    {
    public:

      //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
      static const int emb_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<emb_distype>::numNodePerElement;
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;


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


      virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside );


      virtual void emb_vel(const DRT::Discretization &  embdis,
          const std::string            state,
          const vector<int>&           lm);


      virtual void addembdisp(const DRT::Discretization &  embdis,
          const std::string            state,
          const vector<int>&           lm,
          Epetra_SerialDenseMatrix  &  side_xyze);


      virtual void element_length( double & hk_emb );


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
          double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
          double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
          LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
          LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
          LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
          double &                      press,          // bg p^n
          LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
          LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
      );


      // set prescribed WDBC at Gaussian point
      virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & emb_velint )
      {
        emb_velint.Clear();
        emb_velint.Multiply(emb_vel_,emb_funct_);
      }


      LINALG::Matrix<4*emb_nen_ ,4*nen_>     C_uiu_;    // row: sidenode1:u1,u2,u3,p, sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<4*nen_,4*emb_nen_>      C_uui_;    // includes (ui,pi)
      LINALG::Matrix<4*emb_nen_,1>           rhC_ui_;   // includes (ui,pi)
      LINALG::Matrix<4*emb_nen_,4*emb_nen_>  C_uiui_;   // includes (ui,pi) only for Nitsche coupling
      LINALG::Matrix<emb_nen_,1>             emb_funct_;
      LINALG::Matrix<nsd_,emb_nen_>          emb_deriv_;
      LINALG::Matrix<nsd_,nsd_>              emb_vderxy_;
      LINALG::Matrix<nsd_,emb_nen_>          emb_derxy_;
      LINALG::Matrix<3,emb_nen_>             emb_vel_;
      LINALG::Matrix<emb_nen_,1>             emb_pres_;
      LINALG::Matrix<3,emb_nen_>             emb_disp_;
      LINALG::Matrix<3,emb_nen_>             emb_xyze_;
    };





    } // XFLUID
  } //ELEMENTS
} //DRT




#endif

