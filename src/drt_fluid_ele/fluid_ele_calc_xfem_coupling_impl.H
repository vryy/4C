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

      // for stress-based fluid-fluid coupling
      SideImpl(DRT::Element * side,
          Epetra_SerialDenseMatrix & C_uiu,
          Epetra_SerialDenseMatrix & C_uui,
          Epetra_SerialDenseMatrix & rhC_ui,
          Epetra_SerialDenseMatrix & Gsui,
          Epetra_SerialDenseMatrix & Guis,
          Epetra_SerialDenseMatrix & side_xyze
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        K_sui_(Gsui.A(),true),
        K_uis_(Guis.A(),true),
        xyze_(side_xyze.A(),true)
      {
      }

      // for Nitsche-based fluid-fluid coupling
      SideImpl(DRT::Element * side,
          Epetra_SerialDenseMatrix & C_uiu,
          Epetra_SerialDenseMatrix & C_uui,
          Epetra_SerialDenseMatrix & rhC_ui,
          Epetra_SerialDenseMatrix & C_uiui,
          Epetra_SerialDenseMatrix & side_xyze
      )
      : C_uiu_(C_uiu.A(),true),
        C_uui_(C_uui.A(),true),
        rhC_ui_(rhC_ui.A(),true),
        C_uiui_(C_uiui.A(),true),
        xyze_(side_xyze.A(),true)
      {
      }

      // without any coupling
      SideImpl(DRT::Element * side,
          Epetra_SerialDenseMatrix & side_xyze
      )
      : xyze_(side_xyze.A(),true)
      {
      }

      virtual void Evaluate(const LINALG::Matrix<nsd_-1,1>  & eta,
          LINALG::Matrix<nsd_,1>          & x,
          LINALG::Matrix<nsd_,1>          & normal,
          double                          & drs
      );


      virtual void ProjectOnSide( LINALG::Matrix<3,1> & x_gp_lin,
          LINALG::Matrix<3,1> & x_side,
          LINALG::Matrix<2,1> & xi_side
      );


      virtual void eivel(const DRT::Discretization &  cutdis,
          const std::string            state,
          const vector<int>&           lm);


      virtual void addeidisp(const DRT::Discretization &  cutdis,
          const std::string            state,
          const vector<int>&           lm,
          Epetra_SerialDenseMatrix  &  side_xyze);


      virtual void buildInterfaceForce( const Teuchos::RCP<Epetra_Vector> &   iforcecol,
          const DRT::Discretization &           cutdis,
          const vector<int>&                    lm,
          LINALG::Matrix<nsd_,1> &              traction,
          double &                              fac );


      virtual void buildCouplingMatrices(LINALG::Matrix<nsd_,1>        & normal,
          const double                fac,
          LINALG::Matrix<nen_,1>     & funct,
          LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
      );


      virtual void buildFinalCouplingMatrices(LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  BinvK_ss,
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,
          LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
      );


      virtual void buildCouplingMatricesNitsche(  Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
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
      virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & ivelint )
      {
        ivelint.Clear();
        ivelint.Multiply(eivel_,side_funct_);
      }


      LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,6,3> BK_sui_;
      LINALG::BlockMatrix<LINALG::Matrix<nen_,   1>,6,1>     rhsi_;
      LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,3,6> BK_uis_;

      LINALG::Matrix<numdof*side_nen_,4*nen_>            C_uiu_;    // row: sidenode1:u1,u2,u3(,p), sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<4*nen_,numdof*side_nen_>            C_uui_;    // includes ui(,pi)
      LINALG::Matrix<numdof*side_nen_,1>                 rhC_ui_;   // includes (ui,pi)
      LINALG::Matrix<numdof*side_nen_,numdof*side_nen_>  C_uiui_;   // includes (ui,pi) only for Nitsche coupling
      LINALG::Matrix<side_nen_,1>                        side_funct_;
      LINALG::Matrix<nsd_-1,side_nen_>                   side_deriv_;
      LINALG::Matrix<6*nen_,numdof*side_nen_>            K_sui_;
      LINALG::Matrix<numdof*side_nen_,6*nen_>            K_uis_;    // G_uis
      LINALG::Matrix<3,side_nen_>                        eivel_;
      LINALG::Matrix<side_nen_,1>                        eipres_;
      LINALG::Matrix<3,side_nen_>                        eidisp_;
      LINALG::Matrix<3,side_nen_>                        xyze_;
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


      // for Nitsche-based fluid-fluid coupling
      EmbImpl(DRT::Element * emb_ele,
          Epetra_SerialDenseMatrix & C_uiu,
          Epetra_SerialDenseMatrix & C_uui,
          Epetra_SerialDenseMatrix & rhC_ui,
          Epetra_SerialDenseMatrix & C_uiui,
          Epetra_SerialDenseMatrix & emb_xyze
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


      virtual void buildCouplingMatricesNitscheTwoSided(
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

