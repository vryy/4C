/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling.H

\brief Interface class for coupling of sides and elements from two different meshes

<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_COUPLING_H
#define FLUID_ELE_CALC_XFEM_COUPLING_H


#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../linalg/linalg_fixedsizeblockmatrix.H"




namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {

      template<DRT::Element::DiscretizationType distype>
      class SideInterface
      {
      public:
        SideInterface(){};

        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & C_uiu,
                                                          Epetra_SerialDenseMatrix & C_uui,
                                                          Epetra_SerialDenseMatrix & rhC_ui,
                                                          Epetra_SerialDenseMatrix & Gsui,
                                                          Epetra_SerialDenseMatrix & Guis,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & C_uiu,
                                                          Epetra_SerialDenseMatrix & C_uui,
                                                          Epetra_SerialDenseMatrix & rhC_ui,
                                                          Epetra_SerialDenseMatrix & C_uiui,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        virtual ~SideInterface(){}

        virtual void Evaluate(const LINALG::Matrix<2,1> & eta,
                              LINALG::Matrix<3,1> &       x,
                              LINALG::Matrix<3,1> &       normal,
                              double &                    drs
                              ) = 0;


        virtual void ProjectOnSide(LINALG::Matrix<3,1> & x_gp_lin,
                                   LINALG::Matrix<3,1> & x_side,
                                   LINALG::Matrix<2,1> & xi_side
                                   ) = 0;


        virtual void eivel(const DRT::Discretization &  cutdis,
                           const std::string            state,
                           const vector<int>&           lm
                           ) = 0;


        virtual void addeidisp(const DRT::Discretization &  cutdis,
                               const std::string            state,
                               const vector<int>&           lm,
                               Epetra_SerialDenseMatrix  &  side_xyze
                               ) = 0;


        virtual void buildInterfaceForce( const Teuchos::RCP<Epetra_Vector> &  iforcecol,
                                          const DRT::Discretization &          cutdis,
                                          const vector<int>&                   lm,
                                          LINALG::Matrix<nsd_,1> &             traction,
                                          double &                             fac
                                          ) = 0;


        virtual void buildCouplingMatrices(LINALG::Matrix<3,1> &                                normal,
                                           const double                                         fac,
                                           LINALG::Matrix<nen_,1> &                             funct,
                                           LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
                                           ) = 0;

        virtual void buildFinalCouplingMatrices(LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  invK_ss,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &   rhs
                                                ) = 0;

        virtual void buildCouplingMatricesNitsche(   Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
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
                                                     LINALG::Matrix<nsd_,1> &      velint,          // bg u^n
                                                     LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                                  ) = 0;

        virtual void get_vel_WeakDBC (LINALG::Matrix<nsd_,1> & ivelint) = 0;

      };



      template<DRT::Element::DiscretizationType distype>
      class EmbCoupling
      {
      public:


        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;


        static Teuchos::RCP<EmbCoupling<distype> > TwoSidedImpl(DRT::Element *             emb_ele,
                                                                Epetra_SerialDenseMatrix & C_uiu,
                                                                Epetra_SerialDenseMatrix & C_uui,
                                                                Epetra_SerialDenseMatrix & rhC_ui,
                                                                Epetra_SerialDenseMatrix & C_uiui,
                                                                Epetra_SerialDenseMatrix & emb_xyze
                                                                );


        virtual ~EmbCoupling() {}

        virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside ) = 0;


        virtual void emb_vel(const DRT::Discretization &  cutdis,
                           const std::string              state,
                           const vector<int>&             lm
                           ) = 0;


        virtual void addembdisp(const DRT::Discretization &  cutdis,
                               const std::string             state,
                               const vector<int>&            lm,
                               Epetra_SerialDenseMatrix  &   side_xyze
                               ) = 0;

        virtual void element_length( double & hk_emb ) = 0;


        virtual void buildCouplingMatricesNitscheTwoSided(   Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
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
                                                             LINALG::Matrix<nsd_,1> &      velint,          // bg u^n
                                                             LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                                             ) = 0;

        virtual void get_vel_WeakDBC (LINALG::Matrix<nsd_,1> & emb_velint) = 0;

      };

    } // XFLUID
  } //ELEMENTS
} //DRT




#endif

