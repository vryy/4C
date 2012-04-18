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

        //! Mixed/Stress/Hybrid (MSH) coupling
        static Teuchos::RCP<SideInterface<distype> > Impl(
            DRT::Element *              side,            ///< side element
            Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
            Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
            Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
            Epetra_SerialDenseMatrix &  Gsui,            ///< interface sigma-u_interface coupling
            Epetra_SerialDenseMatrix &  Guis,            ///< interface u_interface-sigma coupling
            Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
        );

        //! Nitsche (NIT) coupling between background element and side element
        static Teuchos::RCP<SideInterface<distype> > Impl(
            DRT::Element *              side,            ///< side element
            Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
            Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
            Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
            Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
            Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
        );

        //! simple side coupling together with embedded Nitsche coupling
        static Teuchos::RCP<SideInterface<distype> > Impl(
            DRT::Element *              side,            ///< side element
            Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
        );

        virtual ~SideInterface(){}


        //! add side's interface displacements and set current side node coordiantes
        virtual void addeidisp(
            const DRT::Discretization &  cutdis,       ///< cut discretization
            const std::string            state,        ///< state
            const vector<int>&           lm,           ///< local map
            Epetra_SerialDenseMatrix  &  side_xyze     ///< side node coordinates
            ) = 0;

        //! extract/set side's interface velocity
        virtual void eivel(
            const DRT::Discretization &  cutdis,       ///< cut discretization
            const std::string            state,        ///< state
            const vector<int>&           lm            ///< local map
            ) = 0;

        //! evaluate shape function, derivatives, normal and transformation w.r.t side element at gaussian point
        virtual void Evaluate(
            const LINALG::Matrix<nsd_-1,1>  & eta,     ///< local coordinates w.r.t side element
            LINALG::Matrix<nsd_,1>          & x,       ///< global coordinates of gaussian point
            LINALG::Matrix<nsd_,1>          & normal,  ///< normal vector
            double                          & drs      ///< transformation factor
        ) = 0;

        //! set prescribed WDBC at Gaussian point
        virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & ivelint ) = 0;

        //! compute interface force for side nodes
        virtual void InterfaceForce(
            Epetra_SerialDenseVector &   iforce,     ///< interface force vector
            LINALG::Matrix<nsd_,1> &     traction,   ///< traction vector at gaussian point
            const double &               fac         ///< integration factor
        ) = 0;

        //! project gaussian point from linearized interfac in normal direction onto corresponding side
        virtual void ProjectOnSide(
            LINALG::Matrix<3,1> & x_gp_lin,            ///< global coordinates of gaussian point w.r.t linearized interface
            LINALG::Matrix<3,1> & x_side,              ///< projected gaussian point on side
            LINALG::Matrix<2,1> & xi_side              ///< local coordinates of projected gaussian point w.r.t side
        ) = 0;

        //! build coupling matrices for Mixed/Stress/Hybrid (MSH) method
        virtual void MSH_buildCouplingMatrices(
            LINALG::Matrix<nsd_,1> &                             normal,     ///< normal vector
            const double                                         fac,        ///< integration factor
            LINALG::Matrix<nen_,1> &                             funct,      ///< shape function
            LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs         ///< rhs block matrix
        ) = 0;

        //! build final coupling matrices for Mixed/Stress/Hybrid (MSH) method
        virtual void MSH_buildFinalCouplingMatrices(
            LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  BinvK_ss,     ///< block inverse K_sigma_sigma matrix
            LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,         ///< block K_iK matrix
            LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,         ///< block K_su matrix
            LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &   rhs           ///< block rhs vector
        ) = 0;

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
        ) = 0;

      };



      template<DRT::Element::DiscretizationType distype>
      class EmbCoupling
      {
      public:


        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

        //!  Nitsche (NIT) coupling between background element and embedded element
        static Teuchos::RCP<EmbCoupling<distype> > TwoSidedImpl(
            DRT::Element *              emb_ele,         ///< side element
            Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
            Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
            Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
            Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
            Epetra_SerialDenseMatrix &  emb_xyze         ///< global node coordinates
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


        virtual void NIT2_buildCouplingMatrices(   Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
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

