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

#include "../drt_inpar/inpar_xfem.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {

    /// Factory class that provides the creation of the Side implementation class (SideImpl) for XFEM interface coupling
    /*!
      This class is a factory class that creates templated Side implementation class (SideImpl).
      Dependent on the template pairs between cut background element and cutting interface Side the SideImpl
      is created, which provides the implementation of coupling matrices for different coupling strategies,
      like Nitsche's method of Mixed/Stress/Hybrid method.

      <h3>Usability</h3>

      The calculations are done by the EvaluateXfemInterface...() methods that call this SideInterface factory class

      \author schott
      \date 04/12
     */
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
            const vector<int>&           lm            ///< local map
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

        virtual void Stab_InflowCoercivity(
            Epetra_SerialDenseMatrix &          C_uu_,            ///< standard bg-bg-matrix
            Epetra_SerialDenseVector &          rhs_Cu_,          ///< standard bg-rhs
            const bool &                        coupling,         ///< assemble coupling terms (yes/no)
            const bool &                        bg_mortaring,     ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
            const LINALG::Matrix<nsd_,1> &      normal,           ///< normal vector
            const double &                      timefacfac,       ///< theta*dt
            const double &                      visceff_1,        ///< viscosity in background fluid
            const double &                      visceff_2,        ///< viscosity in embedded fluid
            const double &                      kappa1,           ///< mortaring weighting
            const double &                      kappa2,           ///< mortaring weighting
            const double &                      stabfac_conv,     ///< Nitsche convective non-dimensionless stabilization factor
            const LINALG::Matrix<nen_,1> &      funct_,            ///< bg shape functions
            const LINALG::Matrix<nsd_,nen_> &   derxy_,           ///< bg deriv
            const LINALG::Matrix<nsd_,nsd_> &   vderxy_,          ///< bg deriv^n
            const LINALG::Matrix<nsd_,1> &      velint,           ///< bg u^n
            const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,///< Dirichlet velocity vector or prescribed jump vector
            INPAR::XFEM::ConvStabScaling        conv_stab_scaling,///< Inflow term strategies
            string                              coupl_method      ///< coupling method (NIT or MSH)
          ) = 0;

        //! build coupling matrices and assemble terms for Nitsche's (NIT) method
        virtual void NIT_buildCouplingMatrices(
            Epetra_SerialDenseMatrix &          C_uu_,            ///< standard bg-bg-matrix
            Epetra_SerialDenseVector &          rhs_Cu_,          ///< standard bg-rhs
            const bool &                        coupling,         ///< assemble coupling terms (yes/no)
            const bool &                        bg_mortaring,     ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
            const LINALG::Matrix<nsd_,1> &      normal,           ///< normal vector
            const double &                      timefacfac,       ///< theta*dt
            const double &                      visceff_1,        ///< viscosity in background fluid
            const double &                      visceff_2,        ///< viscosity in embedded fluid
            const double &                      kappa1,           ///< mortaring weighting
            const double &                      kappa2,           ///< mortaring weighting
            const double &                      stabfac,          ///< Nitsche non-dimensionless stabilization factor
            const double &                      stabfac_conv,     ///< Nitsche convective non-dimensionless stabilization factor
            const LINALG::Matrix<nen_,1> &      funct,            ///< bg shape functions
            const LINALG::Matrix<nsd_,nen_> &   derxy_,           ///< bg deriv
            const LINALG::Matrix<nsd_,nsd_> &   vderxy_,          ///< bg deriv^n
            const double &                      press,            ///< bg p^n
            const LINALG::Matrix<nsd_,1> &      velint,           ///< bg u^n
            const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,///< Dirichlet velocity vector or prescribed jump vector
            INPAR::XFEM::ConvStabScaling        conv_stab_scaling
          ) = 0;

      };


      /// Factory class that provides the creation of the embedded element implementation class (EmbCoupling) for XFEM interface coupling
      /*!
        This class is a factory class that creates templated Element implementation class (EmbImpl).
        Dependent on the template pairs between cut background element and cutting element the EmbImpl
        is created, which provides the implementation of coupling matrices for different coupling strategies,
        like Nitsche's method.

        <h3>Usability</h3>

        The calculations are done by the EvaluateXfemInterface...() methods that call this EmbCoupling factory class

        \author schott
        \date 04/12
       */
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

        //! add embedded element displacements and set current element node coordinates
        virtual void addembdisp(
            const DRT::Discretization &  embdis,       ///< embedded discretization
            const std::string            state,        ///< state
            const vector<int>&           lm            ///< local map
        ) = 0;

        //! extract/set embedded element velocity
        virtual void emb_vel(
            const DRT::Discretization &  embdis,       ///< embedded discretization
            const std::string            state,        ///< state
            const vector<int>&           lm            ///< local map
        ) = 0;

        //! compute embedded element's element length
        virtual void element_length( double & hk_emb ) = 0;

        //!  evaluate shape function, derivatives and transformation w.r.t embedded element at gaussian point
        virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside ) = 0;

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
            double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
            double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
            bool &                        velgrad_interface_stab,  /// penalty term for velocity gradients at the interface
            double &                      velgrad_interface_fac,  /// velgrad_interface_stab stabilization factor
            LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
            LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
            LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
            double &                      press,          // bg p^n
            LINALG::Matrix<nsd_,1> &      velint,          // bg u^n
            LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP, // Dirichlet velocity vector or prescribed jump vector
            INPAR::XFEM::ConvStabScaling  conv_stab_scaling // Inflow term strategies
          ) = 0;


      };

    } // XFLUID
  } //ELEMENTS
} //DRT




#endif

