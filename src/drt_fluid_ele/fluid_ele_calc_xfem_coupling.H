/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class for providing an implementation for coupling with
       Mixed/Stress/Hybrid methods or Nitsche's method to enforce
       interface conditions in the XFEM weakly

\level 2

\maintainer  Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_COUPLING_H
#define FLUID_ELE_CALC_XFEM_COUPLING_H

/*
 * When we enforce coupling conditions or dirichlet conditions on embedded interfaces in the XFEM
 * weakly, we obtain contributions from coupling terms, that have to appear in the element matrix &
 * rhs of the XFEM-element (further denoted as 'master') and the element we want to couple with
 * (named 'slave', can be an XFEM-element itself!). We try to keep this interface as abstract as
 * possible, such that we can use the same routines for xfluid WDBC-problems, xfluid-fluid coupling,
 * XFSI and two-phase flow (and future applications). The classes in here provide an interface to
 * evaluate the slave element and handle the coupling terms, especially all those, that are added to
 * the slave element system. There are two coupling methods available, which are inherently similar:
 * Nitsche's method and a mixed/hybrid stress-based Lagrange multiplier ("HybridLM") approach.
 * kruse, 08/14
 */

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../linalg/linalg_fixedsizeblockmatrix.H"

#include "../drt_inpar/inpar_xfem.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterXFEM;

    namespace XFLUID
    {
      /**
       * Interface class for slave element (whole 3D element or just a 2D boundary element).
       * Handles evaluation of slave element's shape functions, provides access to
       * velocities, velocity gradients and pressure at the coupling interface, and more.
       * This class includes all functionalities, that are related to the element to be coupled,
       * not to the coupling method itself.
       * It can be used as a factory for pure coupling slave element representations.
       * These are useful, if you want to access the coupled element
       * without assembling any coupling terms (e.g. in two-sided xfluid-fluid coupling
       * for gauss-point projection or error calculation in case of a given analytical solution).
       */
      template <DRT::Element::DiscretizationType distype>
      class SlaveElementInterface
      {
       public:
        /// number of spatial dimensions
        static const unsigned nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

        //! set names of displacement and velocity states (differ dependent on the slave element)
        static void DefineStateNames(
            DRT::Element::DiscretizationType slave_distype,  ///< coupling slave discretization type
            std::string& disp_statename,  ///< name of displacement state at current step
            std::string& vel_statename,   ///< name of velocity state at current step
            std::string& veln_statename   ///< name of velocity state at previous step
        );

        //! create a pure coupling slave element representation
        static Teuchos::RCP<SlaveElementInterface<distype>> CreateSlaveElementRepresentation(
            DRT::Element* slave_ele,  ///< coupling slave element
            Epetra_SerialDenseMatrix&
                slave_xyz  ///< global node coordinates of coupling slave element
        );

        //! dtor
        virtual ~SlaveElementInterface(){};

        //! add slave element's displacements and set current element's nodal coordinates
        virtual void AddSlaveEleDisp(
            const DRT::Discretization& slavedis,  ///< coupling slave discretization
            const std::vector<int>& lm,           ///< local map
            std::vector<double>& mymatrix         ///< slave element displacement vector
        )
        {
          dserror("There is no concrete slave element available.");
        };

        virtual void AddSlaveEleDisp(
            const DRT::Discretization& slavedis,  ///< coupling slave discretization
            const std::vector<int>& lm            ///< local map
        )
        {
          dserror("There is no concrete slave element available.");
        };

        //! @name accessors
        //@{

        //! set slave element's nodal velocities
        virtual void SetSlaveState(
            const DRT::Discretization& slavedis,  ///< embedded discretization
            const std::vector<int>& lm            ///< local map
        )
        {
          dserror("There is no concrete slave element available.");
        };

        //! set slave element's nodal velocities
        virtual void SetSlaveStaten(
            const DRT::Discretization& slavedis,  ///< embedded discretization
            const std::vector<int>& lm            ///< local map
        )
        {
          dserror("There is no concrete slave element available.");
        };

        /*!
         * get interface velocity - return 0 velocity (if we don't have a concrete slave
         * element (level-set applications), this is correct - otherwise this method is
         * overloaded!!)
         */
        virtual void GetInterfaceVelnp(
            LINALG::Matrix<nsd_, 1>& ivelint  ///< interface velocity at coupling slave side
            ) const
        {
          ivelint = LINALG::Matrix<nsd_, 1>(true);
        };

        /*!
         * get interface velocity - return 0 velocity (if we don't have a concrete slave
         * element (level-set applications), this is correct - otherwise this method is
         * overloaded!!)
         */
        virtual void GetInterfaceVeln(
            LINALG::Matrix<nsd_, 1>& ivelint  ///< interface velocity at coupling slave side
            ) const
        {
          ivelint = LINALG::Matrix<nsd_, 1>(true);
        };

        //! get interface pressure
        virtual void GetInterfacePresnp(double& ipres  ///< interface pressure
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //! get interface pressure
        virtual void GetInterfacePresn(double& ipres  ///< interface pressure
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //! get interface velocity gradient
        virtual void GetInterfaceVelGradnp(LINALG::Matrix<nsd_, nsd_>&
                velgradint  ///< interface velocity gradients at coupling slave side
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //! get interface velocity gradient
        virtual void GetInterfaceVelGradn(LINALG::Matrix<nsd_, nsd_>&
                velgradint  ///< interface velocity gradients at coupling slave side
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //! set state for interface velocity jump
        virtual void SetInterfaceJumpStatenp(
            const DRT::Discretization& cutterdis,  ///< cutter discretization
            const std::string state,               ///< state
            const std::vector<int>& lm             ///< local map
        )
        {
          dserror("There is no concrete slave element available.");
        };

        //! set state for interface velocity jump for previous time step
        virtual void SetInterfaceJumpStaten(
            const DRT::Discretization& cutterdis,  ///< cutter discretization
            const std::string state,               ///< state
            const std::vector<int>& lm             ///< local map
        )
        {
          dserror("There is no concrete slave element available.");
        };

        //! get interface velocity jump at Gaussian point
        virtual void GetInterfaceJumpVelnp(
            LINALG::Matrix<nsd_, 1>& ivelint_jump  ///< cutter element interface velocity jump or
                                                   ///< prescribed DBC at Gaussian point
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //! get interface velocity jump at Gaussian point for previous time step
        virtual void GetInterfaceJumpVeln(
            LINALG::Matrix<nsd_, 1>& ivelintn_jump  ///< cutter element interface velocity jump or
                                                    ///< prescribed DBC at Gaussian point
            ) const
        {
          dserror("There is no concrete slave element available.");
        };

        //@}

        //! evaluate shape function, derivatives and transformation w.r.t coupling slave element at
        //! gaussian point
        virtual void Evaluate(LINALG::Matrix<nsd_, 1>& xside)
        {
          dserror("There is no concrete slave element available.");
        }

        //! evaluate shape function, derivatives and transformation w.r.t coupling slave element at
        //! gaussian point
        virtual void Evaluate(LINALG::Matrix<nsd_, 1>& xside, LINALG::Matrix<nsd_, 1>& rst_slave)
        {
          dserror("There is no concrete slave element available.");
        }

        //! compute interface force for side nodes
        virtual void ComputeInterfaceForce(
            Epetra_SerialDenseVector& iforce,   ///< interface force vector
            LINALG::Matrix<nsd_, 1>& traction,  ///< traction vector at gaussian point
            const double& fac                   ///< integration factor
        )
        {
          dserror("There is no concrete slave element available.");
        }

        //! project gaussian point from linearized interface in normal direction onto corresponding
        //! side (2D slave element)
        virtual void ProjectOnSide(
            LINALG::Matrix<nsd_, 1>&
                x_gp_lin,  ///< global coordinates of gaussian point w.r.t linearized interface
            LINALG::Matrix<nsd_, 1>& x_side,  ///< projected gaussian point on side
            LINALG::Matrix<nsd_, 1>&
                xi_side  ///< local coordinates of projected gaussian point w.r.t side
        )
        {
          dserror("There is no concrete slave element available.");
        }

        //! evaluate element volume
        virtual double EvalElementVolume()
        {
          dserror("There is no concrete slave element available.");
          return 0.0;
        }
      };

      //! factory class for coupling with Nitsche's method
      template <DRT::Element::DiscretizationType distype>
      class NitscheInterface : virtual public SlaveElementInterface<distype>
      {
       public:
        /// number of nodes per master element
        static const unsigned nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        /// number of spatial dimensions
        static const unsigned nsd_ = SlaveElementInterface<distype>::nsd_;
        /// number of nodal DOF for master element (always a xfem-fluid element)
        static const unsigned master_numdof_ = nsd_ + 1;

        static Teuchos::RCP<NitscheInterface<distype>> CreateNitscheCoupling_XFluidWDBC(
            Epetra_SerialDenseMatrix& C_umum,  ///< C_umum coupling matrix
            Epetra_SerialDenseMatrix& rhC_um,  ///< C_um coupling rhs
            const DRT::ELEMENTS::FluidEleParameterXFEM&
                fldparaxfem  ///< specific XFEM based fluid parameters
        );

        //! create a coupling interface for Nitsche's method for xfluid weak dirichlet problems
        static Teuchos::RCP<NitscheInterface<distype>> CreateNitscheCoupling_XFluidWDBC(
            DRT::Element* bele,                  ///< boundary element
            Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
            Epetra_SerialDenseMatrix& C_umum,    ///< C_umum coupling matrix
            Epetra_SerialDenseMatrix& rhC_um,    ///< C_um coupling rhs
            const DRT::ELEMENTS::FluidEleParameterXFEM&
                fldparaxfem  ///< specific XFEM based fluid parameters
        );

        //! create a coupling interface for Nitsche's method for xfluid-sided coupling strategy
        static Teuchos::RCP<NitscheInterface<distype>> CreateNitscheCoupling_XFluidSided(
            DRT::Element* bele,                  ///< boundary element
            Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
            Epetra_SerialDenseMatrix& C_umum,    ///< C_umum coupling matrix
            Epetra_SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
            Epetra_SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
            Epetra_SerialDenseMatrix& C_usus,    ///< C_usus coupling matrix
            Epetra_SerialDenseMatrix& rhC_um,    ///< C_um coupling rhs
            Epetra_SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
            const DRT::ELEMENTS::FluidEleParameterXFEM&
                fldparaxfem  ///< specific XFEM based fluid parameters
        );

        //! create a coupling interface for Nitsche's method for two-sided coupling strategy
        //! (weighted or fully embedded-sided)
        static Teuchos::RCP<NitscheInterface<distype>> CreateNitscheCoupling_TwoSided(
            DRT::Element* vele,                  ///< volumetric element to couple with
            Epetra_SerialDenseMatrix& vele_xyz,  ///< global node coordinates of volumetric element
            Epetra_SerialDenseMatrix& C_umum,    ///< C_umum coupling matrix
            Epetra_SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
            Epetra_SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
            Epetra_SerialDenseMatrix& C_usus,    ///< C_usus coupling matrix
            Epetra_SerialDenseMatrix& rhC_um,    ///< C_um coupling rhs
            Epetra_SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
            const DRT::ELEMENTS::FluidEleParameterXFEM&
                fldparaxfem  ///< specific XFEM based fluid parameters
        );

        //! dtor
        virtual ~NitscheInterface(){};

        //! add contributions from convective stabilization
        //! this method can also be applied in a non-Nitsche context (e.g. MHVS) by
        //! employing shape functions and velocities from another coupling object (HybridLMCoupling)
        //! slave_ele
        virtual void ApplyConvStabTerms(
            const Teuchos::RCP<SlaveElementInterface<distype>>&
                slave_ele,                            ///< associated slave element coupling object
            const LINALG::Matrix<nen_, 1>& funct_m,   ///< master shape functions
            const LINALG::Matrix<nsd_, 1>& velint_m,  ///< vector of slave shape functions
            const LINALG::Matrix<nsd_, 1>& normal,    ///< normal vector n^b
            const double& density_m,                  ///< fluid density (master)
            const double& NIT_stab_fac_conv,          ///< full Nitsche's penalty term scaling
                                                      ///< (viscous+convective part)
            const double& timefacfac,                 ///< theta*dt
            const LINALG::Matrix<nsd_, 1>&
                ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height
                               ///< for coupled problems
            const INPAR::XFEM::EleCouplingCondType& cond_type  ///< condition type
            ) = 0;

        //! build coupling matrices for Nitsche's method (NIT)
        virtual void NIT_evaluateCoupling(
            const LINALG::Matrix<nsd_, 1>&
                normal,  ///< outward pointing normal (defined by the coupling partner, that
                         ///< determines the interface traction)
            const double& timefacfac,                ///< theta*dt*fac
            const double& pres_timefacfac,           ///< scaling for pressure part
            const double& visceff_m,                 ///< viscosity in coupling master fluid
            const double& visceff_s,                 ///< viscosity in coupling slave fluid
            const double& density_m,                 ///< fluid density (master) USED IN XFF
            const LINALG::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
            const LINALG::Matrix<nsd_, nen_>&
                derxy_m,  ///< spatial derivatives of coupling master shape functions
            const LINALG::Matrix<nsd_, nsd_>&
                vderxy_m,          ///< coupling master spatial velocity derivatives
            const double& pres_m,  ///< coupling master pressure
            const LINALG::Matrix<nsd_, 1>& velint_m,  ///< coupling master interface velocity
            const LINALG::Matrix<nsd_, 1>&
                ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height
                               ///< for coupled problems
            const LINALG::Matrix<nsd_, 1>& itraction_jump,  ///< prescribed interface traction, jump
                                                            ///< height for coupled problems
            const LINALG::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
            const LINALG::Matrix<nsd_, nsd_>&
                LB_proj_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
            const std::vector<Epetra_SerialDenseMatrix>&
                solid_stress,  ///< structural cauchy stress and linearization
            std::map<INPAR::XFEM::CoupTerm, std::pair<bool, double>>&
                configmap  ///< Interface Terms configuration map
            ) = 0;

        //! add rhs contributions from old time step in Nitsche's (NIT) method
        virtual void NIT_evaluateCouplingOldState(
            const LINALG::Matrix<nsd_, 1>&
                normal,  ///< outward pointing normal (defined by the coupling partner, that
                         ///< determines the interface traction)
            const double& timefacfac,                ///< dt*(1-theta)*fac
            bool isImplPressure,                     ///< flag for implicit pressure treatment
            const double& visceff_m,                 ///< viscosity in coupling master fluid
            const double& visceff_s,                 ///< viscosity in coupling slave fluid
            const double& density_m,                 ///< fluid density (master) USED IN XFF
            const LINALG::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
            const LINALG::Matrix<nsd_, nen_>&
                derxy_m,  ///< spatial derivatives of coupling master shape functions
            const LINALG::Matrix<nsd_, nsd_>&
                vderxy_m,          ///< coupling master spatial velocity derivatives
            const double& pres_m,  ///< coupling master pressure
            const LINALG::Matrix<nsd_, 1>& velint_m,  ///< coupling master interface velocity
            const LINALG::Matrix<nsd_, 1>&
                ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height
                               ///< for coupled problems
            const LINALG::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
            const LINALG::Matrix<nsd_, 1>& itraction_jump,  ///< prescribed interface traction, jump
                                                            ///< height for coupled problems
            std::map<INPAR::XFEM::CoupTerm, std::pair<bool, double>>&
                configmap  ///< Interface Terms configuration map
            ) = 0;
      };

      //! factory class for coupling with mixed/hybrid stress-based Lagrange multipliers
      template <DRT::Element::DiscretizationType distype>
      class HybridLMInterface : virtual public SlaveElementInterface<distype>
      {
       public:
        /// number of nodes per master (xfluid) element
        static const unsigned nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        /// number of spatial dimensions of the master element (xfem-fluid)
        static const unsigned nsd_ = SlaveElementInterface<distype>::nsd_;
        /// number of nodal dof for master element (coupling master is always a fluid element!)
        static const unsigned master_numdof_ = nsd_ + 1;
        /// number of independent stress-dof
        static const unsigned numstressdof_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

        //! create a coupling interface for mixed/hybrid LM approach for xfluid weak dirichlet
        //! problems
        static Teuchos::RCP<HybridLMInterface<distype>> CreateHybridLMCoupling_XFluidWDBC(
            bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard
                                          ///< & adjoint viscous term
        );

        //! create a coupling interface for mixed/hybrid LM approach for xfluid weak dirichlet
        //! problems
        static Teuchos::RCP<HybridLMInterface<distype>> CreateHybridLMCoupling_XFluidWDBC(
            DRT::Element* bele,                  ///< boundary element
            Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
            bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard
                                          ///< & adjoint viscous term
        );

        //! create a coupling interface for mixed/hybrid LM approach for xfluid-sided coupling
        //! strategy
        static Teuchos::RCP<HybridLMInterface<distype>> CreateHybridLMCoupling_XFluidSided(
            DRT::Element* bele,                  ///< boundary element
            Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
            Epetra_SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
            Epetra_SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
            Epetra_SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
            Epetra_SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
            Epetra_SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
            bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard
                                          ///< & adjoint viscous term
        );

        //! create a coupling interface for mixed/hybrid LM approach for two-sided coupling strategy
        //! (weighted or fully embedded-sided)
        static Teuchos::RCP<HybridLMInterface<distype>> CreateHybridLMCoupling_TwoSided(
            DRT::Element* vele,                  ///< volumetric element to couple with
            Epetra_SerialDenseMatrix& vele_xyz,  ///< global node coordinates of volumetric element
            Epetra_SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
            Epetra_SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
            Epetra_SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
            Epetra_SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
            Epetra_SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
            DRT::Element* bele,                  ///< boundary element
            Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of slave element
            bool is_viscAdjointSymmetric         ///< flag, that indicates equal signs of Nitsche's
                                                 ///< standard & adjoint viscous term
        )
        {
          dserror("Embedded-sided mixed/hybrid stress-based LM is not implemented yet!");
          return Teuchos::null;
        }

        //! dtor
        virtual ~HybridLMInterface(){};

        //! evaluate interface matrices for mixed/hybrid Cauchy stress-based (MHCS) coupling
        virtual void MHCS_buildCouplingMatrices(
            const LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector
            const double& fac,                      ///< integration factor
            const LINALG::Matrix<nen_, 1>& funct,   ///< shape function
            LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
                rhs_s,  ///< block rhs vector \f$ rhs_{\sigma} \f$
            const LINALG::Matrix<nsd_, 1>&
                ivelint_jump,  ///< prescribed interface velocity or interface jump height
            const LINALG::Matrix<nsd_, 1>&
                itraction_jump  ///< prescribed interface traction or interface jump height
            ) = 0;

        //! evaluate interface matrices for mixed/hybrid viscous stress-based (MHVS) coupling
        virtual void MHVS_buildCouplingMatrices(
            const LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector
            const double& fac,                      ///< integration factor
            const LINALG::Matrix<nen_, 1>& funct,   ///< background element shape functions
            LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
                rhs_s,                          ///< block rhs vector \f$ rhs_{\sigma}\f$
            const double& press,                ///< background element pressure
            LINALG::Matrix<nen_, 1>& rhs_pmus,  ///< contribution to block rhs vector \f$rhs_p\f$
                                                ///< (includes interface slave velocity terms)
            const LINALG::Matrix<nsd_, 1>&
                ivelint_jump,  ///< prescribed interface velocity or interface jump height
            const LINALG::Matrix<nsd_, 1>&
                itraction_jump  ///< prescribed interface traction or interface jump height
            ) = 0;

        //! build the final coupling matrices for mixed/hybrid Cauchy or viscous stress-based
        //! coupling (MHCS or MHVS)
        virtual void HybridLM_buildFinalCouplingMatrices(
            LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, numstressdof_, numstressdof_>&
                BinvK_ss,  ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
            LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, master_numdof_, numstressdof_>&
                BKumsInvKss,  ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
            LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, numstressdof_, master_numdof_>&
                BK_sum,  ///< block matrix \f$ K_{\sigma u} \f$
            LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
                rhs_s  ///< block rhs vector \f$ rhs_{\sigma}\f$
            ) = 0;
      };
    }  // end namespace XFLUID
  }    // end namespace ELEMENTS
}  // end namespace DRT

#endif /* FLUID_ELE_CALC_XFEM_COUPLING_H_ */
