/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of XFluid element interface coupling

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_XFEM_HPP
#define FOUR_C_FLUID_ELE_CALC_XFEM_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_fluid_ele_calc_xfem_coupling.hpp"
#include "4C_utils_singleton_owner.hpp"
#include "4C_xfem_condition_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace XFEM
{
  class CouplingBase;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterXFEM;

    /// Fluid element interface coupling implementation with XFEM
    /*!
      This internal class keeps all the working arrays needed to
      calculate the interface stabilization for XFEM with fluid elements.
      The method element_xfem_interface_hybrid_lm() provides a clean and fast element implementation
      for Mixed/Stress/Hybrid interface coupling, using either Cauchy stress-based or viscous
      stress-based Lagrange multipliers.
      The method element_xfem_interface_nit() provides a clean and fast element implementation
      for interface coupling using Nitsche's method (only the cut element sided mortaring
      for xfluid and xfluidfluid applications).
      The method ElementXfemInterfaceNIT2() provides a clean and fast element implementation
      for interface coupling using Nitsche's method also for non-"cut element sided mortaring"
      (can be used for embedded mortaring in case of Xfluidfluid).

      <h3>Purpose</h3>

      The fluid element will allocate exactly one object of this class for all
      fluid elements with the same number of nodes in the mesh. This
      allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>Usability</h3>

      The calculations are done by the EvaluateXfemInterface...() methods.

      \author schott
      \date 04/12
     */
    template <CORE::FE::CellType distype>
    class FluidEleCalcXFEM : public FluidEleCalc<distype>
    {
      /// private Constructor since we are a Singleton.
      FluidEleCalcXFEM();

      /// private copy Constructor since we are a Singleton.
      FluidEleCalcXFEM(FluidEleCalcXFEM const& copy);

      /// private assignment operator since we are a Singleton.
      FluidEleCalcXFEM& operator=(FluidEleCalcXFEM const& copy);

      typedef FluidEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;
      using my::numdofpernode_;

     private:
      /// pointer to the cast object, fluid parameter list for XFEM
      DRT::ELEMENTS::FluidEleParameterXFEM* fldparaxfem_;


      /// number of stress-dof
      static constexpr int numstressdof_ = CORE::FE::DisTypeToNumDeriv2<distype>::numderiv2;

      //! @name useful constants for DOF-index numbering
      //@{
      static constexpr unsigned Velx = 0;
      static constexpr unsigned Vely = 1;
      static constexpr unsigned Velz = 2;
      static constexpr unsigned Pres = 3;

      static constexpr unsigned Sigmaxx = 0;
      static constexpr unsigned Sigmaxy = 1;
      static constexpr unsigned Sigmaxz = 2;
      static constexpr unsigned Sigmayx = 1;
      static constexpr unsigned Sigmayy = 3;
      static constexpr unsigned Sigmayz = 4;
      static constexpr unsigned Sigmazx = 2;
      static constexpr unsigned Sigmazy = 4;
      static constexpr unsigned Sigmazz = 5;
      //@}

      /// get stress dof-index
      unsigned stress_index(unsigned xi, unsigned xj)
      {
        return (xi * xj > 0) ? xi + xj + 1 : xi + xj;
      }

     public:
      /// Singleton access method
      static FluidEleCalcXFEM<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// evaluate the XFEM cut element
      int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<CORE::MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const std::vector<CORE::FE::GaussIntegration>& intpoints,
          const CORE::GEO::CUT::plain_volumecell_set& cells, bool offdiag = false) override;

      /// evaluate the shape functions in the XFEM
      int integrate_shape_function_xfem(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          const std::vector<CORE::FE::GaussIntegration>& intpoints,
          const CORE::GEO::CUT::plain_volumecell_set& cells) override;

      /// error computation
      int compute_error(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<CORE::MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& ele_dom_norms) override;

      int compute_error(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<CORE::MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& ele_dom_norms,
          const CORE::FE::GaussIntegration& intpoints) override;

      int compute_error_interface(DRT::ELEMENTS::Fluid* ele,         ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<CORE::MAT::Material>& mat,                    ///< material
          CORE::LINALG::SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::FE::GaussIntegration>>&
              bintpoints,                                     ///< boundary integration points
          const CORE::GEO::CUT::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                      ///< parameter list
          ) override;

      /// add terms from mixed/hybrid Lagrange multiplier coupling approach to element matrix and
      /// rhs
      void element_xfem_interface_hybrid_lm(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<CORE::FE::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Teuchos::ParameterList& params,          ///< parameter list
          Teuchos::RCP<CORE::MAT::Material>& mat,  ///< material
          CORE::LINALG::SerialDenseMatrix&
              elemat1_epetra,  ///< local system matrix of intersected element
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,                      ///< local element vector of intersected element
          CORE::LINALG::SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const CORE::GEO::CUT::plain_volumecell_set& vcSet  ///< set of plain volume cells
          ) override;

      /// add Nitsche (NIT) interface condition to element matrix and rhs
      void element_xfem_interface_nit(DRT::ELEMENTS::Fluid* ele,     ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,                     ///< parameter list
          Teuchos::RCP<CORE::MAT::Material>& mat_master,      ///< material for the coupled side
          Teuchos::RCP<CORE::MAT::Material>& mat_slave,       ///< material for the coupled side
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,    ///< element matrix
          CORE::LINALG::SerialDenseVector& elevec1_epetra,    ///< element vector
          const CORE::GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          CORE::LINALG::SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
          ) override;

      ///
      void calculate_continuity_xfem(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                              ///< discretization
          const std::vector<int>& lm,                            ///< local map
          CORE::LINALG::SerialDenseVector& elevec1_epetra,       ///< element vector
          const CORE::FE::GaussIntegration& intpoints            ///< integration points
          ) override;

      ///
      void calculate_continuity_xfem(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                              ///< discretization
          const std::vector<int>& lm,                            ///< local map
          CORE::LINALG::SerialDenseVector& elevec1_epetra        ///< element vector
          ) override;

     private:
      //! evaluate analytical reference solution
      void analytical_reference(const int calcerr,   ///< which reference solution
          const int calcerrfunctno,                  ///< error function number
          CORE::LINALG::Matrix<nsd_, 1>& u,          ///< exact jump vector (coupled)
          CORE::LINALG::Matrix<nsd_, nsd_>& grad_u,  ///< exact velocity gradient
          double& p,                                 ///< exact pressure
          CORE::LINALG::Matrix<nsd_, 1>& xyzint,     ///< xyz position of gaussian point
          const double& t,                           ///< time
          Teuchos::RCP<CORE::MAT::Material> mat = Teuchos::null);

      //! get the interface jump vectors for velocity and traction at the Gaussian point
      void get_interface_jump_vectors(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          Teuchos::RCP<XFEM::CouplingBase> coupling,  ///< coupling object
          CORE::LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface jump vector for velocity
          CORE::LINALG::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface jump vector for traction
          CORE::LINALG::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
          CORE::LINALG::Matrix<nsd_, nsd_>&
              LB_proj_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
          const CORE::LINALG::Matrix<nsd_, 1>& x,       ///< global coordinates of Gaussian point
          const CORE::LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>>
              si,                                       ///< side implementation for cutter element
          CORE::LINALG::Matrix<3, 1>& rst,              ///< local coordinates of GP for bg element
          double& kappa_m,                              ///< fluid sided weighting
          double& visc_m,                               ///< fluid sided weighting
          double& visc_s,                               ///< slave sided dynamic viscosity
          CORE::LINALG::Matrix<3, 1>& rst_slave,        ///< local coord of gp in slave element
          std::vector<double>& eledisp,                 ///< slave element displacement vector
          CORE::Elements::Element* coupl_ele = nullptr  ///< slave coupling element
      );

      //! get the interface jump vectors for velocity and traction at the Gaussian point for
      //! previous time step
      void get_interface_jump_vectors_old_state(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          Teuchos::RCP<XFEM::CouplingBase> coupling,  ///< coupling object
          CORE::LINALG::Matrix<nsd_, 1>&
              ivelintn_jump,  ///< prescribed interface jump vector for velocity
          CORE::LINALG::Matrix<nsd_, 1>&
              itractionn_jump,  ///< prescribed interface jump vector for traction
          const CORE::LINALG::Matrix<nsd_, 1>& x,       ///< global coordinates of Gaussian point
          const CORE::LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>>
              si,                          ///< side implementation for cutter element
          const double& presn_m,           ///< coupling master pressure
          CORE::LINALG::Matrix<3, 1>& rst  ///< local coordinates of GP for bg element
      );

      //! build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
      void nit_build_patch_cuiui(
          CORE::LINALG::SerialDenseMatrix&
              Cuiui,  ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              Cuiui_coupling  ///< Cuiui matrices for all cutting sides
      );

      //! compute viscous part of Nitsche's penalty term scaling for Nitsche's method
      double nit_compute_visc_penalty_stabfac(
          const CORE::FE::CellType ele_distype,  ///< the discretization type of the element w.r.t
                                                 ///< which the stabilization factor is computed
          const double inv_hk,                   ///< the inverse characteristic element length
          const double kappa1,                   ///< Weight parameter (parameter +/master side)
          const double kappa2                    ///< Weight parameter (parameter -/slave  side)
      );

      //! get the constant which satisfies the trace inequality depending on the spatial dimension
      //! and polynomial order of the element
      double nit_get_trace_estimate_constant(const CORE::FE::CellType ele_distype);

      //! prepare coupling matrices, that include contributions from convective stabilization
      void hybrid_lm_create_special_contribution_matrices(
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          std::set<int>& begids,  ///< ids of intersecting boundary elements
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              side_coupling_extra  ///< contributions to coupling matrices from convective
                                   ///< stabilizations
      );

      //! evaluate Neumann boundary condition
      void evaluate_neumann(const double& timefacfac,    ///< theta*dt
          const CORE::LINALG::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
          const CORE::LINALG::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface traction, jump height for coupled problems
          CORE::LINALG::SerialDenseMatrix::Base& elevec1_epetra  ///< element vector
      );

      //! build traction vector w.r.t fluid domain
      void build_traction_vector(CORE::LINALG::Matrix<nsd_, 1>& traction,  ///< traction vector
          double& press,                         ///< pressure at gaussian point
          CORE::LINALG::Matrix<nsd_, 1>& normal  ///< normal vector
      );

      //! assemble side's interface force
      void assemble_interface_force(
          Teuchos::RCP<Epetra_Vector> iforcecol,   ///< interface force column vector
          DRT::Discretization& cutdis,             ///< cut discretization
          std::vector<int>& lm,                    ///< local dof map
          CORE::LINALG::SerialDenseVector& iforce  ///< interface force vector
      );

      //! evaluate shape function and derivative at point with local coordinates rst
      void eval_func_and_deriv(CORE::LINALG::Matrix<3, 1>& rst);

      //! build matrices from volume-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM-coupling \author kruse \date 06/14
      void hybrid_lm_build_vol_based(const std::vector<CORE::FE::GaussIntegration>& intpoints,
          const CORE::GEO::CUT::plain_volumecell_set& cells,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          const CORE::LINALG::Matrix<nen_, 1>& epreaf,     ///< element pressure
          CORE::LINALG::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          CORE::LINALG::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, nsd_, nsd_>&
              K_uu,  ///< K_uu matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, nsd_, 1>&
              rhs_uu,                    ///< rhs_u(u) vector
          const bool is_MHVS = false,    ///< viscous (true) or Cauchy (false) stress-based LM
          const double mhvs_param = 1.0  ///< stabilizing parameter for viscous stress-based LM
      );

      //! evaluate matrices from volume-based terms for Cauchy stress-based mixed/hybrid LM coupling
      //! at current Gauss-point
      void mhcs_evaluate_vol_based(
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          const CORE::LINALG::Matrix<nen_, 1>& epreaf,     ///< element pressure
          CORE::LINALG::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          CORE::LINALG::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s  ///< rhs_s vector
      );

      //! evaluate matrices from volume-based terms for viscous stress-based mixed/hybrid LM
      //! coupling at current Gauss-point \author kruse \date 06/14
      void mhvs_evaluate_vol_based(
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          CORE::LINALG::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          CORE::LINALG::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, nsd_, nsd_>&
              K_uu,  ///< K_uu matrix
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, nsd_, 1>&
              rhs_uu,               ///< rhs_uu vector
          const double& mhvs_param  ///< stabilizing parameter
      );

      //! evaluate matrices from surface-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM coupling at current Gauss-point \author kruse \date 06/14
      void hybrid_lm_evaluate_surf_based(
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>>& si,
          const CORE::LINALG::Matrix<nen_, nen_>& bK_ss,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, numstressdof_, 1>& rhs_s,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, nsd_, nsd_>& K_uu,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, nsd_, 1>& rhs_uu,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, nsd_, 1>& G_up,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, nen_>, 1, nsd_>& G_pu,
          CORE::LINALG::BlockMatrix<CORE::LINALG::Matrix<nen_, 1>, nsd_, 1>& rhs_up,
          CORE::LINALG::Matrix<nen_, 1>& rhs_pu, const CORE::LINALG::Matrix<nsd_, 1>& normal,
          const double& timesurffac, const CORE::LINALG::Matrix<nsd_, 1>& ivelint_jump,
          const CORE::LINALG::Matrix<nsd_, 1>& itraction_jump, const bool eval_side_coupling,
          const bool is_MHVS);

      //! Initiates dummy variables and calls FLuidEleCalc get_material_params routine.
      void get_material_parameters_volume_cell(Teuchos::RCP<const CORE::MAT::Material> material,
          double& densaf,  // done
          double& viscaf,  // done
          double& gamma    // done
      );

      // Density on each side of interface
      double densaf_master_;  ///< density at master side
      double densaf_slave_;   ///< density at slave side

      // Viscosity on each side of interface
      double viscaf_master_;  ///< viscosity at master side
      double viscaf_slave_;   ///< viscosity at slave side

      // Surface tension on each side (If it is not the same, error should be thrown.)
      double gamma_m_;  ///< surface tension coefficient at master side
      double gamma_s_;  ///< surface tension coefficient at master side

      CORE::LINALG::Matrix<nsd_, nen_>
          evelaf_;  ///< element velocity at time t^n+af, implemented also in base class, not
                    ///< accessable via extract_values_from_global_vector
      CORE::LINALG::Matrix<nen_, 1>
          epreaf_;  ///< element pressure at time t^n+af, implemented also in base class, not
                    ///< accessable via extract_values_from_global_vector

      CORE::LINALG::Matrix<nsd_, nen_>
          eveln_;  ///< element velocity at time t^n, implemented also in base class, not accessable
                   ///< via extract_values_from_global_vector
      CORE::LINALG::Matrix<nen_, 1>
          epren_;  ///< element velocity at time t^n, implemented also in base
                   ///< class, not accessable via extract_values_from_global_vector

      CORE::LINALG::Matrix<nsd_, 1> ivelint_jump_;        ///< interface velocity jump at t^n+1
      CORE::LINALG::Matrix<nsd_, 1> itraction_jump_;      ///< interface traction jump at t^n+1
      CORE::LINALG::Matrix<nsd_, nsd_> proj_tangential_;  ///< tangential projection matrix at t^n+1
      CORE::LINALG::Matrix<nsd_, nsd_>
          lb_proj_matrix_;  ///< interface matrix jump (for Laplace-Beltrami) at t^n+1

      std::vector<CORE::LINALG::SerialDenseMatrix>
          solid_stress_;  ///< hold information about solid stress ([0]...traction,
                          ///< [1]...dtraction_dv, [2-4]...d2traction_dv2)

      CORE::LINALG::Matrix<nsd_, 1> ivelintn_jump_;    ///< interface velocity jump at t^n
      CORE::LINALG::Matrix<nsd_, 1> itractionn_jump_;  ///< interface traction jump at t^n

      CORE::LINALG::Matrix<nsd_, 1>
          velint_s_;  ///< velocity of slave side at time t^n+1 at integration point
      CORE::LINALG::Matrix<nsd_, 1>
          velintn_s_;  ///< velocity of slave side at time t^n at integration point

      CORE::LINALG::Matrix<nsd_, 1>
          rst_;  ///< local coordinates of Gaussian point w.r.t background element

      CORE::LINALG::Matrix<nsd_, 1> normal_;  ///< normal vector
      CORE::LINALG::Matrix<nsd_, 1> x_side_;  ///< gauss-point coordinates
      CORE::LINALG::Matrix<nsd_, 1>
          x_gp_lin_;  ///< gauss-point in xyz-system on linearized interface
    };


  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
