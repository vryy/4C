/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of XFluid element interface coupling

\level 2

\maintainer  Christoph Ager

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_H
#define FLUID_ELE_CALC_XFEM_H

#include "../drt_lib/drt_utils.H"

#include "fluid_ele_calc.H"
#include "fluid_ele_calc_xfem_coupling.H"
#include "../drt_xfem/xfem_condition_manager.H"

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
      The method ElementXfemInterfaceHybridLM() provides a clean and fast element implementation
      for Mixed/Stress/Hybrid interface coupling, using either Cauchy stress-based or viscous
      stress-based Lagrange multipliers.
      The method ElementXfemInterfaceNIT() provides a clean and fast element implementation
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
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcXFEM : public FluidEleCalc<distype>
    {
      /// private Constructor since we are a Singleton.
      FluidEleCalcXFEM();

      /// private copy Constructor since we are a Singleton.
      FluidEleCalcXFEM(FluidEleCalcXFEM const& copy);

      /// private assignment operator since we are a Singleton.
      FluidEleCalcXFEM& operator=(FluidEleCalcXFEM const& copy);

      typedef FluidEleCalc<distype> my;

     private:
      /// pointer to the cast object, fluid parameter list for XFEM
      DRT::ELEMENTS::FluidEleParameterXFEM* fldparaxfem_;


      /// number of stress-dof
      static const int numstressdof_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      //! @name useful constants for DOF-index numbering
      //@{
      static const unsigned Velx = 0;
      static const unsigned Vely = 1;
      static const unsigned Velz = 2;
      static const unsigned Pres = 3;

      static const unsigned Sigmaxx = 0;
      static const unsigned Sigmaxy = 1;
      static const unsigned Sigmaxz = 2;
      static const unsigned Sigmayx = 1;
      static const unsigned Sigmayy = 3;
      static const unsigned Sigmayz = 4;
      static const unsigned Sigmazx = 2;
      static const unsigned Sigmazy = 4;
      static const unsigned Sigmazz = 5;
      //@}

      /// get stress dof-index
      unsigned stressIndex(unsigned xi, unsigned xj)
      {
        return (xi * xj > 0) ? xi + xj + 1 : xi + xj;
      }

     public:
      virtual ~FluidEleCalcXFEM() {}

      /// Singleton access method
      static FluidEleCalcXFEM<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// evaluate the XFEM cut element
      virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells, bool offdiag = false);

      /// evaluate the shape functions in the XFEM
      virtual int IntegrateShapeFunctionXFEM(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells);

      /// error computation
      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& ele_dom_norms);

      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& ele_dom_norms,
          const DRT::UTILS::GaussIntegration& intpoints);

      virtual int ComputeErrorInterface(DRT::ELEMENTS::Fluid* ele,   ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<MAT::Material>& mat,                          ///< material
          Epetra_SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,                               ///< boundary integration points
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                ///< parameter list
      );

      /// add terms from mixed/hybrid Lagrange multiplier coupling approach to element matrix and
      /// rhs
      virtual void ElementXfemInterfaceHybridLM(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,                         ///< side coupling matrices
          Teuchos::ParameterList& params,            ///< parameter list
          Teuchos::RCP<MAT::Material>& mat,          ///< material
          Epetra_SerialDenseMatrix& elemat1_epetra,  ///< local system matrix of intersected element
          Epetra_SerialDenseVector&
              elevec1_epetra,               ///< local element vector of intersected element
          Epetra_SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const GEO::CUT::plain_volumecell_set& vcSet  ///< set of plain volume cells
      );

      /// add Nitsche (NIT) interface condition to element matrix and rhs
      virtual void ElementXfemInterfaceNIT(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,               ///< parameter list
          Teuchos::RCP<MAT::Material>& mat_master,      ///< material for the coupled side
          Teuchos::RCP<MAT::Material>& mat_slave,       ///< material for the coupled side
          Epetra_SerialDenseMatrix& elemat1_epetra,     ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,     ///< element vector
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,               ///< side coupling matrices
          Epetra_SerialDenseMatrix& Cuiui  ///< ui-ui coupling matrix
      );

      ///
      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                    ///< discretization
          const std::vector<int>& lm,                                  ///< local map
          Epetra_SerialDenseVector& elevec1_epetra,                    ///< element vector
          const DRT::UTILS::GaussIntegration& intpoints                ///< integration points
      );

      ///
      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                    ///< discretization
          const std::vector<int>& lm,                                  ///< local map
          Epetra_SerialDenseVector& elevec1_epetra                     ///< element vector
      );

     private:
      //! evaluate analytical reference solution
      void AnalyticalReference(const int calcerr,      ///< which reference solution
          const int calcerrfunctno,                    ///< error function number
          LINALG::Matrix<my::nsd_, 1>& u,              ///< exact jump vector (coupled)
          LINALG::Matrix<my::nsd_, my::nsd_>& grad_u,  ///< exact velocity gradient
          double& p,                                   ///< exact pressure
          LINALG::Matrix<my::nsd_, 1>& xyzint,         ///< xyz position of gaussian point
          const double& t,                             ///< time
          Teuchos::RCP<MAT::Material> mat = Teuchos::null);

      //! get the interface jump vectors for velocity and traction at the Gaussian point
      void GetInterfaceJumpVectors(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          Teuchos::RCP<XFEM::CouplingBase> coupling,  ///< coupling object
          LINALG::Matrix<my::nsd_, 1>&
              ivelint_jump,  ///< prescribed interface jump vector for velocity
          LINALG::Matrix<my::nsd_, 1>&
              itraction_jump,  ///< prescribed interface jump vector for traction
          LINALG::Matrix<my::nsd_, my::nsd_>& proj_tangential,  ///< tangential projection matrix
          LINALG::Matrix<my::nsd_, my::nsd_>&
              LB_proj_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
          const LINALG::Matrix<my::nsd_, 1>& x,       ///< global coordinates of Gaussian point
          const LINALG::Matrix<my::nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>>
              si,                           ///< side implementation for cutter element
          LINALG::Matrix<3, 1>& rst,        ///< local coordinates of GP for bg element
          double& kappa_m,                  ///< fluid sided weighting
          double& visc_m,                   ///< fluid sided weighting
          double& visc_s,                   ///< slave sided dynamic viscosity
          LINALG::Matrix<3, 1>& rst_slave,  ///< local coord of gp in slave element
          std::vector<double>& eledisp,     ///< slave element displacement vector
          DRT::Element* coupl_ele = NULL    ///< slave coupling element
      );

      //! get the interface jump vectors for velocity and traction at the Gaussian point for
      //! previous time step
      void GetInterfaceJumpVectorsOldState(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          Teuchos::RCP<XFEM::CouplingBase> coupling,  ///< coupling object
          LINALG::Matrix<my::nsd_, 1>&
              ivelintn_jump,  ///< prescribed interface jump vector for velocity
          LINALG::Matrix<my::nsd_, 1>&
              itractionn_jump,                   ///< prescribed interface jump vector for traction
          const LINALG::Matrix<my::nsd_, 1>& x,  ///< global coordinates of Gaussian point
          const LINALG::Matrix<my::nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>>
              si,                    ///< side implementation for cutter element
          const double& presn_m,     ///< coupling master pressure
          LINALG::Matrix<3, 1>& rst  ///< local coordinates of GP for bg element
      );

      //! build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
      void NIT_BuildPatchCuiui(
          Epetra_SerialDenseMatrix&
              Cuiui,  ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              Cuiui_coupling  ///< Cuiui matrices for all cutting sides
      );

      //! compute viscous part of Nitsche's penalty term scaling for Nitsche's method
      double NIT_Compute_ViscPenalty_Stabfac(
          const DRT::Element::DiscretizationType
              ele_distype,      ///< the discretization type of the element w.r.t which the
                                ///< stabilization factor is computed
          const double inv_hk,  ///< the inverse characteristic element length
          const double kappa1,  ///< Weight parameter (parameter +/master side)
          const double kappa2   ///< Weight parameter (parameter -/slave  side)
      );

      //! get the constant which satisfies the trace inequality depending on the spatial dimension
      //! and polynomial order of the element
      double NIT_getTraceEstimateConstant(const DRT::Element::DiscretizationType ele_distype);

      //! prepare coupling matrices, that include contributions from convective stabilization
      void HybridLM_CreateSpecialContributionMatrices(
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          std::set<int>& begids,  ///< ids of intersecting boundary elements
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              conv_side_coupling  ///< contributions to coupling matrices from convective
                                  ///< stabilizations
      );

      //! evaluate Neumann boundary condition
      void EvaluateNeumann(const double& timefacfac,   ///< theta*dt
          const LINALG::Matrix<my::nen_, 1>& funct_m,  ///< coupling master shape functions
          const LINALG::Matrix<my::nsd_, 1>&
              itraction_jump,  ///< prescribed interface traction, jump height for coupled problems
          Epetra_SerialDenseMatrix& elevec1_epetra  ///< element vector
      );

      //! build traction vector w.r.t fluid domain
      void BuildTractionVector(LINALG::Matrix<my::nsd_, 1>& traction,  ///< traction vector
          double& press,                       ///< pressure at gaussian point
          LINALG::Matrix<my::nsd_, 1>& normal  ///< normal vector
      );

      //! assemble side's interface force
      void AssembleInterfaceForce(
          Teuchos::RCP<Epetra_Vector> iforcecol,  ///< interface force column vector
          DRT::Discretization& cutdis,            ///< cut discretization
          std::vector<int>& lm,                   ///< local dof map
          Epetra_SerialDenseVector& iforce        ///< interface force vector
      );

      //! evaluate shape function and derivative at point with local coordinates rst
      void EvalFuncAndDeriv(LINALG::Matrix<3, 1>& rst);

      //! build matrices from volume-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM-coupling \author kruse \date 06/14
      void HybridLM_Build_VolBased(const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,  ///< element velocity
          const LINALG::Matrix<my::nen_, 1>& epreaf,         ///< element pressure
          LINALG::Matrix<my::nen_, my::nen_>& bK_ss,         ///< block K_ss matrix
          LINALG::Matrix<my::nen_, my::nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, numstressdof_,
              my::numdofpernode_>& K_su,  ///< K_su matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::nsd_, my::nsd_>&
              K_uu,  ///< K_uu matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, my::nsd_, 1>&
              rhs_uu,                    ///< rhs_u(u) vector
          const bool is_MHVS = false,    ///< viscous (true) or Cauchy (false) stress-based LM
          const double mhvs_param = 1.0  ///< stabilizing parameter for viscous stress-based LM
      );

      //! evaluate matrices from volume-based terms for Cauchy stress-based mixed/hybrid LM coupling
      //! at current Gauss-point
      void MHCS_Evaluate_VolBased(
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,  ///< element velocity
          const LINALG::Matrix<my::nen_, 1>& epreaf,         ///< element pressure
          LINALG::Matrix<my::nen_, my::nen_>& bK_ss,         ///< block K_ss matrix
          LINALG::Matrix<my::nen_, my::nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, numstressdof_,
              my::numdofpernode_>& K_su,  ///< K_su matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, numstressdof_, 1>&
              rhs_s  ///< rhs_s vector
      );

      //! evaluate matrices from volume-based terms for viscous stress-based mixed/hybrid LM
      //! coupling at current Gauss-point \author kruse \date 06/14
      void MHVS_Evaluate_VolBased(
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,  ///< element velocity
          LINALG::Matrix<my::nen_, my::nen_>& bK_ss,         ///< block K_ss matrix
          LINALG::Matrix<my::nen_, my::nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, numstressdof_,
              my::numdofpernode_>& K_su,  ///< K_su matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::nsd_, my::nsd_>&
              K_uu,                                                               ///< K_uu matrix
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, my::nsd_, 1>& rhs_uu,  ///< rhs_uu vector
          const double& mhvs_param  ///< stabilizing parameter
      );

      //! evaluate matrices from surface-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM coupling at current Gauss-point \author kruse \date 06/14
      void HybridLM_Evaluate_SurfBased(
          Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>>& si,
          const LINALG::Matrix<my::nen_, my::nen_>& bK_ss,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, numstressdof_,
              my::numdofpernode_>& K_su,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::numdofpernode_,
              numstressdof_>& K_us,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, numstressdof_, 1>& rhs_s,
          const LINALG::Matrix<my::nen_, 1>& epreaf,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::nsd_, my::nsd_>& K_uu,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, my::nsd_, 1>& rhs_uu,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::nsd_, 1>& G_up,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, 1, my::nsd_>& G_pu,
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_, 1>, my::nsd_, 1>& rhs_up,
          LINALG::Matrix<my::nen_, 1>& rhs_pu, const LINALG::Matrix<my::nsd_, 1>& normal,
          const double& timesurffac, const LINALG::Matrix<my::nsd_, 1>& ivelint_jump,
          const LINALG::Matrix<my::nsd_, 1>& itraction_jump, const bool eval_side_coupling,
          const bool is_MHVS);

      //! Initiates dummy variables and calls FLuidEleCalc GetMaterialParams routine.
      void GetMaterialParametersVolumeCell(Teuchos::RCP<const MAT::Material> material,
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

      LINALG::Matrix<my::nsd_, my::nen_>
          evelaf_;  ///< element velocity at time t^n+af, implemented also in base class, not
                    ///< accessable via ExtractValuesFromGlobalVector
      LINALG::Matrix<my::nen_, 1>
          epreaf_;  ///< element pressure at time t^n+af, implemented also in base class, not
                    ///< accessable via ExtractValuesFromGlobalVector

      LINALG::Matrix<my::nsd_, my::nen_>
          eveln_;  ///< element velocity at time t^n, implemented also in base class, not accessable
                   ///< via ExtractValuesFromGlobalVector
      LINALG::Matrix<my::nen_, 1>
          epren_;  ///< element velocity at time t^n, implemented also in base class, not accessable
                   ///< via ExtractValuesFromGlobalVector

      LINALG::Matrix<my::nsd_, 1> ivelint_jump_;    ///< interface velocity jump at t^n+1
      LINALG::Matrix<my::nsd_, 1> itraction_jump_;  ///< interface traction jump at t^n+1
      LINALG::Matrix<my::nsd_, my::nsd_>
          proj_tangential_;  ///< tangential projection matrix at t^n+1
      LINALG::Matrix<my::nsd_, my::nsd_>
          LB_proj_matrix_;  ///< interface matrix jump (for Laplace-Beltrami) at t^n+1

      std::vector<Epetra_SerialDenseMatrix>
          solid_stress_;  ///< hold information about solid stress ([0]...traction,
                          ///< [1]...dtraction_dv, [2-4]...d2traction_dv2)

      LINALG::Matrix<my::nsd_, 1> ivelintn_jump_;    ///< interface velocity jump at t^n
      LINALG::Matrix<my::nsd_, 1> itractionn_jump_;  ///< interface traction jump at t^n

      LINALG::Matrix<my::nsd_, 1>
          velint_s_;  ///< velocity of slave side at time t^n+1 at integration point
      LINALG::Matrix<my::nsd_, 1>
          velintn_s_;  ///< velocity of slave side at time t^n at integration point

      LINALG::Matrix<my::nsd_, 1>
          rst_;  ///< local coordinates of Gaussian point w.r.t background element

      LINALG::Matrix<my::nsd_, 1> normal_;    ///< normal vector
      LINALG::Matrix<my::nsd_, 1> x_side_;    ///< gauss-point coordinates
      LINALG::Matrix<my::nsd_, 1> x_gp_lin_;  ///< gauss-point in xyz-system on linearized interface
    };


  }  // namespace ELEMENTS
}  // namespace DRT

#endif
