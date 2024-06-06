/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG fluid element


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_HDG_HPP
#define FOUR_C_FLUID_ELE_CALC_HDG_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_shapevalues_hdg.hpp"
#include "4C_fluid_ele_hdg.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Fluid HDG element implementation
    /*!

      \author kronbichler
      \date 05/13
    */
    template <Core::FE::CellType distype>
    class FluidEleCalcHDG : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr unsigned int nen_ = Core::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr unsigned int nsd_ = Core::FE::dim<distype>;

      ///! number of faces on element
      static constexpr unsigned int nfaces_ = Core::FE::num_faces<distype>;


      int integrate_shape_function(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const Core::FE::GaussIntegration& intpoints) override
      {
        FOUR_C_THROW("Not implemented!");
        return 1;
      }

      int integrate_shape_function_xfem(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells) override
      {
        FOUR_C_THROW("Not implemented!");
        return 1;
      };


      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */
      int EvaluateService(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
        \brief calculate dissipation of various terms (evaluation of turbulence models)
      */
      virtual int calc_dissipation(Fluid* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat)
      {
        FOUR_C_THROW("Not implemented!");
        return 1;
      }

      /// Evaluate element ERROR
      /*!
          general function to compute the error (analytical solution) for particular problem type
       */
      virtual int compute_error(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec);

      int compute_error(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec,
          const Core::FE::GaussIntegration&) override
      {
        return compute_error(ele, params, mat, discretization, lm, elevec);
      }

      /// projection of function field
      virtual int ProjectField(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      /*! \brief Interpolates an HDG solution to the element nodes for output
       */
      virtual int interpolate_solution_to_nodes(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1);

      /*! \brief Interpolates an HDG solution for homogeneous isotropic turbulence postprocessing
       */
      virtual int interpolate_solution_for_hit(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1);

      /*! \brief Project force from equidistant points on interior node dof vector
       */
      virtual int project_force_on_dof_vec_for_hit(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      /*! \brief Project initial field for hit
       */
      virtual int project_initial_field_for_hit(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3);
      /*!
      \brief Initialize the shape functions and solver to the given element (degree is runtime
      parameter)
       */
      void InitializeShapes(const Discret::ELEMENTS::Fluid* ele);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      int Evaluate(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points
      int Evaluate(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints, bool offdiag = false) override;

      int compute_error_interface(Discret::ELEMENTS::Fluid* ele,     ///< fluid element
          Discret::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<Core::Mat::Material>& mat,                    ///< material
          Core::LinAlg::SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,                                     ///< boundary integration points
          const Core::Geo::Cut::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                      ///< parameter list
          ) override
      {
        FOUR_C_THROW("Not implemented!");
        return 1;
      }

      /// Evaluate the XFEM cut element
      int EvaluateXFEM(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells, bool offdiag = false) override
      {
        FOUR_C_THROW("Not implemented!");
        return 1;
      }


      void element_xfem_interface_hybrid_lm(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Discret::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<Core::FE::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Teuchos::ParameterList& params,          ///< parameter list
          Teuchos::RCP<Core::Mat::Material>& mat,  ///< material
          Core::LinAlg::SerialDenseMatrix&
              elemat1_epetra,  ///< local system matrix of intersected element
          Core::LinAlg::SerialDenseVector&
              elevec1_epetra,                      ///< local element vector of intersected element
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const Core::Geo::Cut::plain_volumecell_set& vcSet  ///< set of plain volume cells
          ) override
      {
        FOUR_C_THROW("Not implemented!");
        return;
      }

      void element_xfem_interface_nit(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Discret::Discretization& dis,                               ///< background discretization
          const std::vector<int>& lm,                                 ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,   ///< XFEM condition manager
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,                     ///< parameter list
          Teuchos::RCP<Core::Mat::Material>& mat_master,      ///< material master side
          Teuchos::RCP<Core::Mat::Material>& mat_slave,       ///< material slave side
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,    ///< element matrix
          Core::LinAlg::SerialDenseVector& elevec1_epetra,    ///< element vector
          const Core::Geo::Cut::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
          ) override
      {
        FOUR_C_THROW("Not implemented!");
        return;
      }

      void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& dis,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
          const Core::FE::GaussIntegration& intpoints) override
      {
        FOUR_C_THROW("Not implemented!");
        return;
      }

      void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& dis,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra) override
      {
        FOUR_C_THROW("Not implemented!");
        return;
      }

      /// Evaluate the pressure average inside the element from an analytical expression
      virtual int evaluate_pressure_average(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseVector& elevec);

      /// Singleton access method
      static FluidEleCalcHDG<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      /// Print local residuals
      void PrintLocalResiduals(Discret::ELEMENTS::Fluid* ele);

      /// Print local variables
      void PrintLocalVariables(Discret::ELEMENTS::Fluid* ele);

      /// Print local correction term
      void print_local_correction(
          Discret::ELEMENTS::Fluid* ele, std::vector<double>& interiorecorrectionterm);

      /// Print local body force
      void PrintLocalBodyForce(
          Discret::ELEMENTS::Fluid* ele, std::vector<double>& interiorebodyforce);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcHDG();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        static constexpr unsigned int nsd_ = FluidEleCalcHDG<distype>::nsd_;
        static constexpr unsigned int nfaces_ = FluidEleCalcHDG<distype>::nfaces_;

        LocalSolver(const Discret::ELEMENTS::Fluid* ele,
            const Core::FE::ShapeValues<distype>& shapeValues,
            Core::FE::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly);

        void compute_interior_residual(const Teuchos::RCP<Core::Mat::Material>& mat,
            const std::vector<double>& valnp, const std::vector<double>& accel,
            const double avgPressure, const Core::LinAlg::Matrix<nsd_, nen_>& ebodyforce,
            const std::vector<double>& intebodyforce, Core::LinAlg::SerialDenseVector& eleVec,
            const std::vector<double>& interiorecorrectionterm,
            const std::vector<double>& interiorebodyforce);

        void ComputeFaceResidual(const int face, const Teuchos::RCP<Core::Mat::Material>& mat,
            const std::vector<double>& val, const std::vector<double>& traceval,
            Core::LinAlg::SerialDenseVector& eleVec);

        void compute_interior_matrices(
            const Teuchos::RCP<Core::Mat::Material>& mat, const bool evaluateOnlyNonlinear);

        void ComputeFaceMatrices(const int face, const Teuchos::RCP<Core::Mat::Material>& mat,
            const bool evaluateOnlyNonlinear, Core::LinAlg::SerialDenseMatrix& elemat);

        // inverts the velocity gradient matrix and puts its contribution into the velocity matrix
        // (pre-factorization). Should only be done once per element even if multiple velocities are
        // used
        void eliminate_velocity_gradient(Core::LinAlg::SerialDenseMatrix& elemat);

        // solves the local problem, including factorization of the matrix
        void SolveResidual();

        // condense the local matrix (involving cell velocity gradients, velocities and pressure)
        // into the element matrix for the trace and similarly for the residuals
        void CondenseLocalPart(
            Core::LinAlg::SerialDenseMatrix& elemat, Core::LinAlg::SerialDenseVector& elevec);

        // compute the correction term on the rhs for the weakly compressible benchmark
        void compute_correction_term(
            std::vector<double>& interiorecorrectionterm, int corrtermfuncnum);

        // compute the body force on the rhs for the weakly compressible benchmark
        void ComputeBodyForce(std::vector<double>& interiorebodyforce, int bodyforcefuncnum);

        const unsigned int ndofs_;

        bool stokes;

        bool weaklycompressible;

        // convention: we sort the entries in the matrices the following way:
        // first come the velocity gradients, then the velocities, and finally the pressure
        // we also build the matrix in a block-fashion, keeping the dofs for individual components
        // closest to each other. I.e. the blocks are in 2D for g_00, g_01, g_10, g_11, v_0, v_1, p
        // and similarly for 3D

        const Core::FE::ShapeValues<distype>& shapes_;    /// evaluated shape values
        Core::FE::ShapeValuesFace<distype>& shapesface_;  /// evaluated shape values

        double stabilization[nfaces_];  /// stabilization parameters

        Core::LinAlg::SerialDenseMatrix
            uuMat;  /// terms for block with velocity and pressure (constant ones)
        Core::LinAlg::SerialDenseMatrix uuMatFinal;  /// terms for block with velocity and pressure
                                                     /// (including convection and stabilization)
        Core::LinAlg::SerialDenseMatrix
            ugMat;  /// coupling between velocity and velocity gradient (not fully stored)
        Core::LinAlg::SerialDenseMatrix
            guMat;  /// evaluated divergence of velocity gradient and velocity (not fully stored)

        Core::LinAlg::SerialDenseMatrix
            gfMat;  /// evaluated coupling between velocity gradient and trace
        Core::LinAlg::SerialDenseMatrix
            fgMat;  /// evaluated coupling between trace and velocity gradient
        Core::LinAlg::SerialDenseMatrix ufMat;  /// evaluated coupling between velocity and trace
        Core::LinAlg::SerialDenseMatrix fuMat;  /// evaluated coupling between trace and velocity

        Core::LinAlg::SerialDenseMatrix
            massPart;  /// temporary matrix for mass matrix on all quadrature points
        Core::LinAlg::SerialDenseMatrix massPartW;  /// temporary matrix for mass matrix weighted by
                                                    /// integration factor on all quadrature points
        Core::LinAlg::SerialDenseMatrix
            gradPart;  /// temporary matrix for gradient matrix on all quadrature points
        Core::LinAlg::SerialDenseMatrix uPart;  /// temporary matrix for convection

        Core::LinAlg::SerialDenseMatrix
            massMat;  /// local mass matrix (will be inverted during init)
        Core::LinAlg::SerialDenseMatrix uuconv;      /// convection matrix
        Core::LinAlg::SerialDenseMatrix tmpMat;      /// matrix holding temporary results
        Core::LinAlg::SerialDenseMatrix tmpMatGrad;  /// matrix holding temporary results

        Core::LinAlg::SerialDenseMatrix trMat;     /// temporary matrix for trace assembly
        Core::LinAlg::SerialDenseMatrix trMatAvg;  /// temporary matrix for trace assembly

        Core::LinAlg::SerialDenseMatrix velnp;  /// velocities evaluated on all quadrature points
        Core::LinAlg::SerialDenseMatrix
            fvelnp;  /// trace velocities evaluated on all face quadrature points

        Core::LinAlg::SerialDenseMatrix uucomp;  /// compressibility matrix
        Core::LinAlg::SerialDenseVector presnp;  /// pressure evaluated on all quadrature points
        Core::LinAlg::SerialDenseMatrix
            gradpresnp;  /// pressure gradient evaluated on all quadrature points
        Core::LinAlg::SerialDenseVector
            ifpresnp;  /// pressure evaluated on all face quadrature points

        Core::LinAlg::SerialDenseVector gRes;   /// residual vector on velocity gradients
        Core::LinAlg::SerialDenseVector upRes;  /// residual vector on velocity and pressure
        Core::LinAlg::SerialDenseVector gUpd;   /// update vector for velocity gradients
        Core::LinAlg::SerialDenseVector upUpd;  /// update vector for velocity and pressure

        std::vector<int> pivots;  /// pivots for factorization of matrices

        Teuchos::RCP<Discret::ELEMENTS::FluidEleParameter> fldpara_;  //! pointer to parameter list
        Teuchos::RCP<Discret::ELEMENTS::FluidEleParameterTimInt>
            fldparatimint_;  //! pointer to time parameter list
      };

      /// reads from global vectors
      void read_global_vectors(const Core::Elements::Element& ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          const bool updateLocally);

      // writes the updated solution vector to the secondary vector stored in the discretization
      void update_secondary_solution(const Core::Elements::Element& ele,
          Discret::Discretization& discretization, const Core::LinAlg::SerialDenseVector& updateG,
          const Core::LinAlg::SerialDenseVector& updateUp);

      void evaluate_velocity(const int startfunc, const Inpar::FLUID::InitialField initfield,
          const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::Matrix<nsd_, 1>& u) const;

      void evaluate_all(const int startfunc, const Inpar::FLUID::InitialField initfield,
          const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::Matrix<nsd_, 1>& u,
          Core::LinAlg::Matrix<nsd_, nsd_>& grad, double& p) const;

      /// local data object
      Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes_;
      Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

      /// local solver object
      Teuchos::RCP<LocalSolver> local_solver_;

      Core::LinAlg::Matrix<nsd_, nen_> ebofoaf_;     /// body force (see fluid_ele_calc.cpp)
      Core::LinAlg::Matrix<nsd_, nen_> eprescpgaf_;  /// pressure gradient body force
      Core::LinAlg::Matrix<nen_, 1> escabofoaf_;     /// scalar body force for loma
      std::vector<double> interiorebofoaf_;          /// extracted body force at n+alpha_f

      std::vector<double>
          interiorecorrectionterm_;  /// local correction term for the weakly compressible benchmark
      std::vector<double>
          interiorebodyforce_;  /// local body force for the weakly compressible benchmark

      std::vector<double> trace_val_;  /// extracted values from trace solution vector at n+alpha_f
      std::vector<double> interior_val_;  /// extracted local values (velocity gradients,
                                          /// velocities, pressure) at n+alpha_f
      std::vector<double> interior_acc_;  /// extracted local accelerations at n+alpha_m

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
