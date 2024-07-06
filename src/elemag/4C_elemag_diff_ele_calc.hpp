/*----------------------------------------------------------------------*/
/*! \file

\brief All functionality for electromagnetic diffusionm element evaluations

<pre>
\level 2

</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_DIFF_ELE_CALC_HPP
#define FOUR_C_ELEMAG_DIFF_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_elemag_diff_ele.hpp"
#include "4C_elemag_ele_interface.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_shapevalues_hdg.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Electromagnetic diffusion element implementation

    template <Core::FE::CellType distype>
    class ElemagDiffEleCalc : public ElemagEleInterface
    {
     public:
      /// nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr unsigned int nen_ = Core::FE::num_nodes<distype>;

      /// number of space dimensions
      static constexpr unsigned int nsd_ = Core::FE::dim<distype>;

      /// number of faces on element
      static constexpr unsigned int nfaces_ = Core::FE::num_faces<distype>;

      /// Integrate shape functions (not implemented)
      int integrate_shape_function(Discret::ELEMENTS::Elemag* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1) override
      {
        FOUR_C_THROW("Not implemented");
        return 1;
      };

      /// Zero initialization of elements.
      virtual void element_init(Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params);

      /// Interpolates an HDG solution to the element nodes for output.
      virtual int interpolate_solution_to_nodes(Discret::ELEMENTS::ElemagDiff* ele,
          Core::FE::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1);

      /// Initialize the shape functions and solver to the given element (degree is runtime
      /// parameter).
      void initialize_shapes(const Discret::ELEMENTS::ElemagDiff* ele);

      /// Evaluate the element.
      /// Generic virtual interface function.Called via base pointer.
      int evaluate(Discret::ELEMENTS::Elemag* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points.
      virtual int evaluate(Discret::ELEMENTS::Elemag* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints, bool offdiag = false);

      /// Singleton access method.
      static ElemagDiffEleCalc<distype>* instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      /// Used to print the trace values as debugging utility.
      void print_trace(Core::Elements::Element* ele);


     private:
      /// private Constructor since we are a Singleton.
      ElemagDiffEleCalc();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        /// Number of Spatial Dimensions
        static constexpr unsigned int nsd_ = ElemagDiffEleCalc<distype>::nsd_;
        /// Number of FACES
        static constexpr unsigned int nfaces_ = ElemagDiffEleCalc<distype>::nfaces_;

        /// Init function for the struct members
        LocalSolver(const Discret::ELEMENTS::ElemagDiff* ele,
            Core::FE::ShapeValues<distype>& shapeValues,
            Teuchos::RCP<Core::FE::ShapeValuesFace<distype>>& shapeValuesFace,
            Inpar::EleMag::DynamicType& dyna, Core::FE::ShapeValues<distype>& postproc_shapeValues);

        /// Compute the residual
        void compute_residual(Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& eleVec, double dt, Discret::ELEMENTS::ElemagDiff& ele);

        /// Computes the source term in the element.
        void compute_source(Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& interiorSourcen,
            Core::LinAlg::SerialDenseVector& interiorSourcenp);

        /// Add terms corresponding to the absorbing boundary condition.
        void compute_absorbing_bc(Core::FE::Discretization& discretization,
            Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            Teuchos::RCP<Core::Mat::Material>& mat, int face,
            Core::LinAlg::SerialDenseMatrix& elemat, int indexstart,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Add terms corresponding to the absorbing boundary condition.
        void compute_boundary_integral(
            Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params, int face);

        /// Calls local solver to compute matrices: internal and face
        void compute_matrices(Core::FE::Discretization& discretization,
            const Teuchos::RCP<Core::Mat::Material>& mat, Discret::ELEMENTS::ElemagDiff& ele,
            double dt, Inpar::EleMag::DynamicType dyna, const double tau);

        /// Set up interior matrices
        void compute_interior_matrices(double dt, double sigma, double mu, double epsilon);

        /// Set up face matrices
        void compute_face_matrices(const int face, double dt, int indexstart, int newindex,
            double sigma, double mu, const double tau);

        /// Condense the local matrix into the element matrix for the trace and similarly for the
        /// residuals.
        void condense_local_part(Core::LinAlg::SerialDenseMatrix& elemat);

        /*!
        \brief Make the matrix symmetric in case of dirichlet boundary conditions in the element

        \note So far the function only acts on the matrices and not on the RHS! This is ok for PEC
        boundary conditions because the DBC is anyway zero.
        */
        void symmetrify(Discret::ELEMENTS::ElemagDiff& ele, Core::LinAlg::SerialDenseMatrix& elemat,
            bool dodirich = true);

        /// Projection of function field.
        /// The function is used to project the field in the initialization phase.
        int project_field(Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2);

        /// Projection of scatra computed gradients into elements
        int project_electric_field_from_scatra(Discret::ELEMENTS::ElemagDiff* ele,
            Teuchos::ParameterList& params, const Teuchos::RCP<Core::Mat::Material>& mat,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Compute the error with respect to an analytical field.
        void compute_error(Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Postprocess the solution to obtain a better approximation
        void post_processing(Discret::ELEMENTS::ElemagDiff& ele);

        /*!
         * \brief Projection of a given field on the interior variables for testing purposes.
         *
         * This is a debug utility necessary for code development and prototyping. The selection
         * of the field to project is hard-coded whithin the method itself with the flags
         * do_internal and do_postprocess. The flags are used as follows:
         *  - do_internal projects the field in the internal dofs of the element
         *  - do_postprocess projects the field in the postprocessed internal dofs of the element
         *
         * \note It is possible to project one or both the variables, depending on the debugging
         * needs.
         */
        int project_field_test(Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2);

        /// Projection of a given field on the trace for testing purposes.
        int project_field_test_trace(Discret::ELEMENTS::ElemagDiff* ele,
            Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector& elevec1);

        /// Projection of Dirichlet function field.
        int project_dirich_field(Discret::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Function evaluation routine
        void evaluate_all(const int start_func, const double t,
            const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::SerialDenseVector& v) const;

        /// Function gradient evaluation routine
        void compute_function_gradient(int start_func, double t,
            const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::SerialDenseMatrix& v) const;

        /// Function time derivative evaluation routine
        void compute_function_time_derivative(int start_func, double t, double dt,
            const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::SerialDenseVector& v) const;

        // convention: we sort the entries in the matrices the following way:
        // first comes the electric field then the magnetic one.

        /// Number of Degrees Of Freedom
        const unsigned int ndofs_;

        /// evaluated shape values
        Core::FE::ShapeValues<distype>& shapes_;

        /// evaluated shape values
        Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

        /// evaluated higher degree shape values
        Core::FE::ShapeValues<distype>& postproc_shapes_;

        // System matrices
        /// u matrix
        Core::LinAlg::SerialDenseMatrix Amat;
        /// Inverse of u matrix
        Core::LinAlg::SerialDenseMatrix invAmat;
        /// u electric
        Core::LinAlg::SerialDenseMatrix Bmat;
        /// u trace
        Core::LinAlg::SerialDenseMatrix Cmat;
        Core::LinAlg::SerialDenseMatrix Fmat;
        /// electric second time
        Core::LinAlg::SerialDenseMatrix Dmat;
        /// electric evolution
        Core::LinAlg::SerialDenseMatrix Emat;
        /// electric electric
        Core::LinAlg::SerialDenseMatrix Gmat;
        /// electric trace
        Core::LinAlg::SerialDenseMatrix Hmat;
        /// trace trace
        Core::LinAlg::SerialDenseMatrix Lmat;

        // auxiliary stuff
        /// final mass matrix used for the projection
        Core::LinAlg::SerialDenseMatrix massMat;
        /// part of the mass matrix (only contains the shape functions)
        Core::LinAlg::SerialDenseMatrix massPart;
        /// other part of the mass matrix (with quadrature weights)
        Core::LinAlg::SerialDenseMatrix massPartW;

        /// Chosen dynamics/time integrator
        Inpar::EleMag::DynamicType& dyna_;
      };  // stuct LocalSolver

      /// Updates interior variables and calculates residual.
      void update_interior_variables_and_compute_residual(Teuchos::ParameterList& params,
          Discret::ELEMENTS::ElemagDiff& ele, const Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseVector& elevec, double dt, bool errormaps, bool updateonly);

      /// Reads from global vectors.
      void read_global_vectors(Core::Elements::Element* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm);

      /// Writes internal fields from elements to global vectors.
      void fill_restart_vectors(
          Core::Elements::Element* ele, Core::FE::Discretization& discretization);

      /// Reads internal field from global vectors to element vectors.
      void element_init_from_restart(
          Core::Elements::Element* ele, Core::FE::Discretization& discretization);

      /// Calculate error maps with local postprocessing.
      double estimate_error(Discret::ELEMENTS::ElemagDiff& ele, Core::LinAlg::SerialDenseVector& p);

      /// Local data object for element.
      Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes_;

      /// Higher degree local data object for element.
      Teuchos::RCP<Core::FE::ShapeValues<distype>> postproc_shapes_;

      /// Local data object for face element.
      Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

      /// Local solver object.
      Teuchos::RCP<LocalSolver> local_solver_;

      std::vector<double> localtrace_;  /// extracted values from trace solution vector

      /// Local values from interior solution vector (electric and magnetic fields) at n
      Core::LinAlg::SerialDenseVector interior_electricnp_;
      Core::LinAlg::SerialDenseVector interior_magneticnp_;
      Core::LinAlg::SerialDenseVector interiorauxiliary_pml_;
      /// Local values from interior solution vector (electric and magnetic fields) at n-1
      Core::LinAlg::SerialDenseVector interior_electricnm_;
      Core::LinAlg::SerialDenseVector interior_magneticnm_;


      /// Chosen dynamics/time integrator
      Inpar::EleMag::DynamicType dyna_;

      /// Use complete polynomial or tensor product
      bool usescompletepoly_;
    };  // class ElemagDiffEleCalc
  }     // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
