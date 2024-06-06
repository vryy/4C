/*--------------------------------------------------------------------------*/
/*! \file

\brief All functionality for electromagnetic element evaluations

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_ELE_CALC_HPP
#define FOUR_C_ELEMAG_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_shapevalues_hdg.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_elemag_ele_interface.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Electromagnetic element implementation

    template <Core::FE::CellType distype>
    class ElemagEleCalc : public ElemagEleInterface
    {
     public:
      /// nen_: number of element nodes (T. Hughes: The Finite Element Method).
      static constexpr unsigned int nen_ = Core::FE::num_nodes<distype>;

      /// Number of space dimensions.
      static constexpr unsigned int nsd_ = Core::FE::dim<distype>;

      /// Number of faces on element.
      static constexpr unsigned int nfaces_ = Core::FE::num_faces<distype>;

      int integrate_shape_function(Discret::ELEMENTS::Elemag* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1) override
      {
        FOUR_C_THROW("Not implemented");
        return 1;
      };

      /// Zero initialization of elements.
      virtual void ElementInit(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params);

      /// Interpolates an HDG solution to the element nodes for output.
      virtual int interpolate_solution_to_nodes(Discret::ELEMENTS::Elemag* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1);

      /// Initialize the shape functions and solver to the given element (degree is runtime
      /// parameter).
      void InitializeShapes(const Discret::ELEMENTS::Elemag* ele);

      /// Evaluate the element.
      /// Generic virtual interface function.Called via base pointer.
      int Evaluate(Discret::ELEMENTS::Elemag* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points.
      virtual int Evaluate(Discret::ELEMENTS::Elemag* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints, bool offdiag = false);



      /// Singleton access method.
      static ElemagEleCalc<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      /// Used to print the trace values as debugging utility.
      void PrintTrace(Core::Elements::Element* ele);


     private:
      /// Private Constructor since we are a Singleton.
      ElemagEleCalc();

      /// Local solver that inverts local problem on an element and can solve with various vectors.
      struct LocalSolver
      {
        /// Number of Spatial Dimensions
        static constexpr unsigned int nsd_ = ElemagEleCalc<distype>::nsd_;
        /// Number of FACES
        static constexpr unsigned int nfaces_ = ElemagEleCalc<distype>::nfaces_;

        /// Init function for the struct members
        LocalSolver(const Discret::ELEMENTS::Elemag* ele,
            Core::FE::ShapeValues<distype>& shapeValues,
            Teuchos::RCP<Core::FE::ShapeValuesFace<distype>>& shapeValuesFace,
            Inpar::EleMag::DynamicType& dyna);

        /// Compute the residual
        void ComputeResidual(Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& eleVec, Discret::ELEMENTS::Elemag& ele);

        /// Computes the source term in the element.
        void ComputeSource(Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& interiorSourcen,
            Core::LinAlg::SerialDenseVector& interiorSourcenp);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeAbsorbingBC(Discret::Discretization& discretization,
            Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Teuchos::RCP<Core::Mat::Material>& mat, int face,
            Core::LinAlg::SerialDenseMatrix& elemat, int indexstart,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Add terms corresponding to the absorbing boundary condition.
        void compute_boundary_integral(
            Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params, int face);

        /// Calls local solver to compute matrices: internal and face
        void ComputeMatrices(Discret::Discretization& discretization,
            const Teuchos::RCP<Core::Mat::Material>& mat, Discret::ELEMENTS::Elemag& ele, double dt,
            Inpar::EleMag::DynamicType dyna, const double tau);

        /// Set up interior matrices
        void compute_interior_matrices(double dt, double sigma, double mu, double epsilon);

        /// Set up face matrices
        void ComputeFaceMatrices(const int face, double dt, int indexstart, int newindex,
            double sigma, double mu, const double tau);

        /// Condense the local matrx into the element matrix for the trace and similarly for the
        /// residuals.
        void CondenseLocalPart(Core::LinAlg::SerialDenseMatrix& elemat);

        /// Projection of function field.
        /// The function is used to project the field in the initialization phase.
        int ProjectField(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2);

        /// Compute the error with respect to an analytical field.
        void compute_error(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Projection of a given field on the interior variables for testing purposes.
        int ProjectFieldTest(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2);

        /// Projection of a given field on the trace for testing purposes.
        int project_field_test_trace(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Projection of Dirichlet function field.
        int ProjectDirichField(Discret::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Core::LinAlg::SerialDenseVector& elevec1);

        /// Function evaluation routine
        void evaluate_all(const int start_func, const double t,
            const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::SerialDenseVector& v) const;

        // convention: we sort the entries in the matrices the following way:
        // first comes the electric field then the magnetic one.

        /// Number of Degrees Of Freedom
        const unsigned int ndofs_;

        /// evaluated shape values
        Core::FE::ShapeValues<distype>& shapes_;
        /// evaluated shape values
        Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

        // System matrices
        /// magnetic evolution matrix
        Core::LinAlg::SerialDenseMatrix Amat;
        /// Inverse of magnetic evolution matrix
        Core::LinAlg::SerialDenseMatrix invAmat;
        /// magnetic electric
        Core::LinAlg::SerialDenseMatrix Cmat;
        /// magnetic trace
        Core::LinAlg::SerialDenseMatrix Dmat;
        /// electric evolution
        Core::LinAlg::SerialDenseMatrix Emat;
        /// electric magnetic
        Core::LinAlg::SerialDenseMatrix Fmat;
        /// electric electric
        Core::LinAlg::SerialDenseMatrix Gmat;
        /// electric trace
        Core::LinAlg::SerialDenseMatrix Hmat;
        /// trace magnetic
        Core::LinAlg::SerialDenseMatrix Imat;
        /// trace electric
        Core::LinAlg::SerialDenseMatrix Jmat;
        /// trace trace
        Core::LinAlg::SerialDenseMatrix Lmat;

        // auxiliary stuff
        Core::LinAlg::SerialDenseMatrix massMat;  // final mass matrix used for the projection
        Core::LinAlg::SerialDenseMatrix
            massPart;  // part of the mass matrix (only contains the shape functions)
        Core::LinAlg::SerialDenseMatrix
            massPartW;  // other part of the mass matrix (with quadrature weights)

        /// Chosen dynamics/time integrator
        Inpar::EleMag::DynamicType& dyna_;
      };

      /// Updates interior variables and calculates residual.
      void update_interior_variables_and_compute_residual(Teuchos::ParameterList& params,
          Discret::ELEMENTS::Elemag& ele, const Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseVector& elevec, double dt, bool errormaps, bool updateonly);

      /// Reads from global vectors.
      void read_global_vectors(Core::Elements::Element* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm);

      /// Writes internal fields from elements to global vectors.
      void fill_restart_vectors(
          Core::Elements::Element* ele, Discret::Discretization& discretization);

      /// Reads internal field from global vectors to element vectors.
      void element_init_from_restart(
          Core::Elements::Element* ele, Discret::Discretization& discretization);

      /// Calculate error maps with local postprocessing.
      double estimate_error(Discret::ELEMENTS::Elemag& ele, Core::LinAlg::SerialDenseVector& p);

      /// Local data object for element.
      Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes_;
      /// Local data object for face element.
      Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

      /// Local solver object.
      Teuchos::RCP<LocalSolver> local_solver_;

      std::vector<double> localtrace_;  /// extracted values from trace solution vector

      /// Local values from interior solution vector (electric and magnetic fields) at n
      Core::LinAlg::SerialDenseVector interior_electricnp_;
      Core::LinAlg::SerialDenseVector interior_magneticnp_;
      /// Local values from interior solution vector (electric and magnetic fields) at n-1
      Core::LinAlg::SerialDenseVector interior_electricnm_;
      Core::LinAlg::SerialDenseVector interior_magneticnm_;
      Core::LinAlg::SerialDenseVector interiorauxiliary_pml_;

      /// Chosen dynamics/time integrator
      Inpar::EleMag::DynamicType dyna_;

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
