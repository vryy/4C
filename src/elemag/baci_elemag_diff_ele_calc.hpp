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

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_discretization_fem_general_utils_shapevalues_hdg.hpp"
#include "baci_elemag_diff_ele.hpp"
#include "baci_elemag_ele_interface.hpp"
#include "baci_inpar_elemag.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    /// Electromagnetic diffusion element implementation

    template <CORE::FE::CellType distype>
    class ElemagDiffEleCalc : public ElemagEleInterface
    {
     public:
      /// nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr unsigned int nen_ = CORE::FE::num_nodes<distype>;

      /// number of space dimensions
      static constexpr unsigned int nsd_ = CORE::FE::dim<distype>;

      /// number of faces on element
      static constexpr unsigned int nfaces_ = CORE::FE::num_faces<distype>;

      /// Integrate shape functions (not implemented)
      int IntegrateShapeFunction(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1) override
      {
        dserror("Not implemented");
        return 1;
      };

      /// Zero initialization of elements.
      virtual void ElementInit(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params);

      /// Interpolates an HDG solution to the element nodes for output.
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::ElemagDiff* ele,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& elevec1);

      /// Initialize the shape functions and solver to the given element (degree is runtime
      /// parameter).
      void InitializeShapes(const DRT::ELEMENTS::ElemagDiff* ele);

      /// Evaluate the element.
      /// Generic virtual interface function.Called via base pointer.
      int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points.
      virtual int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const CORE::FE::GaussIntegration& intpoints, bool offdiag = false);

      /// Singleton access method.
      static ElemagDiffEleCalc<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Used to print the trace values as debugging utility.
      void PrintTrace(DRT::Element* ele);


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
        LocalSolver(const DRT::ELEMENTS::ElemagDiff* ele,
            CORE::FE::ShapeValues<distype>& shapeValues,
            Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>>& shapeValuesFace,
            INPAR::ELEMAG::DynamicType& dyna, CORE::FE::ShapeValues<distype>& postproc_shapeValues);

        /// Compute the residual
        void ComputeResidual(Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& eleVec, double dt, DRT::ELEMENTS::ElemagDiff& ele);

        /// Computes the source term in the element.
        void ComputeSource(Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& interiorSourcen,
            CORE::LINALG::SerialDenseVector& interiorSourcenp);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeAbsorbingBC(DRT::Discretization& discretization, DRT::ELEMENTS::ElemagDiff* ele,
            Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat, int face,
            CORE::LINALG::SerialDenseMatrix& elemat, int indexstart,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeBoundaryIntegral(
            DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params, int face);

        /// Calls local solver to compute matrices: internal and face
        void ComputeMatrices(DRT::Discretization& discretization,
            const Teuchos::RCP<MAT::Material>& mat, DRT::ELEMENTS::ElemagDiff& ele, double dt,
            INPAR::ELEMAG::DynamicType dyna, const double tau);

        /// Set up interior matrices
        void ComputeInteriorMatrices(double dt, double sigma, double mu, double epsilon);

        /// Set up face matrices
        void ComputeFaceMatrices(const int face, double dt, int indexstart, int newindex,
            double sigma, double mu, const double tau);

        /// Condense the local matrix into the element matrix for the trace and similarly for the
        /// residuals.
        void CondenseLocalPart(CORE::LINALG::SerialDenseMatrix& elemat);

        /*!
        \brief Make the matrix symmetric in case of dirichlet boundary conditions in the element

        \note So far the function only acts on the matrices and not on the RHS! This is ok for PEC
        boundary conditions because the DBC is anyway zero.
        */
        void Symmetrify(DRT::ELEMENTS::ElemagDiff& ele, CORE::LINALG::SerialDenseMatrix& elemat,
            bool dodirich = true);

        /// Projection of function field.
        /// The function is used to project the field in the initialization phase.
        int ProjectField(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2);

        /// Projection of scatra computed gradients into elements
        int ProjectElectricFieldFromScatra(DRT::ELEMENTS::ElemagDiff* ele,
            Teuchos::ParameterList& params, const Teuchos::RCP<MAT::Material>& mat,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Compute the error with respect to an analytical field.
        void ComputeError(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Postprocess the solution to obtain a better approximation
        void PostProcessing(DRT::ELEMENTS::ElemagDiff& ele);

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
        int ProjectFieldTest(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2);

        /// Projection of a given field on the trace for testing purposes.
        int ProjectFieldTestTrace(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Projection of Dirichlet function field.
        int ProjectDirichField(DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Function evaluation routine
        void EvaluateAll(const int start_func, const double t,
            const CORE::LINALG::Matrix<nsd_, 1>& xyz, CORE::LINALG::SerialDenseVector& v) const;

        /// Function gradient evaluation routine
        void ComputeFunctionGradient(int start_func, double t,
            const CORE::LINALG::Matrix<nsd_, 1>& xyz, CORE::LINALG::SerialDenseMatrix& v) const;

        /// Function time derivative evaluation routine
        void ComputeFunctionTimeDerivative(int start_func, double t, double dt,
            const CORE::LINALG::Matrix<nsd_, 1>& xyz, CORE::LINALG::SerialDenseVector& v) const;

        // convention: we sort the entries in the matrices the following way:
        // first comes the electric field then the magnetic one.

        /// Number of Degrees Of Freedom
        const unsigned int ndofs_;

        /// evaluated shape values
        CORE::FE::ShapeValues<distype>& shapes_;

        /// evaluated shape values
        Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;

        /// evaluated higher degree shape values
        CORE::FE::ShapeValues<distype>& postproc_shapes_;

        // System matrices
        /// u matrix
        CORE::LINALG::SerialDenseMatrix Amat;
        /// Inverse of u matrix
        CORE::LINALG::SerialDenseMatrix invAmat;
        /// u electric
        CORE::LINALG::SerialDenseMatrix Bmat;
        /// u trace
        CORE::LINALG::SerialDenseMatrix Cmat;
        CORE::LINALG::SerialDenseMatrix Fmat;
        /// electric second time
        CORE::LINALG::SerialDenseMatrix Dmat;
        /// electric evolution
        CORE::LINALG::SerialDenseMatrix Emat;
        /// electric electric
        CORE::LINALG::SerialDenseMatrix Gmat;
        /// electric trace
        CORE::LINALG::SerialDenseMatrix Hmat;
        /// trace trace
        CORE::LINALG::SerialDenseMatrix Lmat;

        // auxiliary stuff
        /// final mass matrix used for the projection
        CORE::LINALG::SerialDenseMatrix massMat;
        /// part of the mass matrix (only contains the shape functions)
        CORE::LINALG::SerialDenseMatrix massPart;
        /// other part of the mass matrix (with quadrature weights)
        CORE::LINALG::SerialDenseMatrix massPartW;

        /// Chosen dynamics/time integrator
        INPAR::ELEMAG::DynamicType& dyna_;
      };  // stuct LocalSolver

      /// Updates interior variables and calculates residual.
      void UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList& params,
          DRT::ELEMENTS::ElemagDiff& ele, const Teuchos::RCP<MAT::Material>& mat,
          CORE::LINALG::SerialDenseVector& elevec, double dt, bool errormaps, bool updateonly);

      /// Reads from global vectors.
      void ReadGlobalVectors(
          DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm);

      /// Writes internal fields from elements to global vectors.
      void FillRestartVectors(DRT::Element* ele, DRT::Discretization& discretization);

      /// Reads internal field from global vectors to element vectors.
      void ElementInitFromRestart(DRT::Element* ele, DRT::Discretization& discretization);

      /// Calculate error maps with local postprocessing.
      double EstimateError(DRT::ELEMENTS::ElemagDiff& ele, CORE::LINALG::SerialDenseVector& p);

      /// Local data object for element.
      Teuchos::RCP<CORE::FE::ShapeValues<distype>> shapes_;

      /// Higher degree local data object for element.
      Teuchos::RCP<CORE::FE::ShapeValues<distype>> postproc_shapes_;

      /// Local data object for face element.
      Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;

      /// Local solver object.
      Teuchos::RCP<LocalSolver> localSolver_;

      std::vector<double> localtrace_;  /// extracted values from trace solution vector

      /// Local values from interior solution vector (electric and magnetic fields) at n
      CORE::LINALG::SerialDenseVector interiorElectricnp_;
      CORE::LINALG::SerialDenseVector interiorMagneticnp_;
      CORE::LINALG::SerialDenseVector interiorauxiliaryPML_;
      /// Local values from interior solution vector (electric and magnetic fields) at n-1
      CORE::LINALG::SerialDenseVector interiorElectricnm_;
      CORE::LINALG::SerialDenseVector interiorMagneticnm_;


      /// Chosen dynamics/time integrator
      INPAR::ELEMAG::DynamicType dyna_;

      /// Use complete polynomial or tensor product
      bool usescompletepoly_;
    };  // class ElemagDiffEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_DIFF_ELE_CALC_H
