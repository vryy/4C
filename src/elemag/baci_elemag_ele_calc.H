/*--------------------------------------------------------------------------*/
/*! \file

\brief All functionality for electromagnetic element evaluations

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_ELEMAG_ELE_CALC_H
#define BACI_ELEMAG_ELE_CALC_H

#include "baci_config.H"

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_shapevalues_hdg.H"
#include "baci_elemag_ele.H"
#include "baci_elemag_ele_interface.H"
#include "baci_inpar_elemag.H"
#include "baci_utils_singleton_owner.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    /// Electromagnetic element implementation

    template <CORE::FE::CellType distype>
    class ElemagEleCalc : public ElemagEleInterface
    {
     public:
      /// nen_: number of element nodes (T. Hughes: The Finite Element Method).
      static const unsigned int nen_ = CORE::FE::num_nodes<distype>;

      /// Number of space dimensions.
      static const unsigned int nsd_ = CORE::FE::dim<distype>;

      /// Number of faces on element.
      static const unsigned int nfaces_ = CORE::FE::num_faces<distype>;

      int IntegrateShapeFunction(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1) override
      {
        dserror("Not implemented");
        return 1;
      };

      /// Zero initialization of elements.
      virtual void ElementInit(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params);

      /// Interpolates an HDG solution to the element nodes for output.
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::Elemag* ele,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& elevec1);

      /// Initialize the shape functions and solver to the given element (degree is runtime
      /// parameter).
      void InitializeShapes(const DRT::ELEMENTS::Elemag* ele);

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
      static ElemagEleCalc<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Used to print the trace values as debugging utility.
      void PrintTrace(DRT::Element* ele);


     private:
      /// Private Constructor since we are a Singleton.
      ElemagEleCalc();

      /// Local solver that inverts local problem on an element and can solve with various vectors.
      struct LocalSolver
      {
        /// Number of Spatial Dimensions
        static const unsigned int nsd_ = ElemagEleCalc<distype>::nsd_;
        /// Number of FACES
        static const unsigned int nfaces_ = ElemagEleCalc<distype>::nfaces_;

        /// Init function for the struct members
        LocalSolver(const DRT::ELEMENTS::Elemag* ele, CORE::FE::ShapeValues<distype>& shapeValues,
            Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>>& shapeValuesFace,
            INPAR::ELEMAG::DynamicType& dyna);

        /// Compute the residual
        void ComputeResidual(Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& eleVec, DRT::ELEMENTS::Elemag& ele);

        /// Computes the source term in the element.
        void ComputeSource(Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& interiorSourcen,
            CORE::LINALG::SerialDenseVector& interiorSourcenp);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeAbsorbingBC(DRT::Discretization& discretization, DRT::ELEMENTS::Elemag* ele,
            Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat, int face,
            CORE::LINALG::SerialDenseMatrix& elemat, int indexstart,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeBoundaryIntegral(
            DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params, int face);

        /// Calls local solver to compute matrices: internal and face
        void ComputeMatrices(DRT::Discretization& discretization,
            const Teuchos::RCP<MAT::Material>& mat, DRT::ELEMENTS::Elemag& ele, double dt,
            INPAR::ELEMAG::DynamicType dyna, const double tau);

        /// Set up interior matrices
        void ComputeInteriorMatrices(double dt, double sigma, double mu, double epsilon);

        /// Set up face matrices
        void ComputeFaceMatrices(const int face, double dt, int indexstart, int newindex,
            double sigma, double mu, const double tau);

        /// Condense the local matrx into the element matrix for the trace and similarly for the
        /// residuals.
        void CondenseLocalPart(CORE::LINALG::SerialDenseMatrix& elemat);

        /// Projection of function field.
        /// The function is used to project the field in the initialization phase.
        int ProjectField(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2);

        /// Compute the error with respect to an analytical field.
        void ComputeError(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Projection of a given field on the interior variables for testing purposes.
        int ProjectFieldTest(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2);

        /// Projection of a given field on the trace for testing purposes.
        int ProjectFieldTestTrace(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Projection of Dirichlet function field.
        int ProjectDirichField(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            CORE::LINALG::SerialDenseVector& elevec1);

        /// Function evaluation routine
        void EvaluateAll(const int start_func, const double t,
            const CORE::LINALG::Matrix<nsd_, 1>& xyz, CORE::LINALG::SerialDenseVector& v) const;

        // convention: we sort the entries in the matrices the following way:
        // first comes the electric field then the magnetic one.

        /// Number of Degrees Of Freedom
        const unsigned int ndofs_;

        /// evaluated shape values
        CORE::FE::ShapeValues<distype>& shapes_;
        /// evaluated shape values
        Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;

        // System matrices
        /// magnetic evolution matrix
        CORE::LINALG::SerialDenseMatrix Amat;
        /// Inverse of magnetic evolution matrix
        CORE::LINALG::SerialDenseMatrix invAmat;
        /// magnetic electric
        CORE::LINALG::SerialDenseMatrix Cmat;
        /// magnetic trace
        CORE::LINALG::SerialDenseMatrix Dmat;
        /// electric evolution
        CORE::LINALG::SerialDenseMatrix Emat;
        /// electric magnetic
        CORE::LINALG::SerialDenseMatrix Fmat;
        /// electric electric
        CORE::LINALG::SerialDenseMatrix Gmat;
        /// electric trace
        CORE::LINALG::SerialDenseMatrix Hmat;
        /// trace magnetic
        CORE::LINALG::SerialDenseMatrix Imat;
        /// trace electric
        CORE::LINALG::SerialDenseMatrix Jmat;
        /// trace trace
        CORE::LINALG::SerialDenseMatrix Lmat;

        // auxiliary stuff
        CORE::LINALG::SerialDenseMatrix massMat;  // final mass matrix used for the projection
        CORE::LINALG::SerialDenseMatrix
            massPart;  // part of the mass matrix (only contains the shape functions)
        CORE::LINALG::SerialDenseMatrix
            massPartW;  // other part of the mass matrix (with quadrature weights)

        /// Chosen dynamics/time integrator
        INPAR::ELEMAG::DynamicType& dyna_;
      };

      /// Updates interior variables and calculates residual.
      void UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList& params,
          DRT::ELEMENTS::Elemag& ele, const Teuchos::RCP<MAT::Material>& mat,
          CORE::LINALG::SerialDenseVector& elevec, double dt, bool errormaps, bool updateonly);

      /// Reads from global vectors.
      void ReadGlobalVectors(
          DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm);

      /// Writes internal fields from elements to global vectors.
      void FillRestartVectors(DRT::Element* ele, DRT::Discretization& discretization);

      /// Reads internal field from global vectors to element vectors.
      void ElementInitFromRestart(DRT::Element* ele, DRT::Discretization& discretization);

      /// Calculate error maps with local postprocessing.
      double EstimateError(DRT::ELEMENTS::Elemag& ele, CORE::LINALG::SerialDenseVector& p);

      /// Local data object for element.
      Teuchos::RCP<CORE::FE::ShapeValues<distype>> shapes_;
      /// Local data object for face element.
      Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_;

      /// Local solver object.
      Teuchos::RCP<LocalSolver> localSolver_;

      std::vector<double> localtrace_;  /// extracted values from trace solution vector

      /// Local values from interior solution vector (electric and magnetic fields) at n
      CORE::LINALG::SerialDenseVector interiorElectricnp_;
      CORE::LINALG::SerialDenseVector interiorMagneticnp_;
      /// Local values from interior solution vector (electric and magnetic fields) at n-1
      CORE::LINALG::SerialDenseVector interiorElectricnm_;
      CORE::LINALG::SerialDenseVector interiorMagneticnm_;
      CORE::LINALG::SerialDenseVector interiorauxiliaryPML_;

      /// Chosen dynamics/time integrator
      INPAR::ELEMAG::DynamicType dyna_;

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_ELE_CALC_H
