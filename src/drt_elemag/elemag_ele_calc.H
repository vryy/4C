/*--------------------------------------------------------------------------*/
/*! \file

\brief All functionality for electromagnetic element evaluations

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
*/
/*--------------------------------------------------------------------------*/

#ifndef ELEMAG_ELE_CALC_H
#define ELEMAG_ELE_CALC_H

#include "elemag_ele_interface.H"
#include "elemag_ele.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_inpar/inpar_elemag.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Electromagnetic element implementation

    template <DRT::Element::DiscretizationType distype>
    class ElemagEleCalc : public ElemagEleInterface
    {
     public:
      /// nen_: number of element nodes (T. Hughes: The Finite Element Method).
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      /// Number of space dimensions.
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      /// Number of faces on element.
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Elemag* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1)
      {
        dserror("Not implemented");
        return 1;
      };

      /// Zero initialization of elements.
      virtual void ElementInit(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params);

      /// Interpolates an HDG solution to the element nodes for output.
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::Elemag* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1);

      /// Initialize the shape functions and solver to the given element (degree is runtime
      /// parameter).
      void InitializeShapes(const DRT::ELEMENTS::Elemag* ele);

      /// Evaluate the element.
      /// Generic virtual interface function.Called via base pointer.
      virtual int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

      /// Evaluate the element at specified gauss points.
      virtual int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints, bool offdiag = false);


      virtual ~ElemagEleCalc() {}

      /// Singleton access method.
      static ElemagEleCalc<distype>* Instance(bool create = true);

      /// Called upon destruction.
      virtual void Done();

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
        LocalSolver(const DRT::ELEMENTS::Elemag* ele, DRT::UTILS::ShapeValues<distype>& shapeValues,
            Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>>& shapeValuesFace,
            INPAR::ELEMAG::DynamicType& dyna);

        /// Compute the residual
        void ComputeResidual(Teuchos::ParameterList& params, Epetra_SerialDenseVector& eleVec,
            DRT::ELEMENTS::Elemag& ele);

        /// Computes the source term in the element.
        void ComputeSource(Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& interiorSourcen, Epetra_SerialDenseVector& interiorSourcenp);

        /// Add terms corresponding to the absorbing boundary condition.
        void ComputeAbsorbingBC(DRT::Discretization& discretization, DRT::ELEMENTS::Elemag* ele,
            Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat, int face,
            Epetra_SerialDenseMatrix& elemat, int indexstart, Epetra_SerialDenseVector& elevec1);

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
        void CondenseLocalPart(Epetra_SerialDenseMatrix& elemat);

        /// Projection of function field.
        /// The function is used to project the field in the initialization phase.
        int ProjectField(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2);

        /// Compute the error with respect to an analytical field.
        void ComputeError(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& elevec1);

        /// Projection of a given field on the interior variables for testing purposes.
        int ProjectFieldTest(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2);

        /// Projection of a given field on the trace for testing purposes.
        int ProjectFieldTestTrace(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& elevec1);

        /// Projection of Dirichlet function field.
        int ProjectDirichField(DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
            Epetra_SerialDenseVector& elevec1);

        /// Function evaluation routine
        void EvaluateAll(const int start_func, const double t, const LINALG::Matrix<nsd_, 1>& xyz,
            Epetra_SerialDenseVector& v) const;

        // convention: we sort the entries in the matrices the following way:
        // first comes the electric field then the magnetic one.

        /// Number of Degrees Of Freedom
        const unsigned int ndofs_;

        /// evaluated shape values
        DRT::UTILS::ShapeValues<distype>& shapes_;
        /// evaluated shape values
        Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

        // System matrices
        /// magnetic evolution matrix
        Epetra_SerialDenseMatrix Amat;
        /// Inverse of magnetic evolution matrix
        Epetra_SerialDenseMatrix invAmat;
        /// magnetic electric
        Epetra_SerialDenseMatrix Cmat;
        /// magnetic trace
        Epetra_SerialDenseMatrix Dmat;
        /// electric evolution
        Epetra_SerialDenseMatrix Emat;
        /// electric magnetic
        Epetra_SerialDenseMatrix Fmat;
        /// electric electric
        Epetra_SerialDenseMatrix Gmat;
        /// electric trace
        Epetra_SerialDenseMatrix Hmat;
        /// trace magnetic
        Epetra_SerialDenseMatrix Imat;
        /// trace electric
        Epetra_SerialDenseMatrix Jmat;
        /// trace trace
        Epetra_SerialDenseMatrix Lmat;

        // auxiliary stuff
        Epetra_SerialDenseMatrix massMat;  // final mass matrix used for the projection
        Epetra_SerialDenseMatrix
            massPart;  // part of the mass matrix (only contains the shape functions)
        Epetra_SerialDenseMatrix
            massPartW;  // other part of the mass matrix (with quadrature weights)

        /// Chosen dynamics/time integrator
        INPAR::ELEMAG::DynamicType& dyna_;
      };

      /// Updates interior variables and calculates residual.
      void UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList& params,
          DRT::ELEMENTS::Elemag& ele, const Teuchos::RCP<MAT::Material>& mat,
          Epetra_SerialDenseVector& elevec, double dt, bool errormaps, bool updateonly);

      /// Reads from global vectors.
      void ReadGlobalVectors(
          DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm);

      /// Writes internal fields from elements to global vectors.
      void FillRestartVectors(DRT::Element* ele, DRT::Discretization& discretization);

      /// Reads internal field from global vectors to element vectors.
      void ElementInitFromRestart(DRT::Element* ele, DRT::Discretization& discretization);

      /// Calculate error maps with local postprocessing.
      double EstimateError(DRT::ELEMENTS::Elemag& ele, Epetra_SerialDenseVector& p);

      /// Local data object for element.
      Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_;
      /// Local data object for face element.
      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

      /// Local solver object.
      Teuchos::RCP<LocalSolver> localSolver_;

      std::vector<double> localtrace_;  /// extracted values from trace solution vector

      /// Local values from interior solution vector (electric and magnetic fields) at n
      Epetra_SerialDenseVector interiorElectricnp_;
      Epetra_SerialDenseVector interiorMagneticnp_;
      /// Local values from interior solution vector (electric and magnetic fields) at n-1
      Epetra_SerialDenseVector interiorElectricnm_;
      Epetra_SerialDenseVector interiorMagneticnm_;
      Epetra_SerialDenseVector interiorauxiliaryPML_;

      /// Chosen dynamics/time integrator
      INPAR::ELEMAG::DynamicType dyna_;

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* ELEMAG_ELE_CALC_H */
