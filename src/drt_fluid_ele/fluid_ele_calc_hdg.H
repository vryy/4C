/*----------------------------------------------------------------------*/
/*!

\brief main file containing routines for calculation of HDG fluid element

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_HDG_H
#define FLUID_ELE_CALC_HDG_H


#include "fluid_ele_interface.H"
#include "fluid_ele_hdg.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"



namespace DRT
{
  namespace ELEMENTS
  {
    /// Fluid HDG element implementation
    /*!

      \author kronbichler
      \date 05/13
    */
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcHDG : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      ///! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;


      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, const DRT::UTILS::GaussIntegration& intpoints)
      {
        dserror("Not implemented!");
        return 1;
      }

      virtual int IntegrateShapeFunctionXFEM(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells)
      {
        dserror("Not implemented!");
        return 1;
      };


      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */
      virtual int EvaluateService(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      /*!
        \brief calculate dissipation of various terms (evaluation of turbulence models)
      */
      virtual int CalcDissipation(Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat)
      {
        dserror("Not implemented!");
        return 1;
      }

      /// Evaluate element ERROR
      /*!
          general function to compute the error (analytical solution) for particular problem type
       */
      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec);

      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec,
          const DRT::UTILS::GaussIntegration&)
      {
        return ComputeError(ele, params, mat, discretization, lm, elevec);
      }

      /// projection of function field
      virtual int ProjectField(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2);

      /*! \brief Interpolates an HDG solution to the element nodes for output
       */
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1);

      /*! \brief Interpolates an HDG solution for homogeneous isotropic turbulence postprocessing
       */
      virtual int InterpolateSolutionForHIT(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1);

      /*! \brief Project force from equidistant points on interior node dof vector
       */
      virtual int ProjectForceOnDofVecForHIT(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2);

      /*! \brief Project initial field for hit
       */
      virtual int ProjectInitialFieldForHIT(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);
      /*!
      \brief Initialize the shape functions and solver to the given element (degree is runtime
      parameter)
       */
      void InitializeShapes(const DRT::ELEMENTS::Fluid* ele);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

      /// Evaluate the element at specified gauss points
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints, bool offdiag = false);

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
      )
      {
        dserror("Not implemented!");
        return 1;
      }

      /// Evaluate the XFEM cut element
      virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells, bool offdiag = false)
      {
        dserror("Not implemented!");
        return 1;
      }


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
      )
      {
        dserror("Not implemented!");
        return;
      }

      virtual void ElementXfemInterfaceNIT(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,               ///< parameter list
          Teuchos::RCP<MAT::Material>& mat_master,      ///< material master side
          Teuchos::RCP<MAT::Material>& mat_slave,       ///< material slave side
          Epetra_SerialDenseMatrix& elemat1_epetra,     ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,     ///< element vector
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,               ///< side coupling matrices
          Epetra_SerialDenseMatrix& Cuiui  ///< ui-ui coupling matrix
      )
      {
        dserror("Not implemented!");
        return;
      }

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra,
          const DRT::UTILS::GaussIntegration& intpoints)
      {
        dserror("Not implemented!");
        return;
      }

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra)
      {
        dserror("Not implemented!");
        return;
      }

      /// Evaluate the pressure average inside the element from an analytical expression
      virtual int EvaluatePressureAverage(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseVector& elevec);

      virtual ~FluidEleCalcHDG() {}

      /// Singleton access method
      static FluidEleCalcHDG<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// Print local residuals
      void PrintLocalResiduals(DRT::ELEMENTS::Fluid* ele);

      /// Print local variables
      void PrintLocalVariables(DRT::ELEMENTS::Fluid* ele);

      /// Print local correction term
      void PrintLocalCorrection(
          DRT::ELEMENTS::Fluid* ele, std::vector<double>& interiorecorrectionterm);

      /// Print local body force
      void PrintLocalBodyForce(DRT::ELEMENTS::Fluid* ele, std::vector<double>& interiorebodyforce);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcHDG();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        static const unsigned int nsd_ = FluidEleCalcHDG<distype>::nsd_;
        static const unsigned int nfaces_ = FluidEleCalcHDG<distype>::nfaces_;

        LocalSolver(const DRT::ELEMENTS::Fluid* ele,
            const DRT::UTILS::ShapeValues<distype>& shapeValues,
            DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly);

        void ComputeInteriorResidual(const Teuchos::RCP<MAT::Material>& mat,
            const std::vector<double>& valnp, const std::vector<double>& valn,
            const double avgPressure, const LINALG::Matrix<nsd_, nen_>& ebodyforce,
            const std::vector<double>& intebodyforce, Epetra_SerialDenseVector& eleVec,
            const std::vector<double>& interiorecorrectionterm,
            const std::vector<double>& interiorebodyforce);

        void ComputeFaceResidual(const int face, const Teuchos::RCP<MAT::Material>& mat,
            const std::vector<double>& val, const std::vector<double>& traceval,
            Epetra_SerialDenseVector& eleVec);

        void ComputeInteriorMatrices(
            const Teuchos::RCP<MAT::Material>& mat, const bool evaluateOnlyNonlinear);

        void ComputeFaceMatrices(const int face, const Teuchos::RCP<MAT::Material>& mat,
            const bool evaluateOnlyNonlinear, Epetra_SerialDenseMatrix& elemat);

        // inverts the velocity gradient matrix and puts its contribution into the velocity matrix
        // (pre-factorization). Should only be done once per element even if multiple velocities are
        // used
        void EliminateVelocityGradient(Epetra_SerialDenseMatrix& elemat);

        // solves the local problem, including factorization of the matrix
        void SolveResidual();

        // condense the local matrix (involving cell velocity gradients, velocities and pressure)
        // into the element matrix for the trace and similarly for the residuals
        void CondenseLocalPart(Epetra_SerialDenseMatrix& elemat, Epetra_SerialDenseVector& elevec);

        // compute the correction term on the rhs for the weakly compressible benchmark
        void ComputeCorrectionTerm(
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

        const DRT::UTILS::ShapeValues<distype>& shapes_;    /// evaluated shape values
        DRT::UTILS::ShapeValuesFace<distype>& shapesface_;  /// evaluated shape values

        double stabilization[nfaces_];  /// stabilization parameters

        Epetra_SerialDenseMatrix
            uuMat;  /// terms for block with velocity and pressure (constant ones)
        Epetra_SerialDenseMatrix uuMatFinal;  /// terms for block with velocity and pressure
                                              /// (including convection and stabilization)
        Epetra_SerialDenseMatrix
            ugMat;  /// coupling between velocity and velocity gradient (not fully stored)
        Epetra_SerialDenseMatrix
            guMat;  /// evaluated divergence of velocity gradient and velocity (not fully stored)

        Epetra_SerialDenseMatrix gfMat;  /// evaluated coupling between velocity gradient and trace
        Epetra_SerialDenseMatrix fgMat;  /// evaluated coupling between trace and velocity gradient
        Epetra_SerialDenseMatrix ufMat;  /// evaluated coupling between velocity and trace
        Epetra_SerialDenseMatrix fuMat;  /// evaluated coupling between trace and velocity

        Epetra_SerialDenseMatrix
            massPart;  /// temporary matrix for mass matrix on all quadrature points
        Epetra_SerialDenseMatrix massPartW;  /// temporary matrix for mass matrix weighted by
                                             /// integration factor on all quadrature points
        Epetra_SerialDenseMatrix
            gradPart;  /// temporary matrix for gradient matrix on all quadrature points
        Epetra_SerialDenseMatrix uPart;  /// temporary matrix for convection

        Epetra_SerialDenseMatrix massMat;     /// local mass matrix (will be inverted during init)
        Epetra_SerialDenseMatrix uuconv;      /// convection matrix
        Epetra_SerialDenseMatrix tmpMat;      /// matrix holding temporary results
        Epetra_SerialDenseMatrix tmpMatGrad;  /// matrix holding temporary results

        Epetra_SerialDenseMatrix trMat;     /// temporary matrix for trace assembly
        Epetra_SerialDenseMatrix trMatAvg;  /// temporary matrix for trace assembly

        Epetra_SerialDenseMatrix velnp;  /// velocities evaluated on all quadrature points
        Epetra_SerialDenseMatrix
            fvelnp;  /// trace velocities evaluated on all face quadrature points

        Epetra_SerialDenseMatrix uucomp;  /// compressibility matrix
        Epetra_SerialDenseVector presnp;  /// pressure evaluated on all quadrature points
        Epetra_SerialDenseMatrix
            gradpresnp;                     /// pressure gradient evaluated on all quadrature points
        Epetra_SerialDenseVector ifpresnp;  /// pressure evaluated on all face quadrature points

        Epetra_SerialDenseVector gRes;   /// residual vector on velocity gradients
        Epetra_SerialDenseVector upRes;  /// residual vector on velocity and pressure
        Epetra_SerialDenseVector gUpd;   /// update vector for velocity gradients
        Epetra_SerialDenseVector upUpd;  /// update vector for velocity and pressure

        std::vector<int> pivots;  /// pivots for factorization of matrices

        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;  //! pointer to parameter list
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameterTimInt>
            fldparatimint_;  //! pointer to time parameter list
      };

      /// reads from global vectors
      void ReadGlobalVectors(const DRT::Element& ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, const bool updateLocally);

      // writes the updated solution vector to the secondary vector stored in the discretization
      void UpdateSecondarySolution(const DRT::Element& ele, DRT::Discretization& discretization,
          const Epetra_SerialDenseVector& updateG, const Epetra_SerialDenseVector& updateUp);

      void EvaluateVelocity(const int start_func, const INPAR::FLUID::InitialField initfield,
          const LINALG::Matrix<nsd_, 1>& xyz, LINALG::Matrix<nsd_, 1>& u) const;

      void EvaluateAll(const int start_func, const INPAR::FLUID::InitialField initfield,
          const LINALG::Matrix<nsd_, 1>& xyz, LINALG::Matrix<nsd_, 1>& u,
          LINALG::Matrix<nsd_, nsd_>& grad, double& p) const;

      /// local data object
      Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_;
      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

      /// local solver object
      Teuchos::RCP<LocalSolver> localSolver_;

      LINALG::Matrix<nsd_, nen_> ebofoaf_;     /// body force (see fluid_ele_calc.cpp)
      LINALG::Matrix<nsd_, nen_> eprescpgaf_;  /// pressure gradient body force
      LINALG::Matrix<nen_, 1> escabofoaf_;     /// scalar body force for loma
      std::vector<double> interiorebofoaf_;    /// extracted body force at n+alpha_f

      std::vector<double>
          interiorecorrectionterm_;  /// local correction term for the weakly compressible benchmark
      std::vector<double>
          interiorebodyforce_;  /// local body force for the weakly compressible benchmark

      std::vector<double> traceVal_;  /// extracted values from trace solution vector at n+alpha_f
      std::vector<double> interiorVal_;  /// extracted local values (velocity gradients, velocities,
                                         /// pressure) at n+alpha_f
      std::vector<double> interiorAcc_;  /// extracted local accelerations at n+alpha_m

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
