/*----------------------------------------------------------------------------*/
/*! \file
\brief Routines for calculation of HDG weakly compressible fluid element

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_HDG_WEAK_COMP_H
#define FLUID_ELE_CALC_HDG_WEAK_COMP_H


#include "fluid_ele_interface.H"
#include "fluid_ele_hdg_weak_comp.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fem_general/drt_utils_shapevalues_hdg.H"



namespace DRT
{
  namespace ELEMENTS
  {
    /// Weakly compressible fluid HDG element implementation
    /*!

      \author laspina
      \date 08/19
    */
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcHDGWeakComp : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! mixed variable dimension according to Voigt notation
      static const unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;

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

      /// Update local solution
      virtual int UpdateLocalSolution(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& interiorinc);

      /// projection of function field
      virtual int ProjectField(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2);

      /*! \brief Interpolates an HDG solution to the element nodes for output
       */
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1);

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
              side_coupling,                ///< side coupling matrices
          Epetra_SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
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

      virtual ~FluidEleCalcHDGWeakComp() {}

      /// Singleton access method
      static FluidEleCalcHDGWeakComp<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcHDGWeakComp();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        static const unsigned int nsd_ = FluidEleCalcHDGWeakComp<distype>::nsd_;
        static const unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;
        static const unsigned int nfaces_ = FluidEleCalcHDGWeakComp<distype>::nfaces_;

        /// local solver
        LocalSolver(const DRT::ELEMENTS::Fluid* ele,
            const DRT::UTILS::ShapeValues<distype>& shapeValues,
            DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly);

        /// initialize all
        void InitializeAll();

        /// compute material matrix
        void ComputeMaterialMatrix(const Teuchos::RCP<MAT::Material>& mat,
            const LINALG::Matrix<nsd_, 1>& xyz, Epetra_SerialDenseMatrix& D_L,
            Epetra_SerialDenseMatrix& D_w);

        /// compute interior residual
        void ComputeInteriorResidual(const Teuchos::RCP<MAT::Material>& mat,
            const std::vector<double>& valnp, const std::vector<double>& valn,
            const std::vector<double>& alevel);

        /// compute face residual
        void ComputeFaceResidual(const int f, const Teuchos::RCP<MAT::Material>& mat,
            const std::vector<double>& val, const std::vector<double>& traceval,
            const std::vector<double>& alevel);

        /// compute interior matrices
        void ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material>& mat);

        /// compute face matrices
        void ComputeFaceMatrices(const int f, const Teuchos::RCP<MAT::Material>& mat);

        /// compute local residual
        void ComputeLocalResidual();

        /// compute global residual
        void ComputeGlobalResidual(DRT::ELEMENTS::Fluid& ele);

        /// compute local-local matrix
        void ComputeLocalLocalMatrix();

        /// compute local-global matrix
        void ComputeLocalGlobalMatrix(DRT::ELEMENTS::Fluid& ele);

        /// compute global-local matrix
        void ComputeGlobalLocalMatrix(DRT::ELEMENTS::Fluid& ele);

        /// compute global-global matrix
        void ComputeGlobalGlobalMatrix(DRT::ELEMENTS::Fluid& ele);

        /// invert local-local matrix
        void InvertLocalLocalMatrix();

        /// condense local residual
        void CondenseLocalResidual(Epetra_SerialDenseVector& eleVec);

        /// condense local matrix
        void CondenseLocalMatrix(Epetra_SerialDenseMatrix& eleMat);

        // print matrices and residuals
        void PrintMatricesAndResiduals(DRT::ELEMENTS::Fluid& ele, Epetra_SerialDenseVector& eleVec,
            Epetra_SerialDenseMatrix& eleMat);

        // number of degrees of freedom
        const unsigned int ndofs_;

        // total number of degrees of freedom in the faces
        unsigned int ndofsfaces_;

        // flag for convective flow
        bool convective;

        // flag for unsteady flow
        bool unsteady;

        // flag for ALE approach
        bool ale;

        // convention: we sort the entries in the matrices the following way:
        // first come the mixed variable, then the density, and finally the momentum
        // we also build the matrix in a block-fashion, keeping the dofs for individual components
        // closest to each other. I.e. the blocks are in 2D for L_0, L_1, L_2, r, w_0, w_1
        // and similarly for 3D

        const DRT::UTILS::ShapeValues<distype>& shapes_;    /// evaluated shape values
        DRT::UTILS::ShapeValuesFace<distype>& shapesface_;  /// evaluated shape values

        // Stabilization parameters
        double tau_r;  /// stabilization of density
        double tau_w;  /// stabilization of momentum

        // Auxiliary matrices
        Epetra_SerialDenseMatrix massPart;   /// temporary matrix for mass matrix
        Epetra_SerialDenseMatrix massPartW;  /// temporary matrix for mass matrix with weights
        Epetra_SerialDenseMatrix massMat;    /// local mass matrix

        // Unknown variables
        Epetra_SerialDenseMatrix Leg;  /// mixed variable evaluated on interior quadrature points
        Epetra_SerialDenseVector reg;  /// density evaluated on interior quadrature points
        Epetra_SerialDenseMatrix weg;  /// momentum evaluated on interior quadrature points
        Epetra_SerialDenseVector rhatefg;  /// trace of density evaluated on face quadrature points
        Epetra_SerialDenseMatrix whatefg;  /// trace of momentum evaluated on face quadrature points

        // ALE variables
        Epetra_SerialDenseMatrix aeg;   /// ALE velocity evaluated on interior quadrature points
        Epetra_SerialDenseMatrix aefg;  /// ALE velocity evaluated on face quadrature points
        Epetra_SerialDenseMatrix
            dadxyzeg;  /// derivatives of ALE velocity evaluated on interior quadrature points

        // Matrices
        Epetra_SerialDenseMatrix ALL;  /// matrix mixed variable - mixed variable
        Epetra_SerialDenseMatrix ALr;  /// matrix mixed variable - density
        Epetra_SerialDenseMatrix ALw;  /// matrix mixed variable - momentum
        Epetra_SerialDenseMatrix ALR;  /// matrix mixed variable - trace of density
        Epetra_SerialDenseMatrix ALW;  /// matrix mixed variable - trace of momentum
        Epetra_SerialDenseMatrix Arr;  /// matrix density - density
        Epetra_SerialDenseMatrix Arw;  /// matrix density - momentum
        Epetra_SerialDenseMatrix ArR;  /// matrix density - trace of density
        Epetra_SerialDenseMatrix ArW;  /// matrix density - trace of momentum
        Epetra_SerialDenseMatrix AwL;  /// matrix momentum - mixed variable
        Epetra_SerialDenseMatrix Awr;  /// matrix momentum - density
        Epetra_SerialDenseMatrix Aww;  /// matrix momentum - momentum
        Epetra_SerialDenseMatrix AwR;  /// matrix momentum - trace of density
        Epetra_SerialDenseMatrix AwW;  /// matrix momentum - trace of momentum
        Epetra_SerialDenseMatrix ARr;  /// matrix trace of density - densitym
        Epetra_SerialDenseMatrix ARR;  /// matrix trace of density - trace of density
        Epetra_SerialDenseMatrix AWL;  /// matrix trace of momentum - mixed variable
        Epetra_SerialDenseMatrix AWw;  /// matrix trace of momentum - momentum
        Epetra_SerialDenseMatrix AWR;  /// matrix trace of momentum - trace of density
        Epetra_SerialDenseMatrix AWW;  /// matrix trace of momentum - trace of momentum

        // Residuals
        Epetra_SerialDenseVector RL;  /// residual vector for mixed variable
        Epetra_SerialDenseVector Rr;  /// residual vector for density
        Epetra_SerialDenseVector Rw;  /// residual vector for momentu
        Epetra_SerialDenseVector RR;  /// residual vector for trace of density
        Epetra_SerialDenseVector RW;  /// residual vector for trace of momentum

        // Local/Global matrices/vectors
        Epetra_SerialDenseMatrix Klocallocal;           /// local-local matrix
        Epetra_SerialDenseMatrix Klocalglobal;          /// local-global matrix
        Epetra_SerialDenseMatrix Kgloballocal;          /// global-local matrix
        Epetra_SerialDenseMatrix Kglobalglobal;         /// global-global matrix
        Epetra_SerialDenseVector Rlocal;                /// local residual vector
        Epetra_SerialDenseVector Rglobal;               /// global residual vector
        Epetra_SerialDenseMatrix KlocallocalInv;        /// inverse local-local matrix
        Epetra_SerialDenseSolver KlocallocalInvSolver;  /// solver for inverse local-local matrix

        // Voigt related quantities
        int VoigtP[msd_ - nsd_][2];  /// pair of indices in Voigt notation

        std::vector<int> pivots;  /// pivots for factorization of matrices

        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;  //! pointer to parameter list
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameterTimInt>
            fldparatimint_;  //! pointer to time parameter list
      };

      /// reads from global vectors
      void ReadGlobalVectors(
          const DRT::Element& ele, DRT::Discretization& discretization, const std::vector<int>& lm);

      /// reads ale vectors
      void ReadAleVectors(const DRT::Element& ele, DRT::Discretization& discretization);

      /// evaluate mixed variable, density and momentum
      void EvaluateAll(const int funcnum, const LINALG::Matrix<nsd_, 1>& xyz, const double t,
          LINALG::Matrix<msd_, 1>& L, double& r, LINALG::Matrix<nsd_, 1>& w) const;

      /// evaluate density and momentum
      void EvaluateDensityMomentum(const int funcnum, const LINALG::Matrix<nsd_, 1>& xyz,
          const double t, double& r, LINALG::Matrix<nsd_, 1>& u) const;

      /// local data object
      Teuchos::RCP<DRT::UTILS::ShapeValues<distype>> shapes_;
      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

      /// local solver object
      Teuchos::RCP<LocalSolver> localSolver_;

      std::vector<double> traceVal_;  /// extracted values from trace solution vector at n+alpha_f
      std::vector<double> interiorVal_;  /// extracted local values at n+alpha_f
      std::vector<double> interiorAcc_;  /// extracted local accelerations at n+alpha_m
      std::vector<double> aleDis_;       /// extracted ale mesh displacement
      std::vector<double> aleVel_;       /// extracted ale mesh velocity

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
