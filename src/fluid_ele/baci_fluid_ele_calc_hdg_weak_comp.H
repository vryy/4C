/*----------------------------------------------------------------------------*/
/*! \file
\brief Routines for calculation of HDG weakly compressible fluid element

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_CALC_HDG_WEAK_COMP_H
#define BACI_FLUID_ELE_CALC_HDG_WEAK_COMP_H


#include "baci_discretization_fem_general_utils_shapevalues_hdg.H"
#include "baci_fluid_ele_hdg_weak_comp.H"
#include "baci_fluid_ele_interface.H"
#include "baci_inpar_fluid.H"
#include "baci_utils_singleton_owner.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Weakly compressible fluid HDG element implementation
    /*!

      \author laspina
      \date 08/19
    */
    template <CORE::FE::CellType distype>
    class FluidEleCalcHDGWeakComp : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions
      static const unsigned int nsd_ = CORE::FE::dim<distype>;

      //! mixed variable dimension according to Voigt notation
      static const unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;

      ///! number of faces on element
      static const unsigned int nfaces_ = CORE::DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;


      int IntegrateShapeFunction(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          const CORE::DRT::UTILS::GaussIntegration& intpoints) override
      {
        dserror("Not implemented!");
        return 1;
      }

      int IntegrateShapeFunctionXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          const std::vector<CORE::DRT::UTILS::GaussIntegration>& intpoints,
          const CORE::GEO::CUT::plain_volumecell_set& cells) override
      {
        dserror("Not implemented!");
        return 1;
      };


      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */
      int EvaluateService(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

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
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec);

      int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec,
          const CORE::DRT::UTILS::GaussIntegration&) override
      {
        return ComputeError(ele, params, mat, discretization, lm, elevec);
      }

      /// Update local solution
      virtual int UpdateLocalSolution(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& interiorinc);

      /// projection of function field
      virtual int ProjectField(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2);

      /*! \brief Interpolates an HDG solution to the element nodes for output
       */
      virtual int InterpolateSolutionToNodes(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& elevec1);

      /*!
      \brief Initialize the shape functions and solver to the given element (degree is runtime
      parameter)
       */
      void InitializeShapes(const DRT::ELEMENTS::Fluid* ele);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points
      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const CORE::DRT::UTILS::GaussIntegration& intpoints, bool offdiag = false) override;

      int ComputeErrorInterface(DRT::ELEMENTS::Fluid* ele,           ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<MAT::Material>& mat,                          ///< material
          CORE::LINALG::SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::DRT::UTILS::GaussIntegration>>&
              bintpoints,                                     ///< boundary integration points
          const CORE::GEO::CUT::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                      ///< parameter list
          ) override
      {
        dserror("Not implemented!");
        return 1;
      }

      /// Evaluate the XFEM cut element
      int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const std::vector<CORE::DRT::UTILS::GaussIntegration>& intpoints,
          const CORE::GEO::CUT::plain_volumecell_set& cells, bool offdiag = false) override
      {
        dserror("Not implemented!");
        return 1;
      }


      void ElementXfemInterfaceHybridLM(DRT::ELEMENTS::Fluid* ele,   ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<CORE::DRT::UTILS::GaussIntegration>&
              intpoints,  ///< element gauss points
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              side_coupling,                 ///< side coupling matrices
          Teuchos::ParameterList& params,    ///< parameter list
          Teuchos::RCP<MAT::Material>& mat,  ///< material
          CORE::LINALG::SerialDenseMatrix&
              elemat1_epetra,  ///< local system matrix of intersected element
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,                      ///< local element vector of intersected element
          CORE::LINALG::SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const CORE::GEO::CUT::plain_volumecell_set& vcSet  ///< set of plain volume cells
          ) override
      {
        dserror("Not implemented!");
        return;
      }

      void ElementXfemInterfaceNIT(DRT::ELEMENTS::Fluid* ele,        ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<CORE::DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,                     ///< parameter list
          Teuchos::RCP<MAT::Material>& mat_master,            ///< material master side
          Teuchos::RCP<MAT::Material>& mat_slave,             ///< material slave side
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,    ///< element matrix
          CORE::LINALG::SerialDenseVector& elevec1_epetra,    ///< element vector
          const CORE::GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<CORE::LINALG::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          CORE::LINALG::SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
          ) override
      {
        dserror("Not implemented!");
        return;
      }

      void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1_epetra,
          const CORE::DRT::UTILS::GaussIntegration& intpoints) override
      {
        dserror("Not implemented!");
        return;
      }

      void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1_epetra) override
      {
        dserror("Not implemented!");
        return;
      }

      /// Singleton access method
      static FluidEleCalcHDGWeakComp<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcHDGWeakComp();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
        using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;

        static const unsigned int nsd_ = FluidEleCalcHDGWeakComp<distype>::nsd_;
        static const unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;
        static const unsigned int nfaces_ = FluidEleCalcHDGWeakComp<distype>::nfaces_;

        /// local solver
        LocalSolver(const DRT::ELEMENTS::Fluid* ele,
            const CORE::DRT::UTILS::ShapeValues<distype>& shapeValues,
            CORE::DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly);

        /// initialize all
        void InitializeAll();

        /// compute material matrix
        void ComputeMaterialMatrix(const Teuchos::RCP<MAT::Material>& mat,
            const CORE::LINALG::Matrix<nsd_, 1>& xyz, CORE::LINALG::SerialDenseMatrix& DL,
            CORE::LINALG::SerialDenseMatrix& Dw);

        /// compute interior residual
        void ComputeInteriorResidual(const Teuchos::RCP<MAT::Material>& mat,
            const std::vector<double>& valnp, const std::vector<double>& accel,
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
        void CondenseLocalResidual(CORE::LINALG::SerialDenseVector& eleVec);

        /// condense local matrix
        void CondenseLocalMatrix(CORE::LINALG::SerialDenseMatrix& eleMat);

        // print matrices and residuals
        void PrintMatricesAndResiduals(DRT::ELEMENTS::Fluid& ele,
            CORE::LINALG::SerialDenseVector& eleVec, CORE::LINALG::SerialDenseMatrix& eleMat);

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

        const CORE::DRT::UTILS::ShapeValues<distype>& shapes_;    /// evaluated shape values
        CORE::DRT::UTILS::ShapeValuesFace<distype>& shapesface_;  /// evaluated shape values

        // Stabilization parameters
        double tau_r;  /// stabilization of density
        double tau_w;  /// stabilization of momentum

        // Auxiliary matrices
        CORE::LINALG::SerialDenseMatrix massPart;  /// temporary matrix for mass matrix
        CORE::LINALG::SerialDenseMatrix
            massPartW;                            /// temporary matrix for mass matrix with weights
        CORE::LINALG::SerialDenseMatrix massMat;  /// local mass matrix

        // Unknown variables
        CORE::LINALG::SerialDenseMatrix
            Leg;  /// mixed variable evaluated on interior quadrature points
        CORE::LINALG::SerialDenseVector reg;  /// density evaluated on interior quadrature points
        CORE::LINALG::SerialDenseMatrix weg;  /// momentum evaluated on interior quadrature points
        CORE::LINALG::SerialDenseVector
            rhatefg;  /// trace of density evaluated on face quadrature points
        CORE::LINALG::SerialDenseMatrix
            whatefg;  /// trace of momentum evaluated on face quadrature points

        // ALE variables
        CORE::LINALG::SerialDenseMatrix
            aeg;  /// ALE velocity evaluated on interior quadrature points
        CORE::LINALG::SerialDenseMatrix aefg;  /// ALE velocity evaluated on face quadrature points
        CORE::LINALG::SerialDenseMatrix
            dadxyzeg;  /// derivatives of ALE velocity evaluated on interior quadrature points

        // Matrices
        CORE::LINALG::SerialDenseMatrix ALL;  /// matrix mixed variable - mixed variable
        CORE::LINALG::SerialDenseMatrix ALr;  /// matrix mixed variable - density
        CORE::LINALG::SerialDenseMatrix ALw;  /// matrix mixed variable - momentum
        CORE::LINALG::SerialDenseMatrix ALR;  /// matrix mixed variable - trace of density
        CORE::LINALG::SerialDenseMatrix ALW;  /// matrix mixed variable - trace of momentum
        CORE::LINALG::SerialDenseMatrix Arr;  /// matrix density - density
        CORE::LINALG::SerialDenseMatrix Arw;  /// matrix density - momentum
        CORE::LINALG::SerialDenseMatrix ArR;  /// matrix density - trace of density
        CORE::LINALG::SerialDenseMatrix ArW;  /// matrix density - trace of momentum
        CORE::LINALG::SerialDenseMatrix AwL;  /// matrix momentum - mixed variable
        CORE::LINALG::SerialDenseMatrix Awr;  /// matrix momentum - density
        CORE::LINALG::SerialDenseMatrix Aww;  /// matrix momentum - momentum
        CORE::LINALG::SerialDenseMatrix AwR;  /// matrix momentum - trace of density
        CORE::LINALG::SerialDenseMatrix AwW;  /// matrix momentum - trace of momentum
        CORE::LINALG::SerialDenseMatrix ARr;  /// matrix trace of density - densitym
        CORE::LINALG::SerialDenseMatrix ARR;  /// matrix trace of density - trace of density
        CORE::LINALG::SerialDenseMatrix AWL;  /// matrix trace of momentum - mixed variable
        CORE::LINALG::SerialDenseMatrix AWw;  /// matrix trace of momentum - momentum
        CORE::LINALG::SerialDenseMatrix AWR;  /// matrix trace of momentum - trace of density
        CORE::LINALG::SerialDenseMatrix AWW;  /// matrix trace of momentum - trace of momentum

        // Residuals
        CORE::LINALG::SerialDenseVector RL;  /// residual vector for mixed variable
        CORE::LINALG::SerialDenseVector Rr;  /// residual vector for density
        CORE::LINALG::SerialDenseVector Rw;  /// residual vector for momentu
        CORE::LINALG::SerialDenseVector RR;  /// residual vector for trace of density
        CORE::LINALG::SerialDenseVector RW;  /// residual vector for trace of momentum

        // Local/Global matrices/vectors
        CORE::LINALG::SerialDenseMatrix Klocallocal;     /// local-local matrix
        CORE::LINALG::SerialDenseMatrix Klocalglobal;    /// local-global matrix
        CORE::LINALG::SerialDenseMatrix Kgloballocal;    /// global-local matrix
        CORE::LINALG::SerialDenseMatrix Kglobalglobal;   /// global-global matrix
        CORE::LINALG::SerialDenseVector Rlocal;          /// local residual vector
        CORE::LINALG::SerialDenseVector Rglobal;         /// global residual vector
        CORE::LINALG::SerialDenseMatrix KlocallocalInv;  /// inverse local-local matrix
        Teuchos::SerialDenseSolver<ordinalType, scalarType>
            KlocallocalInvSolver;  /// solver for inverse local-local matrix

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
      void EvaluateAll(const int funcnum, const CORE::LINALG::Matrix<nsd_, 1>& xyz, const double t,
          CORE::LINALG::Matrix<msd_, 1>& L, double& r, CORE::LINALG::Matrix<nsd_, 1>& w) const;

      /// evaluate density and momentum
      void EvaluateDensityMomentum(const int funcnum, const CORE::LINALG::Matrix<nsd_, 1>& xyz,
          const double t, double& r, CORE::LINALG::Matrix<nsd_, 1>& w) const;

      /// local data object
      Teuchos::RCP<CORE::DRT::UTILS::ShapeValues<distype>> shapes_;
      Teuchos::RCP<CORE::DRT::UTILS::ShapeValuesFace<distype>> shapesface_;

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
