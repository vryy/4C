/*----------------------------------------------------------------------------*/
/*! \file
\brief Routines for calculation of HDG weakly compressible fluid element

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_HDG_WEAK_COMP_HPP
#define FOUR_C_FLUID_ELE_CALC_HDG_WEAK_COMP_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_shapevalues_hdg.hpp"
#include "4C_fluid_ele_hdg_weak_comp.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Weakly compressible fluid HDG element implementation
    /*!

      \author laspina
      \date 08/19
    */
    template <Core::FE::CellType distype>
    class FluidEleCalcHDGWeakComp : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr unsigned int nen_ = Core::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr unsigned int nsd_ = Core::FE::dim<distype>;

      //! mixed variable dimension according to Voigt notation
      static constexpr unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;

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

      /// Update local solution
      virtual int UpdateLocalSolution(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& interiorinc);

      /// projection of function field
      virtual int ProjectField(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      /*! \brief Interpolates an HDG solution to the element nodes for output
       */
      virtual int interpolate_solution_to_nodes(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& elevec1);

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

      /// Singleton access method
      static FluidEleCalcHDGWeakComp<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcHDGWeakComp();

      /// local solver that inverts local problem on an element and can solve with various vectors
      struct LocalSolver
      {
        using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
        using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;

        static constexpr unsigned int nsd_ = FluidEleCalcHDGWeakComp<distype>::nsd_;
        static constexpr unsigned int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;
        static constexpr unsigned int nfaces_ = FluidEleCalcHDGWeakComp<distype>::nfaces_;

        /// local solver
        LocalSolver(const Discret::ELEMENTS::Fluid* ele,
            const Core::FE::ShapeValues<distype>& shapeValues,
            Core::FE::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly);

        /// initialize all
        void InitializeAll();

        /// compute material matrix
        void compute_material_matrix(const Teuchos::RCP<Core::Mat::Material>& mat,
            const Core::LinAlg::Matrix<nsd_, 1>& xyz, Core::LinAlg::SerialDenseMatrix& DL,
            Core::LinAlg::SerialDenseMatrix& Dw);

        /// compute interior residual
        void compute_interior_residual(const Teuchos::RCP<Core::Mat::Material>& mat,
            const std::vector<double>& valnp, const std::vector<double>& accel,
            const std::vector<double>& alevel);

        /// compute face residual
        void ComputeFaceResidual(const int f, const Teuchos::RCP<Core::Mat::Material>& mat,
            const std::vector<double>& val, const std::vector<double>& traceval,
            const std::vector<double>& alevel);

        /// compute interior matrices
        void compute_interior_matrices(const Teuchos::RCP<Core::Mat::Material>& mat);

        /// compute face matrices
        void ComputeFaceMatrices(const int f, const Teuchos::RCP<Core::Mat::Material>& mat);

        /// compute local residual
        void compute_local_residual();

        /// compute global residual
        void compute_global_residual(Discret::ELEMENTS::Fluid& ele);

        /// compute local-local matrix
        void compute_local_local_matrix();

        /// compute local-global matrix
        void compute_local_global_matrix(Discret::ELEMENTS::Fluid& ele);

        /// compute global-local matrix
        void compute_global_local_matrix(Discret::ELEMENTS::Fluid& ele);

        /// compute global-global matrix
        void compute_global_global_matrix(Discret::ELEMENTS::Fluid& ele);

        /// invert local-local matrix
        void invert_local_local_matrix();

        /// condense local residual
        void condense_local_residual(Core::LinAlg::SerialDenseVector& eleVec);

        /// condense local matrix
        void CondenseLocalMatrix(Core::LinAlg::SerialDenseMatrix& eleMat);

        // print matrices and residuals
        void print_matrices_and_residuals(Discret::ELEMENTS::Fluid& ele,
            Core::LinAlg::SerialDenseVector& eleVec, Core::LinAlg::SerialDenseMatrix& eleMat);

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

        const Core::FE::ShapeValues<distype>& shapes_;    /// evaluated shape values
        Core::FE::ShapeValuesFace<distype>& shapesface_;  /// evaluated shape values

        // Stabilization parameters
        double tau_r;  /// stabilization of density
        double tau_w;  /// stabilization of momentum

        // Auxiliary matrices
        Core::LinAlg::SerialDenseMatrix massPart;  /// temporary matrix for mass matrix
        Core::LinAlg::SerialDenseMatrix
            massPartW;                            /// temporary matrix for mass matrix with weights
        Core::LinAlg::SerialDenseMatrix massMat;  /// local mass matrix

        // Unknown variables
        Core::LinAlg::SerialDenseMatrix
            Leg;  /// mixed variable evaluated on interior quadrature points
        Core::LinAlg::SerialDenseVector reg;  /// density evaluated on interior quadrature points
        Core::LinAlg::SerialDenseMatrix weg;  /// momentum evaluated on interior quadrature points
        Core::LinAlg::SerialDenseVector
            rhatefg;  /// trace of density evaluated on face quadrature points
        Core::LinAlg::SerialDenseMatrix
            whatefg;  /// trace of momentum evaluated on face quadrature points

        // ALE variables
        Core::LinAlg::SerialDenseMatrix
            aeg;  /// ALE velocity evaluated on interior quadrature points
        Core::LinAlg::SerialDenseMatrix aefg;  /// ALE velocity evaluated on face quadrature points
        Core::LinAlg::SerialDenseMatrix
            dadxyzeg;  /// derivatives of ALE velocity evaluated on interior quadrature points

        // Matrices
        Core::LinAlg::SerialDenseMatrix ALL;  /// matrix mixed variable - mixed variable
        Core::LinAlg::SerialDenseMatrix ALr;  /// matrix mixed variable - density
        Core::LinAlg::SerialDenseMatrix ALw;  /// matrix mixed variable - momentum
        Core::LinAlg::SerialDenseMatrix ALR;  /// matrix mixed variable - trace of density
        Core::LinAlg::SerialDenseMatrix ALW;  /// matrix mixed variable - trace of momentum
        Core::LinAlg::SerialDenseMatrix Arr;  /// matrix density - density
        Core::LinAlg::SerialDenseMatrix Arw;  /// matrix density - momentum
        Core::LinAlg::SerialDenseMatrix ArR;  /// matrix density - trace of density
        Core::LinAlg::SerialDenseMatrix ArW;  /// matrix density - trace of momentum
        Core::LinAlg::SerialDenseMatrix AwL;  /// matrix momentum - mixed variable
        Core::LinAlg::SerialDenseMatrix Awr;  /// matrix momentum - density
        Core::LinAlg::SerialDenseMatrix Aww;  /// matrix momentum - momentum
        Core::LinAlg::SerialDenseMatrix AwR;  /// matrix momentum - trace of density
        Core::LinAlg::SerialDenseMatrix AwW;  /// matrix momentum - trace of momentum
        Core::LinAlg::SerialDenseMatrix ARr;  /// matrix trace of density - densitym
        Core::LinAlg::SerialDenseMatrix ARR;  /// matrix trace of density - trace of density
        Core::LinAlg::SerialDenseMatrix AWL;  /// matrix trace of momentum - mixed variable
        Core::LinAlg::SerialDenseMatrix AWw;  /// matrix trace of momentum - momentum
        Core::LinAlg::SerialDenseMatrix AWR;  /// matrix trace of momentum - trace of density
        Core::LinAlg::SerialDenseMatrix AWW;  /// matrix trace of momentum - trace of momentum

        // Residuals
        Core::LinAlg::SerialDenseVector RL;  /// residual vector for mixed variable
        Core::LinAlg::SerialDenseVector Rr;  /// residual vector for density
        Core::LinAlg::SerialDenseVector Rw;  /// residual vector for momentu
        Core::LinAlg::SerialDenseVector RR;  /// residual vector for trace of density
        Core::LinAlg::SerialDenseVector RW;  /// residual vector for trace of momentum

        // Local/Global matrices/vectors
        Core::LinAlg::SerialDenseMatrix Klocallocal;     /// local-local matrix
        Core::LinAlg::SerialDenseMatrix Klocalglobal;    /// local-global matrix
        Core::LinAlg::SerialDenseMatrix Kgloballocal;    /// global-local matrix
        Core::LinAlg::SerialDenseMatrix Kglobalglobal;   /// global-global matrix
        Core::LinAlg::SerialDenseVector Rlocal;          /// local residual vector
        Core::LinAlg::SerialDenseVector Rglobal;         /// global residual vector
        Core::LinAlg::SerialDenseMatrix KlocallocalInv;  /// inverse local-local matrix
        Teuchos::SerialDenseSolver<ordinalType, scalarType>
            KlocallocalInvSolver;  /// solver for inverse local-local matrix

        // Voigt related quantities
        int VoigtP[msd_ - nsd_][2];  /// pair of indices in Voigt notation

        std::vector<int> pivots;  /// pivots for factorization of matrices

        Teuchos::RCP<Discret::ELEMENTS::FluidEleParameter> fldpara_;  //! pointer to parameter list
        Teuchos::RCP<Discret::ELEMENTS::FluidEleParameterTimInt>
            fldparatimint_;  //! pointer to time parameter list
      };

      /// reads from global vectors
      void read_global_vectors(const Core::Elements::Element& ele,
          Discret::Discretization& discretization, const std::vector<int>& lm);

      /// reads ale vectors
      void read_ale_vectors(
          const Core::Elements::Element& ele, Discret::Discretization& discretization);

      /// evaluate mixed variable, density and momentum
      void evaluate_all(const int funcnum, const Core::LinAlg::Matrix<nsd_, 1>& xyz, const double t,
          Core::LinAlg::Matrix<msd_, 1>& L, double& r, Core::LinAlg::Matrix<nsd_, 1>& w) const;

      /// evaluate density and momentum
      void evaluate_density_momentum(const int funcnum, const Core::LinAlg::Matrix<nsd_, 1>& xyz,
          const double t, double& r, Core::LinAlg::Matrix<nsd_, 1>& w) const;

      /// local data object
      Teuchos::RCP<Core::FE::ShapeValues<distype>> shapes_;
      Teuchos::RCP<Core::FE::ShapeValuesFace<distype>> shapesface_;

      /// local solver object
      Teuchos::RCP<LocalSolver> local_solver_;

      std::vector<double> trace_val_;  /// extracted values from trace solution vector at n+alpha_f
      std::vector<double> interior_val_;  /// extracted local values at n+alpha_f
      std::vector<double> interior_acc_;  /// extracted local accelerations at n+alpha_m
      std::vector<double> ale_dis_;       /// extracted ale mesh displacement
      std::vector<double> ale_vel_;       /// extracted ale mesh velocity

      bool usescompletepoly_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
