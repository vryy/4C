/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

\level 2


*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
#define BACI_SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H

#include "baci_config.H"

#include "baci_coupling_adapter.H"
#include "baci_coupling_adapter_mortar.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_inpar_s2i.H"
#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"
#include "baci_scatra_timint_meshtying_strategy_base.H"

#include <Epetra_FEVector.h>
#include <Epetra_IntVector.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingMortar;
}  // namespace ADAPTER

namespace DRT
{
  class Condition;

  namespace ELEMENTS
  {
    class ScaTraEleParameterBoundary;
  }

}  // namespace DRT

namespace CORE::FE
{
  template <const int NSD>
  class IntPointsAndWeights;
}

namespace CORE::LINALG
{
  class MatrixColTransform;
  class MatrixRowTransform;
  class MatrixRowColTransform;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class Equilibration;
  enum class MatrixType;
}  // namespace CORE::LINALG

namespace MORTAR
{
  class IntCell;
  class Element;
  class Node;
}  // namespace MORTAR

namespace SCATRA
{
  // forward declaration
  class MortarCellAssemblyStrategy;

  /*!
  \brief Scatra-scatra interface coupling strategy for standard scalar transport problems

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the scatra-scatra interface coupling strategy for
  standard scalar transport problems.

  */

  class MeshtyingStrategyS2I : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyS2I(
        SCATRA::ScaTraTimIntImpl* scatratimint,  //!< scalar transport time integrator
        const Teuchos::ParameterList&
            parameters  //!< input parameters for scatra-scatra interface coupling
    );

    //! provide global state vectors for element evaluation
    void AddTimeIntegrationSpecificVectors() const override;

    //! compute time step size
    void ComputeTimeStepSize(double& dt) override;

    //! return map extractor associated with blocks of auxiliary system matrix for master side
    const CORE::LINALG::MultiMapExtractor& BlockMapsMaster() const { return *blockmaps_master_; };

    //! return map extractor associated with blocks of auxiliary system matrix for slave side
    const CORE::LINALG::MultiMapExtractor& BlockMapsSlave() const { return *blockmaps_slave_; };

    //! compute time derivatives of discrete state variables
    void ComputeTimeDerivative() const override;

    void CondenseMatAndRHS(const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const bool calcinittimederiv = false) const override;

    //! return interface coupling adapter
    Teuchos::RCP<const CORE::ADAPTER::Coupling> CouplingAdapter() const { return icoup_; };

    //! return flag for meshtying method
    const INPAR::S2I::CouplingType& CouplingType() const { return couplingtype_; }

    //! return global map of degrees of freedom
    const Epetra_Map& DofRowMap() const override;

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying() override;

    /*!
     * @brief  evaluate mortar integration cells
     *
     * @param idiscret  interface discretization
     * @param params    parameter list for evaluation of mortar integration cells
     * @param strategy  assembly strategy for mortar integration cells
     */
    void EvaluateMortarCells(const DRT::Discretization& idiscret,
        const Teuchos::ParameterList& params, SCATRA::MortarCellAssemblyStrategy& strategy) const;

    //! explicit predictor step to obtain better starting value for Newton-Raphson iteration
    void ExplicitPredictor() const override;

    //! extract selected rows from a sparse matrix
    static void ExtractMatrixRows(const CORE::LINALG::SparseMatrix& matrix,  //!< source matrix
        CORE::LINALG::SparseMatrix& rows,                                    //!< destination matrix
        const Epetra_Map& rowmap  //!< map of matrix rows to be extracted
    );

    /*!
     * @brief finite difference check for extended system matrix involving scatra-scatra interface
     * layer growth (for debugging only)
     *
     * @param extendedsystemmatrix  global system matrix
     * @param extendedresidual      global residual vector
     */
    void FDCheck(const CORE::LINALG::BlockSparseMatrixBase& extendedsystemmatrix,
        const Teuchos::RCP<Epetra_Vector>& extendedresidual) const;

    //! return state vector of discrete scatra-scatra interface layer thicknesses at time n
    const Teuchos::RCP<Epetra_Vector>& GrowthVarN() const { return growthn_; };

    //! return state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    const Teuchos::RCP<Epetra_Vector>& GrowthVarNp() const { return growthnp_; };

    //! perform initialization of scatra-scatra interface coupling
    void InitMeshtying() override;

    bool SystemMatrixInitializationNeeded() const override { return false; }

    Teuchos::RCP<CORE::LINALG::SparseOperator> InitSystemMatrix() const override
    {
      dserror(
          "This meshtying strategy does not need to initialize the system matrix, but relies "
          "instead on the initialization of the field. If this changes, you also need to change "
          "'SystemMatrixInitializationNeeded()' to return true");
      // dummy return
      return Teuchos::null;
    }

    //! return interface map extractor
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> InterfaceMaps() const override
    {
      return interfacemaps_;
    }

    //! return flag for evaluation of scatra-scatra interface coupling involving interface layer
    //! growth
    const INPAR::S2I::GrowthEvaluation& IntLayerGrowthEvaluation() const
    {
      return intlayergrowth_evaluation_;
    };

    //! return the slave-side scatra-scatra interface kinetics conditions applied to a mesh tying
    //! interface
    const std::map<const int, DRT::Condition* const>& KineticsConditionsMeshtyingSlaveSide() const
    {
      return kinetics_conditions_meshtying_slaveside_;
    }

    //! corresponding master conditions to kinetics condiditions
    std::map<const int, DRT::Condition* const>& MasterConditions() { return master_conditions_; }

    //! return vector of Lagrange multiplier dofs
    Teuchos::RCP<const Epetra_Vector> LM() const { return lm_; };

    //! return constraint residual vector associated with Lagrange multiplier dofs
    Teuchos::RCP<const Epetra_Vector> LMResidual() const { return lmresidual_; };

    //! return constraint increment vector associated with Lagrange multiplier dofs
    Teuchos::RCP<const Epetra_Vector> LMIncrement() const { return lmincrement_; };

    //! return auxiliary system matrix for linearizations of slave fluxes w.r.t. master dofs
    const Teuchos::RCP<CORE::LINALG::SparseMatrix>& MasterMatrix() const { return imastermatrix_; };

    //! return type of global system matrix in global system of equations
    const CORE::LINALG::MatrixType& MatrixType() const { return matrixtype_; };

    //! return mortar interface discretization associated with particular condition ID
    DRT::Discretization& MortarDiscretization(const int& condid) const;

    //! output solution for post-processing
    void Output() const override;

    void WriteRestart() const override;

    //! return mortar projector P
    const Teuchos::RCP<CORE::LINALG::SparseMatrix>& P() const { return P_; };

    void ReadRestart(
        const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) const override;

    //! set general parameters for element evaluation
    void SetElementGeneralParameters(Teuchos::ParameterList& parameters) const override;

    /*!
     * \brief Method sets the scatra-scatra interface condition specific values to the scatra
     * element interface condition.
     *
     * \note Parameters are stored to the parameter class using the evaluate call at the end of this
     * method.
     *
     * @param[in] s2icondition Scatra-scatra interface condition of which parameters are read and
     * stored to the parameter class
     */
    void SetConditionSpecificScaTraParameters(DRT::Condition& s2icondition) const;

    /*!
     * \brief Writes S2IKinetics condition specific parameters to parameter list that is stored to
     * the boundary parameter class afterwards
     *
     * @param[in]  s2ikinetics_cond       ScaTra-ScaTra interface condition whose parameters are
     *                                    stored to the parameter list
     * @param[out] s2icouplingparameters  parameter list filled with condition specific parameters
     */
    static void WriteS2IKineticsSpecificScaTraParametersToParameterList(
        DRT::Condition& s2ikinetics_cond, Teuchos::ParameterList& s2icouplingparameters);

    //! compute history vector, i.e., the history part of the right-hand side vector with all
    //! contributions from the previous time step
    void SetOldPartOfRHS() const override;

    //! perform setup of scatra-scatra interface coupling
    void SetupMeshtying() override;

    //! return auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs
    //! (non-mortar case) or slave and master dofs (mortar case)
    const Teuchos::RCP<CORE::LINALG::SparseMatrix>& SlaveMatrix() const { return islavematrix_; };

    void Solve(const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& increment, const Teuchos::RCP<Epetra_Vector>& residual,
        const Teuchos::RCP<Epetra_Vector>& phinp, const int iteration,
        CORE::LINALG::SolverParams& solver_params) const override;

    //! return linear solver for global system of linear equations
    const CORE::LINALG::Solver& Solver() const override;

    //! update solution after convergence of the nonlinear Newton-Raphson iteration
    void Update() const override;

   protected:
    void EquipExtendedSolverWithNullSpaceInfo() const override;

    //! instantiate strategy for Newton-Raphson convergence check
    void InitConvCheckStrategy() override;

    //! interface map extractor (0: other, 1: slave, 2: master)
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> interfacemaps_;

    //! map extractor associated with scatra-scatra interface slave-side blocks of global system
    //! matrix
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockmaps_slave_;
    //! map extractor associated with scatra-scatra interface master-side blocks of global system
    //! matrix
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockmaps_master_;

    //! non-mortar interface coupling adapter
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoup_;

    //! mortar interface coupling adapters
    std::map<int, Teuchos::RCP<CORE::ADAPTER::CouplingMortar>> icoupmortar_;

    //! mortar integration cells
    std::map<int, std::vector<std::pair<Teuchos::RCP<MORTAR::IntCell>, INPAR::SCATRA::ImplType>>>
        imortarcells_;

    //! flag for parallel redistribution of mortar interfaces
    const bool imortarredistribution_;

    //! map of all slave-side degrees of freedom before parallel redistribution
    Teuchos::RCP<Epetra_Map> islavemap_;

    //! map of all master-side degrees of freedom before parallel redistribution
    Teuchos::RCP<Epetra_Map> imastermap_;

    //! vectors for node-to-segment connectivity, i.e., for pairings between slave nodes and master
    //! elements
    std::map<int, Teuchos::RCP<Epetra_IntVector>> islavenodestomasterelements_;

    //! vectors for physical implementation types of slave-side nodes
    std::map<int, Teuchos::RCP<Epetra_IntVector>> islavenodesimpltypes_;

    //! vectors for lumped interface area fractions associated with slave-side nodes
    std::map<int, Teuchos::RCP<Epetra_Vector>> islavenodeslumpedareas_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs (non-mortar
    //! case) or for linearizations of slave fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> islavematrix_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. master dofs (non-mortar
    //! case) or for linearizations of master fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> imastermatrix_;

    //! auxiliary system matrix for linearizations of master fluxes w.r.t. slave dofs
    Teuchos::RCP<CORE::LINALG::SparseMatrix> imasterslavematrix_;

    //! flag for meshtying method
    const INPAR::S2I::CouplingType couplingtype_;

    //! mortar matrix D
    Teuchos::RCP<CORE::LINALG::SparseMatrix> D_;

    //! mortar matrix M
    Teuchos::RCP<CORE::LINALG::SparseMatrix> M_;

    //! mortar matrix E
    Teuchos::RCP<CORE::LINALG::SparseMatrix> E_;

    //! mortar projector P
    Teuchos::RCP<CORE::LINALG::SparseMatrix> P_;

    //! mortar projector Q
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Q_;

    //! vector of Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lm_;

    //! extended map extractor (0: standard dofs, 1: Lagrange multiplier dofs or scatra-scatra
    //! interface layer thickness variables)
    Teuchos::RCP<CORE::LINALG::MapExtractor> extendedmaps_;

    //! constraint residual vector associated with Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lmresidual_;

    //! constraint increment vector associated with Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lmincrement_;

    //! transformation operators for auxiliary system matrices
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> islavetomastercoltransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> islavetomasterrowtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> islavetomasterrowcoltransform_;

    //! auxiliary residual vector for slave residuals
    Teuchos::RCP<Epetra_Vector> islaveresidual_;

    //! auxiliary residual vector for master residuals
    Teuchos::RCP<Epetra_FEVector> imasterresidual_;

    //! time derivative of slave dofs of scatra-scatra interface
    Teuchos::RCP<Epetra_Vector> islavephidtnp_;

    //! time derivative of master dofs transformed to slave side of scatra-scatra interface
    Teuchos::RCP<Epetra_Vector> imasterphidt_on_slave_side_np_;

    //! master dofs transformed to slave side of scatra-scatra interface
    Teuchos::RCP<Epetra_Vector> imasterphi_on_slave_side_np_;

    //! flag for interface side underlying Lagrange multiplier definition
    const INPAR::S2I::InterfaceSides lmside_;

    //! type of global system matrix in global system of equations
    const CORE::LINALG::MatrixType matrixtype_;

    //! node-to-segment projection tolerance
    const double ntsprojtol_;

    //! flag for evaluation of scatra-scatra interface coupling involving interface layer growth
    const INPAR::S2I::GrowthEvaluation intlayergrowth_evaluation_;

    //! local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
    //! interface layer growth
    const double intlayergrowth_convtol_;

    //! maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
    //! involving interface layer growth
    const unsigned intlayergrowth_itemax_;

    //! modified time step size for scatra-scatra interface coupling involving interface layer
    //! growth
    const double intlayergrowth_timestep_;

    //! map extractor associated with all degrees of freedom for scatra-scatra interface layer
    //! growth
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmapgrowth_;

    //! extended map extractor associated with blocks of global system matrix for scatra-scatra
    //! interface coupling involving interface layer growth
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> extendedblockmaps_;

    //! extended system matrix including rows and columns associated with scatra-scatra interface
    //! layer thickness variables
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> extendedsystemmatrix_;

    //! linear solver for monolithic scatra-scatra interface coupling involving interface layer
    //! growth
    Teuchos::RCP<CORE::LINALG::Solver> extendedsolver_;

    //! state vector of discrete scatra-scatra interface layer thicknesses at time n
    Teuchos::RCP<Epetra_Vector> growthn_;

    //! state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    Teuchos::RCP<Epetra_Vector> growthnp_;

    //! state vector of time derivatives of discrete scatra-scatra interface layer thicknesses at
    //! time n
    Teuchos::RCP<Epetra_Vector> growthdtn_;

    //! state vector of time derivatives of discrete scatra-scatra interface layer thicknesses at
    //! time n+1
    Teuchos::RCP<Epetra_Vector> growthdtnp_;

    //! state vector of history values associated with discrete scatra-scatra interface layer
    //! thicknesses
    Teuchos::RCP<Epetra_Vector> growthhist_;

    //! state vector of residual values associated with discrete scatra-scatra interface layer
    //! thicknesses
    Teuchos::RCP<Epetra_Vector> growthresidual_;

    //! state vector of Newton-Raphson increment values associated with discrete scatra-scatra
    //! interface layer thicknesses
    Teuchos::RCP<Epetra_Vector> growthincrement_;

    //! scatra-growth block of extended global system matrix (derivatives of discrete scatra
    //! residuals w.r.t. discrete scatra-scatra interface layer thicknesses)
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatragrowthblock_;

    //! growth-scatra block of extended global system matrix (derivatives of discrete scatra-scatra
    //! interface layer growth residuals w.r.t. discrete scatra degrees of freedom)
    Teuchos::RCP<CORE::LINALG::SparseOperator> growthscatrablock_;

    //! growth-growth block of extended global system matrix (derivatives of discrete scatra-scatra
    //! interface layer growth residuals w.r.t. discrete scatra-scatra interface layer thicknesses)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> growthgrowthblock_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<CORE::LINALG::Equilibration> equilibration_;

   private:
    //! copy constructor
    MeshtyingStrategyS2I(const MeshtyingStrategyS2I& old);

    //! build map extractors associated with blocks of global system matrix
    void BuildBlockMapExtractors();

    //! evaluate and assemble all contributions due to capacitive fluxes at the scatra-scatra
    //! interface
    void EvaluateAndAssembleCapacitiveContributions();

    /*!
     * @brief  evaluate single mortar integration cell
     *
     * @param idiscret       interface discretization
     * @param cell           mortar integration cell
     * @param impltype       physical implementation type of mortar integration cell
     * @param slaveelement   slave-side mortar element
     * @param masterelement  master-side mortar element
     * @param la_slave       slave-side location array
     * @param la_master      master-side location array
     * @param params         parameter list
     * @param cellmatrix1    cell matrix 1
     * @param cellmatrix2    cell matrix 2
     * @param cellmatrix3    cell matrix 3
     * @param cellmatrix4    cell matrix 4
     * @param cellvector1    cell vector 1
     * @param cellvector2    cell vector 2
     */
    void EvaluateMortarCell(const DRT::Discretization& idiscret, MORTAR::IntCell& cell,
        const INPAR::SCATRA::ImplType& impltype, MORTAR::Element& slaveelement,
        MORTAR::Element& masterelement, DRT::Element::LocationArray& la_slave,
        DRT::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
        CORE::LINALG::SerialDenseMatrix& cellmatrix1, CORE::LINALG::SerialDenseMatrix& cellmatrix2,
        CORE::LINALG::SerialDenseMatrix& cellmatrix3, CORE::LINALG::SerialDenseMatrix& cellmatrix4,
        CORE::LINALG::SerialDenseVector& cellvector1,
        CORE::LINALG::SerialDenseVector& cellvector2) const;

    /*!
     * @brief  evaluate single slave-side node for node-to-segment coupling
     *
     * @param idiscret       interface discretization
     * @param slavenode      slave-side node
     * @param lumpedarea     lumped interface area fraction associated with slave-side node
     * @param impltype       physical implementation type of mortar integration cell
     * @param slaveelement   slave-side mortar element
     * @param masterelement  master-side mortar element
     * @param la_slave       slave-side location array
     * @param la_master      master-side location array
     * @param params         parameter list
     * @param ntsmatrix1     node-to-segment matrix 1
     * @param ntsmatrix2     node-to-segment matrix 2
     * @param ntsmatrix3     node-to-segment matrix 3
     * @param ntsmatrix4     node-to-segment matrix 4
     * @param ntsvector1     node-to-segment vector 1
     * @param ntsvector2     node-to-segment vector 2
     */
    void EvaluateSlaveNode(const DRT::Discretization& idiscret, const MORTAR::Node& slavenode,
        const double& lumpedarea, const INPAR::SCATRA::ImplType& impltype,
        MORTAR::Element& slaveelement, MORTAR::Element& masterelement,
        DRT::Element::LocationArray& la_slave, DRT::Element::LocationArray& la_master,
        const Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& ntsmatrix1,
        CORE::LINALG::SerialDenseMatrix& ntsmatrix2, CORE::LINALG::SerialDenseMatrix& ntsmatrix3,
        CORE::LINALG::SerialDenseMatrix& ntsmatrix4, CORE::LINALG::SerialDenseVector& ntsvector1,
        CORE::LINALG::SerialDenseVector& ntsvector2) const;

    /*!
     * @brief  evaluate single mortar element
     *
     * @param idiscret   interface discretization
     * @param element    mortar element
     * @param impltype   physical implementation type of mortar element
     * @param la         location array
     * @param params     parameter list
     * @param elematrix1 element matrix 1
     * @param elematrix2 element matrix 2
     * @param elematrix3 element matrix 3
     * @param elematrix4 element matrix 4
     * @param elevector1 element vector 1
     * @param elevector2 element vector 2
     */
    void EvaluateMortarElement(const DRT::Discretization& idiscret, MORTAR::Element& element,
        const INPAR::SCATRA::ImplType& impltype, DRT::Element::LocationArray& la,
        const Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& elematrix1,
        CORE::LINALG::SerialDenseMatrix& elematrix2, CORE::LINALG::SerialDenseMatrix& elematrix3,
        CORE::LINALG::SerialDenseMatrix& elematrix4, CORE::LINALG::SerialDenseVector& elevector1,
        CORE::LINALG::SerialDenseVector& elevector2) const;

    /*!
     * @brief  evaluate mortar integration cells
     *
     * @param idiscret           interface discretization
     * @param params             parameter list for evaluation of mortar integration cells
     * @param systemmatrix1      system matrix 1
     * @param matrix1_side_rows  interface side associated with rows of system matrix 1
     * @param matrix1_side_cols  interface side associated with columns of system matrix 1
     * @param systemmatrix2      system matrix 2
     * @param matrix2_side_rows  interface side associated with rows of system matrix 2
     * @param matrix2_side_cols  interface side associated with columns of system matrix 2
     * @param systemmatrix3      system matrix 3
     * @param matrix3_side_rows  interface side associated with rows of system matrix 3
     * @param matrix3_side_cols  interface side associated with columns of system matrix 3
     * @param systemmatrix4      system matrix 4
     * @param matrix4_side_rows  interface side associated with rows of system matrix 4
     * @param matrix4_side_cols  interface side associated with columns of system matrix 4
     * @param systemvector1      system vector 1
     * @param vector1_side       interface side associated with system vector 1
     * @param systemvector2      system vector 2
     * @param vector2_side       interface side associated with system vector 2
     */
    void EvaluateMortarCells(const DRT::Discretization& idiscret,
        const Teuchos::ParameterList& params,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix1,
        const INPAR::S2I::InterfaceSides matrix1_side_rows,
        const INPAR::S2I::InterfaceSides matrix1_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix2,
        const INPAR::S2I::InterfaceSides matrix2_side_rows,
        const INPAR::S2I::InterfaceSides matrix2_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix3,
        const INPAR::S2I::InterfaceSides matrix3_side_rows,
        const INPAR::S2I::InterfaceSides matrix3_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix4,
        const INPAR::S2I::InterfaceSides matrix4_side_rows,
        const INPAR::S2I::InterfaceSides matrix4_side_cols,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
        const INPAR::S2I::InterfaceSides vector1_side,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
        const INPAR::S2I::InterfaceSides vector2_side) const;

    /*!
     * @brief  evaluate node-to-segment coupling
     *
     * @param islavenodestomasterelements  vector for node-to-segment connectivity
     * @param islavenodeslumpedareas       vector for lumped interface area fractions associated
     *                                     with slave-side nodes
     * @param islavenodesimpltypes         vector for physical implementation types of slave-side
     *                                     nodes
     * @param idiscret                     interface discretization
     * @param params                       parameter list for evaluation of mortar integration cells
     * @param systemmatrix1                system matrix 1
     * @param matrix1_side_rows            interface side associated with rows of system matrix 1
     * @param matrix1_side_cols            interface side associated with columns of system matrix 1
     * @param systemmatrix2                system matrix 2
     * @param matrix2_side_rows            interface side associated with rows of system matrix 2
     * @param matrix2_side_cols            interface side associated with columns of system matrix 2
     * @param systemmatrix3                system matrix 3
     * @param matrix3_side_rows            interface side associated with rows of system matrix 3
     * @param matrix3_side_cols            interface side associated with columns of system matrix 3
     * @param systemmatrix4                system matrix 4
     * @param matrix4_side_rows     interface side associated with rows of system matrix 4
     * @param matrix4_side_cols     interface side associated with columns of system matrix 4
     * @param systemvector1         system vector 1
     * @param vector1_side          interface side associated with system vector 1
     * @param systemvector2         system vector 2
     * @param vector2_side          interface side associated with system vector 2
     */
    void EvaluateNTS(const Epetra_IntVector& islavenodestomasterelements,
        const Epetra_Vector& islavenodeslumpedareas, const Epetra_IntVector& islavenodesimpltypes,
        const DRT::Discretization& idiscret, const Teuchos::ParameterList& params,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix1,
        const INPAR::S2I::InterfaceSides matrix1_side_rows,
        const INPAR::S2I::InterfaceSides matrix1_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix2,
        const INPAR::S2I::InterfaceSides matrix2_side_rows,
        const INPAR::S2I::InterfaceSides matrix2_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix3,
        const INPAR::S2I::InterfaceSides matrix3_side_rows,
        const INPAR::S2I::InterfaceSides matrix3_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix4,
        const INPAR::S2I::InterfaceSides matrix4_side_rows,
        const INPAR::S2I::InterfaceSides matrix4_side_cols,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
        const INPAR::S2I::InterfaceSides vector1_side,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
        const INPAR::S2I::InterfaceSides vector2_side) const;

    /*!
     * @brief  evaluate mortar elements
     *
     * @param ielecolmap         column map of mortar elements
     * @param ieleimpltypes      vector for physical implementation types of mortar elements
     * @param idiscret           interface discretization
     * @param params             parameter list for evaluation of mortar integration cells
     * @param systemmatrix1      system matrix 1
     * @param matrix1_side_rows  interface side associated with rows of system matrix 1
     * @param matrix1_side_cols  interface side associated with columns of system matrix 1
     * @param systemmatrix2      system matrix 2
     * @param matrix2_side_rows  interface side associated with rows of system matrix 2
     * @param matrix2_side_cols  interface side associated with columns of system matrix 2
     * @param systemmatrix3      system matrix 3
     * @param matrix3_side_rows  interface side associated with rows of system matrix 3
     * @param matrix3_side_cols  interface side associated with columns of system matrix 3
     * @param systemmatrix4      system matrix 4
     * @param matrix4_side_rows  interface side associated with rows of system matrix 4
     * @param matrix4_side_cols  interface side associated with columns of system matrix 4
     * @param systemvector1      system vector 1
     * @param vector1_side       interface side associated with system vector 1
     * @param systemvector2      system vector 2
     * @param vector2_side       interface side associated with system vector 2
     */
    void EvaluateMortarElements(const Epetra_Map& ielecolmap, const Epetra_IntVector& ieleimpltypes,
        const DRT::Discretization& idiscret, const Teuchos::ParameterList& params,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix1,
        const INPAR::S2I::InterfaceSides matrix1_side_rows,
        const INPAR::S2I::InterfaceSides matrix1_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix2,
        const INPAR::S2I::InterfaceSides matrix2_side_rows,
        const INPAR::S2I::InterfaceSides matrix2_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix3,
        const INPAR::S2I::InterfaceSides matrix3_side_rows,
        const INPAR::S2I::InterfaceSides matrix3_side_cols,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix4,
        const INPAR::S2I::InterfaceSides matrix4_side_rows,
        const INPAR::S2I::InterfaceSides matrix4_side_cols,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
        const INPAR::S2I::InterfaceSides vector1_side,
        const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
        const INPAR::S2I::InterfaceSides vector2_side) const;

    //! flag indicating if we have capacitive interface flux contributions
    bool has_capacitive_contributions_;

    //! slave-side scatra-scatra interface kinetics conditions applied to a mesh tying interface
    std::map<const int, DRT::Condition* const> kinetics_conditions_meshtying_slaveside_;

    //! corresponding master conditions to kinetics condiditions
    std::map<const int, DRT::Condition* const> master_conditions_;

    //! flag for evaluation of interface linearizations and residuals on slave side only
    bool slaveonly_;

    //! flag indicating that mesh tying for different conditions should be setup independently
    const bool indepedent_setup_of_conditions_;

  };  // class MeshtyingStrategyS2I


  class MortarCellInterface
  {
   public:
    //! Virtual destructor.
    virtual ~MortarCellInterface() = default;

    //! evaluate single mortar integration cell of particular slave-side and master-side
    //! discretization types
    virtual void Evaluate(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                                  //!< mortar integration cell
        MORTAR::Element& slaveelement,                          //!< slave-side mortar element
        MORTAR::Element& masterelement,                         //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,                  //!< slave-side location array
        DRT::Element::LocationArray& la_master,                 //!< master-side location array
        const Teuchos::ParameterList& params,                   //!< parameter list
        CORE::LINALG::SerialDenseMatrix& cellmatrix1,           //!< cell matrix 1
        CORE::LINALG::SerialDenseMatrix& cellmatrix2,           //!< cell matrix 2
        CORE::LINALG::SerialDenseMatrix& cellmatrix3,           //!< cell matrix 3
        CORE::LINALG::SerialDenseMatrix& cellmatrix4,           //!< cell matrix 4
        CORE::LINALG::SerialDenseVector& cellvector1,           //!< cell vector 1
        CORE::LINALG::SerialDenseVector& cellvector2            //!< cell vector 2
        ) = 0;

    //! evaluate single slave-side node for node-to-segment coupling
    virtual void EvaluateNTS(const DRT::Discretization& idiscret,  //!< interface discretization
        const MORTAR::Node& slavenode,                             //!< slave-side node
        const double&
            lumpedarea,  //!< lumped interface area fraction associated with slave-side node
        MORTAR::Element& slaveelement,                //!< slave-side mortar element
        MORTAR::Element& masterelement,               //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,        //!< slave-side location array
        DRT::Element::LocationArray& la_master,       //!< master-side location array
        const Teuchos::ParameterList& params,         //!< parameter list
        CORE::LINALG::SerialDenseMatrix& ntsmatrix1,  //!< node-to-segment matrix 1
        CORE::LINALG::SerialDenseMatrix& ntsmatrix2,  //!< node-to-segment matrix 2
        CORE::LINALG::SerialDenseMatrix& ntsmatrix3,  //!< node-to-segment matrix 3
        CORE::LINALG::SerialDenseMatrix& ntsmatrix4,  //!< node-to-segment matrix 4
        CORE::LINALG::SerialDenseVector& ntsvector1,  //!< node-to-segment vector 1
        CORE::LINALG::SerialDenseVector& ntsvector2   //!< node-to-segment vector 2
        ) = 0;

    //! evaluate single mortar element
    virtual void EvaluateMortarElement(
        const DRT::Discretization& idiscret,          //!< interface discretization
        MORTAR::Element& element,                     //!< mortar element
        DRT::Element::LocationArray& la,              //!< location array
        const Teuchos::ParameterList& params,         //!< parameter list
        CORE::LINALG::SerialDenseMatrix& elematrix1,  //!< element matrix 1
        CORE::LINALG::SerialDenseMatrix& elematrix2,  //!< element matrix 2
        CORE::LINALG::SerialDenseMatrix& elematrix3,  //!< element matrix 3
        CORE::LINALG::SerialDenseMatrix& elematrix4,  //!< element matrix 4
        CORE::LINALG::SerialDenseVector& elevector1,  //!< element vector 1
        CORE::LINALG::SerialDenseVector& elevector2   //!< element vector 2
        ) = 0;

   protected:
    //! protected constructor for singletons
    MortarCellInterface(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    );

    //! flag for interface side underlying Lagrange multiplier definition
    const INPAR::S2I::InterfaceSides lmside_;

    //! flag for meshtying method
    const INPAR::S2I::CouplingType couplingtype_;

    //! number of slave-side degrees of freedom per node
    const int numdofpernode_slave_;

    //! number of master-side degrees of freedom per node
    const int numdofpernode_master_;
  };


  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class MortarCellCalc : public MortarCellInterface
  {
   public:
    //! singleton access method
    static MortarCellCalc<distypeS, distypeM>* Instance(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
        const std::string& disname        //!< name of mortar discretization
    );

    //! evaluate single mortar integration cell of particular slave-side and master-side
    //! discretization types
    void Evaluate(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                          //!< mortar integration cell
        MORTAR::Element& slaveelement,                  //!< slave-side mortar element
        MORTAR::Element& masterelement,                 //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,          //!< slave-side location array
        DRT::Element::LocationArray& la_master,         //!< master-side location array
        const Teuchos::ParameterList& params,           //!< parameter list
        CORE::LINALG::SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
        CORE::LINALG::SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
        CORE::LINALG::SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
        CORE::LINALG::SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
        CORE::LINALG::SerialDenseVector& cellvector1,   //!< cell vector 1
        CORE::LINALG::SerialDenseVector& cellvector2    //!< cell vector 2
        ) override;

    //! evaluate single slave-side node for node-to-segment coupling
    void EvaluateNTS(const DRT::Discretization& idiscret,  //!< interface discretization
        const MORTAR::Node& slavenode,                     //!< slave-side node
        const double&
            lumpedarea,  //!< lumped interface area fraction associated with slave-side node
        MORTAR::Element& slaveelement,                //!< slave-side mortar element
        MORTAR::Element& masterelement,               //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,        //!< slave-side location array
        DRT::Element::LocationArray& la_master,       //!< master-side location array
        const Teuchos::ParameterList& params,         //!< parameter list
        CORE::LINALG::SerialDenseMatrix& ntsmatrix1,  //!< node-to-segment matrix 1
        CORE::LINALG::SerialDenseMatrix& ntsmatrix2,  //!< node-to-segment matrix 2
        CORE::LINALG::SerialDenseMatrix& ntsmatrix3,  //!< node-to-segment matrix 3
        CORE::LINALG::SerialDenseMatrix& ntsmatrix4,  //!< node-to-segment matrix 4
        CORE::LINALG::SerialDenseVector& ntsvector1,  //!< node-to-segment vector 1
        CORE::LINALG::SerialDenseVector& ntsvector2   //!< node-to-segment vector 2
        ) override;

    //! evaluate single mortar element
    void EvaluateMortarElement(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::Element& element,                                    //!< mortar element
        DRT::Element::LocationArray& la,                             //!< location array
        const Teuchos::ParameterList& params,                        //!< parameter list
        CORE::LINALG::SerialDenseMatrix& elematrix1,                 //!< element matrix 1
        CORE::LINALG::SerialDenseMatrix& elematrix2,                 //!< element matrix 2
        CORE::LINALG::SerialDenseMatrix& elematrix3,                 //!< element matrix 3
        CORE::LINALG::SerialDenseMatrix& elematrix4,                 //!< element matrix 4
        CORE::LINALG::SerialDenseVector& elevector1,                 //!< element vector 1
        CORE::LINALG::SerialDenseVector& elevector2                  //!< element vector 2
        ) override;

   protected:
    //! number of slave element nodes
    static const int nen_slave_ = CORE::FE::num_nodes<distypeS>;

    //! number of master element nodes
    static const int nen_master_ = CORE::FE::num_nodes<distypeM>;

    //! spatial dimensionality of slave elements
    static const int nsd_slave_ = CORE::FE::dim<distypeS>;

    //! spatial dimensionality of master elements
    static const int nsd_master_ = CORE::FE::dim<distypeM>;

    //! protected constructor for singletons
    MortarCellCalc(const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    );

    //! evaluate mortar matrices
    void EvaluateMortarMatrices(MORTAR::IntCell& cell,  //!< mortar integration cell
        MORTAR::Element& slaveelement,                  //!< slave-side mortar element
        MORTAR::Element& masterelement,                 //!< master-side mortar element
        CORE::LINALG::SerialDenseMatrix& D,             //!< mortar matrix D
        CORE::LINALG::SerialDenseMatrix& M,             //!< mortar matrix M
        CORE::LINALG::SerialDenseMatrix& E              //!< mortar matrix E
    );

    //! evaluate and assemble interface linearizations and residuals
    virtual void EvaluateCondition(
        const DRT::Discretization& idiscret,     //!< interface discretization
        MORTAR::IntCell& cell,                   //!< mortar integration cell
        MORTAR::Element& slaveelement,           //!< slave-side mortar element
        MORTAR::Element& masterelement,          //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,   //!< slave-side location array
        DRT::Element::LocationArray& la_master,  //!< master-side location array
        const Teuchos::ParameterList& params,    //!< parameter list
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseVector& r_s,  //!< slave-side residual vector
        CORE::LINALG::SerialDenseVector& r_m   //!< master-side residual vector
    );

    //! evaluate and assemble interface linearizations and residuals for node-to-segment coupling
    virtual void EvaluateConditionNTS(
        DRT::Condition& condition,      //!< scatra-scatra interface coupling condition
        const MORTAR::Node& slavenode,  //!< slave-side node
        const double&
            lumpedarea,  //!< lumped interface area fraction associated with slave-side node
        MORTAR::Element& slaveelement,   //!< slave-side mortar element
        MORTAR::Element& masterelement,  //!< master-side mortar element
        const std::vector<CORE::LINALG::Matrix<nen_slave_, 1>>&
            ephinp_slave,  //!< state variables at slave-side nodes
        const std::vector<CORE::LINALG::Matrix<nen_master_, 1>>&
            ephinp_master,  //!< state variables at master-side nodes
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseVector& r_s,  //!< slave-side residual vector
        CORE::LINALG::SerialDenseVector& r_m   //!< master-side residual vector
    );

    //! evaluate and assemble lumped interface area fractions associated with slave-side element
    //! nodes
    void EvaluateNodalAreaFractions(MORTAR::Element& slaveelement,  //!< slave-side mortar element
        CORE::LINALG::SerialDenseVector&
            areafractions  //!< lumped interface area fractions associated
                           //!< with slave-side element nodes
    );

    //! extract nodal state variables associated with mortar integration cell
    virtual void ExtractNodeValues(
        const DRT::Discretization& idiscret,    //!< interface discretization
        DRT::Element::LocationArray& la_slave,  //!< slave-side location array
        DRT::Element::LocationArray& la_master  //!< master-side location array
    );

    /*!
     * @brief  extract nodal state variables associated with slave element
     *
     * @param estate_slave  state variables at slave-side nodes
     * @param idiscret      interface discretization
     * @param la_slave      slave-side location array
     * @param statename     name of relevant state
     * @param nds          number of relevant dofset
     */
    void ExtractNodeValues(CORE::LINALG::Matrix<nen_slave_, 1>& estate_slave,
        const DRT::Discretization& idiscret, DRT::Element::LocationArray& la_slave,
        const std::string& statename = "iphinp", const int& nds = 0) const;

    /*!
     * @brief extract nodal state variables associated with slave and master elements
     *
     * @param estate_slave   state variables at slave-side nodes
     * @param estate_master  state variables at master-side nodes
     * @param idiscret       interface discretization
     * @param la_slave      slave-side location array
     * @param la_master     master-side location array
     * @param statename      name of relevant state
     * @param nds            number of relevant dofset
     */
    void ExtractNodeValues(std::vector<CORE::LINALG::Matrix<nen_slave_, 1>>& estate_slave,
        std::vector<CORE::LINALG::Matrix<nen_master_, 1>>& estate_master,
        const DRT::Discretization& idiscret, DRT::Element::LocationArray& la_slave,
        DRT::Element::LocationArray& la_master, const std::string& statename = "iphinp",
        const int& nds = 0) const;

    //! evaluate slave-side and master-side shape functions and domain integration factor at cell
    //! integration point
    double EvalShapeFuncAndDomIntFacAtIntPoint(
        MORTAR::Element& slaveelement,                               //!< slave-side mortar element
        MORTAR::Element& masterelement,                              //!< master-side mortar element
        MORTAR::IntCell& cell,                                       //!< mortar integration cell
        const CORE::FE::IntPointsAndWeights<nsd_slave_>& intpoints,  //!< quadrature rule
        const int iquad                                              //!< ID of integration point
    );

    //! evaluate slave-side shape functions and domain integration factor at element integration
    //! point
    double EvalShapeFuncAndDomIntFacAtIntPoint(MORTAR::Element& element,  //!< mortar element
        const CORE::FE::IntPointsAndWeights<nsd_slave_>& intpoints,       //!< quadrature rule
        const int iquad  //!< ID of integration point
    );

    //! evaluate shape functions at position of slave-side node
    void EvalShapeFuncAtSlaveNode(const MORTAR::Node& slavenode,  //!< slave-side node
        MORTAR::Element& slaveelement,                            //!< slave-side element
        MORTAR::Element& masterelement                            //!< master-side element
    );

    //! pointer to scatra boundary parameter list
    DRT::ELEMENTS::ScaTraEleParameterBoundary* scatraparamsboundary_;

    //! nodal, slave-side state variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<CORE::LINALG::Matrix<nen_slave_, 1>> ephinp_slave_;

    //! nodal, master-side state variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<CORE::LINALG::Matrix<nen_master_, 1>> ephinp_master_;

    //! shape and test function values associated with slave-side dofs at integration point
    CORE::LINALG::Matrix<nen_slave_, 1> funct_slave_;

    //! shape and test function values associated with master-side dofs at integration point
    CORE::LINALG::Matrix<nen_master_, 1> funct_master_;

    //! shape function values associated with slave-side Lagrange multipliers at integration point
    CORE::LINALG::Matrix<nen_slave_, 1> shape_lm_slave_;

    //! shape function values associated with master-side Lagrange multipliers at integration point
    CORE::LINALG::Matrix<nen_master_, 1> shape_lm_master_;

    //! test function values associated with slave-side Lagrange multipliers at integration point
    CORE::LINALG::Matrix<nen_slave_, 1> test_lm_slave_;

    //! test function values associated with master-side Lagrange multipliers at integration point
    CORE::LINALG::Matrix<nen_master_, 1> test_lm_master_;
  };  // class MortarCellCalc


  class MortarCellFactory
  {
   public:
    //! provide instance of mortar cell evaluation class of particular slave-side discretization
    //! type
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&
            impltype,  //!< physical implementation type of mortar integration cell
        const MORTAR::Element& slaveelement,           //!< slave-side mortar element
        const MORTAR::Element& masterelement,          //!< master-side mortar element
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const std::string& disname  //!< name of interface discretization
    );

   private:
    //! provide instance of mortar cell evaluation class of particular slave-side and master-side
    //! discretization types
    template <CORE::FE::CellType distypeS>
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&
            impltype,  //!< physical implementation type of mortar integration cell
        const MORTAR::Element& masterelement,          //!< master-side mortar element
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const std::string& disname       //!< name of interface discretization
    );

    //! provide specific instance of mortar cell evaluation class
    template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&
            impltype,  //!< physical implementation type of mortar integration cell
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
        const std::string& disname        //!< name of interface discretization
    );
  };  // class MortarCellFactory


  class MortarCellAssemblyStrategy
  {
   public:
    //! constructor
    MortarCellAssemblyStrategy(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,  //!< system matrix 1
        const INPAR::S2I::InterfaceSides
            matrix1_side_rows,  //!< interface side associated with rows of system matrix 1
        const INPAR::S2I::InterfaceSides
            matrix1_side_cols,  //!< interface side associated with columns of system matrix 1
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,  //!< system matrix 2
        const INPAR::S2I::InterfaceSides
            matrix2_side_rows,  //!< interface side associated with rows of system matrix 2
        const INPAR::S2I::InterfaceSides
            matrix2_side_cols,  //!< interface side associated with columns of system matrix 2
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix3,  //!< system matrix 3
        const INPAR::S2I::InterfaceSides
            matrix3_side_rows,  //!< interface side associated with rows of system matrix 3
        const INPAR::S2I::InterfaceSides
            matrix3_side_cols,  //!< interface side associated with columns of system matrix 3
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix4,  //!< system matrix 4
        const INPAR::S2I::InterfaceSides
            matrix4_side_rows,  //!< interface side associated with rows of system matrix 4
        const INPAR::S2I::InterfaceSides
            matrix4_side_cols,  //!< interface side associated with columns of system matrix 4
        Teuchos::RCP<Epetra_MultiVector> systemvector1,  //!< system vector 1
        const INPAR::S2I::InterfaceSides
            vector1_side,  //!< interface side associated with system vector 1
        Teuchos::RCP<Epetra_MultiVector> systemvector2,  //!< system vector 2
        const INPAR::S2I::InterfaceSides
            vector2_side,        //!< interface side associated with system vector 2
        const int nds_rows = 0,  //!< number of dofset associated with matrix rows
        const int nds_cols = 0   //!< number of dofset associated with matrix columns
    );

    /*!
     * @brief assemble cell matrices and vectors into system matrices and vectors
     *
     * @param la_slave   slave-side location array
     * @param la_master  master-side location array
     * @param assembler_pid_master  ID of processor performing master-side matrix and vector
     *                              assembly
     */
    void AssembleCellMatricesAndVectors(DRT::Element::LocationArray& la_slave,
        DRT::Element::LocationArray& la_master, const int assembler_pid_master) const;

    //! bool flag for assembly of system matrix 1
    bool AssembleMatrix1() const { return systemmatrix1_ != Teuchos::null; };

    //! bool flag for assembly of system matrix 2
    bool AssembleMatrix2() const { return systemmatrix2_ != Teuchos::null; };

    //! bool flag for assembly of system matrix 3
    bool AssembleMatrix3() const { return systemmatrix3_ != Teuchos::null; };

    //! bool flag for assembly of system matrix 4
    bool AssembleMatrix4() const { return systemmatrix4_ != Teuchos::null; };

    //! bool flag for assembly of system vector 1
    bool AssembleVector1() const { return systemvector1_ != Teuchos::null; };

    //! bool flag for assembly of system vector 2
    bool AssembleVector2() const { return systemvector2_ != Teuchos::null; };

    //! return cell matrix 1
    CORE::LINALG::SerialDenseMatrix& CellMatrix1() { return cellmatrix1_; };

    //! return cell matrix 2
    CORE::LINALG::SerialDenseMatrix& CellMatrix2() { return cellmatrix2_; };

    //! return cell matrix 3
    CORE::LINALG::SerialDenseMatrix& CellMatrix3() { return cellmatrix3_; };

    //! return cell matrix 4
    CORE::LINALG::SerialDenseMatrix& CellMatrix4() { return cellmatrix4_; };

    //! return cell vector 1
    CORE::LINALG::SerialDenseVector& CellVector1() { return cellvector1_; };

    //! return cell vector 2
    CORE::LINALG::SerialDenseVector& CellVector2() { return cellvector2_; };

    //! initialize cell matrices and vectors
    void InitCellMatricesAndVectors(
        DRT::Element::LocationArray& la_slave,  //!< slave-side location array
        DRT::Element::LocationArray& la_master  //!< master-side location array
    );

   private:
    /*!
     * @brief  assemble cell matrix into system matrix
     *
     * @param systemmatrix   system matrix
     * @param cellmatrix     cell matrix
     * @param side_rows      interface side associated with matrix rows
     * @param side_cols      interface side associated with matrix columns
     * @param la_slave       slave-side location array
     * @param la_master      master-side location array
     * @param assembler_pid_master   ID of processor performing master-side matrix assembly
     */
    void AssembleCellMatrix(const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,
        const CORE::LINALG::SerialDenseMatrix& cellmatrix,
        const INPAR::S2I::InterfaceSides side_rows, const INPAR::S2I::InterfaceSides side_cols,
        DRT::Element::LocationArray& la_slave, DRT::Element::LocationArray& la_master,
        const int assembler_pid_master) const;

    /*!
     * @brief  assemble cell vector into system vector
     *
     * @param systemvector  system vector
     * @param cellvector    cell vector
     * @param side          interface side associated with system and cell vectors
     * @param la_slave      slave-side location array
     * @param la_master     master-side location array
     * @param assembler_pid_master  ID of processor performing master-side vector assembly
     */
    void AssembleCellVector(const Teuchos::RCP<Epetra_MultiVector>& systemvector,
        const CORE::LINALG::SerialDenseVector& cellvector, const INPAR::S2I::InterfaceSides side,
        DRT::Element::LocationArray& la_slave, DRT::Element::LocationArray& la_master,
        const int assembler_pid_master) const;

    /*!
     * @brief  initialize cell matrix
     *
     * @param cellmatrix  cell matrix
     * @param side_rows   interface side associated with rows of cell matrix
     * @param side_cols   interface side associated with columns of cell matrix
     * @param la_slave    slave-side location array
     * @param la_master  master-side location array
     */
    void InitCellMatrix(CORE::LINALG::SerialDenseMatrix& cellmatrix,
        const INPAR::S2I::InterfaceSides side_rows, const INPAR::S2I::InterfaceSides side_cols,
        DRT::Element::LocationArray& la_slave, DRT::Element::LocationArray& la_master) const;

    /*!
     * @brief  initialize cell vector
     *
     * @param cellvector  cell vector
     * @param side        interface side associated with cell vector
     * @param la_slave    slave-side location array
     * @param la_master   master-side location array
     */
    void InitCellVector(CORE::LINALG::SerialDenseVector& cellvector,
        const INPAR::S2I::InterfaceSides side, DRT::Element::LocationArray& la_slave,
        DRT::Element::LocationArray& la_master) const;

    //! cell matrix 1
    CORE::LINALG::SerialDenseMatrix cellmatrix1_;

    //! cell matrix 2
    CORE::LINALG::SerialDenseMatrix cellmatrix2_;

    //! cell matrix 3
    CORE::LINALG::SerialDenseMatrix cellmatrix3_;

    //! cell matrix 4
    CORE::LINALG::SerialDenseMatrix cellmatrix4_;

    //! cell vector 1
    CORE::LINALG::SerialDenseVector cellvector1_;

    //! cell vector 2
    CORE::LINALG::SerialDenseVector cellvector2_;

    //! interface side associated with rows of system matrix 1
    const INPAR::S2I::InterfaceSides matrix1_side_rows_;

    //! interface side associated with columns of system matrix 1
    const INPAR::S2I::InterfaceSides matrix1_side_cols_;

    //! interface side associated with rows of system matrix 2
    const INPAR::S2I::InterfaceSides matrix2_side_rows_;

    //! interface side associated with columns of system matrix 2
    const INPAR::S2I::InterfaceSides matrix2_side_cols_;

    //! interface side associated with rows of system matrix 3
    const INPAR::S2I::InterfaceSides matrix3_side_rows_;

    //! interface side associated with columns of system matrix 3
    const INPAR::S2I::InterfaceSides matrix3_side_cols_;

    //! interface side associated with rows of system matrix 4
    const INPAR::S2I::InterfaceSides matrix4_side_rows_;

    //! interface side associated with columns of system matrix 4
    const INPAR::S2I::InterfaceSides matrix4_side_cols_;

    //! system matrix 1
    const Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1_;

    //! system matrix 2
    const Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2_;

    //! system matrix 3
    const Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix3_;

    //! system matrix 4
    const Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix4_;

    //! system vector 1
    const Teuchos::RCP<Epetra_MultiVector> systemvector1_;

    //! system vector 2
    const Teuchos::RCP<Epetra_MultiVector> systemvector2_;

    //! interface side associated with system vector 1
    const INPAR::S2I::InterfaceSides vector1_side_;

    //! interface side associated with system vector 2
    const INPAR::S2I::InterfaceSides vector2_side_;

    //! number of dofset associated with matrix rows
    const int nds_rows_;

    //! number of dofset associated with matrix columns
    const int nds_cols_;
  };  // class MortarCellAssembleStrategy
}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif  // SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
