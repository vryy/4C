/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i.H

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H

#include "scatra_timint_meshtying_strategy_base.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_inpar/inpar_scatra.H"

#include <Epetra_FEVector.h>
#include <Epetra_IntVector.h>

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingMortar;
}

namespace DRT
{
  class Condition;

  namespace UTILS
  {
    template<const int NSD> class IntPointsAndWeights;
  }
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixColTransform;
    class MatrixRowTransform;
    class MatrixRowColTransform;
  }
}

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
}

namespace MORTAR
{
  class IntCell;
  class MortarElement;
  class MortarNode;
}

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

\author fang (fang@lnm.mw.tum.de)
*/

class MeshtyingStrategyS2I : public MeshtyingStrategyBase
{
  public:

    //! constructor
    explicit MeshtyingStrategyS2I(
        SCATRA::ScaTraTimIntImpl*       scatratimint,   //!< scalar transport time integrator
        const Teuchos::ParameterList&   parameters      //!< input parameters for scatra-scatra interface coupling
        );

    //! destructor
    virtual ~MeshtyingStrategyS2I(){return;};

    //! provide global state vectors for element evaluation
    void AddTimeIntegrationSpecificVectors() const;

    //! compute time derivatives of discrete state variables
    void ComputeTimeDerivative() const;

    //! condense global system of equations
    void CondenseMatAndRHS(
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,              //!< system matrix
        const Teuchos::RCP<Epetra_Vector>&            residual,                  //!< residual vector
        const bool                                    calcinittimederiv = false  //!< flag for calculation of initial time derivative
        ) const;

    //! return interface coupling adapter
    const Teuchos::RCP<const ADAPTER::Coupling> CouplingAdapter() const {return icoup_;};

    //! return flag for meshtying method
    const INPAR::S2I::CouplingType& CouplingType() const { return couplingtype_;}

    //! compute meshtying residual terms and their linearizations
    virtual void EvaluateMeshtying();

    //! evaluate mortar integration cells
    void EvaluateMortarCells(
        const DRT::Discretization&            idiscret,   //!< interface discretization
        const Teuchos::ParameterList&         params,     //!< parameter list for evaluation of mortar integration cells
        SCATRA::MortarCellAssemblyStrategy&   strategy    //!< assembly strategy for mortar integration cells
        ) const;

    //! explicit predictor step to obtain better starting value for Newton-Raphson iteration
    void ExplicitPredictor() const;

    //! finite difference check for extended system matrix involving scatra-scatra interface layer growth (for debugging only)
    void FDCheck(
        const LINALG::BlockSparseMatrixBase&   extendedsystemmatrix,   //!< global system matrix
        const Teuchos::RCP<Epetra_Vector>&     extendedresidual        //!< global residual vector
        ) const;

    //! return state vector of discrete scatra-scatra interface layer thicknesses at time n
    const Teuchos::RCP<Epetra_Vector>& GrowthVarN() const {return growthn_;};

    //! return state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    const Teuchos::RCP<Epetra_Vector>& GrowthVarNp() const {return growthnp_;};

    //! perform initialization of scatra-scatra interface coupling
    void InitMeshtying();

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

    //! return interface map extractor
    const Teuchos::RCP<LINALG::MultiMapExtractor>& InterfaceMaps() const {return interfacemaps_;};

    //! return flag for evaluation of scatra-scatra interface coupling involving interface layer growth
    const INPAR::S2I::GrowthEvaluation& IntLayerGrowthEvaluation() const {return intlayergrowth_evaluation_;};

    //! return vector of Lagrange multiplier dofs
    const Teuchos::RCP<const Epetra_Vector> LM() const {return lm_;};

    //! return constraint residual vector associated with Lagrange multiplier dofs
    const Teuchos::RCP<const Epetra_Vector> LMResidual() const {return lmresidual_;};

    //! return constraint increment vector associated with Lagrange multiplier dofs
    const Teuchos::RCP<const Epetra_Vector> LMIncrement() const {return lmincrement_;};

    //! return auxiliary system matrix for linearizations of slave fluxes w.r.t. master dofs
    const Teuchos::RCP<LINALG::SparseMatrix>& MasterMatrix() const {return imastermatrix_;};

    //! return master dofs transformed to slave side of scatra-scatra interface
    const Teuchos::RCP<const Epetra_Vector> MasterPhinp() const {return imasterphinp_;};

    //! return mortar interface discretization associated with particular condition ID
    DRT::Discretization& MortarDiscretization(const int& condid) const;

    //! output solution for post-processing
    void Output() const;

    //! output restart information
    void OutputRestart() const;

    //! return mortar projector P
    const Teuchos::RCP<LINALG::SparseMatrix>& P() const {return P_;};

    //! read restart data
    void ReadRestart(
        const int                        step,                   //!< restart step
        Teuchos::RCP<IO::InputControl>   input = Teuchos::null   //!< control file manager
        ) const;

    //! set general parameters for element evaluation
    void SetElementGeneralParameters(Teuchos::ParameterList& parameters) const;

    //! compute history vector, i.e., the history part of the right-hand side vector with all contributions from the previous time step
    void SetOldPartOfRHS() const;

    //! perform setup of scatra-scatra interface coupling
    void SetupMeshtying();

    //! return auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs (non-mortar case) or slave and master dofs (mortar case)
    const Teuchos::RCP<LINALG::SparseMatrix>& SlaveMatrix() const {return islavematrix_;};

    //! solve resulting linear system of equations
    void Solve(
        const Teuchos::RCP<LINALG::Solver>&            solver,         //!< solver
        const Teuchos::RCP<LINALG::SparseOperator>&    systemmatrix,   //!< system matrix
        const Teuchos::RCP<Epetra_Vector>&             increment,      //!< increment vector
        const Teuchos::RCP<Epetra_Vector>&             residual,       //!< residual vector
        const Teuchos::RCP<Epetra_Vector>&             phinp,          //!< state vector at time n+1
        const int&                                     iteration,      //!< number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>&   projector       //!< Krylov projector
        ) const;

    //! update solution after convergence of the nonlinear Newton-Raphson iteration
    virtual void Update() const;

  protected:

    //! build maps associated with blocks of global system matrix
    virtual void BuildBlockMaps(
        const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //!< domain partitioning conditions
        std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //!< empty vector for maps to be built
        ) const;

    //! build null spaces associated with blocks of global system matrix
    virtual void BuildBlockNullSpaces() const;

    //! instantiate strategy for Newton-Raphson convergence check
    virtual void InitConvCheckStrategy();

    //! interface map extractor (0: other, 1: slave, 2: master)
    Teuchos::RCP<LINALG::MultiMapExtractor> interfacemaps_;

    //! map extractors associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_;
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_slave_;
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_master_;

    //! non-mortar interface coupling adapter
    Teuchos::RCP<ADAPTER::Coupling> icoup_;

    //! mortar interface coupling adapters
    std::map<int,Teuchos::RCP<ADAPTER::CouplingMortar> > icoupmortar_;

    //! mortar integration cells
    std::map<int,std::vector<std::pair<Teuchos::RCP<MORTAR::IntCell>,INPAR::SCATRA::ImplType> > > imortarcells_;

    //! vectors for node-to-segment connectivity, i.e., for pairings between slave nodes and master elements
    std::map<int,Teuchos::RCP<Epetra_IntVector> > islavenodestomasterelements_;

    //! vectors for physical implementation types of slave-side nodes
    std::map<int,Teuchos::RCP<Epetra_IntVector> > islavenodesimpltypes_;

    //! vectors for lumped interface area fractions associated with slave-side nodes
    std::map<int,Teuchos::RCP<Epetra_Vector> > islavenodeslumpedareas_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs (non-mortar case)
    //! or for linearizations of slave fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<LINALG::SparseMatrix> islavematrix_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. master dofs (non-mortar case)
    //! or for linearizations of master fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<LINALG::SparseMatrix> imastermatrix_;

    //! flag for meshtying method
    const INPAR::S2I::CouplingType couplingtype_;

    //! mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> D_;

    //! mortar matrix M
    Teuchos::RCP<LINALG::SparseMatrix> M_;

    //! mortar matrix E
    Teuchos::RCP<LINALG::SparseMatrix> E_;

    //! mortar projector P
    Teuchos::RCP<LINALG::SparseMatrix> P_;

    //! mortar projector Q
    Teuchos::RCP<LINALG::SparseMatrix> Q_;

    //! vector of Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lm_;

    //! extended map extractor (0: standard dofs, 1: Lagrange multiplier dofs or scatra-scatra interface layer thickness variables)
    Teuchos::RCP<LINALG::MapExtractor> extendedmaps_;

    //! constraint residual vector associated with Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lmresidual_;

    //! constraint increment vector associated with Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Vector> lmincrement_;

    //! transformation operators for auxiliary system matrices
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> imastertoslaverowtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> islavetomastercoltransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> islavetomasterrowtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> islavetomasterrowcoltransform_;

    //! auxiliary residual vector for slave residuals
    Teuchos::RCP<Epetra_Vector> islaveresidual_;

    //! auxiliary residual vector for master residuals
    Teuchos::RCP<Epetra_FEVector> imasterresidual_;

    //! master dofs transformed to slave side of scatra-scatra interface
    Teuchos::RCP<Epetra_Vector> imasterphinp_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;

    //! inverse sums of absolute values of column entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invcolsums_;

    //! flag for interface side underlying Lagrange multiplier definition
    const INPAR::S2I::InterfaceSides lmside_;

    //! type of global system matrix in global system of equations
    const INPAR::S2I::MatrixType matrixtype_;

    //! node-to-segment projection tolerance
    const double ntsprojtol_;

    //! flag for evaluation of scatra-scatra interface coupling involving interface layer growth
    const INPAR::S2I::GrowthEvaluation intlayergrowth_evaluation_;

    //! local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving interface layer growth
    const double intlayergrowth_convtol_;

    //! maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling involving interface layer growth
    const unsigned intlayergrowth_itemax_;

    //! state vector of discrete scatra-scatra interface layer thicknesses at time n
    Teuchos::RCP<Epetra_Vector> growthn_;

    //! state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    Teuchos::RCP<Epetra_Vector> growthnp_;

    //! state vector of time derivatives of discrete scatra-scatra interface layer thicknesses at time n
    Teuchos::RCP<Epetra_Vector> growthdtn_;

    //! state vector of time derivatives of discrete scatra-scatra interface layer thicknesses at time n+1
    Teuchos::RCP<Epetra_Vector> growthdtnp_;

    //! state vector of history values associated with discrete scatra-scatra interface layer thicknesses
    Teuchos::RCP<Epetra_Vector> growthhist_;

    //! state vector of residual values associated with discrete scatra-scatra interface layer thicknesses
    Teuchos::RCP<Epetra_Vector> growthresidual_;

    //! state vector of Newton-Raphson increment values associated with discrete scatra-scatra interface layer thicknesses
    Teuchos::RCP<Epetra_Vector> growthincrement_;

    //! scatra-growth block of extended global system matrix (derivatives of discrete scatra residuals w.r.t. discrete scatra-scatra interface layer thicknesses)
    Teuchos::RCP<LINALG::SparseMatrix> scatragrowthblock_;

    //! growth-scatra block of extended global system matrix (derivatives of discrete scatra-scatra interface layer growth residuals w.r.t. discrete scatra degrees of freedom)
    Teuchos::RCP<LINALG::SparseMatrix> growthscatrablock_;

    //! growth-growth block of extended global system matrix (derivatives of discrete scatra-scatra interface layer growth residuals w.r.t. discrete scatra-scatra interface layer thicknesses)
    Teuchos::RCP<LINALG::SparseMatrix> growthgrowthblock_;

  private:

    //! copy constructor
    MeshtyingStrategyS2I(const MeshtyingStrategyS2I& old);

    //! build map extractors associated with blocks of global system matrix
    void BuildBlockMapExtractors();

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,   //!< system matrix
        const Teuchos::RCP<Epetra_Vector>&            residual        //!< residual vector
        ) const;

    //! evaluate single mortar integration cell
    void EvaluateMortarCell(
        const DRT::Discretization&       idiscret,        //!< interface discretization
        MORTAR::IntCell&                 cell,            //!< mortar integration cell
        const INPAR::SCATRA::ImplType&   impltype,        //!< physical implementation type of mortar integration cell
        MORTAR::MortarElement&           slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&           masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&     la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&     la_master,       //!< master-side location array
        const Teuchos::ParameterList&    params,          //!< parameter list
        Epetra_SerialDenseMatrix&        cellmatrix1,     //!< cell matrix 1
        Epetra_SerialDenseMatrix&        cellmatrix2,     //!< cell matrix 2
        Epetra_SerialDenseMatrix&        cellmatrix3,     //!< cell matrix 3
        Epetra_SerialDenseMatrix&        cellmatrix4,     //!< cell matrix 4
        Epetra_SerialDenseVector&        cellvector1,     //!< cell vector 1
        Epetra_SerialDenseVector&        cellvector2      //!< cell vector 2
        ) const;

    //! evaluate single slave-side node for node-to-segment coupling
    void EvaluateSlaveNode(
        const DRT::Discretization&       idiscret,        //!< interface discretization
        const MORTAR::MortarNode&        slavenode,       //!< slave-side node
        const double&                    lumpedarea,      //!< lumped interface area fraction associated with slave-side node
        const INPAR::SCATRA::ImplType&   impltype,        //!< physical implementation type of mortar integration cell
        MORTAR::MortarElement&           slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&           masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&     la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&     la_master,       //!< master-side location array
        const Teuchos::ParameterList&    params,          //!< parameter list
        Epetra_SerialDenseMatrix&        ntsmatrix1,      //!< node-to-segment matrix 1
        Epetra_SerialDenseMatrix&        ntsmatrix2,      //!< node-to-segment matrix 2
        Epetra_SerialDenseMatrix&        ntsmatrix3,      //!< node-to-segment matrix 3
        Epetra_SerialDenseMatrix&        ntsmatrix4,      //!< node-to-segment matrix 4
        Epetra_SerialDenseVector&        ntsvector1,      //!< node-to-segment vector 1
        Epetra_SerialDenseVector&        ntsvector2       //!< node-to-segment vector 2
        ) const;

    //! evaluate single mortar element
    void EvaluateMortarElement(
        const DRT::Discretization&       idiscret,        //!< interface discretization
        MORTAR::MortarElement&           element,         //!< mortar element
        const INPAR::SCATRA::ImplType&   impltype,        //!< physical implementation type of mortar element
        DRT::Element::LocationArray&     la,              //!< location array
        const Teuchos::ParameterList&    params,          //!< parameter list
        Epetra_SerialDenseMatrix&        elematrix1,      //!< element matrix 1
        Epetra_SerialDenseMatrix&        elematrix2,      //!< element matrix 2
        Epetra_SerialDenseMatrix&        elematrix3,      //!< element matrix 3
        Epetra_SerialDenseMatrix&        elematrix4,      //!< element matrix 4
        Epetra_SerialDenseVector&        elevector1,      //!< element vector 1
        Epetra_SerialDenseVector&        elevector2       //!< element vector 2
        ) const;

    //! evaluate mortar integration cells
    void EvaluateMortarCells(
        const DRT::Discretization&                    idiscret,            //!< interface discretization
        const Teuchos::ParameterList&                 params,              //!< parameter list for evaluation of mortar integration cells
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix1,       //!< system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_rows,   //!< interface side associated with rows of system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_cols,   //!< interface side associated with columns of system matrix 1
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix2,       //!< system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_rows,   //!< interface side associated with rows of system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_cols,   //!< interface side associated with columns of system matrix 2
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix3,       //!< system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_rows,   //!< interface side associated with rows of system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_cols,   //!< interface side associated with columns of system matrix 3
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix4,       //!< system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_rows,   //!< interface side associated with rows of system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_cols,   //!< interface side associated with columns of system matrix 4
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector1,       //!< system vector 1
        const INPAR::S2I::InterfaceSides              vector1_side,        //!< interface side associated with system vector 1
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector2,       //!< system vector 2
        const INPAR::S2I::InterfaceSides              vector2_side         //!< interface side associated with system vector 2
        ) const;

    //! evaluate node-to-segment coupling
    void EvaluateNTS(
        const Epetra_IntVector&                       islavenodestomasterelements,   //!< vector for node-to-segment connectivity
        const Epetra_Vector&                          islavenodeslumpedareas,        //!< vector for lumped interface area fractions associated with slave-side nodes
        const Epetra_IntVector&                       islavenodesimpltypes,          //!< vector for physical implementation types of slave-side nodes
        const DRT::Discretization&                    idiscret,                      //!< interface discretization
        const Teuchos::ParameterList&                 params,                        //!< parameter list for evaluation of mortar integration cells
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix1,                 //!< system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_rows,             //!< interface side associated with rows of system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_cols,             //!< interface side associated with columns of system matrix 1
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix2,                 //!< system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_rows,             //!< interface side associated with rows of system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_cols,             //!< interface side associated with columns of system matrix 2
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix3,                 //!< system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_rows,             //!< interface side associated with rows of system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_cols,             //!< interface side associated with columns of system matrix 3
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix4,                 //!< system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_rows,             //!< interface side associated with rows of system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_cols,             //!< interface side associated with columns of system matrix 4
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector1,                 //!< system vector 1
        const INPAR::S2I::InterfaceSides              vector1_side,                  //!< interface side associated with system vector 1
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector2,                 //!< system vector 2
        const INPAR::S2I::InterfaceSides              vector2_side                   //!< interface side associated with system vector 2
        ) const;

    //! evaluate mortar elements
    void EvaluateMortarElements(
        const Epetra_Map&                             ielecolmap,          //!< column map of mortar elements
        const Epetra_IntVector&                       ieleimpltypes,       //!< vector for physical implementation types of mortar elements
        const DRT::Discretization&                    idiscret,            //!< interface discretization
        const Teuchos::ParameterList&                 params,              //!< parameter list for evaluation of mortar integration cells
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix1,       //!< system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_rows,   //!< interface side associated with rows of system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_cols,   //!< interface side associated with columns of system matrix 1
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix2,       //!< system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_rows,   //!< interface side associated with rows of system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_cols,   //!< interface side associated with columns of system matrix 2
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix3,       //!< system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_rows,   //!< interface side associated with rows of system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_cols,   //!< interface side associated with columns of system matrix 3
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix4,       //!< system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_rows,   //!< interface side associated with rows of system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_cols,   //!< interface side associated with columns of system matrix 4
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector1,       //!< system vector 1
        const INPAR::S2I::InterfaceSides              vector1_side,        //!< interface side associated with system vector 1
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector2,       //!< system vector 2
        const INPAR::S2I::InterfaceSides              vector2_side         //!< interface side associated with system vector 2
        ) const;

    //! compute inverse sums of absolute values of matrix row entries
    void ComputeInvRowSums(
        const LINALG::SparseMatrix&          matrix,      //!< matrix
        const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< inverse sums of absolute values of row entries in matrix
        ) const;

    //! compute inverse sums of absolute values of matrix column entries
    void ComputeInvColSums(
        const LINALG::SparseMatrix&          matrix,      //!< matrix
        const Teuchos::RCP<Epetra_Vector>&   invcolsums   //!< inverse sums of absolute values of column entries in matrix
        ) const;

    //! equilibrate matrix rows
    void EquilibrateMatrixRows(
        LINALG::SparseMatrix&                matrix,      //!< matrix
        const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< sums of absolute values of row entries in matrix
        ) const;

    //! equilibrate matrix columns
    void EquilibrateMatrixColumns(
        LINALG::SparseMatrix&                matrix,      //!< matrix
        const Teuchos::RCP<Epetra_Vector>&   invcolsums   //!< sums of absolute values of column entries in matrix
        ) const;

    //! unequilibrate global increment vector if necessary
    void UnequilibrateIncrement(
        const Teuchos::RCP<Epetra_Vector>&   increment   //!< increment vector
        ) const;

    //! slave-side scatra-scatra interface coupling conditions
    std::map<const int,DRT::Condition* const> slaveconditions_;

    //! flag for row equilibration of global system of equations
    const bool rowequilibration_;

    //! flag for column equilibration of global system of equations
    const bool colequilibration_;

    //! flag for evaluation of interface linearizations and residuals on slave side only
    bool slaveonly_;
}; // class MeshtyingStrategyS2I


class MortarCellInterface : public DRT::SingletonDestruction
{
  public:

    //! evaluate single mortar integration cell of particular slave-side and master-side discretization types
    virtual void Evaluate(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        MORTAR::IntCell&                cell,            //!< mortar integration cell
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&    la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&    la_master,       //!< master-side location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       cellmatrix1,     //!< cell matrix 1
        Epetra_SerialDenseMatrix&       cellmatrix2,     //!< cell matrix 2
        Epetra_SerialDenseMatrix&       cellmatrix3,     //!< cell matrix 3
        Epetra_SerialDenseMatrix&       cellmatrix4,     //!< cell matrix 4
        Epetra_SerialDenseVector&       cellvector1,     //!< cell vector 1
        Epetra_SerialDenseVector&       cellvector2      //!< cell vector 2
        ) = 0;

    //! evaluate single slave-side node for node-to-segment coupling
    virtual void EvaluateNTS(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        const MORTAR::MortarNode&       slavenode,       //!< slave-side node
        const double&                   lumpedarea,      //!< lumped interface area fraction associated with slave-side node
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&    la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&    la_master,       //!< master-side location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       ntsmatrix1,      //!< node-to-segment matrix 1
        Epetra_SerialDenseMatrix&       ntsmatrix2,      //!< node-to-segment matrix 2
        Epetra_SerialDenseMatrix&       ntsmatrix3,      //!< node-to-segment matrix 3
        Epetra_SerialDenseMatrix&       ntsmatrix4,      //!< node-to-segment matrix 4
        Epetra_SerialDenseVector&       ntsvector1,      //!< node-to-segment vector 1
        Epetra_SerialDenseVector&       ntsvector2       //!< node-to-segment vector 2
        ) = 0;

    //! evaluate single mortar element
    virtual void EvaluateMortarElement(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        MORTAR::MortarElement&          element,         //!< mortar element
        DRT::Element::LocationArray&    la,              //!< location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       elematrix1,      //!< element matrix 1
        Epetra_SerialDenseMatrix&       elematrix2,      //!< element matrix 2
        Epetra_SerialDenseMatrix&       elematrix3,      //!< element matrix 3
        Epetra_SerialDenseMatrix&       elematrix4,      //!< element matrix 4
        Epetra_SerialDenseVector&       elevector1,      //!< element vector 1
        Epetra_SerialDenseVector&       elevector2       //!< element vector 2
        ) = 0;

  protected:

    //! protected constructor for singletons
    MortarCellInterface(
        const INPAR::S2I::CouplingType&     couplingtype,          //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );

    //! protected destructor for singletons
    virtual ~MortarCellInterface(){return;};

    //! flag for interface side underlying Lagrange multiplier definition
    const INPAR::S2I::InterfaceSides lmside_;

    //! flag for meshtying method
    const INPAR::S2I::CouplingType couplingtype_;

    //! number of slave-side degrees of freedom per node
    const int numdofpernode_slave_;

    //! number of master-side degrees of freedom per node
    const int numdofpernode_master_;
};


template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
class MortarCellCalc : public MortarCellInterface
{
  public:

    //! singleton access method
    static MortarCellCalc<distypeS,distypeM>* Instance(
        const INPAR::S2I::CouplingType&     couplingtype,           //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                 //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,    //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master,   //!< number of master-side degrees of freedom per node
        const std::string&                  disname,                //!< name of mortar discretization
        const MortarCellCalc*               delete_me = NULL        //!< pointer to instance to be deleted
        );

    //! singleton destruction
    virtual void Done();

    //! evaluate single mortar integration cell of particular slave-side and master-side discretization types
    virtual void Evaluate(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        MORTAR::IntCell&                cell,            //!< mortar integration cell
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&    la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&    la_master,       //!< master-side location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       cellmatrix1,     //!< cell matrix 1
        Epetra_SerialDenseMatrix&       cellmatrix2,     //!< cell matrix 2
        Epetra_SerialDenseMatrix&       cellmatrix3,     //!< cell matrix 3
        Epetra_SerialDenseMatrix&       cellmatrix4,     //!< cell matrix 4
        Epetra_SerialDenseVector&       cellvector1,     //!< cell vector 1
        Epetra_SerialDenseVector&       cellvector2      //!< cell vector 2
        );

    //! evaluate single slave-side node for node-to-segment coupling
    void EvaluateNTS(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        const MORTAR::MortarNode&       slavenode,       //!< slave-side node
        const double&                   lumpedarea,      //!< lumped interface area fraction associated with slave-side node
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&    la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&    la_master,       //!< master-side location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       ntsmatrix1,      //!< node-to-segment matrix 1
        Epetra_SerialDenseMatrix&       ntsmatrix2,      //!< node-to-segment matrix 2
        Epetra_SerialDenseMatrix&       ntsmatrix3,      //!< node-to-segment matrix 3
        Epetra_SerialDenseMatrix&       ntsmatrix4,      //!< node-to-segment matrix 4
        Epetra_SerialDenseVector&       ntsvector1,      //!< node-to-segment vector 1
        Epetra_SerialDenseVector&       ntsvector2       //!< node-to-segment vector 2
        );

    //! evaluate single mortar element
    void EvaluateMortarElement(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        MORTAR::MortarElement&          element,         //!< mortar element
        DRT::Element::LocationArray&    la,              //!< location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       elematrix1,      //!< element matrix 1
        Epetra_SerialDenseMatrix&       elematrix2,      //!< element matrix 2
        Epetra_SerialDenseMatrix&       elematrix3,      //!< element matrix 3
        Epetra_SerialDenseMatrix&       elematrix4,      //!< element matrix 4
        Epetra_SerialDenseVector&       elevector1,      //!< element vector 1
        Epetra_SerialDenseVector&       elevector2       //!< element vector 2
        );

  protected:

    //! number of slave element nodes
    static const int nen_slave_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;

    //! number of master element nodes
    static const int nen_master_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

    //! spatial dimensionality of slave elements
    static const int nsd_slave_ = DRT::UTILS::DisTypeToDim<distypeS>::dim;

    //! spatial dimensionality of master elements
    static const int nsd_master_ = DRT::UTILS::DisTypeToDim<distypeM>::dim;

    //! protected constructor for singletons
    MortarCellCalc(
        const INPAR::S2I::CouplingType&     couplingtype,          //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );

    //! protected destructor for singletons
    virtual ~MortarCellCalc(){return;};

    //! evaluate mortar matrices
    void EvaluateMortarMatrices(
        MORTAR::IntCell&                cell,            //!< mortar integration cell
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        Epetra_SerialDenseMatrix&       D,               //!< mortar matrix D
        Epetra_SerialDenseMatrix&       M,               //!< mortar matrix M
        Epetra_SerialDenseMatrix&       E                //!< mortar matrix E
        );

    //! evaluate and assemble interface linearizations and residuals
    virtual void EvaluateCondition(
        const DRT::Discretization&      idiscret,        //!< interface discretization
        MORTAR::IntCell&                cell,            //!< mortar integration cell
        MORTAR::MortarElement&          slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&          masterelement,   //!< master-side mortar element
        DRT::Element::LocationArray&    la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&    la_master,       //!< master-side location array
        const Teuchos::ParameterList&   params,          //!< parameter list
        Epetra_SerialDenseMatrix&       k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&       k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseMatrix&       k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&       k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseVector&       r_s,             //!< slave-side residual vector
        Epetra_SerialDenseVector&       r_m              //!< master-side residual vector
        );

    //! evaluate and assemble interface linearizations and residuals for node-to-segment coupling
    virtual void EvaluateConditionNTS(
        DRT::Condition&                                      condition,       //!< scatra-scatra interface coupling condition
        const MORTAR::MortarNode&                            slavenode,       //!< slave-side node
        const double&                                        lumpedarea,      //!< lumped interface area fraction associated with slave-side node
        MORTAR::MortarElement&                               slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&                               masterelement,   //!< master-side mortar element
        const std::vector<LINALG::Matrix<nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
        const std::vector<LINALG::Matrix<nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
        Epetra_SerialDenseMatrix&                            k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                            k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseMatrix&                            k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                            k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseVector&                            r_s,             //!< slave-side residual vector
        Epetra_SerialDenseVector&                            r_m              //!< master-side residual vector
        );

    //! evaluate and assemble lumped interface area fractions associated with slave-side element nodes
    void EvaluateNodalAreaFractions(
        MORTAR::MortarElement&      slaveelement,   //!< slave-side mortar element
        Epetra_SerialDenseVector&   areafractions   //!< lumped interface area fractions associated with slave-side element nodes
        );

    //! extract nodal state variables associated with mortar integration cell
    virtual void ExtractNodeValues(
        const DRT::Discretization&     idiscret,   //!< interface discretization
        DRT::Element::LocationArray&   la_slave,   //!< slave-side location array
        DRT::Element::LocationArray&   la_master   //!< master-side location array
        );

    //! extract nodal state variables associated with slave element
    void ExtractNodeValues(
        LINALG::Matrix<nen_slave_,1>&   estate_slave,           //!< state variables at slave-side nodes
        const DRT::Discretization&      idiscret,               //!< interface discretization
        DRT::Element::LocationArray&    la_slave,               //!< slave-side location array
        const std::string&              statename = "iphinp",   //!< name of relevant state
        const int&                      nds = 0                 //!< number of relevant dofset
        ) const;

    //! extract nodal state variables associated with slave and master elements
    void ExtractNodeValues(
        std::vector<LINALG::Matrix<nen_slave_,1> >&    estate_slave,           //!< state variables at slave-side nodes
        std::vector<LINALG::Matrix<nen_master_,1> >&   estate_master,          //!< state variables at master-side nodes
        const DRT::Discretization&                     idiscret,               //!< interface discretization
        DRT::Element::LocationArray&                   la_slave,               //!< slave-side location array
        DRT::Element::LocationArray&                   la_master,              //!< master-side location array
        const std::string&                             statename = "iphinp",   //!< name of relevant state
        const int&                                     nds = 0                 //!< number of relevant dofset
        ) const;

    //! evaluate slave-side and master-side shape functions and domain integration factor at cell integration point
    double EvalShapeFuncAndDomIntFacAtIntPoint(
        MORTAR::MortarElement&                               slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&                               masterelement,   //!< master-side mortar element
        MORTAR::IntCell&                                     cell,            //!< mortar integration cell
        const DRT::UTILS::IntPointsAndWeights<nsd_slave_>&   intpoints,       //!< quadrature rule
        const int                                            iquad            //!< ID of integration point
        );

    //! evaluate slave-side shape functions and domain integration factor at element integration point
    double EvalShapeFuncAndDomIntFacAtIntPoint(
        MORTAR::MortarElement&                               element,     //!< mortar element
        const DRT::UTILS::IntPointsAndWeights<nsd_slave_>&   intpoints,   //!< quadrature rule
        const int                                            iquad        //!< ID of integration point
        );

    //! evaluate shape functions at position of slave-side node
    void EvalShapeFuncAtSlaveNode(
        const MORTAR::MortarNode&   slavenode,      //!< slave-side node
        MORTAR::MortarElement&      slaveelement,   //!< slave-side element
        MORTAR::MortarElement&      masterelement   //!< master-side element
        );

    //! nodal, slave-side state variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<LINALG::Matrix<nen_slave_,1> > ephinp_slave_;

    //! nodal, master-side state variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<LINALG::Matrix<nen_master_,1> > ephinp_master_;

    //! shape and test function values associated with slave-side dofs at integration point
    LINALG::Matrix<nen_slave_,1> funct_slave_;

    //! shape and test function values associated with master-side dofs at integration point
    LINALG::Matrix<nen_master_,1> funct_master_;

    //! shape function values associated with slave-side Lagrange multipliers at integration point
    LINALG::Matrix<nen_slave_,1> shape_lm_slave_;

    //! shape function values associated with master-side Lagrange multipliers at integration point
    LINALG::Matrix<nen_master_,1> shape_lm_master_;

    //! test function values associated with slave-side Lagrange multipliers at integration point
    LINALG::Matrix<nen_slave_,1> test_lm_slave_;

    //! test function values associated with master-side Lagrange multipliers at integration point
    LINALG::Matrix<nen_master_,1> test_lm_master_;
}; // class MortarCellCalc


class MortarCellFactory
{
  public:

    //! constructor
    MortarCellFactory(){return;};

    //! destructor
    virtual ~MortarCellFactory(){return;};

    //! provide instance of mortar cell evaluation class of particular slave-side discretization type
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&      impltype,        //!< physical implementation type of mortar integration cell
        const MORTAR::MortarElement&        slaveelement,    //!< slave-side mortar element
        const MORTAR::MortarElement&        masterelement,   //!< master-side mortar element
        const INPAR::S2I::CouplingType&     couplingtype,    //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,          //!< flag for interface side underlying Lagrange multiplier definition
        const std::string&                  disname          //!< name of interface discretization
        );

  private:

    //! provide instance of mortar cell evaluation class of particular slave-side and master-side discretization types
    template<DRT::Element::DiscretizationType distypeS>
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&      impltype,              //!< physical implementation type of mortar integration cell
        const MORTAR::MortarElement&        masterelement,         //!< master-side mortar element
        const INPAR::S2I::CouplingType&     couplingtype,          //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const std::string&                  disname                //!< name of interface discretization
        );

    //! provide specific instance of mortar cell evaluation class
    template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
    static MortarCellInterface* MortarCellCalc(
        const INPAR::SCATRA::ImplType&      impltype,               //!< physical implementation type of mortar integration cell
        const INPAR::S2I::CouplingType&     couplingtype,           //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                 //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,    //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master,   //!< number of master-side degrees of freedom per node
        const std::string&                  disname                 //!< name of interface discretization
        );
}; // class MortarCellFactory


class MortarCellAssemblyStrategy
{
  public:

    //! constructor
    MortarCellAssemblyStrategy(
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix1,       //!< system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_rows,   //!< interface side associated with rows of system matrix 1
        const INPAR::S2I::InterfaceSides              matrix1_side_cols,   //!< interface side associated with columns of system matrix 1
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix2,       //!< system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_rows,   //!< interface side associated with rows of system matrix 2
        const INPAR::S2I::InterfaceSides              matrix2_side_cols,   //!< interface side associated with columns of system matrix 2
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix3,       //!< system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_rows,   //!< interface side associated with rows of system matrix 3
        const INPAR::S2I::InterfaceSides              matrix3_side_cols,   //!< interface side associated with columns of system matrix 3
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix4,       //!< system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_rows,   //!< interface side associated with rows of system matrix 4
        const INPAR::S2I::InterfaceSides              matrix4_side_cols,   //!< interface side associated with columns of system matrix 4
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector1,       //!< system vector 1
        const INPAR::S2I::InterfaceSides              vector1_side,        //!< interface side associated with system vector 1
        const Teuchos::RCP<Epetra_MultiVector>&       systemvector2,       //!< system vector 2
        const INPAR::S2I::InterfaceSides              vector2_side,        //!< interface side associated with system vector 2
        const unsigned                                nds_rows = 0,        //!< number of dofset associated with matrix rows
        const unsigned                                nds_cols = 0         //!< number of dofset associated with matrix columns
        );

    //! destructor
    virtual ~MortarCellAssemblyStrategy(){return;};

    //! assemble cell matrices and vectors into system matrices and vectors
    void AssembleCellMatricesAndVectors(
        DRT::Element::LocationArray&   la_slave,              //!< slave-side location array
        DRT::Element::LocationArray&   la_master,             //!< master-side location array
        const int                      assembler_pid_master   //!< ID of processor performing master-side matrix and vector assembly
        ) const;

    //! bool flag for assembly of system matrix 1
    bool AssembleMatrix1() const {return systemmatrix1_ != Teuchos::null;};

    //! bool flag for assembly of system matrix 2
    bool AssembleMatrix2() const {return systemmatrix2_ != Teuchos::null;};

    //! bool flag for assembly of system matrix 3
    bool AssembleMatrix3() const {return systemmatrix3_ != Teuchos::null;};

    //! bool flag for assembly of system matrix 4
    bool AssembleMatrix4() const {return systemmatrix4_ != Teuchos::null;};

    //! bool flag for assembly of system vector 1
    bool AssembleVector1() const {return systemvector1_ != Teuchos::null;};

    //! bool flag for assembly of system vector 2
    bool AssembleVector2() const {return systemvector2_ != Teuchos::null;};

    //! return cell matrix 1
    Epetra_SerialDenseMatrix& CellMatrix1(){return cellmatrix1_;};

    //! return cell matrix 2
    Epetra_SerialDenseMatrix& CellMatrix2(){return cellmatrix2_;};

    //! return cell matrix 3
    Epetra_SerialDenseMatrix& CellMatrix3(){return cellmatrix3_;};

    //! return cell matrix 4
    Epetra_SerialDenseMatrix& CellMatrix4(){return cellmatrix4_;};

    //! return cell vector 1
    Epetra_SerialDenseVector& CellVector1(){return cellvector1_;};

    //! return cell vector 2
    Epetra_SerialDenseVector& CellVector2(){return cellvector2_;};

    //! initialize cell matrices and vectors
    void InitCellMatricesAndVectors(
        DRT::Element::LocationArray&   la_slave,   //!< slave-side location array
        DRT::Element::LocationArray&   la_master   //!< master-side location array
        );

    //! return system matrix 1
    const Teuchos::RCP<LINALG::SparseOperator>& SystemMatrix1() const {return systemmatrix1_;};

    //! return system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator>& SystemMatrix2() const {return systemmatrix2_;};

    //! return system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator>& SystemMatrix3() const {return systemmatrix3_;};

    //! return system matrix 4
    const Teuchos::RCP<LINALG::SparseOperator>& SystemMatrix4() const {return systemmatrix4_;};

    //! return system vector 1
    const Teuchos::RCP<Epetra_MultiVector>& SystemVector1() const {return systemvector1_;};

    //! return system vector 2
    const Teuchos::RCP<Epetra_MultiVector>& SystemVector2() const {return systemvector2_;};

  private:

    //! assemble cell matrix into system matrix
    void AssembleCellMatrix(
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,          //!< system matrix
        const Epetra_SerialDenseMatrix&               cellmatrix,            //!< cell matrix
        const INPAR::S2I::InterfaceSides              side_rows,             //!< interface side associated with matrix rows
        const INPAR::S2I::InterfaceSides              side_cols,             //!< interface side associated with matrix columns
        DRT::Element::LocationArray&                  la_slave,              //!< slave-side location array
        DRT::Element::LocationArray&                  la_master,             //!< master-side location array
        const int                                     assembler_pid_master   //!< ID of processor performing master-side matrix assembly
        ) const;

    //! assemble cell vector into system vector
    void AssembleCellVector(
        const Teuchos::RCP<Epetra_MultiVector>&   systemvector,          //!< system vector
        const Epetra_SerialDenseVector&           cellvector,            //!< cell vector
        const INPAR::S2I::InterfaceSides          side,                  //!< interface side associated with system and cell vectors
        DRT::Element::LocationArray&              la_slave,              //!< slave-side location array
        DRT::Element::LocationArray&              la_master,             //!< master-side location array
        const int                                 assembler_pid_master   //!< ID of processor performing master-side vector assembly
        ) const;

    //! initialize cell matrix
    void InitCellMatrix(
        Epetra_SerialDenseMatrix&          cellmatrix,   //!< cell matrix
        const INPAR::S2I::InterfaceSides   side_rows,    //!< interface side associated with rows of cell matrix
        const INPAR::S2I::InterfaceSides   side_cols,    //!< interface side associated with columns of cell matrix
        DRT::Element::LocationArray&       la_slave,     //!< slave-side location array
        DRT::Element::LocationArray&       la_master     //!< master-side location array
        ) const;

    //! initialize cell vector
    void InitCellVector(
        Epetra_SerialDenseVector&          cellvector,   //!< cell vector
        const INPAR::S2I::InterfaceSides   side,         //!< interface side associated with cell vector
        DRT::Element::LocationArray&       la_slave,     //!< slave-side location array
        DRT::Element::LocationArray&       la_master     //!< master-side location array
        ) const;

    //! cell matrix 1
    Epetra_SerialDenseMatrix cellmatrix1_;

    //! cell matrix 2
    Epetra_SerialDenseMatrix cellmatrix2_;

    //! cell matrix 3
    Epetra_SerialDenseMatrix cellmatrix3_;

    //! cell matrix 4
    Epetra_SerialDenseMatrix cellmatrix4_;

    //! cell vector 1
    Epetra_SerialDenseVector cellvector1_;

    //! cell vector 2
    Epetra_SerialDenseVector cellvector2_;

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
    const Teuchos::RCP<LINALG::SparseOperator> systemmatrix1_;

    //! system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator> systemmatrix2_;

    //! system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator> systemmatrix3_;

    //! system matrix 4
    const Teuchos::RCP<LINALG::SparseOperator> systemmatrix4_;

    //! system vector 1
    const Teuchos::RCP<Epetra_MultiVector> systemvector1_;

    //! system vector 2
    const Teuchos::RCP<Epetra_MultiVector> systemvector2_;

    //! interface side associated with system vector 1
    const INPAR::S2I::InterfaceSides vector1_side_;

    //! interface side associated with system vector 2
    const INPAR::S2I::InterfaceSides vector2_side_;

    //! number of dofset associated with matrix rows
    const unsigned nds_rows_;

    //! number of dofset associated with matrix columns
    const unsigned nds_cols_;
}; // class MortarCellAssembleStrategy
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
