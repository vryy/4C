/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i.H

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

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

#include <Epetra_FEVector.h>

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
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
}

namespace MORTAR
{
  class IntCell;
  class MortarElement;
}

namespace SCATRA
{

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

    //! condense global system of equations
    void CondenseMatAndRHS(
        const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,              //!< system matrix
        const Teuchos::RCP<Epetra_Vector>&            residual,                  //!< residual vector
        const bool                                    calcinittimederiv = false  //!< flag for calculation of initial time derivative
        ) const;

    //! return interface coupling adapter
    const Teuchos::RCP<const ADAPTER::Coupling> CouplingAdapter() const {return icoup_;};

    //! compute meshtying residual terms and their linearizations
    virtual void EvaluateMeshtying();

    //! initialize meshtying objects
    void InitMeshtying();

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

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

    //! return auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs
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

  protected:

    //! build maps associated with blocks of global system matrix
    virtual void BuildBlockMaps(
        const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //!< domain partitioning conditions
        std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //!< empty vector for maps to be built
        ) const;

    //! build null spaces associated with blocks of global system matrix
    virtual void BuildBlockNullSpaces() const;

    //! evaluate single mortar integration cell
    void EvaluateMortarCell(
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix1,          //!< interface matrix 1
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix2,          //!< interface matrix 2
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix3,          //!< interface matrix 3
        Teuchos::RCP<Epetra_Vector>          islaveresidual,    //!< slave-side residual vector
        Teuchos::RCP<Epetra_FEVector>        imasterresidual,   //!< master-side residual vector
        MORTAR::IntCell&                     cell,              //!< mortar integration cell
        const Teuchos::ParameterList&        params,            //!< parameter list
        const DRT::Discretization&           idiscret           //!< interface discretization
        ) const;

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
    std::map<int,std::vector<Teuchos::RCP<MORTAR::IntCell> > > imortarcells_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. slave dofs (non-mortar case)
    //! or for linearizations of slave fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<LINALG::SparseMatrix> islavematrix_;

    //! auxiliary system matrix for linearizations of slave fluxes w.r.t. master dofs (non-mortar case)
    //! or for linearizations of master fluxes w.r.t. slave and master dofs (mortar case)
    Teuchos::RCP<LINALG::SparseMatrix> imastermatrix_;

    //! flag for meshtying method
    const INPAR::S2I::MortarType mortartype_;

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

    //! extended map extractor (0: standard dofs, 1: Lagrange multiplier dofs)
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
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix1,          //!< interface matrix 1
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix2,          //!< interface matrix 2
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix3,          //!< interface matrix 3
        Teuchos::RCP<Epetra_Vector>          islaveresidual,    //!< slave-side residual vector
        Teuchos::RCP<Epetra_FEVector>        imasterresidual,   //!< master-side residual vector
        const DRT::Discretization&           idiscret,          //!< interface discretization
        const Teuchos::ParameterList&        params,            //!< parameter list
        MORTAR::MortarElement&               slaveelement,      //!< slave-side mortar element
        MORTAR::MortarElement&               masterelement,     //!< master-side mortar element
        MORTAR::IntCell&                     cell,              //!< mortar integration cell
        DRT::Element::LocationArray&         la_slave,          //!< slave-side location array
        DRT::Element::LocationArray&         la_master          //!< master-side location array
        ) = 0;

  protected:

    //! protected constructor for singletons
    MortarCellInterface(
        const INPAR::S2I::MortarType&       mortartype,            //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );

    //! protected destructor for singletons
    virtual ~MortarCellInterface(){return;};

    //! flag for interface side underlying Lagrange multiplier definition
    const INPAR::S2I::InterfaceSides lmside_;

    //! flag for meshtying method
    const INPAR::S2I::MortarType mortartype_;

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
        const INPAR::S2I::MortarType&       mortartype,             //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                 //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,    //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master,   //!< number of master-side degrees of freedom per node
        bool                                create = true           //!< creation flag
        );

    //! singleton destruction
    virtual void Done();

    //! evaluate single mortar integration cell of particular slave-side and master-side discretization types
    void Evaluate(
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix1,          //!< interface matrix 1
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix2,          //!< interface matrix 2
        Teuchos::RCP<LINALG::SparseMatrix>   imatrix3,          //!< interface matrix 3
        Teuchos::RCP<Epetra_Vector>          islaveresidual,    //!< slave-side residual vector
        Teuchos::RCP<Epetra_FEVector>        imasterresidual,   //!< master-side residual vector
        const DRT::Discretization&           idiscret,          //!< interface discretization
        const Teuchos::ParameterList&        params,            //!< parameter list
        MORTAR::MortarElement&               slaveelement,      //!< slave-side mortar element
        MORTAR::MortarElement&               masterelement,     //!< master-side mortar element
        MORTAR::IntCell&                     cell,              //!< mortar integration cell
        DRT::Element::LocationArray&         la_slave,          //!< slave-side location array
        DRT::Element::LocationArray&         la_master          //!< master-side location array
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
        const INPAR::S2I::MortarType&       mortartype,            //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );

    //! protected destructor for singletons
    virtual ~MortarCellCalc(){return;};

    //! evaluate mortar matrices
    void EvaluateMortarMatrices(
        Teuchos::RCP<LINALG::SparseMatrix>   D,               //!< mortar matrix D
        Teuchos::RCP<LINALG::SparseMatrix>   M,               //!< mortar matrix M
        Teuchos::RCP<LINALG::SparseMatrix>   E,               //!< mortar matrix E
        const DRT::Discretization&           idiscret,        //!< interface discretization
        MORTAR::MortarElement&               slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&               masterelement,   //!< master-side mortar element
        MORTAR::IntCell&                     cell,            //!< mortar integration cell
        DRT::Element::LocationArray&         la_slave,        //!< slave-side location array
        DRT::Element::LocationArray&         la_master        //!< master-side location array
        );

    //! evaluate and assemble interface linearizations and residuals
    virtual void EvaluateCondition(
        Teuchos::RCP<LINALG::SparseMatrix>             islavematrix,      //!< linearizations of slave-side residuals
        Teuchos::RCP<LINALG::SparseMatrix>             imastermatrix,     //!< linearizations of master-side residuals
        Teuchos::RCP<Epetra_Vector>                    islaveresidual,    //!< slave-side residual vector
        Teuchos::RCP<Epetra_FEVector>                  imasterresidual,   //!< master-side residual vector
        const DRT::Discretization&                     idiscret,          //!< interface discretization
        DRT::Condition&                                condition,         //!< scatra-scatra interface coupling condition
        MORTAR::MortarElement&                         slaveelement,      //!< slave-side mortar element
        MORTAR::MortarElement&                         masterelement,     //!< master-side mortar element
        MORTAR::IntCell&                               cell,              //!< mortar integration cell
        std::vector<LINALG::Matrix<nen_slave_,1> >&    ephinp_slave,      //!< state variables at slave-side nodes
        std::vector<LINALG::Matrix<nen_master_,1> >&   ephinp_master,     //!< state variables at master-side nodes
        DRT::Element::LocationArray&                   la_slave,          //!< slave-side location array
        DRT::Element::LocationArray&                   la_master          //!< master-side location array
        );

    //! extract nodal state variables associated with slave and master elements
    void ExtractNodeValues(
        std::vector<LINALG::Matrix<nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
        std::vector<LINALG::Matrix<nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
        const DRT::Discretization&                     idiscret,         //!< interface discretization
        DRT::Element::LocationArray&                   la_slave,         //!< slave-side location array
        DRT::Element::LocationArray&                   la_master         //!< master-side location array
        ) const;

    //! evaluate shape functions and domain integration factor at integration point
    double EvalShapeFuncAndDomIntFacAtIntPoint(
        MORTAR::MortarElement&                               slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&                               masterelement,   //!< master-side mortar element
        MORTAR::IntCell&                                     cell,            //!< mortar integration cell
        const DRT::UTILS::IntPointsAndWeights<nsd_slave_>&   intpoints,       //!< quadrature rule
        const int                                            iquad            //!< ID of integration point
        );

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
        const MORTAR::MortarElement&        slaveelement,    //!< slave-side mortar element
        const MORTAR::MortarElement&        masterelement,   //!< master-side mortar element
        const MeshtyingStrategyS2I* const   strategy,        //!< scatra-scatra interface coupling strategy
        const INPAR::S2I::MortarType&       mortartype,      //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside           //!< flag for interface side underlying Lagrange multiplier definition
        );

  private:

    //! provide instance of mortar cell evaluation class of particular slave-side and master-side discretization types
    template<DRT::Element::DiscretizationType distypeS>
    static MortarCellInterface* MortarCellCalc(
        const MORTAR::MortarElement&        masterelement,        //!< master-side mortar element
        const MeshtyingStrategyS2I* const   strategy,             //!< scatra-scatra interface coupling strategy
        const INPAR::S2I::MortarType&       mortartype,           //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,               //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave   //!< number of slave-side degrees of freedom per node
        );

    //! provide specific instance of mortar cell evaluation class
    template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
    static MortarCellInterface* MortarCellCalc(
        const MeshtyingStrategyS2I* const   strategy,              //!< scatra-scatra interface coupling strategy
        const INPAR::S2I::MortarType&       mortartype,            //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );
}; // class MortarCellFactory
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
