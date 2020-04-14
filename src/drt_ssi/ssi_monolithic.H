/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

\maintainer Christoph Schmidt

*/
/*--------------------------------------------------------------------------*/

#ifndef SSI_MONOLITHIC_H_
#define SSI_MONOLITHIC_H_

#include "ssi_base.H"

#include "../drt_inpar/inpar_s2i.H"

#include "ssi_monolithic_assemble_strategy.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace LINALG
{
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
}

namespace SSI
{
  class SSI_Mono : public SSI_Base
  {
   public:
    //! constructor
    explicit SSI_Mono(const Epetra_Comm& comm,          //!< communicator
        const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    );

    //! destructor
    virtual ~SSI_Mono() { return; };

    //! initialize monolithic algorithm
    int Init(const Epetra_Comm& comm,                    //!< communicator
        const Teuchos::ParameterList& globaltimeparams,  //!< parameter list for time integration
        const Teuchos::ParameterList& scatraparams,      //!< parameter list for scalar transport
        const Teuchos::ParameterList& structparams,      //!< parameter list for structure
        const std::string struct_disname,                //!< name of structural discretization
        const std::string scatra_disname,  //!< name of scalar transport discretization
        bool isAle                         //!< flag for ALE
    );

    //! return counter for Newton-Raphson iterations
    const unsigned& Iter() const { return iter_; };

    //! return global map extractor (0: scalar transport, 1: structure)
    const Teuchos::RCP<LINALG::MapExtractor>& Maps() const { return maps_; };

    //! retrun map extractor associated with all degrees of freedom inside structural field
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapStructure() const
    {
      return map_structure_;
    };

    //! return meshtying strategy
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> S2IMeshtyingStrategy() const
    {
      return strategy_scatra_;
    };

    //! setup monolithic algorithm
    void Setup();

    //! setup global system of equations
    void SetupSystem();

    //! return algebraic solver for global system of equations
    const LINALG::Solver& Solver() const { return *solver_; };

    //! time loop
    void Timeloop();

   private:
    // strategies for Newton-Raphson convergence check
    class ConvCheckStrategyBase;
    class ConvCheckStrategyStd;
    class ConvCheckStrategyElch;

    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void EvaluateODBlockScatraStructureDomain();

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void EvaluateODBlockScatraStructureInterface();

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    void EvaluateODBlockStructureScatraDomain() const;

    //! build null spaces associated with blocks of global system matrix
    void BuildNullSpaces() const;

    //! compute inverse sums of absolute values of matrix row entries
    void ComputeInvRowSums(const LINALG::SparseMatrix& matrix,  //!< matrix
        const Teuchos::RCP<Epetra_Vector>&
            invrowsums  //!< inverse sums of absolute values of row entries in matrix
        ) const;

    //! return global map of degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

    //! equilibrate matrix rows
    void EquilibrateMatrixRows(LINALG::SparseMatrix& matrix,  //!< matrix
        const Teuchos::RCP<Epetra_Vector>&
            invrowsums  //!< sums of absolute values of row entries in matrix
        ) const;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& residual                //!< residual vector
        ) const;

    //! finite difference check for global system matrix
    void FDCheck();

    //! output solution to screen and files
    void Output();

    //! prepare time step
    void PrepareTimeStep();

    //! set up structural model evaluator for scalar-structure interaction
    void SetupModelEvaluator() const;

    //! evaluate time step using Newton-Raphson iteration
    void Solve();

    //! update scalar transport and structure fields after time step evaluation
    void Update();

    //! time for element evaluation and assembly of global system of equations
    double dtele_;

    //! time for solution of global system of equations
    double dtsolve_;

    //! global increment vector for Newton-Raphson iteration
    Teuchos::RCP<Epetra_Vector> increment_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const LINALG::MultiMapExtractor> map_structure_;

    //! global map extractor (0: scalar transport, 1: structure)
    Teuchos::RCP<LINALG::MapExtractor> maps_;

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_systemmatrix_;

    //! equilibration method applied to system matrix
    const INPAR::SSI::EquilibrationMethod equilibration_method_;

    //! type of global system matrix in global system of equations
    const INPAR::SSI::MatrixType matrixtype_;

    //! type of scalar transport system matrix
    INPAR::S2I::MatrixType matrixtype_scatra_;

    //! global residual vector on right-hand side of global system of equations
    Teuchos::RCP<Epetra_Vector> residual_;

    //! scatra-structure block of global system matrix (derivatives of scalar transport residuals
    //! w.r.t. structural degrees of freedom) containing domain contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;

    //! scatra-structure block of global system matrix (derivatives of scalar transport residuals
    //! w.r.t. structural degrees of freedom) containing interface contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;

    //! algebraic solver for global system of equations
    Teuchos::RCP<LINALG::Solver> solver_;

    //! strategy for Newton-Raphson convergence check
    Teuchos::RCP<const SSI::SSI_Mono::ConvCheckStrategyBase> strategy_convcheck_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSI::AssembleStrategyBase> strategy_assemble_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> strategy_scatra_;

    //! structure-scatra block of global system matrix (derivatives of structural residuals w.r.t.
    //! scalar transport degrees of freedom)
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;

    //! global system matrix
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;

    //! timer for Newton-Raphson iteration
    Teuchos::RCP<Epetra_Time> timer_;
  };
}  // namespace SSI
#endif
