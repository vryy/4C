/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2


*/
/*--------------------------------------------------------------------------*/

#ifndef SSI_MONOLITHIC_H_
#define SSI_MONOLITHIC_H_

#include "ssi_base.H"

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
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
}

namespace SSI
{
  class ScatraStructureOffDiagCoupling;
  class AssembleStrategyBase;

  class SSI_Mono : public SSI_Base
  {
   public:
    //! constructor
    explicit SSI_Mono(const Epetra_Comm& comm,          //!< communicator
        const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    );

    ~SSI_Mono() override = default;

    int Init(const Epetra_Comm& comm,                    //!< communicator
        const Teuchos::ParameterList& globaltimeparams,  //!< parameter list for time integration
        const Teuchos::ParameterList& scatraparams,      //!< parameter list for scalar transport
        const Teuchos::ParameterList& structparams,      //!< parameter list for structure
        const std::string struct_disname,                //!< name of structural discretization
        const std::string scatra_disname,  //!< name of scalar transport discretization
        bool isAle                         //!< flag for ALE
        ) override;

    //! return global map extractor (0: scalar transport, 1: structure)
    const Teuchos::RCP<LINALG::MapExtractor>& MapsSubProblems() const { return maps_sub_problems_; }

    //! return map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsScatra() const { return maps_scatra_; }

    //! return map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapStructure() const { return map_structure_; }

    //! return map extractor associated with blocks of global system matrix
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsSystemMatrix() const
    {
      return maps_systemmatrix_;
    }

    //! return the scatra-scatra interface meshtying strategy
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingStrategyS2I() const
    {
      return meshtying_strategy_s2i_;
    }

    void Setup() override;

    void SetupSystem() override;

    //! return algebraic solver for global system of equations
    const LINALG::Solver& Solver() const { return *solver_; };

    void Timeloop() override;

   private:
    // strategies for Newton-Raphson convergence check
    class ConvCheckStrategyBase;
    class ConvCheckStrategyStd;
    class ConvCheckStrategyElch;

    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! build null spaces associated with blocks of global system matrix
    void BuildNullSpaces() const;

    //! return global map of degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

    //! finite difference check for global system matrix
    void FDCheck();

    void Output() override;

    void PrepareTimeStep() override;

    void SetupModelEvaluator() const override;

    //! evaluate time step using Newton-Raphson iteration
    void NewtonLoop();

    void Update() override;

    //! update ScaTra state within Newton iteration
    void UpdateIterScaTra();

    //! update structure state within Newton iteration
    void UpdateIterStructure();

    //! time for element evaluation and assembly of global system of equations
    double dtele_;

    //! time for solution of global system of equations
    double dtsolve_;

    //! global increment vector for Newton-Raphson iteration
    Teuchos::RCP<Epetra_Vector> increment_;

    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const LINALG::MultiMapExtractor> map_structure_;

    //! global map extractor (0: scalar transport, 1: structure)
    Teuchos::RCP<LINALG::MapExtractor> maps_sub_problems_;

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_systemmatrix_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<LINALG::Equilibration> equilibration_;

    //! equilibration method applied to system matrix
    const LINALG::EquilibrationMethod equilibration_method_;

    //! type of global system matrix in global system of equations
    const LINALG::MatrixType matrixtype_;

    //! global residual vector on right-hand side of global system of equations
    Teuchos::RCP<Epetra_Vector> residual_;

    //! scatra-structure block of global system matrix (derivatives of scalar transport residuals
    //! w.r.t. structural degrees of freedom) containing domain contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;

    //! scatra-structure block of global system matrix (derivatives of scalar transport residuals
    //! w.r.t. structural degrees of freedom) containing interface contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;

    //! all OD evaluation is in here
    Teuchos::RCP<SSI::ScatraStructureOffDiagCoupling> scatrastructureOffDiagcoupling_;

    //! algebraic solver for global system of equations
    Teuchos::RCP<LINALG::Solver> solver_;

    //! strategy for Newton-Raphson convergence check
    Teuchos::RCP<const SSI::SSI_Mono::ConvCheckStrategyBase> strategy_convcheck_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSI::AssembleStrategyBase> strategy_assemble_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

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
