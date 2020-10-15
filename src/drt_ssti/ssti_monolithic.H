/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2


*/
/*--------------------------------------------------------------------------*/

#ifndef SSTI_MONOLITHIC_H_
#define SSTI_MONOLITHIC_H_

#include "ssti_algorithm.H"

#include <Epetra_Time.h>

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace LINALG
{
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
}

namespace STI
{
  class ScatraThermoOffDiagCoupling;
}

namespace SSI
{
  class ScatraStructureOffDiagCoupling;
}  // namespace SSI

namespace SSTI
{
  class AssembleStrategyBase;
  class ConvCheckMono;
  class SSTIMapsMono;
  class ThermoStructureOffDiagCoupling;

  enum class Subproblem
  {
    structure,
    scalar_transport,
    thermo
  };

  class SSTIMono : public SSTIAlgorithm
  {
   public:
    explicit SSTIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    ~SSTIMono() override = default;

    //! Setup of algorithm
    //@{
    void Init(const Epetra_Comm& comm, const Teuchos::ParameterList& sstitimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams,
        const Teuchos::ParameterList& structparams) override;
    void Setup() override;
    void SetupSystem() override;
    //@}

    //! Loop over all time steps
    void Timeloop() override;

    //! return all maps
    const Teuchos::RCP<SSTI::SSTIMapsMono> AllMaps() const { return ssti_maps_mono_; };

    //! number of current Newton Iteration
    unsigned int NewtonIteration() const { return iter_; };

    //! state vectors
    //@{
    const Teuchos::RCP<Epetra_Vector> Increment() const { return increment_; };
    const Teuchos::RCP<Epetra_Vector> Residual() const { return residual_; };
    //}

    //! statistics for evaluation and solving
    std::vector<double> Statistics() const
    {
      return {dtevaluate_ + dtassemble_, dtsolve_, dtnewton_};
    };

   private:
    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global
    //! system matrix
    //! build null spaces associated with blocks of global system matrix
    void BuildNullSpaces();

    //! Get Matrix and Right-Hand-Side for all subproblems incl. coupling
    void EvaluateSubproblems();

    //! get solution increment for given subproblem
    Teuchos::RCP<Epetra_Vector> ExtractSubIncrement(Subproblem sub);

    //! evaluate time step using Newton-Raphson iteration
    void NewtonLoop();

    //! output solution to screen and files
    void Output() override;

    void PrepareNewtonStep();

    //! prepare time step
    void PrepareTimeStep() override;

    //! solve linear system of equations
    void LinearSolve();

    //! update scalar transport and structure fields after time step evaluation
    void Update() override;

    //! update routine after newton iteration
    void UpdateIterStates();

    //! Newton Raphson loop
    //@{
    Teuchos::RCP<Epetra_Vector> increment_;
    unsigned int itermax_;
    double itertol_;
    Teuchos::RCP<Epetra_Vector> residual_;
    double restol_;
    Teuchos::RCP<LINALG::Solver> solver_;
    //@}

    //! subblocks of system matrix
    //@{
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;
    Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain_;
    Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface_;
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface_;
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface_;
    //@}

    //! evaluation of off-diagonal blocks
    //@{
    Teuchos::RCP<SSI::ScatraStructureOffDiagCoupling> scatrastructureoffdiagcoupling_;
    Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> scatrathermooffdiagcoupling_;
    Teuchos::RCP<SSTI::ThermoStructureOffDiagCoupling> thermostructureoffdiagcoupling_;
    //@}

    //! time monitor
    //@{
    double dtassemble_;
    double dtevaluate_;
    double dtnewton_;
    double dtsolve_;
    Teuchos::RCP<Epetra_Time> timer_;
    //@}

    //! control parameters
    //@{
    const LINALG::EquilibrationMethod equilibration_method_;
    const LINALG::MatrixType matrixtype_;
    //@}

    //! all maps
    Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSTI::AssembleStrategyBase> strategy_assemble_;

    //! convergence check of Newton iteration
    Teuchos::RCP<SSTI::ConvCheckMono> convcheck_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<LINALG::Equilibration> equilibration_;
  };
}  // namespace SSTI
#endif
