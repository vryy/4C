/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*--------------------------------------------------------------------------*/

#ifndef SSTI_MONOLITHIC_H
#define SSTI_MONOLITHIC_H

#include "ssti_algorithm.H"

#include <Teuchos_Time.hpp>

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace CORE::LINALG
{
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}  // namespace CORE::LINALG

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
  class SSTIMatrices;
  class ThermoStructureOffDiagCoupling;

  //! equilibration methods applied to system matrix
  struct SSTIMonoEquilibrationMethod
  {
    const CORE::LINALG::EquilibrationMethod global;     //! unique equilibration
    const CORE::LINALG::EquilibrationMethod scatra;     //! equilibration for scatra block
    const CORE::LINALG::EquilibrationMethod structure;  //! equilibration for structure block
    const CORE::LINALG::EquilibrationMethod thermo;     //! equilibration for thermo block
  };

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

    //! get vector containing positions within system matrix for specific subproblem
    std::vector<int> GetBlockPositions(Subproblem subproblem) const;

    //! get position within global dof map for specific subproblem
    int GetProblemPosition(Subproblem subproblem) const;

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
    unsigned int NewtonIteration() const { return Iter(); };

    //! state vectors
    //@{
    const Teuchos::RCP<Epetra_Vector> Increment() const { return increment_; };
    const Teuchos::RCP<Epetra_Vector> Residual() const { return residual_; };
    //}

    //! statistics for evaluation and solving
    std::vector<double> TimeStatistics() const
    {
      return {dtevaluate_ + dtassemble_, dtsolve_, dtnewton_};
    };

   private:
    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! build null spaces associated with blocks of global system matrix
    void BuildNullSpaces();

    //! Get Matrix and Right-Hand-Side for all subproblems incl. coupling
    void EvaluateSubproblems();

    //! get solution increment for given subproblem
    Teuchos::RCP<Epetra_Vector> ExtractSubIncrement(Subproblem sub);

    // build and return vector of equilibration methods for each block of system matrix
    std::vector<CORE::LINALG::EquilibrationMethod> GetBlockEquilibration();

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
    Teuchos::RCP<Epetra_Vector> residual_;
    Teuchos::RCP<CORE::LINALG::Solver> solver_;
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
    Teuchos::RCP<Teuchos::Time> timer_;
    //@}

    //! control parameters
    //@{
    //! equilibration method applied to system matrix
    const struct SSTIMonoEquilibrationMethod equilibration_method_;
    const CORE::LINALG::MatrixType matrixtype_;
    //@}

    //! convergence check of Newton iteration
    Teuchos::RCP<SSTI::ConvCheckMono> convcheck_;

    //! all maps
    Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono_;

    //! system matrix and submatrices
    Teuchos::RCP<SSTI::SSTIMatrices> ssti_matrices_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSTI::AssembleStrategyBase> strategy_assemble_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<CORE::LINALG::Equilibration> strategy_equilibration_;
  };
}  // namespace SSTI
#endif
