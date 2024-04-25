/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_MONOLITHIC_HPP
#define FOUR_C_SSI_MONOLITHIC_HPP

#include "4C_config.hpp"

#include "4C_ssi_base.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CONTACT
{
  class NitscheStrategySsi;
}

namespace INPAR::SSI
{
  enum class ScaTraTimIntType;
}  // namespace INPAR::SSI

namespace CORE::LINALG
{
  class Solver;
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
}  // namespace CORE::LINALG

namespace SSI
{
  namespace UTILS
  {
    class SSIMaps;
    class SSIMatrices;
    class SSIVectors;
  }  // namespace UTILS

  class AssembleStrategyBase;
  class ContactStrategyBase;
  class DBCHandlerBase;
  class ManifoldMeshTyingStrategyBase;
  class MeshtyingStrategyBase;
  class ScatraStructureOffDiagCoupling;
  class ScaTraManifoldScaTraFluxEvaluator;

  //! equilibration methods applied to system matrix
  struct SSIMonoEquilibrationMethod
  {
    const CORE::LINALG::EquilibrationMethod global;     //! unique equilibration
    const CORE::LINALG::EquilibrationMethod scatra;     //! equilibration for scatra block
    const CORE::LINALG::EquilibrationMethod structure;  //! equilibration for structure block
  };

  enum class Subproblem : int
  {
    scalar_transport,
    structure,
    manifold
  };

  class SSIMono : public SSIBase
  {
   public:
    //! constructor
    explicit SSIMono(const Epetra_Comm& comm,           //!< communicator
        const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    );

    //! return global map of degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

    void Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, bool isAle) override;

    //! return contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::NitscheStrategySsi> NitscheStrategySsi() const
    {
      return contact_strategy_nitsche_;
    }

    //! return global map extractor (0: scalar transport, 1: structure, [2: scatra manifold])
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> MapsSubProblems() const;

    //! return map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapScaTra() const;

    //! return map extractor associated with all degrees of freedom inside scatra manifold field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapScaTraManifold() const;

    //! return map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapStructure() const;

    //! return map extractor associated with blocks of global system matrix
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapSystemMatrix() const;

    //! Return matrix type of global system matrix
    CORE::LINALG::MatrixType MatrixType() const { return matrixtype_; };

    void ReadRestart(int restart) override;

    void Setup() override;

    void SetupSystem() override;

    /*!
     * @brief solves the linear system
     *
     * @note in case an equilibration method (scaling of rows and columns) is defined this is also
     * performed within this call
     */
    void SolveLinearSystem();

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

    //! return algebraic solver for global system of equations
    const CORE::LINALG::Solver& Solver() const { return *solver_; };

    void Timeloop() override;

   private:
    //! strategies for Newton-Raphson convergence check
    class ConvCheckStrategyBase;
    class ConvCheckStrategyElch;
    class ConvCheckStrategyElchScaTraManifold;
    class ConvCheckStrategyStd;

    //! apply the contact contributions to matrices and residuals of the sub problems
    void ApplyContactToSubProblems();

    //! apply the Dirichlet boundary conditions to the ssi system, i.e. matrices and residuals
    void ApplyDBCToSystem();

    //! apply mesh tying between manifold domains on matrices and residuals
    void ApplyManifoldMeshtying();

    //! perform mesh tying on matrices and residuals as obtained from sub problems
    void ApplyMeshtyingToSubProblems();

    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! assemble linearization of scatra residuals to system matrix
    void AssembleMatScaTra();

    //! assemble linearization of scatra on manifold residuals to system matrix
    void AssembleMatScaTraManifold();

    //! assemble linearization of structural residuals to system matrix
    void AssembleMatStructure();

    //! build null spaces associated with blocks of global system matrix
    void BuildNullSpaces() const;

    //! calc initial potential field for monolithic SSI problem including scatra and scatra manifold
    //! fields
    void CalcInitialPotentialField();

    //! calc initial time derivative of transported scalars for monolithic SSI problem including
    //! scatra and scatra manifold fields
    void CalcInitialTimeDerivative();

    //! call complete on the sub problem matrices
    void CompleteSubproblemMatrices();

    //! distribute solution to all other fields
    //! \param restore_velocity   restore velocity when StructureField()->SetState() is called
    void DistributeSolutionAllFields(bool restore_velocity = false);

    //! evaluate all off-diagonal matrix contributions
    void EvaluateOffDiagContributions();

    //! Evaluate ScaTra including copy to corresponding ssi matrix
    void EvaluateScaTra();

    //! Evaluate ScaTra on manifold incl. coupling with scatra
    void EvaluateScaTraManifold();

    //! get matrix and right-hand-side for all subproblems incl. coupling
    void EvaluateSubproblems();

    //! build and return vector of equilibration methods for each block of system matrix
    std::vector<CORE::LINALG::EquilibrationMethod> GetBlockEquilibration();

    /*!
     * @note This is only necessary in the first iteration of the simulation, since only there the
     * graph of the matrix changes
     *
     * @return flag indicating if we need to uncomplete the matrices before adding the mesh tying
     * contributions.
     */
    bool IsUncompleteOfMatricesNecessaryForMeshTying() const;

    void Output() override;

    //! do everything, that has to be done once before first time step
    void PrepareTimeLoop();

    void PrepareTimeStep() override;

    //! prepare output for subproblems if needed
    void PrepareOutput();

    //! print system matrix, rhs, and map of system matrix to file
    void PrintSystemMatrixRHSToMatLabFormat();

    //! print time step size, time, and number of time step
    void PrintTimeStepInfo();

    //! set up a pointer to the contact strategy of the structural field and store it
    void SetupContactStrategy();

    //! set scatra manifold solution on scatra field
    void SetScatraManifoldSolution(Teuchos::RCP<const Epetra_Vector> phi);

    void SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const override;

    /*!
     * @brief set contact states needed for evaluation of ssi contact
     *
     * @param[in] phi  scatra state to be set to contact nitsche strategy
     */
    void SetSSIContactStates(Teuchos::RCP<const Epetra_Vector> phi) const;

    //! evaluate time step using Newton-Raphson iteration
    void NewtonLoop();

    void Update() override;

    //! update ScaTra state within Newton iteration
    void UpdateIterScaTra();

    //! update structure state within Newton iteration
    void UpdateIterStructure();

    //! store contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_strategy_nitsche_;

    //! Dirichlet boundary condition handler
    Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler_;

    //! time for element evaluation and assembly of global system of equations
    double dt_eval_ = 0.0;

    //! time for solution of global system of equations
    double dt_solve_ = 0.0;

    //! equilibration method applied to system matrix
    const struct SSIMonoEquilibrationMethod equilibration_method_;

    //! Evaluation of coupling flux between scatra and manifold on scatra
    Teuchos::RCP<SSI::ScaTraManifoldScaTraFluxEvaluator> manifoldscatraflux_;

    //! type of global system matrix in global system of equations
    const CORE::LINALG::MatrixType matrixtype_;

    //! print system matrix, rhs, and map of system matrix to file
    const bool print_matlab_;

    //! relax the tolerance of the linear solver in case it is an iterative solver by scaling the
    //! convergence tolerance with factor @p relax_lin_solver_tolerance_
    const double relax_lin_solver_tolerance_;

    //! relax the tolerance of the linear solver within the first @p relax_lin_solver_step_ steps
    const int relax_lin_solver_iter_step_;

    //! all OD evaluation is in here
    Teuchos::RCP<SSI::ScatraStructureOffDiagCoupling> scatrastructure_off_diagcoupling_;

    //! algebraic solver for global system of equations
    Teuchos::RCP<CORE::LINALG::Solver> solver_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps_;

    //! this object holds the system matrix and all sub blocks
    Teuchos::RCP<SSI::UTILS::SSIMatrices> ssi_matrices_;

    //! this object holds the system residuals and increment
    Teuchos::RCP<SSI::UTILS::SSIVectors> ssi_vectors_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSI::AssembleStrategyBase> strategy_assemble_;

    //! strategy how to apply contact contributions to sub matrices and rhs
    Teuchos::RCP<SSI::ContactStrategyBase> strategy_contact_;

    //! strategy for Newton-Raphson convergence check
    Teuchos::RCP<SSI::SSIMono::ConvCheckStrategyBase> strategy_convcheck_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<CORE::LINALG::Equilibration> strategy_equilibration_;

    //! strategy how to apply mesh tying on manifold domains
    Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> strategy_manifold_meshtying_;

    //! strategy how to apply mesh tying to system matrix and rhs
    Teuchos::RCP<SSI::MeshtyingStrategyBase> strategy_meshtying_;

    //! timer for Newton-Raphson iteration
    Teuchos::RCP<Teuchos::Time> timer_;
  };
}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
