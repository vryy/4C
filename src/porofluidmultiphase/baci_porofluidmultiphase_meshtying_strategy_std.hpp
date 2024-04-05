/*----------------------------------------------------------------------*/
/*! \file
\brief standard case without mesh tying

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_HPP

#include "baci_config.hpp"

#include "baci_porofluidmultiphase_meshtying_strategy_base.hpp"

BACI_NAMESPACE_OPEN

namespace POROFLUIDMULTIPHASE
{
  class MeshtyingStrategyStd : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyStd(POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams);


    //! prepare time loop
    void PrepareTimeLoop() override;

    //! prepare time step
    void PrepareTimeStep() override;

    //! update
    void Update() override;

    //! output
    void Output() override;

    //! Initialize the linear solver
    void InitializeLinearSolver(Teuchos::RCP<CORE::LINALG::Solver> solver) override;

    //! solve linear system of equations
    void LinearSolve(Teuchos::RCP<CORE::LINALG::Solver> solver,
        Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector> increment,
        Teuchos::RCP<Epetra_Vector> residual, CORE::LINALG::SolverParams& solver_params) override;

    //! calculate norms for convergence checks
    void CalculateNorms(std::vector<double>& preresnorm, std::vector<double>& incprenorm,
        std::vector<double>& prenorm, const Teuchos::RCP<const Epetra_Vector> increment) override;

    //! create the field test
    void CreateFieldTest() override;

    //! restart
    void ReadRestart(const int step) override;

    //! evaluate mesh tying
    void Evaluate() override;

    //! extract increments and update mesh tying
    Teuchos::RCP<const Epetra_Vector> ExtractAndUpdateIter(
        const Teuchos::RCP<const Epetra_Vector> inc) override;

    //! access to global (combined) increment of coupled problem
    Teuchos::RCP<const Epetra_Vector> CombinedIncrement(
        const Teuchos::RCP<const Epetra_Vector> inc) const override;

    //! check if initial fields on coupled DOFs are equal
    void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont) const override;

    //! set the element pairs that are close as found by search algorithm
    void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs) override;

    //! setup the strategy
    void Setup() override;

    //! apply the mesh movement
    void ApplyMeshMovement() const override;
  };

}  // namespace POROFLUIDMULTIPHASE



BACI_NAMESPACE_CLOSE

#endif  // POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_H
