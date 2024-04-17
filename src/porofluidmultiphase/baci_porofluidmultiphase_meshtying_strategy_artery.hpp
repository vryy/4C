/*----------------------------------------------------------------------*/
/*! \file
\brief routines for coupling with artery network

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_ARTERY_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_ARTERY_HPP

#include "baci_config.hpp"

#include "baci_porofluidmultiphase_meshtying_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScaTraArtCouplBase;
}

namespace POROFLUIDMULTIPHASE
{
  class MeshtyingStrategyArtery : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyArtery(POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint,
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

    // return arterial network time integrator
    Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt() override { return artnettimint_; }

    //! access dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    //! right-hand side alias the dynamic force residual for coupled system
    Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const override;

    //! access to block system matrix of artery poro problem
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const override;

    //! get global (combined) increment of coupled problem
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

    //! return blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction() override;

   protected:
    //! artery time integration
    Teuchos::RCP<ADAPTER::ArtNet> artnettimint_;

    //! artery discretization
    Teuchos::RCP<DRT::Discretization> arterydis_;

    //! the mesh tying object
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase> arttoporofluidcoupling_;

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> comb_systemmatrix_;

    //! global rhs
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! global increment
    Teuchos::RCP<Epetra_Vector> comb_increment_;

    //! global solution at time n+1
    Teuchos::RCP<Epetra_Vector> comb_phinp_;
  };

}  // namespace POROFLUIDMULTIPHASE



FOUR_C_NAMESPACE_CLOSE

#endif
