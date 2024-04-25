/*----------------------------------------------------------------------*/
/*! \file
\brief Abstract interface for meshtying strategies in porofluidmultiphase problems
       (including standard strategy without meshtying)

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_BASE_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_BASE_HPP

#include "4C_config.hpp"

#include "4C_inpar_porofluidmultiphase.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class ArtNet;
}

namespace CORE::LINALG
{
  struct SolverParams;
}

namespace POROFLUIDMULTIPHASE
{
  class MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBase(POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams)
        : porofluidmultitimint_(porofluidmultitimint),
          params_(probparams),
          poroparams_(poroparams),
          vectornormfres_(CORE::UTILS::IntegralValue<INPAR::POROFLUIDMULTIPHASE::VectorNorm>(
              poroparams_, "VECTORNORM_RESF")),
          vectornorminc_(CORE::UTILS::IntegralValue<INPAR::POROFLUIDMULTIPHASE::VectorNorm>(
              poroparams_, "VECTORNORM_INC"))
    {
      return;
    }

    //! destructor
    virtual ~MeshtyingStrategyBase() = default;

    //! prepare time loop
    virtual void PrepareTimeLoop() = 0;

    //! prepare time step
    virtual void PrepareTimeStep() = 0;

    //! update
    virtual void Update() = 0;

    //! output
    virtual void Output() = 0;

    //! Initialize the linear solver
    virtual void InitializeLinearSolver(Teuchos::RCP<CORE::LINALG::Solver> solver) = 0;

    //! solve linear system of equations
    virtual void LinearSolve(Teuchos::RCP<CORE::LINALG::Solver> solver,
        Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector> increment,
        Teuchos::RCP<Epetra_Vector> residual, CORE::LINALG::SolverParams& solver_params) = 0;

    //! calculate norms for convergence checks
    virtual void CalculateNorms(std::vector<double>& preresnorm, std::vector<double>& incprenorm,
        std::vector<double>& prenorm, const Teuchos::RCP<const Epetra_Vector> increment) = 0;

    //! create the field test
    virtual void CreateFieldTest() = 0;

    //! restart
    virtual void ReadRestart(const int step) = 0;

    //! evaluate mesh tying
    virtual void Evaluate() = 0;

    //! extract increments and update mesh tying
    virtual Teuchos::RCP<const Epetra_Vector> ExtractAndUpdateIter(
        const Teuchos::RCP<const Epetra_Vector> inc) = 0;

    // return arterial network time integrator
    virtual Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt()
    {
      FOUR_C_THROW("ArtNetTimInt() not implemented in base class, wrong mesh tying object?");
      return Teuchos::null;
    }

    //! access dof row map
    virtual Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const
    {
      FOUR_C_THROW("ArteryDofRowMap() not implemented in base class, wrong mesh tying object?");
      return Teuchos::null;
    }

    //! access to block system matrix of artery poro problem
    virtual Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const
    {
      FOUR_C_THROW(
          "ArteryPorofluidSysmat() not implemented in base class, wrong mesh tying object?");
      return Teuchos::null;
    }

    //! right-hand side alias the dynamic force residual for coupled system
    virtual Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const
    {
      FOUR_C_THROW("ArteryPorofluidRHS() not implemented in base class, wrong mesh tying object?");
      return Teuchos::null;
    }

    //! access to global (combined) increment of coupled problem
    virtual Teuchos::RCP<const Epetra_Vector> CombinedIncrement(
        Teuchos::RCP<const Epetra_Vector> inc) const = 0;

    //! check if initial fields on coupled DOFs are equal (only for node-based coupling)
    virtual void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont) const = 0;

    //! set the element pairs that are close as found by search algorithm
    virtual void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs) = 0;

    //! setup the strategy
    virtual void Setup() = 0;

    //! apply the mesh movement
    virtual void ApplyMeshMovement() const = 0;

    //! return blood vessel volume fraction
    virtual Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction()
    {
      FOUR_C_THROW(
          "BloodVesselVolumeFraction() not implemented in base class, wrong mesh tying object?");
      return Teuchos::null;
    }

   protected:
    //! porofluid multi time integrator
    POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint_;

    //! parameter list of global control problem
    const Teuchos::ParameterList& params_;

    //! parameter list of poro fluid multiphase problem
    const Teuchos::ParameterList& poroparams_;

    // vector norm for residuals
    enum INPAR::POROFLUIDMULTIPHASE::VectorNorm vectornormfres_;

    // vector norm for increments
    enum INPAR::POROFLUIDMULTIPHASE::VectorNorm vectornorminc_;
  };

}  // namespace POROFLUIDMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif
