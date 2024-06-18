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

namespace Adapter
{
  class ArtNet;
}

namespace Core::LinAlg
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
          vectornormfres_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::VectorNorm>(
              poroparams_, "VECTORNORM_RESF")),
          vectornorminc_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::VectorNorm>(
              poroparams_, "VECTORNORM_INC"))
    {
      return;
    }

    //! destructor
    virtual ~MeshtyingStrategyBase() = default;

    //! prepare time loop
    virtual void prepare_time_loop() = 0;

    //! prepare time step
    virtual void prepare_time_step() = 0;

    //! update
    virtual void Update() = 0;

    //! output
    virtual void Output() = 0;

    //! Initialize the linear solver
    virtual void initialize_linear_solver(Teuchos::RCP<Core::LinAlg::Solver> solver) = 0;

    //! solve linear system of equations
    virtual void linear_solve(Teuchos::RCP<Core::LinAlg::Solver> solver,
        Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector> increment,
        Teuchos::RCP<Epetra_Vector> residual, Core::LinAlg::SolverParams& solver_params) = 0;

    //! calculate norms for convergence checks
    virtual void CalculateNorms(std::vector<double>& preresnorm, std::vector<double>& incprenorm,
        std::vector<double>& prenorm, const Teuchos::RCP<const Epetra_Vector> increment) = 0;

    //! create the field test
    virtual void CreateFieldTest() = 0;

    //! restart
    virtual void read_restart(const int step) = 0;

    //! evaluate mesh tying
    virtual void evaluate() = 0;

    //! extract increments and update mesh tying
    virtual Teuchos::RCP<const Epetra_Vector> extract_and_update_iter(
        const Teuchos::RCP<const Epetra_Vector> inc) = 0;

    // return arterial network time integrator
    virtual Teuchos::RCP<Adapter::ArtNet> ArtNetTimInt()
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
    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat() const
    {
      FOUR_C_THROW(
          "artery_porofluid_sysmat() not implemented in base class, wrong mesh tying object?");
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
    virtual Teuchos::RCP<const Epetra_Vector> blood_vessel_volume_fraction()
    {
      FOUR_C_THROW(
          "blood_vessel_volume_fraction() not implemented in base class, wrong mesh tying object?");
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
    enum Inpar::POROFLUIDMULTIPHASE::VectorNorm vectornormfres_;

    // vector norm for increments
    enum Inpar::POROFLUIDMULTIPHASE::VectorNorm vectornorminc_;
  };

}  // namespace POROFLUIDMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif
