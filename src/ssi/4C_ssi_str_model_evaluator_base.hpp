/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for scalar-structure interaction

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SSI_STR_MODEL_EVALUATOR_BASE_HPP
#define FOUR_C_SSI_STR_MODEL_EVALUATOR_BASE_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR::MODELEVALUATOR
{
  class BaseSSI : public Generic
  {
   public:
    bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override { return true; }

    bool AssembleJacobian(
        CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override
    {
      return true;
    }

    void DetermineEnergy() override {}

    void determine_optional_quantity() override {}

    void determine_stress_strain() override;

    bool EvaluateForce() override { return true; }

    bool EvaluateForceStiff() override { return true; }

    bool EvaluateStiff() override { return true; }

    [[nodiscard]] Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override
    {
      FOUR_C_THROW("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override
    {
      FOUR_C_THROW("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_mechanical_stress_state() const override
    {
      return mechanical_stress_state_;
    }

    void OutputStepState(IO::DiscretizationWriter& iowriter) const override {}

    void PostEvaluate() override {}

    void PostOutput() override {}

    void Predict(const INPAR::STR::PredEnum& pred_type) override {}

    void PreEvaluate() override {}

    void ReadRestart(IO::DiscretizationReader& ioreader) override {}

    void Reset(const Epetra_Vector& x) override {}

    void ResetStepState() override {}

    void RunPostComputeX(
        const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
    {
    }

    void RunPostIterate(const ::NOX::Solver::Generic& solver) override {}

    void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
        const NOX::NLN::Group& curr_grp) override
    {
    }

    void Setup() override;

    [[nodiscard]] INPAR::STR::ModelType Type() const override
    {
      return INPAR::STR::model_basic_coupling;
    }

    void UpdateStepElement() override {}

    void UpdateStepState(const double& timefac_n) override {}

    void WriteRestart(
        IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override
    {
    }

   private:
    //! mechanical stress state
    Teuchos::RCP<Epetra_Vector> mechanical_stress_state_;
  };
}  // namespace STR::MODELEVALUATOR
FOUR_C_NAMESPACE_CLOSE

#endif