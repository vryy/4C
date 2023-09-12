/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for scalar-structure interaction

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef BACI_SSI_STR_MODEL_EVALUATOR_BASE_H
#define BACI_SSI_STR_MODEL_EVALUATOR_BASE_H

#include "baci_structure_new_model_evaluator_generic.H"

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

    void DetermineOptionalQuantity() override {}

    void DetermineStressStrain() override;

    bool EvaluateForce() override { return true; }

    bool EvaluateForceStiff() override { return true; }

    bool EvaluateStiff() override { return true; }

    [[nodiscard]] Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override
    {
      dserror("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override
    {
      dserror("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> GetMechanicalStressState() const override
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

    void RunPostIterate(const NOX::Solver::Generic& solver) override {}

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
#endif