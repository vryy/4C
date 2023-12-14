/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of constraint terms.


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef BACI_CONSTRAINT_FRAMEWORK_MODEL_EVALUATOR_H
#define BACI_CONSTRAINT_FRAMEWORK_MODEL_EVALUATOR_H


#include "baci_constraint_framework_submodelevaluator_base.H"
#include "baci_inpar_constraint_framework.H"
#include "baci_structure_new_model_evaluator_generic.H"

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace MODELEVALUATOR
  {
    /**
     * \brief This class serves as a model evaluator for different types
     * of constraints applied to structural degrees of freedom.
     *
     * Through submodel evaluators, this class evaluates and assembles
     * the contributions resulting from periodic displacement boundary
     * conditions and coupling terms resulting from the constraint
     * enforcement for embedded mesh methods. The implementation of
     * these applications and their submodel evaluators is still a
     * work in progress.
     */
    class Constraints : public Generic
    {
     public:
      using Vector = std::vector<Teuchos::RCP<CONSTRAINTS::SUBMODELEVALUATOR::Constraint>>;

      /**
       * \brief Constructor.
       */
      Constraints() = default;

      /*! \brief Setup of the model evaluator and submodel evaluator
       *
       */
      void Setup() override;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_constraints; }

      void Reset(const Epetra_Vector& x) override;

      bool EvaluateForce() override;

      bool EvaluateStiff() override;

      bool EvaluateForceStiff() override;

      void PreEvaluate() override {}

      void PostEvaluate() override {}

      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      bool AssembleJacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
      }

      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
      }

      void RunPostIterate(const ::NOX::Solver::Generic& solver) override {}

      void Predict(const INPAR::STR::PredEnum& pred_type) override;

      void UpdateStepState(const double& timefac_n) override;

      void UpdateStepElement() override;

      void DetermineStressStrain() override;

      void DetermineEnergy() override;

      void DetermineOptionalQuantity() override;

      void ResetStepState() override;

      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      void RuntimePreOutputStepState() override;

      void RuntimeOutputStepState() const override;

      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      void PostOutput() override;

      void EvaluateJacobianContributionsFromElementLevelForPTC() override;

      void AssembleJacobianContributionsFromElementLevelForPTC(
          Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n) override;

      void CreateBackupState(const Epetra_Vector& dir) override;

      void RecoverFromBackupState() override;

      //! @}

     private:
      //! Vector of submodel evaluators
      Teuchos::RCP<STR::MODELEVALUATOR::Constraints::Vector> me_vec_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr_;

      //! Enforcement strategy
      INPAR::CONSTRAINTS::Strategy strategy_ = INPAR::CONSTRAINTS::Strategy::penalty_regularization;
    };
  }  // namespace MODELEVALUATOR
}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif