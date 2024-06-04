/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of constraint terms.


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_MODEL_EVALUATOR_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_MODEL_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_submodelevaluator_base.hpp"
#include "4C_inpar_constraint_framework.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  class RveMultiPointConstraintManager;

}
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
      using SubmodelevaluatorVector =
          std::vector<Teuchos::RCP<CONSTRAINTS::SUBMODELEVALUATOR::ConstraintBase>>;

      Constraints() = default;

      /*! \brief Setup of the model evaluator and submodel evaluator
       *
       */
      void Setup() override;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_constraints; }

      void Reset(const Epetra_Vector& x) override;

      bool evaluate_force() override;

      bool evaluate_stiff() override;

      bool evaluate_force_stiff() override;

      void pre_evaluate() override;

      void post_evaluate() override {}

      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      void write_restart(
          CORE::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      void read_restart(CORE::IO::DiscretizationReader& ioreader) override;

      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
      }

      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
      }

      void run_post_iterate(const ::NOX::Solver::Generic& solver) override {}

      void Predict(const INPAR::STR::PredEnum& pred_type) override;

      void UpdateStepState(const double& timefac_n) override;

      void UpdateStepElement() override;

      void determine_stress_strain() override;

      void DetermineEnergy() override;

      void determine_optional_quantity() override;

      void ResetStepState() override;

      void OutputStepState(CORE::IO::DiscretizationWriter& iowriter) const override;

      void runtime_pre_output_step_state() override;

      void runtime_output_step_state() const override;

      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      void PostOutput() override;

      void evaluate_jacobian_contributions_from_element_level_for_ptc() override;

      void assemble_jacobian_contributions_from_element_level_for_ptc(
          Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n) override;

      void CreateBackupState(const Epetra_Vector& dir) override;

      void recover_from_backup_state() override;
      //! @}

     private:
      //!@name routines for submodel management
      //! @{

      //! Set Submodeltypes depending on input
      void set_sub_model_types();

      //! build, init and setup submodel evaluator
      void create_sub_model_evaluators();

      //! @}

     private:
      //!@name data for submodel management
      //! @{
      /// active model types for the model evaluator
      std::set<enum INPAR::CONSTRAINTS::SubModelType> submodeltypes_;

      //! vector of submodelevaluators
      SubmodelevaluatorVector sub_model_vec_ptr_;

      //! constraint stiffness matrix
      Teuchos::RCP<CORE::LINALG::SparseMatrix> constraint_stiff_ptr_;

      //! constraint force vector
      Teuchos::RCP<Epetra_Vector> constraint_force_ptr_;
      //! @}
    };
  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
