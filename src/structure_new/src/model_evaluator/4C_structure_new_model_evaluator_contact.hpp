/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all contact terms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_CONTACT_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_CONTACT_HPP

#include "4C_config.hpp"

#include "4C_structure_new_enum_lists.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CONTACT
{
  class Manager;
  class AbstractStrategy;
}  // namespace CONTACT

namespace MORTAR
{
  class StrategyBase;
}  // namespace MORTAR

namespace STR
{
  namespace MODELEVALUATOR
  {
    class ContactData;

    class Contact : public Generic
    {
     public:
      //! setup class variables [derived]
      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //!@{

      //! [derived]
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_contact; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override;

      //! [derived]
      bool evaluate_force() override;

      //! [derived]
      bool evaluate_stiff() override;

      //! [derived]
      bool evaluate_force_stiff() override;

      //! [derived]
      void pre_evaluate() override;

      //! [derived]
      void post_evaluate() override;

      //! [derived]
      void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) override;

      //! [derived]
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! Perform a correction of adaptive parameters
      bool CorrectParameters(NOX::NLN::CorrectionType type) override;

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! [derived]
      void read_restart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override{};

      //! recover condensed Lagrange multipliers
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! [derived]
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;

      //! [derived]
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override;

      /// [derived]
      void RunPreSolve(const ::NOX::Solver::Generic& solver) override;

      //! [derived]
      void run_post_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void UpdateStepState(const double& timefac_n) override;

      //! [derived]
      void UpdateStepElement() override;

      //! [derived]
      void determine_stress_strain() override;

      //! [derived]
      void DetermineEnergy() override;

      //! [derived]
      void determine_optional_quantity() override;

      //! [derived]
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! [derived]
      void ResetStepState() override;

      //! [derived]
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! [derived]
      void PostOutput() override;

      //! [derived]
      bool evaluate_cheap_soc_rhs() override;

      //! [derived]
      bool assemble_cheap_soc_rhs(Epetra_Vector& f, const double& timefac_np) const override;

      //! @}

      //! @name Call-back routines
      //!@{

      Teuchos::RCP<const CORE::LINALG::SparseMatrix> GetJacobianBlock(const MatBlockType bt) const;

      /** \brief Assemble the structural right-hand side vector
       *
       *  \param[in] without_these_models  Exclude all models defined in this vector
       *                                   during the assembly
       *  \param[in] apply_dbc             Apply Dirichlet boundary conditions
       *
       *  \author hiermeier \date 08/17 */
      Teuchos::RCP<Epetra_Vector> assemble_force_of_models(
          const std::vector<INPAR::STR::ModelType>* without_these_models = nullptr,
          const bool apply_dbc = false) const;

      virtual Teuchos::RCP<CORE::LINALG::SparseOperator> GetAuxDisplJacobian() const;

      void evaluate_weighted_gap_gradient_error();

      //!@}

      //! @name Accessors
      //!@{

      //! Returns a pointer to the underlying contact strategy object
      const Teuchos::RCP<CONTACT::AbstractStrategy>& StrategyPtr();

      //! Returns the underlying contact strategy object
      CONTACT::AbstractStrategy& Strategy();
      const CONTACT::AbstractStrategy& Strategy() const;

      //!@}

     protected:
      STR::MODELEVALUATOR::ContactData& EvalContact();
      const STR::MODELEVALUATOR::ContactData& EvalContact() const;

      virtual void CheckPseudo2D() const;

     private:
      void post_setup(Teuchos::ParameterList& cparams);

      /// Set the correct time integration parameters within the contact strategy
      void set_time_integration_info(CONTACT::AbstractStrategy& strategy) const;

      void post_update_step_state();

      void extend_lagrange_multiplier_domain(Teuchos::RCP<Epetra_Vector>& lm_vec) const;

      //! contact evaluation data container
      Teuchos::RCP<STR::MODELEVALUATOR::ContactData> eval_contact_ptr_;

      //! contact strategy
      Teuchos::RCP<CONTACT::AbstractStrategy> strategy_ptr_;

    };  // class Contact

  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
