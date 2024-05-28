/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of constraint terms.


\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_constraint_framework_model_evaluator.hpp"

#include "4C_constraint_framework_submodelevaluator_mpc.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/


void STR::MODELEVALUATOR::Constraints::Setup()
{
  check_init();

  constraint_stiff_ptr_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*g_state().DofRowMapView(), 81, true, true));

  constraint_force_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().DofRowMapView(), true));

  set_sub_model_types();
  create_sub_model_evaluators();


  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::set_sub_model_types()
{
  check_init();

  submodeltypes_ = std::set<enum INPAR::CONSTRAINTS::SubModelType>();

  // ---------------------------------------------------------------------------
  // check for multi point constraints
  // ---------------------------------------------------------------------------
  std::vector<Teuchos::RCP<CORE::Conditions::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation;

  discret_ptr()->GetCondition("LinePeriodicRve", linePeriodicRve);
  discret_ptr()->GetCondition("SurfacePeriodicRve", surfPeriodicRve);
  discret_ptr()->GetCondition("PointLinearCoupledEquation", pointLinearCoupledEquation);

  if (linePeriodicRve.size() > 0 || surfPeriodicRve.size() > 0 ||
      pointLinearCoupledEquation.size() > 0)
  {
    submodeltypes_.insert(INPAR::CONSTRAINTS::SubModelType::submodel_pbc_rve);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::create_sub_model_evaluators()
{
  // Create vector with the Sub-model-evaluators
  sub_model_vec_ptr_ = STR::MODELEVALUATOR::Constraints::SubmodelevaluatorVector(0);

  for (const auto& mt : submodeltypes_)
  {
    switch (mt)
    {
      case INPAR::CONSTRAINTS::SubModelType::submodel_pbc_rve:
      {
        sub_model_vec_ptr_.emplace_back(
            Teuchos::rcp(new CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager(
                discret_ptr(), constraint_stiff_ptr_.get())));

        break;
      }

      default:
      {
        FOUR_C_THROW(
            "Something went wrong: Apparently a Constraint ME was created that is not "
            "required. Check the Adapter");
      }
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::Reset(const Epetra_Vector& x)
{
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->Reset();
  }
  constraint_stiff_ptr_->Zero();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::evaluate_force()
{
  pre_evaluate();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->evaluate_force_stiff(Teuchos::null, constraint_force_ptr_);
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::evaluate_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->UnComplete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->evaluate_force_stiff(constraint_stiff_ptr_, Teuchos::null);
  }
  if (not constraint_stiff_ptr_->Filled()) constraint_stiff_ptr_->Complete();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::evaluate_force_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->UnComplete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->evaluate_force_stiff(constraint_stiff_ptr_, constraint_force_ptr_);
  }
  if (not constraint_stiff_ptr_->Filled()) constraint_stiff_ptr_->Complete();

  return true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::pre_evaluate()
{
  for (auto& sme : sub_model_vec_ptr_)
  {
    sme->evaluate_coupling_terms(*g_state_ptr());
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *constraint_force_ptr_);
  constraint_force_ptr_->PutScalar(0.0);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::assemble_jacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);

  jac_dd_ptr->Add(*constraint_stiff_ptr_, false, timefac_np, 1.0);

  constraint_stiff_ptr_->Zero();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::write_restart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // There is nothing to write for now
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::read_restart(IO::DiscretizationReader& ioreader)
{
  // There is nothing to read for now
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::Predict(const INPAR::STR::PredEnum& pred_type) {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::UpdateStepState(const double& timefac_n)
{
  if (not constraint_force_ptr_.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstruct_ptr = g_state().GetFstructureOld();
    fstruct_ptr->Update(timefac_n, *constraint_force_ptr_, 1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::UpdateStepElement() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::determine_stress_strain() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::DetermineEnergy()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::determine_optional_quantity() {}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::ResetStepState()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::OutputStepState(IO::DiscretizationWriter& iowriter) const {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::runtime_pre_output_step_state() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::runtime_output_step_state() const {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Constraints::get_block_dof_row_map_ptr() const
{
  return GState().dof_row_map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Constraints::get_current_solution_ptr() const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::Constraints::get_last_time_step_solution_ptr() const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::PostOutput() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::evaluate_jacobian_contributions_from_element_level_for_ptc()
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::assemble_jacobian_contributions_from_element_level_for_ptc(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::CreateBackupState(const Epetra_Vector& dir)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::recover_from_backup_state()
{
  FOUR_C_THROW("This function is not yet implemented");
}

FOUR_C_NAMESPACE_CLOSE
