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
#include "4C_fem_general_utils_createdis.hpp"
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
void Solid::ModelEvaluator::Constraints::setup()
{
  check_init();

  constraint_stiff_ptr_ = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*global_state().dof_row_map_view(), 81, true, true));

  constraint_force_ptr_ = Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map_view(), true));

  set_sub_model_types();
  create_sub_model_evaluators();


  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::set_sub_model_types()
{
  check_init();

  submodeltypes_ = std::set<enum Inpar::CONSTRAINTS::SubModelType>();

  // ---------------------------------------------------------------------------
  // check for multi point constraints
  // ---------------------------------------------------------------------------
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation;

  discret_ptr()->get_condition("LinePeriodicRve", linePeriodicRve);
  discret_ptr()->get_condition("SurfacePeriodicRve", surfPeriodicRve);
  discret_ptr()->get_condition("PointLinearCoupledEquation", pointLinearCoupledEquation);

  if (linePeriodicRve.size() > 0 || surfPeriodicRve.size() > 0 ||
      pointLinearCoupledEquation.size() > 0)
  {
    submodeltypes_.insert(Inpar::CONSTRAINTS::SubModelType::submodel_pbc_rve);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::create_sub_model_evaluators()
{
  // Create vector with the Sub-model-evaluators
  sub_model_vec_ptr_ = Solid::ModelEvaluator::Constraints::SubmodelevaluatorVector(0);

  for (const auto& mt : submodeltypes_)
  {
    switch (mt)
    {
      case Inpar::CONSTRAINTS::SubModelType::submodel_pbc_rve:
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
void Solid::ModelEvaluator::Constraints::reset(const Epetra_Vector& x)
{
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->reset();
  }
  constraint_stiff_ptr_->zero();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::evaluate_force()
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
bool Solid::ModelEvaluator::Constraints::evaluate_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->un_complete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->evaluate_force_stiff(constraint_stiff_ptr_, Teuchos::null);
  }
  if (not constraint_stiff_ptr_->filled()) constraint_stiff_ptr_->complete();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::evaluate_force_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->un_complete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->evaluate_force_stiff(constraint_stiff_ptr_, constraint_force_ptr_);
  }
  if (not constraint_stiff_ptr_->filled()) constraint_stiff_ptr_->complete();

  return true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::pre_evaluate()
{
  for (auto& sme : sub_model_vec_ptr_)
  {
    sme->evaluate_coupling_terms(*global_state_ptr());
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  Core::LinAlg::AssembleMyVector(1.0, f, timefac_np, *constraint_force_ptr_);
  constraint_force_ptr_->PutScalar(0.0);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> jac_dd_ptr = global_state().extract_displ_block(jac);

  jac_dd_ptr->add(*constraint_stiff_ptr_, false, timefac_np, 1.0);

  constraint_stiff_ptr_->zero();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // There is nothing to write for now
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  // There is nothing to read for now
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::predict(const Inpar::Solid::PredEnum& pred_type) {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::update_step_state(const double& timefac_n)
{
  if (not constraint_force_ptr_.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstruct_ptr = global_state().get_fstructure_old();
    fstruct_ptr->Update(timefac_n, *constraint_force_ptr_, 1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::update_step_element() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_stress_strain() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_energy()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_optional_quantity() {}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::reset_step_state()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::runtime_pre_output_step_state() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::runtime_output_step_state() const {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Solid::ModelEvaluator::Constraints::get_block_dof_row_map_ptr() const
{
  return global_state().dof_row_map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Solid::ModelEvaluator::Constraints::get_current_solution_ptr()
    const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
Solid::ModelEvaluator::Constraints::get_last_time_step_solution_ptr() const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::post_output() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::
    evaluate_jacobian_contributions_from_element_level_for_ptc()
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::assemble_jacobian_contributions_from_element_level_for_ptc(
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& modjac, const double& timefac_n)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::create_backup_state(const Epetra_Vector& dir)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::recover_from_backup_state()
{
  FOUR_C_THROW("This function is not yet implemented");
}

FOUR_C_NAMESPACE_CLOSE
