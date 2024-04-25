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
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/


void STR::MODELEVALUATOR::Constraints::Setup()
{
  CheckInit();

  constraint_stiff_ptr_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

  constraint_force_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView(), true));

  SetSubModelTypes();
  CreateSubModelEvaluators();


  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::SetSubModelTypes()
{
  CheckInit();

  submodeltypes_ = std::set<enum INPAR::CONSTRAINTS::SubModelType>();

  // ---------------------------------------------------------------------------
  // check for multi point constraints
  // ---------------------------------------------------------------------------
  std::vector<Teuchos::RCP<DRT::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation;

  DiscretPtr()->GetCondition("LinePeriodicRve", linePeriodicRve);
  DiscretPtr()->GetCondition("SurfacePeriodicRve", surfPeriodicRve);
  DiscretPtr()->GetCondition("PointLinearCoupledEquation", pointLinearCoupledEquation);

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
void STR::MODELEVALUATOR::Constraints::CreateSubModelEvaluators()
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
                DiscretPtr(), constraint_stiff_ptr_.get())));

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
bool STR::MODELEVALUATOR::Constraints::EvaluateForce()
{
  PreEvaluate();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->EvaluateForceStiff(Teuchos::null, constraint_force_ptr_);
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::EvaluateStiff()
{
  PreEvaluate();

  constraint_stiff_ptr_->UnComplete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->EvaluateForceStiff(constraint_stiff_ptr_, Teuchos::null);
  }
  if (not constraint_stiff_ptr_->Filled()) constraint_stiff_ptr_->Complete();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::EvaluateForceStiff()
{
  PreEvaluate();

  constraint_stiff_ptr_->UnComplete();
  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    sme_iter->EvaluateForceStiff(constraint_stiff_ptr_, constraint_force_ptr_);
  }
  if (not constraint_stiff_ptr_->Filled()) constraint_stiff_ptr_->Complete();

  return true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::PreEvaluate()
{
  for (auto& sme : sub_model_vec_ptr_)
  {
    sme->EvaluateCouplingTerms(*GStatePtr());
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *constraint_force_ptr_);
  constraint_force_ptr_->PutScalar(0.0);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraints::AssembleJacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);

  jac_dd_ptr->Add(*constraint_stiff_ptr_, false, timefac_np, 1.0);

  constraint_stiff_ptr_->Zero();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // There is nothing to write for now
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::ReadRestart(IO::DiscretizationReader& ioreader)
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
    Teuchos::RCP<Epetra_Vector>& fstruct_ptr = GState().GetFstructureOld();
    fstruct_ptr->Update(timefac_n, *constraint_force_ptr_, 1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::UpdateStepElement() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::DetermineStressStrain() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::DetermineEnergy()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::DetermineOptionalQuantity() {}
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
void STR::MODELEVALUATOR::Constraints::RuntimePreOutputStepState() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::RuntimeOutputStepState() const {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Constraints::GetBlockDofRowMapPtr() const
{
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Constraints::GetCurrentSolutionPtr() const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Constraints::GetLastTimeStepSolutionPtr()
    const
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::PostOutput() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::EvaluateJacobianContributionsFromElementLevelForPTC()
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraints::AssembleJacobianContributionsFromElementLevelForPTC(
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
void STR::MODELEVALUATOR::Constraints::RecoverFromBackupState()
{
  FOUR_C_THROW("This function is not yet implemented");
}

FOUR_C_NAMESPACE_CLOSE
