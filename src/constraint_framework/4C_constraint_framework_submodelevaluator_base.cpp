/*-----------------------------------------------------------*/
/*! \file
\brief Generic class for all constraint submodel evaluators.
\level 3
*/
#include "4C_config.hpp"

#include "4C_constraint_framework_submodelevaluator_base.hpp"

#include "4C_io_pstream.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN


bool CONSTRAINTS::SUBMODELEVALUATOR::ConstraintBase::evaluate_force_stiff(
    Teuchos::RCP<CORE::LINALG::SparseMatrix> me_stiff_ptr, Teuchos::RCP<Epetra_Vector> me_force_ptr)
{
  if (me_stiff_ptr == Teuchos::null && me_force_ptr == Teuchos::null)
    FOUR_C_THROW("Both stiffness and force point are null");

  if (me_stiff_ptr != Teuchos::null)
  {
    if (!(Q_Ld_->Filled() && Q_dd_->Filled() && Q_dL_->Filled()))
      FOUR_C_THROW("Call evaluate_coupling_terms() first.");

    // evaluate the stiffness contribution of this sme:
    auto sme_stiff_ptr = CORE::LINALG::Multiply(*Q_dL_, false, *Q_Ld_, false, false);
    sme_stiff_ptr->Scale(penalty_parameter_);
    sme_stiff_ptr->Add(*Q_dd_, false, 1.0, 1.0);
    sme_stiff_ptr->Complete();

    // add it to the modelevaluator stiffness
    me_stiff_ptr->Add(*sme_stiff_ptr->EpetraMatrix(), false, 1., 1.);
  }

  if (me_force_ptr != Teuchos::null)
  {
    //  Calculate force contribution
    Teuchos::RCP<Epetra_Vector> r_pen = Teuchos::rcp(new Epetra_Vector(stiff_ptr_->RowMap(), true));
    Q_Ld_->Multiply(true, *constraint_vector_, *r_pen);
    CORE::LINALG::AssembleMyVector(1.0, *me_force_ptr, penalty_parameter_, *r_pen);
  }
  return true;
}

void CONSTRAINTS::SUBMODELEVALUATOR::ConstraintBase::evaluate_coupling_terms(
    STR::TIMINT::BaseDataGlobalState& gstate)
{
  // Get the number of multipoint equations
  int ncon_ = 0;
  for (const auto& mpc : listMPCs_) ncon_ += mpc->GetNumberOfMPCs();

  // ToDo: Add an offset to the contraint dof map.
  n_condition_map_ = Teuchos::rcp(new Epetra_Map(ncon_, 0, stiff_ptr_->Comm()));

  // initialise all global coupling objects
  constraint_vector_ = Teuchos::rcp(new Epetra_Vector(*n_condition_map_, true));
  Q_Ld_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*n_condition_map_, 4));
  Q_dL_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(stiff_ptr_->RowMap(), 4));
  Q_dd_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(stiff_ptr_->RowMap(), 0));

  // set Q_dd to zero as default
  Q_dd_->Zero();
  // Evaluate the Constraint Pairs / equations objects
  Teuchos::RCP<const Epetra_Vector> dis_np = gstate.GetDisNp();
  for (const auto& obj : listMPCs_)
  {
    obj->EvaluateEquation(*Q_dd_, *Q_dL_, *Q_Ld_, *constraint_vector_, *dis_np);
  }
  IO::cout(IO::verbose) << "Evaluated all constraint objects" << IO::endl;

  // Complete
  Q_dd_->Complete();
  Q_Ld_->Complete(stiff_ptr_->DomainMap(), *n_condition_map_);
  Q_dL_->Complete(*n_condition_map_, stiff_ptr_->DomainMap());
}
FOUR_C_NAMESPACE_CLOSE