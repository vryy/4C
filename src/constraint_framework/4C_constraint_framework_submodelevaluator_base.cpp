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
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN


bool CONSTRAINTS::SUBMODELEVALUATOR::ConstraintBase::evaluate_force_stiff(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> me_stiff_ptr, Teuchos::RCP<Epetra_Vector> me_force_ptr)
{
  if (me_stiff_ptr == Teuchos::null && me_force_ptr == Teuchos::null)
    FOUR_C_THROW("Both stiffness and force point are null");

  if (me_stiff_ptr != Teuchos::null)
  {
    if (!(Q_Ld_->filled() && Q_dd_->filled() && Q_dL_->filled()))
      FOUR_C_THROW("Call evaluate_coupling_terms() first.");

    // evaluate the stiffness contribution of this sme:
    auto sme_stiff_ptr = Core::LinAlg::MatrixMultiply(*Q_dL_, false, *Q_Ld_, false, false);
    sme_stiff_ptr->scale(penalty_parameter_);
    sme_stiff_ptr->add(*Q_dd_, false, 1.0, 1.0);
    sme_stiff_ptr->complete();

    // add it to the modelevaluator stiffness
    me_stiff_ptr->add(*sme_stiff_ptr, false, 1., 1.);
  }

  if (me_force_ptr != Teuchos::null)
  {
    //  Calculate force contribution
    Teuchos::RCP<Epetra_Vector> r_pen =
        Teuchos::rcp(new Epetra_Vector(stiff_ptr_->row_map(), true));
    Q_Ld_->multiply(true, *constraint_vector_, *r_pen);
    Core::LinAlg::AssembleMyVector(1.0, *me_force_ptr, penalty_parameter_, *r_pen);
  }
  return true;
}

void CONSTRAINTS::SUBMODELEVALUATOR::ConstraintBase::evaluate_coupling_terms(
    Solid::TimeInt::BaseDataGlobalState& gstate)
{
  // Get the number of multipoint equations
  int ncon_ = 0;
  for (const auto& mpc : listMPCs_) ncon_ += mpc->get_number_of_mp_cs();

  // ToDo: Add an offset to the contraint dof map.
  n_condition_map_ = Teuchos::rcp(new Epetra_Map(ncon_, 0, stiff_ptr_->Comm()));

  // initialise all global coupling objects
  constraint_vector_ = Teuchos::rcp(new Epetra_Vector(*n_condition_map_, true));
  Q_Ld_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*n_condition_map_, 4));
  Q_dL_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(stiff_ptr_->row_map(), 4));
  Q_dd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(stiff_ptr_->row_map(), 0));

  // set Q_dd to zero as default
  Q_dd_->zero();
  // Evaluate the Constraint Pairs / equations objects
  Teuchos::RCP<const Epetra_Vector> dis_np = gstate.get_dis_np();
  for (const auto& obj : listMPCs_)
  {
    obj->evaluate_equation(*Q_dd_, *Q_dL_, *Q_Ld_, *constraint_vector_, *dis_np);
  }
  Core::IO::cout(Core::IO::verbose) << "Evaluated all constraint objects" << Core::IO::endl;

  // Complete
  Q_dd_->complete();
  Q_Ld_->complete(stiff_ptr_->domain_map(), *n_condition_map_);
  Q_dL_->complete(*n_condition_map_, stiff_ptr_->domain_map());
}
FOUR_C_NAMESPACE_CLOSE