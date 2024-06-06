/*-----------------------------------------------------------*/
/*! \file

\brief Tangential displacement predictor.



\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_predict_tangdis.hpp"

#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Predict::TangDis::TangDis()
    : dbc_incr_ptr_(Teuchos::null), apply_linear_reaction_forces_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::TangDis::Setup()
{
  check_init();
  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln group in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_grp_opt = nox_params().sublist("Group Options");
  // Get the current map. If there is no map, return a new empty one. (reference)
  NOX::Nln::GROUP::PrePostOperator::Map& prepostgroup_map =
      NOX::Nln::GROUP::PrePostOp::GetMap(p_grp_opt);
  // create the new tangdis pre/post operator
  Teuchos::RCP<NOX::Nln::Abstract::PrePostOperator> preposttangdis_ptr =
      Teuchos::rcp(new NOX::Nln::GROUP::PrePostOp::TangDis(Teuchos::rcp(this, false)));
  // insert/replace the old pointer in the map
  prepostgroup_map[NOX::Nln::GROUP::prepost_tangdis] = preposttangdis_ptr;

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::TangDis::Compute(::NOX::Abstract::Group& grp)
{
  check_init_setup();
  NOX::Nln::Group* grp_ptr = dynamic_cast<NOX::Nln::Group*>(&grp);
  FOUR_C_ASSERT(grp_ptr != nullptr, "Dynamic cast failed!");
  grp_ptr->reset_pre_post_operator(nox_params().sublist("Group Options"));

  impl_int().EvalData().SetPredictorType(Inpar::STR::pred_tangdis);

  // ---------------------------------------------------------------------------
  // calculate the dbc increment on the dirichlet boundary
  // ---------------------------------------------------------------------------
  dbc_incr_ptr_ = dbc().get_dirichlet_increment();

  // ---------------------------------------------------------------------------
  // We create at this point a new solution vector and initialize it
  // with the values of the last converged time step.
  // ---------------------------------------------------------------------------
  Teuchos::RCP<::NOX::Epetra::Vector> x_ptr = global_state().CreateGlobalVector(
      TimeInt::BaseDataGlobalState::VecInitType::last_time_step, impl_int().ModelEvalPtr());
  // Set the solution vector in the nox group. This will reset all isValid
  // flags.
  grp.setX(*x_ptr);

  // ---------------------------------------------------------------------------
  // Compute F and jacobian and apply the linear reaction forces due to changing
  // Dirichlet boundary conditions.
  // ---------------------------------------------------------------------------
  apply_linear_reaction_forces_ = true;
  grp_ptr->computeFandJacobian();
  apply_linear_reaction_forces_ = false;

  // ---------------------------------------------------------------------------
  // Check if we are using a Newton direction
  // ---------------------------------------------------------------------------
  std::string dir_str = nox_params().sublist("Direction").get<std::string>("Method");
  if (dir_str == "User Defined")
    dir_str = nox_params().sublist("Direction").get<std::string>("User Defined Method");
  if (dir_str != "Newton" and dir_str != "Modified Newton")
    FOUR_C_THROW(
        "The TangDis predictor is currently only working for the direction-"
        "methods \"Newton\" and \"Modified Newton\".");

  // ---------------------------------------------------------------------------
  // (re)set the linear solver parameters
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p =
      nox_params().sublist("Direction").sublist("Newton").sublist("Linear Solver");
  p.set<int>("Number of Nonlinear Iterations", 0);
  p.set<int>("Current Time Step", global_state().GetStepNp());
  // ToDo Get the actual tolerance value
  p.set<double>("Wanted Tolerance", 1.0e-6);

  // ---------------------------------------------------------------------------
  // solve the linear system of equations and update the current state
  // ---------------------------------------------------------------------------
  // compute the Newton direction
  grp_ptr->computeNewton(p);
  // reset isValid flags
  grp_ptr->computeX(*grp_ptr, grp_ptr->getNewton(), 1.0);
  // add the DBC values to the current state vector
  Teuchos::RCP<Epetra_Vector> dbc_incr_exp_ptr =
      Teuchos::rcp(new Epetra_Vector(global_state().GlobalProblemMap(), true));
  Core::LinAlg::Export(*dbc_incr_ptr_, *dbc_incr_exp_ptr);
  grp_ptr->computeX(*grp_ptr, *dbc_incr_exp_ptr, 1.0);
  // Reset the state variables
  const ::NOX::Epetra::Vector& x_eptra =
      dynamic_cast<const ::NOX::Epetra::Vector&>(grp_ptr->getX());
  // set the consistent state in the models (e.g. structure and contact models)
  impl_int().ResetModelStates(x_eptra.getEpetraVector());

  // For safety purposes, we set the dbc_incr vector to zero
  dbc_incr_ptr_->PutScalar(0.0);

  impl_int().ModelEval().Predict(GetType());

  impl_int().EvalData().SetPredictorType(Inpar::STR::pred_vague);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::Predict::TangDis::GetDbcIncr() const
{
  FOUR_C_ASSERT(!dbc_incr_ptr_.is_null(), "The dbc increment is not initialized!");
  return *dbc_incr_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& STR::Predict::TangDis::is_apply_linear_reaction_forces() const
{
  return apply_linear_reaction_forces_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Predict::TangDis::pre_apply_force_external(Epetra_Vector& fextnp) const
{
  check_init_setup();

  if (GetType() != Inpar::STR::pred_tangdis_constfext) return false;

  if (apply_linear_reaction_forces_)
  {
    fextnp.Scale(1.0, *GlobalState().GetFextN());
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOp::TangDis::TangDis(
    const Teuchos::RCP<const STR::Predict::TangDis>& tang_predict_ptr)
    : tang_predict_ptr_(tang_predict_ptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::GROUP::PrePostOp::TangDis::runPostComputeF(
    Epetra_Vector& F, const NOX::Nln::Group& grp)
{
  // If we do not want to apply linear reaction forces due to changing Dirichlet
  // boundary conditions, we just return.
  if (not tang_predict_ptr_->is_apply_linear_reaction_forces()) return;

  // get the new dirichlet boundary increment
  const Epetra_Vector& dbc_incr = tang_predict_ptr_->GetDbcIncr();

  double dbc_incr_nrm2 = 0.0;
  dbc_incr.Norm2(&dbc_incr_nrm2);

  // If there are only Neumann loads, do a direct return.
  if (dbc_incr_nrm2 == 0.0) return;

  /* Alternatively, it's also possible to get a const pointer on the jacobian
   * by calling grp.getLinearSystem()->getJacobianOperator()... */
  Teuchos::RCP<const Core::LinAlg::SparseMatrix> stiff_ptr =
      tang_predict_ptr_->GlobalState().get_jacobian_displ_block();

  // check if the jacobian is filled
  if (not stiff_ptr->Filled()) FOUR_C_THROW("The jacobian is not yet filled!");

  Teuchos::RCP<Epetra_Vector> freact_ptr =
      Teuchos::rcp(new Epetra_Vector(*tang_predict_ptr_->GlobalState().DofRowMapView()));
  if (stiff_ptr->Multiply(false, dbc_incr, *freact_ptr)) FOUR_C_THROW("Multiply failed!");

  // finally add the linear reaction forces to the current rhs
  Core::LinAlg::AssembleMyVector(1.0, F, 1.0, *freact_ptr);
}

FOUR_C_NAMESPACE_CLOSE
