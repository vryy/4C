/*-----------------------------------------------------------*/
/*!
\file str_predict_tangdis.cpp

\brief Tangential displacemnt predictor.

\maintainer Michael Hiermeier

\date Sep 1, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_predict_tangdis.H"
#include "str_timint_base.H"
#include "str_dbc.H"
#include "str_utils.H"
#include "str_impl_generic.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dserror.H"

#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"

#include <NOX_Epetra_Vector.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::TangDis::TangDis()
    : dbc_incr_ptr_(Teuchos::null),
      applyLinearReactionForces_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  STR::PREDICT::TangDis::Setup()
{
  CheckInit();
  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln group in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_grp_opt = NoxParams().sublist("Group Options");
  // Get the current map. If there is no map, return a new empty one. (reference)
  NOX::NLN::GROUP::PrePostOperator::Map& prepostgroup_map =
      NOX::NLN::GROUP::PrePostOp::GetMutableMap(p_grp_opt);
  // create the new tangdis pre/post operator
  Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> preposttangdis_ptr =
      Teuchos::rcp(new NOX::NLN::GROUP::PrePostOp::TangDis(Teuchos::rcp(this,false)));
  // insert/replace the old pointer in the map
  prepostgroup_map[NOX::NLN::GROUP::prepost_tangdis] = preposttangdis_ptr;

  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::TangDis::Compute(NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  NOX::NLN::Group* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (grp_ptr == NULL)
    dserror("Dynamic cast failed!");
  grp_ptr->ResetPrePostOperator(NoxParams().sublist("Group Options"));

  // ---------------------------------------------------------------------------
  // calculate the dbc increment on the dirichlet boundary
  // ---------------------------------------------------------------------------
  dbc_incr_ptr_ = Dbc().GetDirichletIncrement();

  // ---------------------------------------------------------------------------
  // We create at this point a new solution vector and initialize it
  // with the values of the last converged time step.
  // ---------------------------------------------------------------------------
  Teuchos::RCP<NOX::Epetra::Vector> x_ptr =
      GlobalState().CreateGlobalVector(STR::vec_init_last_time_step,
          ImplInt().ModelEvalPtr());
  // Set the solution vector in the nox group. This will reset all isValid
  // flags.
  grp.setX(*x_ptr);

  // ---------------------------------------------------------------------------
  // Compute F and jacobian and apply the linear reaction forces due to changing
  // Dirichlet boundary conditions.
  // ---------------------------------------------------------------------------
  applyLinearReactionForces_ = true;
  grp_ptr->computeFandJacobian();
  applyLinearReactionForces_ = false;

  // ---------------------------------------------------------------------------
  // Check if we are using a Newton direction
  // ---------------------------------------------------------------------------
  const std::string& dir_str =
      NoxParams().sublist("Direction").get<std::string>("Method");
  if (dir_str != "Newton")
    dserror("The TangDis predictor is currently only working for the direction-"
        "method \"Newton\".");

  // ---------------------------------------------------------------------------
  // (re)set the linear solver parameters
  // solve the linear system of equations and update the current state
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p =
      NoxParams().sublist("Direction").sublist("Newton").sublist("Linear Solver");
  p.set<int>("Number of Nonlinear Iterations",0);
  p.set<int>("Current Time Step",GlobalState().GetStepNp());
  // ToDo Get the actual tolerance value
  p.set<double>("Wanted Tolerance",1.0e-6);
  // compute the Newton direction
  grp_ptr->computeNewton(p);
  // reset isValid flags
  grp_ptr->computeX(*grp_ptr,grp_ptr->getNewton(),1.0);
  // add the DBC values to the current state vector
  Teuchos::RCP<Epetra_Vector> dbc_incr_exp_ptr =
      Teuchos::rcp(new Epetra_Vector(GlobalState().GlobalProblemMap(),true));
  LINALG::Export(*dbc_incr_ptr_,*dbc_incr_exp_ptr);
  grp_ptr->computeX(*grp_ptr,*dbc_incr_exp_ptr,1.0);
  // Reset the state variables
  const NOX::Epetra::Vector& x_eptra =
      dynamic_cast<const NOX::Epetra::Vector&>(grp_ptr->getX());
  // set the consistent state in the active implicit time integrator
  ImplInt().SetState(x_eptra.getEpetraVector());

  // For safety purposes, we set the dbc_incr vector to zero
  dbc_incr_ptr_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::PREDICT::TangDis::GetDbcIncr() const
{
  if (dbc_incr_ptr_.is_null())
    dserror("The dbc increment is not initialized!");
  return *dbc_incr_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& STR::PREDICT::TangDis::IsApplyLinearReactionForces() const
{
  return applyLinearReactionForces_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOp::TangDis::TangDis(
    const Teuchos::RCP<const ::STR::PREDICT::TangDis>& tang_predict_ptr)
    : tang_predict_ptr_(tang_predict_ptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOp::TangDis::runPostComputeF(
    Epetra_Vector& F, const NOX::NLN::Group& grp)
{
  // If we do not want to apply linear reaction forces due to changing Dirichlet
  // boundary conditions, we just return.
  if (not tang_predict_ptr_->IsApplyLinearReactionForces())
    return;

  // get the new dirichlet boundary increment
  const Epetra_Vector& dbc_incr = tang_predict_ptr_->GetDbcIncr();

  double dbc_incr_nrm2 = 0.0;
  dbc_incr.Norm2(&dbc_incr_nrm2);

  // If there are only Neumann loads, do a direct return.
  if (dbc_incr_nrm2 == 0.0)
    return;

  /* Alternatively, it's also possible to get a const pointer on the jacobian
   * by calling grp.getLinearSystem()->getJacobianOperator()... */
  Teuchos::RCP<const LINALG::SparseMatrix> stiff_ptr =
      tang_predict_ptr_->GlobalState().GetJacobianDisplBlock();

  // check if the jacobian is filled
  if (not stiff_ptr->Filled())
    dserror("The jacobian is not yet filled!");

  Teuchos::RCP<Epetra_Vector> freact_ptr =
      Teuchos::rcp(new Epetra_Vector(*tang_predict_ptr_->GlobalState().DofRowMapView()));
  if (stiff_ptr->Multiply(false,dbc_incr,*freact_ptr))
    dserror("Multiply failed!");

  // finally add the linear reaction forces to the current rhs
  ::STR::AssembleVector(1.0,F,1.0,*freact_ptr);

  return;
}
