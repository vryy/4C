/*-----------------------------------------------------------*/
/*! \file

\brief Manager of the model evaluator calls.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_model_evaluator.hpp"

#include "4C_linalg_blocksparsematrix.hpp"  // debugging
#include "4C_linalg_sparseoperator.hpp"
#include "4C_structure_new_integrator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_model_evaluator_factory.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator::ModelEvaluator()
    : isinit_(false),
      issetup_(false),
      me_map_ptr_(Teuchos::null),
      me_vec_ptr_(Teuchos::null),
      eval_data_ptr_(Teuchos::null),
      sdyn_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      gio_ptr_(Teuchos::null),
      int_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::check_init() const { FOUR_C_ASSERT(is_init(), "Call Init() first!"); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  issetup_ = false;

  eval_data_ptr_ = eval_data_ptr;
  sdyn_ptr_ = sdyn_ptr;
  gstate_ptr_ = gstate_ptr;
  gio_ptr_ = gio_ptr;
  int_ptr_ = int_ptr;
  timint_ptr_ = timint_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Setup()
{
  check_init();

  me_map_ptr_ = STR::MODELEVALUATOR::build_model_evaluators(
      sdyn_ptr_->GetModelTypes(), sdyn_ptr_->CouplingModelPtr());

  me_vec_ptr_ = transform_to_vector(*me_map_ptr_);

  Map::iterator me_iter;
  int dof_offset = 0;
  for (me_iter = me_map_ptr_->begin(); me_iter != me_map_ptr_->end(); ++me_iter)
  {
    me_iter->second->Init(eval_data_ptr_, gstate_ptr_, gio_ptr_, int_ptr_, timint_ptr_, dof_offset);
    me_iter->second->Setup();
    // setup the block information for saddle point problems
    dof_offset = gstate_ptr_->setup_block_information(*(me_iter->second), me_iter->first);
  }
  gstate_ptr_->setup_multi_map_extractor();
  gstate_ptr_->setup_element_technology_map_extractors();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::setup_multi_map_extractor()
{
  // setup the block information for saddle point problems
  for (Vector::iterator me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    gstate_ptr_->setup_block_information(**me_iter, (**me_iter).Type());

  gstate_ptr_->setup_multi_map_extractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::initialize_inertia_and_damping(
    const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac)
{
  check_init_setup();

  // initialize stiffness matrix to zero
  jac.Zero();
  // get structural model evaluator
  STR::MODELEVALUATOR::Structure& str_model =
      dynamic_cast<STR::MODELEVALUATOR::Structure&>(Evaluator(INPAR::STR::model_structure));

  str_model.Reset(x);

  return str_model.initialize_inertia_and_damping();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::assemble_force(const double timefac_np, Epetra_Vector& f,
    const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  if (not without_these_models) return assemble_force(timefac_np, f);

  Vector partial_me_vec;
  partial_me_vec.reserve(me_vec_ptr_->size());
  split_model_vector(partial_me_vec, *without_these_models);

  return assemble_force(partial_me_vec, timefac_np, f);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::assemble_force(
    bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->assemble_force(f, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::assemble_jacobian(const double timefac_np,
    CORE::LINALG::SparseOperator& jac,
    const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  if (not without_these_models) return assemble_jacobian(timefac_np, jac);

  Vector partial_me_vec;
  split_model_vector(partial_me_vec, *without_these_models);

  return assemble_jacobian(partial_me_vec, timefac_np, jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::assemble_jacobian(bool& ok, const Vector& me_vec, const double timefac_np,
    CORE::LINALG::SparseOperator& jac) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->assemble_jacobian(jac, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::assemble_jacobian_contributions_from_element_level_for_ptc(
    const Vector& me_vec, const double timefac_np, Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac)
{
  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    (*cit)->assemble_jacobian_contributions_from_element_level_for_ptc(modjac, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::evaluate_force(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  pre_evaluate(ok, me_vec);

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->evaluate_force() : false);

  post_evaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::evaluate_stiff(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  pre_evaluate(ok, me_vec);

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->evaluate_stiff() : false);

  post_evaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::evaluate_force_stiff(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  pre_evaluate(ok, me_vec);


  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->evaluate_force_stiff() : false);

  post_evaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::pre_evaluate(bool ok, const Vector& me_vec) const
{
  for (auto& me : me_vec) me->pre_evaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::post_evaluate(bool ok, const Vector& me_vec) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    (*cit)->post_evaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyInitialForce(const Epetra_Vector& x, Epetra_Vector& f)
{
  check_init_setup();

  Vector::iterator me_iter;
  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x, false);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->evaluate_initial_force() : false);

  post_evaluate(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together, including mass and viscous contributions
  // ---------------------------------------------------------------------------
  assemble_force(ok, *me_vec_ptr_, 1.0, f);

  // ---------------------------------------------------------------------------
  // subtract mass and viscous contributions from initial force vector
  // ---------------------------------------------------------------------------
  f.Scale(-1.);
  int_ptr_->add_visco_mass_contributions(f);
  f.Scale(-1.);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x) const
{
  // default ResetStates call
  ResetStates(x, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x, bool setstate) const
{
  ResetStates(x, setstate, *me_vec_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x, bool setstate, Vector& me_vec) const
{
  if (setstate) int_ptr_->set_state(x);
  for (auto& me_iter : me_vec) me_iter->Reset(x);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForce(
    const Epetra_Vector& x, Epetra_Vector& f, const double& timefac_np) const
{
  check_init_setup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // and also update the state variables of the current time integrator
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  evaluate_force(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  assemble_force(ok, *me_vec_ptr_, timefac_np, f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(
    const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  check_init_setup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix to zero
  jac.Zero();

  // ---------------------------------------------------------------------------
  // update the state variables of the current time integrator
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  evaluate_stiff(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  assemble_jacobian(ok, *me_vec_ptr_, timefac_np, jac);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(const INPAR::STR::ModelType& mt, const Epetra_Vector& x,
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  check_init_setup();
  bool ok = true;
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> model_ptr = me_map_ptr_->at(mt);
  const Vector me_vec(1, model_ptr);

  // initialize stiffness matrix to zero
  jac.Zero();

  // update the state variables of the current time integrator
  int_ptr_->set_state(x);
  model_ptr->Reset(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  evaluate_stiff(ok, me_vec);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  assemble_jacobian(ok, me_vec, timefac_np, jac);

  return ok;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForceStiff(const Epetra_Vector& x, Epetra_Vector& f,
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  check_init_setup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix and right hand side to zero
  f.PutScalar(0.0);
  jac.Zero();

  // ---------------------------------------------------------------------------
  // update the state variables of the current time integrator
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  evaluate_force_stiff(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  assemble_force(ok, *me_vec_ptr_, timefac_np, f);
  assemble_jacobian(ok, *me_vec_ptr_, timefac_np, jac);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyCheapSOCRhs(const enum NOX::NLN::CorrectionType type,
    const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
    Epetra_Vector& f, const double& timefac_np) const
{
  check_init_setup();

  Vector constraint_me_vec;
  extract_model_vector(constraint_me_vec, constraint_models);

  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);

  // ---------------------------------------------------------------------------
  // update the state variables of the current time integrator
  // reset model specific variables of the constraint models
  // ---------------------------------------------------------------------------
  ResetStates(x, true, constraint_me_vec);

  // ---------------------------------------------------------------------------
  // evaluate all rhs terms of the constraint models
  // ---------------------------------------------------------------------------
  evaluate_cheap_soc_rhs(ok, constraint_me_vec);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  assemble_cheap_soc_rhs(ok, constraint_me_vec, timefac_np, f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::evaluate_cheap_soc_rhs(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->evaluate_cheap_soc_rhs() : false);

  post_evaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::assemble_cheap_soc_rhs(
    bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->assemble_cheap_soc_rhs(f, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::correct_parameters(const enum NOX::NLN::CorrectionType type) const
{
  bool ok = true;
  for (auto& cit : *me_vec_ptr_)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? cit->correct_parameters(type) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Predict(const INPAR::STR::PredEnum& pred_type) const
{
  check_init_setup();
  for (Vector::iterator me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->Predict(pred_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::write_restart(
    CORE::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::read_restart(CORE::IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->read_restart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunRecover() const
{
  check_init_setup();
  // set some parameters for the element evaluation
  eval_data_ptr_->SetIsDefaultStep(true);
  eval_data_ptr_->SetStepLength(1.0);
  eval_data_ptr_->ResetMyNorms(true);
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::run_post_compute_x(const Epetra_Vector& xold, const Epetra_Vector& dir,
    const double& step, const Epetra_Vector& xnew, const bool isdefaultstep) const
{
  check_init_setup();
  // set some parameters for the element evaluation
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);
  eval_data_ptr_->ResetMyNorms(isdefaultstep);
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->run_post_compute_x(xold, dir, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
    const double& step, const NOX::NLN::Group& curr_grp, const bool isdefaultstep) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->run_pre_compute_x(xold, dir_mutable, curr_grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::run_post_iterate(const ::NOX::Solver::Generic& solver, const double step,
    const bool isdefaultstep, const int num_corrs) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);
  eval_data_ptr_->set_number_of_modified_newton_corrections(num_corrs);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->run_post_iterate(solver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPreSolve(
    const ::NOX::Solver::Generic& solver, const double step, const bool isdefaultstep) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPreSolve(solver);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::run_post_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp) const
{
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->run_post_apply_jacobian_inverse(rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp) const
{
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->run_pre_apply_jacobian_inverse(rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::ModelEvaluator::GetGlobalState() const
{
  check_init_setup();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::ModelEvaluator::global_state_ptr()
{
  check_init_setup();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<const STR::TIMINT::Base>& STR::ModelEvaluator::GetTimIntPtr() const
{
  check_init_setup();
  return timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::ModelEvaluator::Evaluator(const enum INPAR::STR::ModelType& mt)
{
  check_init_setup();
  // sanity check, if there is a model evaluator for the given model type
  STR::ModelEvaluator::Map::const_iterator me_iter = me_map_ptr_->find(mt);
  if (me_iter == me_map_ptr_->end())
    FOUR_C_THROW("There is no model evaluator for the model type %s",
        INPAR::STR::ModelTypeString(mt).c_str());

  return *(me_iter->second);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::ModelEvaluator::Evaluator(
    const enum INPAR::STR::ModelType& mt) const
{
  check_init_setup();
  // sanity check, if there is a model evaluator for the given model type
  STR::ModelEvaluator::Map::const_iterator me_iter = me_map_ptr_->find(mt);
  if (me_iter == me_map_ptr_->end())
    FOUR_C_THROW("There is no model evaluator for the model type %s",
        INPAR::STR::ModelTypeString(mt).c_str());

  return *(me_iter->second);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepState(const double& timefac_n)
{
  check_init_setup();
  /* Reset old structural right hand side.
   * It will be filled within the model evaluators */
  gstate_ptr_->GetFstructureOld()->PutScalar(0.0);
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->UpdateStepState(timefac_n);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::compute_jacobian_contributions_from_element_level_for_ptc(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& scalingMatrixOpPtr)
{
  // evaluate ptc contributions at t^n+1
  double timefac_np = 1.0;
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->evaluate_jacobian_contributions_from_element_level_for_ptc();

  assemble_jacobian_contributions_from_element_level_for_ptc(
      *me_vec_ptr_, timefac_np, scalingMatrixOpPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const
{
  check_init_setup();
  for (auto& me : *me_vec_ptr_) me->remove_condensed_contributions_from_rhs(rhs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepElement()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::update_residual()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->update_residual();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::determine_stress_strain()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->determine_stress_strain();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineEnergy()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->DetermineEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::determine_optional_quantity()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->determine_optional_quantity();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::OutputStepState(CORE::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();
  Vector::const_iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::runtime_pre_output_step_state()
{
  check_init_setup();
  for (auto me : *me_vec_ptr_) me->runtime_pre_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::runtime_output_step_state() const
{
  check_init_setup();
  Vector::const_iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->runtime_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::PostOutput()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStepState()
{
  check_init_setup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Vector> STR::ModelEvaluator::transform_to_vector(
    const STR::ModelEvaluator::Map& model_map) const
{
  Teuchos::RCP<STR::ModelEvaluator::Vector> me_vec_ptr =
      Teuchos::rcp(new STR::ModelEvaluator::Vector(0));
  me_vec_ptr->reserve(model_map.size());

  // --------------------------------------------------------------------------
  // There must be a structural model evaluator at the first position
  // --------------------------------------------------------------------------
  if (model_map.begin()->first != INPAR::STR::model_structure)
    FOUR_C_THROW(
        "The first model evaluator in the model_map must be a "
        "structural model evaluator!");

  // --------------------------------------------------------------------------
  // Fill the vector
  // --------------------------------------------------------------------------
  for (Map::const_iterator cit = model_map.begin(); cit != model_map.end(); ++cit)
    me_vec_ptr->push_back(cit->second);

  return me_vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::split_model_vector(STR::ModelEvaluator::Vector& partial_me_vec,
    const std::vector<INPAR::STR::ModelType>& without_these_models) const
{
  partial_me_vec.reserve(me_vec_ptr_->size());
  for (Vector::const_iterator cit = me_vec_ptr_->begin(); cit != me_vec_ptr_->end(); ++cit)
  {
    const STR::MODELEVALUATOR::Generic& model = **cit;

    auto citer = without_these_models.cbegin();

    while (citer != without_these_models.cend())
    {
      if (*citer == model.Type()) break;

      ++citer;
    }

    if (citer == without_these_models.cend()) partial_me_vec.push_back(*cit);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::extract_model_vector(STR::ModelEvaluator::Vector& partial_me_vec,
    const std::vector<INPAR::STR::ModelType>& only_these_models) const
{
  partial_me_vec.reserve(only_these_models.size());
  for (const auto mtype : only_these_models)
  {
    auto cit = me_vec_ptr_->cbegin();
    while (cit != me_vec_ptr_->cend())
    {
      if ((*cit)->Type() == mtype) break;

      ++cit;
    }

    FOUR_C_ASSERT(cit != me_vec_ptr_->end(), "Couldn't find the model type in me_vec_ptr_.");

    partial_me_vec.push_back(*cit);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::CreateBackupState(const Epetra_Vector& dir)
{
  check_init_setup();
  for (const auto& me_iter : *me_vec_ptr_) me_iter->CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::recover_from_backup_state()
{
  check_init_setup();
  for (const auto& me_iter : *me_vec_ptr_) me_iter->recover_from_backup_state();
}

FOUR_C_NAMESPACE_CLOSE
