/*-----------------------------------------------------------*/
/*!

\brief Manager of the model evaluator calls.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator.H"
#include "str_model_evaluator_factory.H"
#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"
#include "str_integrator.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_blocksparsematrix.H"  // debugging

#include <Epetra_Vector.h>

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
void STR::ModelEvaluator::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::CheckInit() const
{
  if (not IsInit()) dserror("Call Init() first!");
}

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
  CheckInit();

  me_map_ptr_ = STR::MODELEVALUATOR::BuildModelEvaluators(
      sdyn_ptr_->GetModelTypes(), sdyn_ptr_->CouplingModelPtr());

  me_vec_ptr_ = TransformToVector(*me_map_ptr_);

  Map::iterator me_iter;
  int dof_offset = 0;
  unsigned int i = 0;
  for (me_iter = me_map_ptr_->begin(); me_iter != me_map_ptr_->end(); ++me_iter)
  {
    me_iter->second->Init(eval_data_ptr_, gstate_ptr_, gio_ptr_, int_ptr_, timint_ptr_, dof_offset);
    me_iter->second->Setup();
    // setup the block information for saddle point problems
    dof_offset = gstate_ptr_->SetupBlockInformation(*(me_iter->second), me_iter->first);
    ++i;
  }
  gstate_ptr_->SetupMultiMapExtractor();
  gstate_ptr_->SetupElementTechnologyMapExtractors();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::SetupMultiMapExtractor()
{
  // setup the block information for saddle point problems
  for (Vector::iterator me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    gstate_ptr_->SetupBlockInformation(**me_iter, (**me_iter).Type());

  gstate_ptr_->SetupMultiMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::InitializeInertiaAndDamping(
    const Epetra_Vector& x, LINALG::SparseOperator& jac)
{
  CheckInitSetup();

  // initialize stiffness matrix to zero
  jac.Zero();
  // get structural model evaluator
  STR::MODELEVALUATOR::Structure& str_model =
      dynamic_cast<STR::MODELEVALUATOR::Structure&>(Evaluator(INPAR::STR::model_structure));

  str_model.Reset(x);

  return str_model.InitializeInertiaAndDamping();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::AssembleForce(const double timefac_np, Epetra_Vector& f,
    const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  if (not without_these_models) return AssembleForce(timefac_np, f);

  Vector partial_me_vec;
  partial_me_vec.reserve(me_vec_ptr_->size());
  SplitModelVector(partial_me_vec, *without_these_models);

  return AssembleForce(partial_me_vec, timefac_np, f);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::AssembleForce(
    bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->AssembleForce(f, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::AssembleJacobian(const double timefac_np, LINALG::SparseOperator& jac,
    const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  if (not without_these_models) return AssembleJacobian(timefac_np, jac);

  Vector partial_me_vec;
  SplitModelVector(partial_me_vec, *without_these_models);

  return AssembleJacobian(partial_me_vec, timefac_np, jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::AssembleJacobian(
    bool& ok, const Vector& me_vec, const double timefac_np, LINALG::SparseOperator& jac) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->AssembleJacobian(jac, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::AssembleJacobianContributionsFromElementLevelForPTC(
    const Vector& me_vec, const double timefac_np, Teuchos::RCP<LINALG::SparseMatrix>& modjac)
{
  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    (*cit)->AssembleJacobianContributionsFromElementLevelForPTC(modjac, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::EvaluateForce(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  PreEvaluate(ok, me_vec);

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->EvaluateForce() : false);

  PostEvaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::EvaluateStiff(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  PreEvaluate(ok, me_vec);

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->EvaluateStiff() : false);

  PostEvaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::EvaluateForceStiff(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  PreEvaluate(ok, me_vec);


  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->EvaluateForceStiff() : false);

  PostEvaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::PreEvaluate(bool ok, const Vector& me_vec) const
{
  for (auto& me : me_vec) me->PreEvaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::PostEvaluate(bool ok, const Vector& me_vec) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    (*cit)->PostEvaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyInitialForce(const Epetra_Vector& x, Epetra_Vector& f)
{
  CheckInitSetup();

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
    ok = (ok ? (*me_iter)->EvaluateInitialForce() : false);

  PostEvaluate(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together, including mass and viscous contributions
  // ---------------------------------------------------------------------------
  AssembleForce(ok, *me_vec_ptr_, 1.0, f);

  // ---------------------------------------------------------------------------
  // subtract mass and viscous contributions from initial force vector
  // ---------------------------------------------------------------------------
  f.Scale(-1.);
  int_ptr_->AddViscoMassContributions(f);
  f.Scale(-1.);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x) const
{
  // default ResetStates call
  ResetStates(x, true);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x, bool setstate) const
{
  ResetStates(x, setstate, *me_vec_ptr_);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x, bool setstate, Vector& me_vec) const
{
  if (setstate) int_ptr_->SetState(x);
  for (auto& me_iter : me_vec) me_iter->Reset(x);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForce(
    const Epetra_Vector& x, Epetra_Vector& f, const double& timefac_np) const
{
  CheckInitSetup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);
  // update the state variables of the current time integrator
  int_ptr_->SetState(x);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  EvaluateForce(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  AssembleForce(ok, *me_vec_ptr_, timefac_np, f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(
    const Epetra_Vector& x, LINALG::SparseOperator& jac, const double& timefac_np) const
{
  CheckInitSetup();
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
  EvaluateStiff(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  AssembleJacobian(ok, *me_vec_ptr_, timefac_np, jac);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(const INPAR::STR::ModelType& mt, const Epetra_Vector& x,
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  CheckInitSetup();
  bool ok = true;
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> model_ptr = me_map_ptr_->at(mt);
  const Vector me_vec(1, model_ptr);

  // initialize stiffness matrix to zero
  jac.Zero();

  // update the state variables of the current time integrator
  int_ptr_->SetState(x);
  model_ptr->Reset(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  EvaluateStiff(ok, me_vec);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  AssembleJacobian(ok, me_vec, timefac_np, jac);

  return ok;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForceStiff(const Epetra_Vector& x, Epetra_Vector& f,
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  CheckInitSetup();
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
  EvaluateForceStiff(ok, *me_vec_ptr_);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  AssembleForce(ok, *me_vec_ptr_, timefac_np, f);
  AssembleJacobian(ok, *me_vec_ptr_, timefac_np, jac);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyCheapSOCRhs(const enum NOX::NLN::CorrectionType type,
    const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
    Epetra_Vector& f, const double& timefac_np) const
{
  CheckInitSetup();

  Vector constraint_me_vec;
  ExtractModelVector(constraint_me_vec, constraint_models);

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
  EvaluateCheapSOCRhs(ok, constraint_me_vec);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  AssembleCheapSOCRhs(ok, constraint_me_vec, timefac_np, f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::EvaluateCheapSOCRhs(bool& ok, const Vector& me_vec) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->EvaluateCheapSOCRhs() : false);

  PostEvaluate(ok, me_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::AssembleCheapSOCRhs(
    bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const
{
  if (not ok) return;

  for (Vector::const_iterator cit = me_vec.begin(); cit != me_vec.end(); ++cit)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*cit)->AssembleCheapSOCRhs(f, timefac_np) : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::CorrectParameters(const enum NOX::NLN::CorrectionType type) const
{
  bool ok = true;
  for (auto& cit : *me_vec_ptr_)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? cit->CorrectParameters(type) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Predict(const INPAR::STR::PredEnum& pred_type) const
{
  CheckInitSetup();
  for (Vector::iterator me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->Predict(pred_type);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->WriteRestart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->ReadRestart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPostComputeX(const Epetra_Vector& xold, const Epetra_Vector& dir,
    const double& step, const Epetra_Vector& xnew, const bool isdefaultstep) const
{
  CheckInitSetup();
  // set some parameters for the element evaluation
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);
  eval_data_ptr_->ResetMyNorms(isdefaultstep);
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPostComputeX(xold, dir, xnew);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
    const double& step, const NOX::NLN::Group& curr_grp, const bool isdefaultstep) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPreComputeX(xold, dir_mutable, curr_grp);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPostIterate(const NOX::Solver::Generic& solver, const double step,
    const bool isdefaultstep, const int num_corrs) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);
  eval_data_ptr_->SetNumberOfModifiedNewtonCorrections(num_corrs);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPostIterate(solver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPreSolve(
    const NOX::Solver::Generic& solver, const double step, const bool isdefaultstep) const
{
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);

  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPreSolve(solver);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPostApplyJacobianInverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp) const
{
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPostApplyJacobianInverse(rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RunPreApplyJacobianInverse(const Epetra_Vector& rhs,
    Epetra_Vector& result, const Epetra_Vector& xold, const NOX::NLN::Group& grp) const
{
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RunPreApplyJacobianInverse(rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::ModelEvaluator::GetGlobalState() const
{
  CheckInitSetup();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::ModelEvaluator::GlobalStatePtr()
{
  CheckInitSetup();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<const STR::TIMINT::Base>& STR::ModelEvaluator::GetTimIntPtr() const
{
  CheckInitSetup();
  return timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::ModelEvaluator::Evaluator(const enum INPAR::STR::ModelType& mt)
{
  CheckInitSetup();
  // sanity check, if there is a model evaluator for the given model type
  STR::ModelEvaluator::Map::const_iterator me_iter = me_map_ptr_->find(mt);
  if (me_iter == me_map_ptr_->end())
    dserror("There is no model evaluator for the model type %s",
        INPAR::STR::ModelTypeString(mt).c_str());

  return *(me_iter->second);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::ModelEvaluator::Evaluator(
    const enum INPAR::STR::ModelType& mt) const
{
  CheckInitSetup();
  // sanity check, if there is a model evaluator for the given model type
  STR::ModelEvaluator::Map::const_iterator me_iter = me_map_ptr_->find(mt);
  if (me_iter == me_map_ptr_->end())
    dserror("There is no model evaluator for the model type %s",
        INPAR::STR::ModelTypeString(mt).c_str());

  return *(me_iter->second);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepState(const double& timefac_n)
{
  CheckInitSetup();
  /* Reset old structural right hand side.
   * It will be filled within the model evaluators */
  gstate_ptr_->GetMutableFstructureOld()->Scale(0.0);
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->UpdateStepState(timefac_n);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ComputeJacobianContributionsFromElementLevelForPTC(
    Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr)
{
  // evaluate ptc contributions at t^n+1
  double timefac_np = 1.0;
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->EvaluateJacobianContributionsFromElementLevelForPTC();

  AssembleJacobianContributionsFromElementLevelForPTC(*me_vec_ptr_, timefac_np, scalingMatrixOpPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) const
{
  CheckInitSetup();
  for (auto& me : *me_vec_ptr_) me->RemoveCondensedContributionsFromRhs(rhs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepElement()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineStressStrain()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->DetermineStressStrain();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineEnergy()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->DetermineEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineOptionalQuantity()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->DetermineOptionalQuantity();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();
  Vector::const_iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RuntimeOutputStepState() const
{
  CheckInitSetup();
  Vector::const_iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->RuntimeOutputStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::PostOutput()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStepState()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter)
    (*me_iter)->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Vector> STR::ModelEvaluator::TransformToVector(
    const STR::ModelEvaluator::Map& model_map) const
{
  Teuchos::RCP<STR::ModelEvaluator::Vector> me_vec_ptr =
      Teuchos::rcp(new STR::ModelEvaluator::Vector(0));
  me_vec_ptr->reserve(model_map.size());

  // --------------------------------------------------------------------------
  // There must be a structural model evaluator at the first position
  // --------------------------------------------------------------------------
  if (model_map.begin()->first != INPAR::STR::model_structure)
    dserror(
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
void STR::ModelEvaluator::SplitModelVector(STR::ModelEvaluator::Vector& partial_me_vec,
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
void STR::ModelEvaluator::ExtractModelVector(STR::ModelEvaluator::Vector& partial_me_vec,
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

    if (cit == me_vec_ptr_->end()) dserror("Couldn't find the model type in me_vec_ptr_.");

    partial_me_vec.push_back(*cit);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::CreateBackupState(const Epetra_Vector& dir)
{
  CheckInitSetup();
  for (const auto& me_iter : *me_vec_ptr_) me_iter->CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RecoverFromBackupState()
{
  CheckInitSetup();
  for (const auto& me_iter : *me_vec_ptr_) me_iter->RecoverFromBackupState();
}
