/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator.cpp

\brief Manager of the model evaluator calls.

\maintainer Michael Hiermeier

\date Nov 30, 2015

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
#include "../linalg/linalg_blocksparsematrix.H" // debugging

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
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
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

  me_map_ptr_ =
      STR::MODELEVALUATOR::BuildModelEvaluators(sdyn_ptr_->GetModelTypes(),
          sdyn_ptr_->CouplingModelPtr());
  std::vector<enum INPAR::STR::ModelType> sorted_modeltypes(0);

  me_vec_ptr_ = Sort(*me_map_ptr_,sorted_modeltypes);
  Vector::iterator me_iter;
  int dof_offset = 0;
  unsigned int i = 0;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
  {
    (*me_iter)->Init(eval_data_ptr_,gstate_ptr_,gio_ptr_,int_ptr_,
        timint_ptr_,dof_offset);
    (*me_iter)->Setup();
    // setup the block information for saddle point problems
    dof_offset = gstate_ptr_->SetupBlockInformation(**me_iter,
        sorted_modeltypes[i]);
    ++i;
  }
  gstate_ptr_->SetupMultiMapExtractor();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::SetupMultiMapExtractor()
{
  // setup the block information for saddle point problems
  for (Vector::iterator me_iter=me_vec_ptr_->begin();
      me_iter!=me_vec_ptr_->end();++me_iter)
    gstate_ptr_->SetupBlockInformation(**me_iter,
        (**me_iter).Type());

  gstate_ptr_->SetupMultiMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::InitializeInertiaAndDamping(const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();

  // initialize stiffness matrix to zero
  jac.Zero();
  // get structural model evaluator
  STR::MODELEVALUATOR::Structure& str_model = dynamic_cast<
      STR::MODELEVALUATOR::Structure&>(Evaluator(INPAR::STR::model_structure));

  str_model.Reset(x);
  bool ok = str_model.InitializeInertiaAndDamping();

  // if the model evaluator failed, skip the assembly and return false
  ok = (ok ? str_model.AssembleJacobian(jac,1.0) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyInitialForce(const Epetra_Vector& x,
    Epetra_Vector& f)
{
  CheckInitSetup();

  Vector::iterator me_iter;
  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x,false);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->EvaluateInitialForce() : false);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->AssembleForce(f,1.0) : false);

  return ok;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x) const
{
  // default ResetStates call
  ResetStates(x,true);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStates(const Epetra_Vector& x, bool setstate) const
{
  Vector::iterator me_iter;
  if(setstate)
    int_ptr_->SetState(x);
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->Reset(x);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForce(const Epetra_Vector& x, Epetra_Vector& f,
    const double& timefac_np)
    const
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
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->EvaluateForce() : false);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->AssembleForce(f,timefac_np) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(const Epetra_Vector& x,
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  CheckInitSetup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix to zero
  jac.Zero();
  // update the state variables of the current time integrator
  int_ptr_->SetState(x);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->EvaluateStiff() : false);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->AssembleJacobian(jac,timefac_np) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(
    const INPAR::STR::ModelType& mt,
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac,
    const double& timefac_np) const
{
  CheckInitSetup();
  bool ok = true;
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> model_ptr = me_map_ptr_->at(mt);
  // initialize stiffness matrix to zero
  jac.Zero();
  // update the state variables of the current time integrator
  int_ptr_->SetState(x);
  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  model_ptr->Reset(x);
  ok = model_ptr->EvaluateStiff();

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  ok = (ok ? model_ptr->AssembleJacobian(jac,timefac_np) : false);

  return ok;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForceStiff(const Epetra_Vector& x,
    Epetra_Vector& f, LINALG::SparseOperator& jac,
    const double& timefac_np) const
{
  CheckInitSetup();
  Vector::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix and right hand side to zero
  f.PutScalar(0.0);
  jac.Zero();
  // update the state variables of the current time integrator
  int_ptr_->SetState(x);

  // ---------------------------------------------------------------------------
  // reset model specific variables
  // ---------------------------------------------------------------------------
  ResetStates(x);

  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->EvaluateForceStiff() : false);

  // ---------------------------------------------------------------------------
  // put everything together
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
  {
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? (*me_iter)->AssembleForce(f,timefac_np) : false);
    ok = (ok ? (*me_iter)->AssembleJacobian(jac,timefac_np) : false);
  }

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Predict(const INPAR::STR::PredEnum& pred_type) const
{
  CheckInitSetup();
  for (Vector::iterator me_iter=me_vec_ptr_->begin();
      me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->Predict(pred_type);

  return ;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->WriteRestart(iowriter,forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->ReadRestart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const double& step,
    const Epetra_Vector& xnew,
    const bool& isdefaultstep) const
{
  CheckInitSetup();
  // set some parameters for the element evaluation
  eval_data_ptr_->SetIsDefaultStep(isdefaultstep);
  eval_data_ptr_->SetStepLength(step);
  eval_data_ptr_->ResetMyNorms(isdefaultstep);
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin(); me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->RecoverState(xold,dir,xnew);

  return;
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
const Teuchos::RCP<const STR::TIMINT::Base>& STR::ModelEvaluator::GetTimIntPtr()
const
{
  CheckInitSetup();
  return timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::ModelEvaluator::Evaluator(
    const enum INPAR::STR::ModelType& mt)
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
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->UpdateStepState(timefac_n);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepElement()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineStressStrain()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->DetermineStressStrain();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::DetermineEnergy()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->DetermineEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::OutputStepState(IO::DiscretizationWriter& iowriter)
    const
{
  CheckInitSetup();
  Vector::const_iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::PostOutput()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::ResetStepState()
{
  CheckInitSetup();
  Vector::iterator me_iter;
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
    (*me_iter)->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Vector> STR::ModelEvaluator::Sort(
    STR::ModelEvaluator::Map model_map,
    std::vector<INPAR::STR::ModelType>& sorted_model_types
    ) const
{
  Teuchos::RCP<STR::ModelEvaluator::Vector> me_vec_ptr =
      Teuchos::rcp(new STR::ModelEvaluator::Vector(0));

  STR::ModelEvaluator::Map::iterator miter;
  // --------------------------------------------------------------------------
  // There has to be a structural model evaluator and we put it at first place
  // --------------------------------------------------------------------------
  miter = model_map.find(INPAR::STR::model_structure);
  if (miter==model_map.end())
    dserror("The structural model evaluator could not be found!");
  me_vec_ptr->push_back(miter->second);
  sorted_model_types.push_back(miter->first);

  // erase the structural model evaluator
  model_map.erase(miter);

  // --------------------------------------------------------------------------
  // insert the remaining model evaluators into the model vector
  // --------------------------------------------------------------------------
  for (miter=model_map.begin();miter!=model_map.end();++miter)
  {
    me_vec_ptr->push_back(miter->second);
    sorted_model_types.push_back(miter->first);
  }

  return me_vec_ptr;
}
