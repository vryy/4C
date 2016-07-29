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
#include "str_model_evaluator_generic.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_blocksparsematrix.H" // debugging

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator::ModelEvaluator()
    : isinit_(false),
      issetup_(false),
      modeltypes_ptr_(Teuchos::null),
      me_map_ptr_(Teuchos::null),
      me_vec_ptr_(Teuchos::null),
      eval_data_ptr_(Teuchos::null),
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
void STR::ModelEvaluator::Init(const std::set<enum INPAR::STR::ModelType>& mt,
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  issetup_ = false;

  modeltypes_ptr_ = Teuchos::rcp(&mt,false);
  eval_data_ptr_ = eval_data_ptr;
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
      STR::MODELEVALUATOR::BuildModelEvaluators(*modeltypes_ptr_);
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
void STR::ModelEvaluator::Reset(const Epetra_Vector& x)
{
  Vector::iterator me_iter;
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
  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
  {
    // if one model evaluator failed, skip the remaining ones and return false
    if (not ok) break;

    (*me_iter)->Reset(x);
    ok = (*me_iter)->EvaluateForce();
  }

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
  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
  {
    // if one model evaluator failed, skip the remaining ones and return false
    if (not ok) break;

    (*me_iter)->Reset(x);
    ok = (*me_iter)->EvaluateStiff();
  }

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
  // ---------------------------------------------------------------------------
  // evaluate all terms
  // ---------------------------------------------------------------------------
  for (me_iter=me_vec_ptr_->begin();me_iter!=me_vec_ptr_->end();++me_iter)
  {
    // if one model evaluator failed, skip the remaining ones and return false
    if (not ok) break;

    (*me_iter)->Reset(x);
    ok = (*me_iter)->EvaluateForceStiff();
  }

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
