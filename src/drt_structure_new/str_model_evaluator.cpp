/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator.cpp

\maintainer Michael Hiermeier

\date Nov 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator.H"
#include "str_model_evaluator_factory.H"
#include "str_model_evaluator_generic.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator::ModelEvaluator()
    : isinit_(false),
      issetup_(false),
      modeltypes_ptr_(Teuchos::null),
      modelevaluators_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      gio_ptr_(Teuchos::null),
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
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  issetup_ = false;

  modeltypes_ptr_ = Teuchos::rcp(&mt,false);
  gstate_ptr_ = gstate_ptr;
  gio_ptr_ = gio_ptr;
  timint_ptr_ = timint_ptr;


  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::Setup()
{
  CheckInit();

  modelevaluators_ptr_ =
      STR::MODELEVALUATOR::BuildModelEvaluators(*modeltypes_ptr_);

  Map::iterator me_iter;
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
  {
    me_iter->second->Init(gstate_ptr_,gio_ptr_,timint_ptr_);
    me_iter->second->Setup();
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForce(const Epetra_Vector& x, Epetra_Vector& f)
    const
{
  CheckInitSetup();
  Map::iterator me_iter;
  bool ok = true;
  // initialize right hand side to zero
  f.PutScalar(0.0);
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? me_iter->second->ApplyForce(x,f) : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyStiff(const Epetra_Vector& x,
    LINALG::SparseOperator& jac) const
{
  CheckInitSetup();
  Map::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix to zero
  jac.Zero();
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? me_iter->second->ApplyStiff(x,jac) : false);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::ModelEvaluator::ApplyForceStiff(const Epetra_Vector& x,
    Epetra_Vector& f, LINALG::SparseOperator& jac) const
{
  CheckInitSetup();
  Map::iterator me_iter;
  bool ok = true;
  // initialize stiffness matrix and right hand side to zero
  f.PutScalar(0.0);
  jac.Zero();
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    // if one model evaluator failed, skip the remaining ones and return false
    ok = (ok ? me_iter->second->ApplyForceStiff(x,f,jac) : false);

  jac.Complete();

  return ok;
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
void STR::ModelEvaluator::UpdateStepState()
{
  CheckInitSetup();
  Map::iterator me_iter;
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    me_iter->second->UpdateStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::UpdateStepElement()
{
  CheckInitSetup();
  Map::iterator me_iter;
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    me_iter->second->UpdateStepElement();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ModelEvaluator::OutputStepState()
{
  CheckInitSetup();
  Map::iterator me_iter;
  for (me_iter=modelevaluators_ptr_->begin();
      me_iter!=modelevaluators_ptr_->end();++me_iter)
    me_iter->second->OutputStepState();
}
