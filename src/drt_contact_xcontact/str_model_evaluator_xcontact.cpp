/*----------------------------------------------------------------------------*/
/**
\file str_model_evaluator_xcontact.cpp

\brief Model evaluator for the eXtended contact formulation

\maintainer Matthias Mayr

\date Jul 6, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "str_model_evaluator_xcontact.H"
#include "../drt_structure_new/str_model_evaluator_data.H"

#include "xcontact_multi_discretization_wrapper.H"
#include "xcontact_strategy.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::XContact::XContact()
    : xstrategy_ptr_(Teuchos::null), levelset_values_ptr_(Teuchos::null)
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::XContact::Setup()
{
  // call the setup routine of the base class for now
  Contact::Setup();
  xstrategy_ptr_ = Teuchos::rcp_dynamic_cast<XCONTACT::Strategy>(Contact::StrategyPtr(), true);

  /* set the pointer to the created interface discretizations in the
   * discretization wrapper object */
  XCONTACT::MultiDiscretizationWrapper& xstr_discret_wrapper =
      dynamic_cast<XCONTACT::MultiDiscretizationWrapper&>(Discret());
  xstr_discret_wrapper.AddContactIDiscret(xstrategy_ptr_->GetIDiscrets());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::XContact::CheckPseudo2D() const
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::XContact::SetContactStatus(const bool& is_in_contact)
{
  XStrategy().SetContactStatus(is_in_contact);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::XContact::EvaluateWeightedGap()
{
  CheckInitSetup();
  // --- evaluate weighted gap values ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_weighted_gap);
  Strategy().Evaluate(EvalData().Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::XContact::GetWeightedGap() const
{
  return XStrategy().GetWeightedGap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const XCONTACT::Strategy& STR::MODELEVALUATOR::XContact::XStrategy() const
{
  CheckInitSetup();
  return *xstrategy_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::Strategy& STR::MODELEVALUATOR::XContact::XStrategy()
{
  CheckInitSetup();
  return *xstrategy_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::XContact::SetLevelSetValuesPtr(
    const Teuchos::RCP<Epetra_Vector>& levelset_values_ptr)
{
  levelset_values_ptr_ = levelset_values_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::XContact::LevelSetValues() const
{
  if (levelset_values_ptr_.is_null())
    dserror("Set level set values are not set correctly! (NULL ptr)");

  return *levelset_values_ptr_;
}
