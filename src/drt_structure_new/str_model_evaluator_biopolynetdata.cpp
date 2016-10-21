/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_biopolynetdata.cpp

\brief Concrete implementation of the statmech parameter interface

\maintainer Jonas Eichinger

\date Jun 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "str_model_evaluator_data.H"
#include "../drt_lib/drt_globalproblem.H"
#include "str_timint_base.H"
#include "str_timint_databiopolynetdyn.H"
#include "str_timint_databiopolynetdyn.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::StatMechData::StatMechData()
    : isinit_(false),
      issetup_(false),
      smdyn_ptr_(Teuchos::null),
      randomforces_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 |  Init statmech data                         (public)  eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::StatMechData::Init(
    const Teuchos::RCP<const STR::MODELEVALUATOR::Data>& str_data_ptr,
    const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& sdyn_ptr_)
{
  issetup_ = false;

  smdyn_ptr_ = sdyn_ptr_->GetDataSMDynPtr();

  // set flag
  isinit_ = true;

  return;

} // Init()

/*----------------------------------------------------------------------*
 |  Setup                                      (public)  eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::StatMechData::Setup()
{
  CheckInit();

  // set flag
  issetup_ = true;

  return;
} // Setup()

/*----------------------------------------------------------------------*
 |  Reszize Multivector with random forces     (public)  eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::StatMechData::ResizeRandomForceMVector(
    Teuchos::RCP<DRT::Discretization> discret_ptr,
    int                               maxrandnumelement)
{
  CheckInitSetup();

  // resize in case of new crosslinkers that were set and are now part of the discretization
  randomforces_=
      Teuchos::rcp( new Epetra_MultiVector(*(discret_ptr->ElementColMap()),maxrandnumelement,true) );

  return;
} // ResizeRandomForceMVector()
