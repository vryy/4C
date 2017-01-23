/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_browniandyndata.cpp

\brief Concrete implementation of the brownian dynamic parameter interface

\maintainer Jonas Eichinger

\date Jun 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure_new/str_timint_base.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::BrownianDynData::BrownianDynData()
    : isinit_(false),
      issetup_(false),
      str_data_ptr_(Teuchos::null),
      viscosity_(0.0),
      kt_(0.0),
      maxrandforce_(0.0),
      timeintconstrandnumb_(0.0),
      randomforces_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::Init(
    const Teuchos::RCP<const STR::MODELEVALUATOR::Data>& str_data_ptr)
{
  issetup_ = false;

  str_data_ptr_ = str_data_ptr;

  const Teuchos::ParameterList& brwondyn_params_list =
      DRT::Problem::Instance()->BrownianDynamicsParams();

  // viscosity
  viscosity_ =  brwondyn_params_list.get<double> ("VISCOSITY");
  // thermal energy
  kt_ = brwondyn_params_list.get<double> ("KT");
  // maximum random force
  maxrandforce_ = brwondyn_params_list.get<double> ("MAXRANDFORCE");
  // time intervall with constant random forces
  timeintconstrandnumb_ = brwondyn_params_list.get<double> ("TIMEINTCONSTRANDFORCES");

  // set flag
  isinit_ = true;

  return;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::Setup()
{
  CheckInit();

  // set flag
  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BrownianDynData::ResizeRandomForceMVector(
    Teuchos::RCP<DRT::Discretization> discret_ptr,
    int                               maxrandnumelement)
{
  CheckInitSetup();

  // resize in case of new crosslinkers that were set and are now part of the discretization
  randomforces_= Teuchos::rcp( new Epetra_MultiVector(
      *(discret_ptr->ElementColMap()),maxrandnumelement,true) );

  return;
}
