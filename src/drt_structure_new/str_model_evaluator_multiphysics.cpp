/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_multiphysics.cpp

\brief Generic class for all model evaluators.

\maintainer Andreas Rauch

\date Nov 28, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_multiphysics.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Multiphysics::Multiphysics()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr,
    const int& dof_offset)
{
  STR::MODELEVALUATOR::Generic::Init(eval_data_ptr,gstate_ptr,gio_ptr,int_ptr,timint_ptr,dof_offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Setup()
{

}

