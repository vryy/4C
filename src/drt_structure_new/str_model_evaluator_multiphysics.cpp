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
#include "../drt_immersed_problem/immersed_field_exchange_manager.H"
#include "../drt_structure_new/str_timint_base.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Multiphysics::Multiphysics()
:
active_mt_(mt_none)
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
  // so far only one case
  active_mt_ = mt_fsi;

  STR::MODELEVALUATOR::Generic::Init(eval_data_ptr,gstate_ptr,gio_ptr,int_ptr,timint_ptr,dof_offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Setup()
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::
    AssembleForce(Epetra_Vector& f,
      const double & timefac_np) const
{
  if(active_mt_ == mt_none)
    dserror("No active model evaluator set for Multiphysics");

  GetModelEvaluatorFromMap(active_mt_)->AssembleForce(f,timefac_np);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  STR::MODELEVALUATOR::Multiphysics::
    OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  // here you can write data, e.g.:
  // iowriter.WriteVector("cell_adhesion_force",GState().GetFreactNp());

}

