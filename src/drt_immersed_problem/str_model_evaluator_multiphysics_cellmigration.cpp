/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_multiphysics_cellmigration.cpp

\brief Generic class for all model evaluators.

\maintainer Andreas Rauch

\date Nov 28, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_multiphysics_cellmigration.H"

#include "../drt_fsi/fsi_str_model_evaluator_partitioned.H"
#include "../drt_ssi/ssi_str_model_evaluator_partitioned.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::CellMigration::CellMigration()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::CellMigration::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr,
    const int& dof_offset)
{
  STR::MODELEVALUATOR::Multiphysics::Init(eval_data_ptr,gstate_ptr,gio_ptr,int_ptr,timint_ptr,dof_offset);

  // construct model evaluators of sub modules
  Teuchos::RCP<STR::MODELEVALUATOR::PartitionedFSI>
      fsi_model_evaluator = Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedFSI());
  Teuchos::RCP<STR::MODELEVALUATOR::PartitionedSSI>
      ssi_model_evaluator = Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedSSI());

  // initialize model evaluators of sub modules
  fsi_model_evaluator -> Init(eval_data_ptr,gstate_ptr,gio_ptr,int_ptr,timint_ptr,dof_offset);
  ssi_model_evaluator -> Init(eval_data_ptr,gstate_ptr,gio_ptr,int_ptr,timint_ptr,dof_offset);

  // register model evaluators in map
  GetModelEvalutaorMap().
      insert( std::pair<enum STR::MODELEVALUATOR::MultiphysicType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >(
          mt_fsi, fsi_model_evaluator) );
  GetModelEvalutaorMap().
      insert( std::pair<enum STR::MODELEVALUATOR::MultiphysicType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >(
          mt_ssi, ssi_model_evaluator) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::CellMigration::Setup()
{
  STR::MODELEVALUATOR::Multiphysics::Setup();
}
