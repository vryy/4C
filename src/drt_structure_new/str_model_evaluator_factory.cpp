/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_factory.cpp

\brief Factory to create the desired model evaluators.

\maintainer Michael Hiermeier

\date Sep 10, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_factory.H"
#include "../drt_inpar/inpar_structure.H"

// supported model evaluators
#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_cardiovascular0d.H"
#include "str_model_evaluator_springdashpot.H"
#include "str_model_evaluator_contact.H"
#include "str_model_evaluator_lagpenconstraint.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Map> STR::MODELEVALUATOR::Factory::
    BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes,
        const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr
        ) const
{
  // create a new standard map
  Teuchos::RCP<STR::ModelEvaluator::Map> model_map =
      Teuchos::rcp(new STR::ModelEvaluator::Map());

  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter=modeltypes.begin();mt_iter!=modeltypes.end();++mt_iter)
  {
    switch(*mt_iter)
    {
      case INPAR::STR::model_structure:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Structure());
        break;
      case INPAR::STR::model_springdashpot:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::SpringDashpot());
        break;
      case INPAR::STR::model_contact:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Contact());
        break;
      case INPAR::STR::model_lag_pen_constraint:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::LagPenConstraint());
        break;
      case INPAR::STR::model_cardiovascular0d:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Cardiovascular0D());
        break;
      case INPAR::STR::model_partitioned_coupling:
      {
        if (coupling_model_ptr.is_null())
          dserror("The partitioned coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_meshtying:
      default:
        dserror("Not yet implemented!");
        break;
    }
  }

  return model_map;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Map> STR::MODELEVALUATOR::
    BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes,
        const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr)
{
  Factory factory;
  return factory.BuildModelEvaluators(modeltypes,coupling_model_ptr);
}
