/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_factory.cpp

\maintainer Michael Hiermeier

\date Sep 10, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_factory.H"
#include "../drt_inpar/inpar_structure.H"

// supported model evaluators
#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_springdashpot.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Map> STR::MODELEVALUATOR::Factory::
    BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes
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
      case INPAR::STR::model_meshtying:
      case INPAR::STR::model_lag_pen_constraint:
      case INPAR::STR::model_windkessel:
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
    BuildModelEvaluators(const std::set<enum INPAR::STR::ModelType>& modeltypes
    )
{
  Factory factory;
  return factory.BuildModelEvaluators(modeltypes);
}
