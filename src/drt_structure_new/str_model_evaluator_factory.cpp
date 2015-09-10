/*
 * str_model_evaluator_factory.cpp
 *
 *  Created on: Sep 10, 2015
 *      Author: hiermeier
 */


#include "str_model_evaluator_factory.H"

// supported model evaluators
#include "str_model_evaluator_structure.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > >
STR::MODELEVALUATOR::Factory::BuildModelEvaluators(
    const std::vector<const enum INPAR::STR::ModelType>& modeltypes
    ) const
{
  // create a new standard map
  Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > models =
      Teuchos::rcp(new std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >());

  std::vector<const enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter=modeltypes.begin();mt_iter!=modeltypes.end();++mt_iter)
  {
    switch(*mt_iter)
    {
      case INPAR::STR::model_structure:
        (*models)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Structure());
        break;
      case INPAR::STR::model_meshtying_contact:
      case INPAR::STR::model_lag_pen_constraint:
      case INPAR::STR::model_windkessel:
      case INPAR::STR::model_springdashpot:
      default:
        dserror("Not yet implemented!");
        break;
    }
  }


  return models;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > >
STR::MODELEVALUATOR::BuildModelEvaluators(
    const std::vector<const enum INPAR::STR::ModelType>& modeltypes
    )
{
  Factory factory;
  return factory.BuildModelEvaluators(modeltypes);
}
