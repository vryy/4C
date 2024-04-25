/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired model evaluators.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_model_evaluator_factory.hpp"

#include "4C_beamcontact_str_model_evaluator_beaminteraction_old.hpp"
#include "4C_beaminteraction_str_model_evaluator.hpp"
#include "4C_browniandyn_str_model_evaluator.hpp"
#include "4C_cardiovascular0d_structure_new_model_evaluator.hpp"
#include "4C_constraint_framework_model_evaluator.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_structure_new_model_evaluator_lagpenconstraint.hpp"
#include "4C_structure_new_model_evaluator_meshtying.hpp"
#include "4C_structure_new_model_evaluator_springdashpot.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Map> STR::MODELEVALUATOR::Factory::BuildModelEvaluators(
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr) const
{
  // create a new standard map
  Teuchos::RCP<STR::ModelEvaluator::Map> model_map = Teuchos::rcp(new STR::ModelEvaluator::Map());

  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    switch (*mt_iter)
    {
      case INPAR::STR::model_structure:
        (*model_map)[*mt_iter] = BuildStructureModelEvaluator();
        break;
      case INPAR::STR::model_springdashpot:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::SpringDashpot());
        break;
      case INPAR::STR::model_browniandyn:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::BrownianDyn());
        break;
      case INPAR::STR::model_beaminteraction:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteraction());
        break;
      case INPAR::STR::model_contact:
      {
        (*model_map)[*mt_iter] = BuildContactModelEvaluator();
        break;
      }
      case INPAR::STR::model_beam_interaction_old:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteractionOld());
        break;
      case INPAR::STR::model_lag_pen_constraint:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::LagPenConstraint());
        break;
      case INPAR::STR::model_cardiovascular0d:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Cardiovascular0D());
        break;
      case INPAR::STR::model_monolithic_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The monolithic coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_partitioned_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The partitioned coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_basic_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The basic coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_meshtying:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Meshtying());
        break;
      case INPAR::STR::model_constraints:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Constraints());
        break;
      default:
        FOUR_C_THROW("Not yet implemented!");
        break;
    }
  }

  return model_map;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::Generic>
STR::MODELEVALUATOR::Factory::BuildContactModelEvaluator() const
{
  return Teuchos::rcp(new STR::MODELEVALUATOR::Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::Generic>
STR::MODELEVALUATOR::Factory::BuildStructureModelEvaluator() const
{
  return Teuchos::rcp(new STR::MODELEVALUATOR::Structure());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::ModelEvaluator::Map> STR::MODELEVALUATOR::BuildModelEvaluators(
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr)
{
  Factory factory;
  return factory.BuildModelEvaluators(modeltypes, coupling_model_ptr);
}

FOUR_C_NAMESPACE_CLOSE
