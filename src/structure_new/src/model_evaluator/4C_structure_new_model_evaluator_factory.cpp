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
Solid::MODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::ModelEvaluator::Map> Solid::MODELEVALUATOR::Factory::build_model_evaluators(
    const std::set<enum Inpar::Solid::ModelType>& modeltypes,
    const Teuchos::RCP<Solid::MODELEVALUATOR::Generic>& coupling_model_ptr) const
{
  // create a new standard map
  Teuchos::RCP<Solid::ModelEvaluator::Map> model_map =
      Teuchos::rcp(new Solid::ModelEvaluator::Map());

  std::set<enum Inpar::Solid::ModelType>::const_iterator mt_iter;
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    switch (*mt_iter)
    {
      case Inpar::Solid::model_structure:
        (*model_map)[*mt_iter] = build_structure_model_evaluator();
        break;
      case Inpar::Solid::model_springdashpot:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::SpringDashpot());
        break;
      case Inpar::Solid::model_browniandyn:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::BrownianDyn());
        break;
      case Inpar::Solid::model_beaminteraction:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::BeamInteraction());
        break;
      case Inpar::Solid::model_contact:
      {
        (*model_map)[*mt_iter] = build_contact_model_evaluator();
        break;
      }
      case Inpar::Solid::model_beam_interaction_old:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::BeamInteractionOld());
        break;
      case Inpar::Solid::model_lag_pen_constraint:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::LagPenConstraint());
        break;
      case Inpar::Solid::model_cardiovascular0d:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::Cardiovascular0D());
        break;
      case Inpar::Solid::model_monolithic_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The monolithic coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case Inpar::Solid::model_partitioned_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The partitioned coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case Inpar::Solid::model_basic_coupling:
      {
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The basic coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case Inpar::Solid::model_meshtying:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::Meshtying());
        break;
      case Inpar::Solid::model_constraints:
        (*model_map)[*mt_iter] = Teuchos::rcp(new Solid::MODELEVALUATOR::Constraints());
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
Teuchos::RCP<Solid::MODELEVALUATOR::Generic>
Solid::MODELEVALUATOR::Factory::build_contact_model_evaluator() const
{
  return Teuchos::rcp(new Solid::MODELEVALUATOR::Contact());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::MODELEVALUATOR::Generic>
Solid::MODELEVALUATOR::Factory::build_structure_model_evaluator() const
{
  return Teuchos::rcp(new Solid::MODELEVALUATOR::Structure());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::ModelEvaluator::Map> Solid::MODELEVALUATOR::build_model_evaluators(
    const std::set<enum Inpar::Solid::ModelType>& modeltypes,
    const Teuchos::RCP<Solid::MODELEVALUATOR::Generic>& coupling_model_ptr)
{
  Factory factory;
  return factory.build_model_evaluators(modeltypes, coupling_model_ptr);
}

FOUR_C_NAMESPACE_CLOSE
