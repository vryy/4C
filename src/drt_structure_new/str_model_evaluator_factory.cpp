/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired model evaluators.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_factory.H"

#include "../drt_beamcontact/str_model_evaluator_beaminteraction_old.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction.H"
#include "../drt_browniandyn/str_model_evaluator_browniandyn.H"
#include "../drt_inpar/inpar_structure.H"

// supported model evaluators
#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_cardiovascular0d.H"
#include "str_model_evaluator_springdashpot.H"
#include "str_model_evaluator_contact.H"
#include "str_model_evaluator_meshtying.H"
#include "str_model_evaluator_lagpenconstraint.H"
#include "../drt_contact_xcontact/str_model_evaluator_xcontact.H"
#include "../drt_struct_ale/struct_ale_str_model_evaluator.H"

// problem types
#include "../drt_lib/drt_globalproblem.H"


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
          dserror("The monolithic coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_partitioned_coupling:
      {
        if (coupling_model_ptr.is_null())
          dserror("The partitioned coupling model evaluator is not defined.");
        (*model_map)[*mt_iter] = coupling_model_ptr;
        break;
      }
      case INPAR::STR::model_meshtying:
        (*model_map)[*mt_iter] = Teuchos::rcp(new STR::MODELEVALUATOR::Meshtying());
        break;
      default:
        dserror("Not yet implemented!");
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
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> contact_model_ptr = Teuchos::null;
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  switch (probtype)
  {
    case prb_xcontact:
    {
      contact_model_ptr = Teuchos::rcp(new STR::MODELEVALUATOR::XContact());
      break;
    }
    default:
    {
      contact_model_ptr = Teuchos::rcp(new STR::MODELEVALUATOR::Contact());
      break;
    }
  }
  return contact_model_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::Generic>
STR::MODELEVALUATOR::Factory::BuildStructureModelEvaluator() const
{
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> structure_model_ptr = Teuchos::null;
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  switch (probtype)
  {
    case prb_struct_ale:
    case prb_immersed_cell:
    {
      structure_model_ptr = Teuchos::rcp(new STR::MODELEVALUATOR::StructAle());
      break;
    }
    default:
    {
      structure_model_ptr = Teuchos::rcp(new STR::MODELEVALUATOR::Structure());
      break;
    }
  }
  return structure_model_ptr;
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
