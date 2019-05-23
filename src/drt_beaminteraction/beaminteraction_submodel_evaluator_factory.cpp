/*-----------------------------------------------------------*/
/*!

\brief Factory to create the desired subemodel evaluators.

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#include "beaminteraction_submodel_evaluator_factory.H"

#include "../drt_inpar/inpar_beaminteraction.H"

// supported submodel evaluators
#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_beamcontact.H"

// problem types
#include "../drt_lib/drt_globalproblem.H"
#include "beaminteraction_submodel_evaluator_crosslinking.H"
#include "beaminteraction_submodel_evaluator_potential.H"
#include "beaminteraction_submodel_evaluator_spherebeamlinking.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::Factory::Factory()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map>
BEAMINTERACTION::SUBMODELEVALUATOR::Factory::BuildModelEvaluators(
    const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes) const
{
  // create a new standard map
  Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> model_map =
      Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteraction::Map());

  std::set<enum INPAR::BEAMINTERACTION::SubModelType>::const_iterator mt_iter;
  for (mt_iter = submodeltypes.begin(); mt_iter != submodeltypes.end(); ++mt_iter)
  {
    switch (*mt_iter)
    {
      case INPAR::BEAMINTERACTION::submodel_beamcontact:
        (*model_map)[*mt_iter] =
            Teuchos::rcp(new BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact());
        break;
      case INPAR::BEAMINTERACTION::submodel_crosslinking:
        (*model_map)[*mt_iter] =
            Teuchos::rcp(new BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking());
        break;
      case INPAR::BEAMINTERACTION::submodel_spherebeamlink:
        (*model_map)[*mt_iter] =
            Teuchos::rcp(new BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking());
        break;
      case INPAR::BEAMINTERACTION::submodel_potential:
        (*model_map)[*mt_iter] =
            Teuchos::rcp(new BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential());
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
Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map>
BEAMINTERACTION::SUBMODELEVALUATOR::BuildModelEvaluators(
    const std::set<enum INPAR::BEAMINTERACTION::SubModelType>& submodeltypes)
{
  Factory factory;
  return factory.BuildModelEvaluators(submodeltypes);
}
