/*-----------------------------------------------------------*/
/*! \file


\brief Factory class to build predictor objects

\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_predict_factory.hpp"

// supported predictor classes
#include "4C_structure_new_predict_constdisvelaccpress.hpp"
#include "4C_structure_new_predict_tangdis.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Predict::Factory::Factory()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Predict::Generic> Solid::Predict::Factory::build_predictor(
    const enum Inpar::Solid::PredEnum& predType) const
{
  Teuchos::RCP<Solid::Predict::Generic> predictor = Teuchos::null;

  switch (predType)
  {
    case Inpar::Solid::pred_constdis:
    case Inpar::Solid::pred_constvel:
    case Inpar::Solid::pred_constacc:
    case Inpar::Solid::pred_constdisvelacc:
    case Inpar::Solid::pred_constdispres:
    case Inpar::Solid::pred_constdisvelaccpres:
      predictor = Teuchos::rcp(new Solid::Predict::ConstDisVelAccPress());
      break;
    case Inpar::Solid::pred_tangdis:
    case Inpar::Solid::pred_tangdis_constfext:
      predictor = Teuchos::rcp(new Solid::Predict::TangDis());
      break;
    case Inpar::Solid::pred_vague:
    default:
      FOUR_C_THROW("Unknown predictor type!");
      break;
  }

  return predictor;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Predict::Generic> Solid::Predict::build_predictor(
    const enum Inpar::Solid::PredEnum& predType)
{
  Factory factory;
  return factory.build_predictor(predType);
}

FOUR_C_NAMESPACE_CLOSE
