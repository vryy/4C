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
STR::Predict::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Predict::Generic> STR::Predict::Factory::build_predictor(
    const enum Inpar::STR::PredEnum& predType) const
{
  Teuchos::RCP<STR::Predict::Generic> predictor = Teuchos::null;

  switch (predType)
  {
    case Inpar::STR::pred_constdis:
    case Inpar::STR::pred_constvel:
    case Inpar::STR::pred_constacc:
    case Inpar::STR::pred_constdisvelacc:
    case Inpar::STR::pred_constdispres:
    case Inpar::STR::pred_constdisvelaccpres:
      predictor = Teuchos::rcp(new STR::Predict::ConstDisVelAccPress());
      break;
    case Inpar::STR::pred_tangdis:
    case Inpar::STR::pred_tangdis_constfext:
      predictor = Teuchos::rcp(new STR::Predict::TangDis());
      break;
    case Inpar::STR::pred_vague:
    default:
      FOUR_C_THROW("Unknown predictor type!");
      break;
  }

  return predictor;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Predict::Generic> STR::Predict::build_predictor(
    const enum Inpar::STR::PredEnum& predType)
{
  Factory factory;
  return factory.build_predictor(predType);
}

FOUR_C_NAMESPACE_CLOSE
