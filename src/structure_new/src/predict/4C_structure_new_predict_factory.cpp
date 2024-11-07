// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
std::shared_ptr<Solid::Predict::Generic> Solid::Predict::Factory::build_predictor(
    const enum Inpar::Solid::PredEnum& predType) const
{
  std::shared_ptr<Solid::Predict::Generic> predictor = nullptr;

  switch (predType)
  {
    case Inpar::Solid::pred_constdis:
    case Inpar::Solid::pred_constvel:
    case Inpar::Solid::pred_constacc:
    case Inpar::Solid::pred_constdisvelacc:
    case Inpar::Solid::pred_constdispres:
    case Inpar::Solid::pred_constdisvelaccpres:
      predictor = std::make_shared<Solid::Predict::ConstDisVelAccPress>();
      break;
    case Inpar::Solid::pred_tangdis:
    case Inpar::Solid::pred_tangdis_constfext:
      predictor = std::make_shared<Solid::Predict::TangDis>();
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
std::shared_ptr<Solid::Predict::Generic> Solid::Predict::build_predictor(
    const enum Inpar::Solid::PredEnum& predType)
{
  Factory factory;
  return factory.build_predictor(predType);
}

FOUR_C_NAMESPACE_CLOSE
