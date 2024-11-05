// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_factory.hpp"

#include "4C_adapter_str_structure_new.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::StructureFactory::StructureFactory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Adapter::StructureBaseAlgorithmNew>
Adapter::StructureFactory::build_structure_algorithm(const Teuchos::ParameterList& sdyn) const
{
  std::shared_ptr<Adapter::StructureBaseAlgorithmNew> adapterbase = nullptr;

  const auto intstrat =
      Teuchos::getIntegralValue<Inpar::Solid::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case Inpar::Solid::int_standard:
      adapterbase = std::make_shared<Adapter::StructureBaseAlgorithmNew>();
      break;
    default:
      FOUR_C_THROW("Unknown integration strategy!");
      break;
  }

  return adapterbase;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Adapter::StructureBaseAlgorithmNew> Adapter::build_structure_algorithm(
    const Teuchos::ParameterList& sdyn)
{
  StructureFactory factory;
  return factory.build_structure_algorithm(sdyn);
}

FOUR_C_NAMESPACE_CLOSE
