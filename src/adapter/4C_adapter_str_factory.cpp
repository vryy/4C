/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_adapter_str_factory.hpp"

#include "4C_adapter_str_structure_new.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::StructureFactory::StructureFactory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Adapter::StructureBaseAlgorithmNew>
Adapter::StructureFactory::build_structure_algorithm(const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<Adapter::StructureBaseAlgorithmNew> adapterbase = Teuchos::null;

  const enum Inpar::Solid::IntegrationStrategy intstrat =
      Core::UTILS::IntegralValue<Inpar::Solid::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case Inpar::Solid::int_standard:
      adapterbase = Teuchos::rcp(new Adapter::StructureBaseAlgorithmNew());
      break;
    default:
      FOUR_C_THROW("Unknown integration strategy!");
      break;
  }

  return adapterbase;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Adapter::StructureBaseAlgorithmNew> Adapter::build_structure_algorithm(
    const Teuchos::ParameterList& sdyn)
{
  StructureFactory factory;
  return factory.build_structure_algorithm(sdyn);
}

FOUR_C_NAMESPACE_CLOSE
