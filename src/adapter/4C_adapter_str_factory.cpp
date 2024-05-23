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
ADAPTER::StructureFactory::StructureFactory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew>
ADAPTER::StructureFactory::build_structure_algorithm(const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase = Teuchos::null;

  const enum INPAR::STR::IntegrationStrategy intstrat =
      CORE::UTILS::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case INPAR::STR::int_standard:
      adapterbase = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithmNew());
      break;
    default:
      FOUR_C_THROW("Unknown integration strategy!");
      break;
  }

  return adapterbase;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> ADAPTER::build_structure_algorithm(
    const Teuchos::ParameterList& sdyn)
{
  StructureFactory factory;
  return factory.build_structure_algorithm(sdyn);
}

FOUR_C_NAMESPACE_CLOSE
