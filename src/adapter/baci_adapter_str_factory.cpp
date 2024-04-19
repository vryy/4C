/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_adapter_str_factory.hpp"

#include "baci_adapter_str_structure_new.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::StructureFactory::StructureFactory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> ADAPTER::StructureFactory::BuildStructureAlgorithm(
    const Teuchos::ParameterList& sdyn) const
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
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> ADAPTER::BuildStructureAlgorithm(
    const Teuchos::ParameterList& sdyn)
{
  StructureFactory factory;
  return factory.BuildStructureAlgorithm(sdyn);
}

FOUR_C_NAMESPACE_CLOSE
