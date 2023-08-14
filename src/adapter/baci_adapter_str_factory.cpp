/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_adapter_str_factory.H"

#include "baci_adapter_str_structure_new.H"
#include "baci_inpar_structure.H"
#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::STR::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> ADAPTER::STR::Factory::BuildStructureAlgorithm(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase = Teuchos::null;

  const enum INPAR::STR::IntegrationStrategy intstrat =
      DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case INPAR::STR::int_standard:
      adapterbase = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithmNew());
      break;
    default:
      dserror("Unknown integration strategy!");
      break;
  }

  return adapterbase;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> ADAPTER::STR::BuildStructureAlgorithm(
    const Teuchos::ParameterList& sdyn)
{
  STR::Factory factory;
  return factory.BuildStructureAlgorithm(sdyn);
}
