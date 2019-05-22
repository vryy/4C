/*-----------------------------------------------------------*/
/*!

\maintainer Anh-Tu Vuong

\brief  structure adapter using the LOCA library

\level 3

*/
/*-----------------------------------------------------------*/

#include "ad_str_loca.H"
#include "ad_str_loca_wrapper.H"

#include "../drt_structure_new/str_timint_base.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::StructureLocaAlgorithm::StructureLocaAlgorithm()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureLocaAlgorithm::Setup()
{
  if (not IsInit()) dserror("Call Init() first!");

  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(*sdyn_, "DYNAMICTYP") !=
      INPAR::STR::dyna_statics)
    dserror("Currently LOCA supports only the static case!");

  // call the setup routine of the base class
  ADAPTER::StructureBaseAlgorithmNew::SetupTimInt();

  if (not IsSetup()) dserror("The base class Setup() routine should have set this flag!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureLocaAlgorithm::SetParams(Teuchos::ParameterList& ioflags,
    Teuchos::ParameterList& xparams, Teuchos::ParameterList& taflags)
{
  // get the problem instance and the problem type
  DRT::Problem* problem = DRT::Problem::Instance();

  // Get LOCA specific parameters
  Teuchos::RCP<Teuchos::ParameterList> sloca =
      Teuchos::rcp(new Teuchos::ParameterList(problem->LocaParams()));
  Teuchos::ParameterList& loca = xparams.sublist("LOCA");
  loca = *sloca;

  // Set the base class variables (though most are unused in the LOCA case)
  ADAPTER::StructureBaseAlgorithmNew::SetParams(ioflags, xparams, taflags);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureLocaAlgorithm::SetStructureWrapper(const Teuchos::ParameterList& ioflags,
    const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
    const Teuchos::ParameterList& taflags, Teuchos::RCP<STR::TIMINT::Base> ti_strategy)
{
  str_wrapper_ = Teuchos::rcp(new StructureLocaWrapper(ti_strategy));

  return;
}
