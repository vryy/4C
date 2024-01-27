/*-----------------------------------------------------------*/
/*! \file

\brief factory for time integration base strategy and data container


\level 3

*/
/*-----------------------------------------------------------*/


#include "baci_structure_new_timint_factory.H"

#include "baci_global_data.H"
#include "baci_inpar_structure.H"
#include "baci_structure_new_timint_basedatasdyn.H"
#include "baci_structure_new_timint_explicit.H"
#include "baci_structure_new_timint_implicit.H"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Factory::Factory()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildStrategy(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> ti_strategy = Teuchos::null;

  const enum INPAR::STR::IntegrationStrategy intstrat =
      INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case INPAR::STR::int_standard:
    {
      // Check first if a implicit integration strategy is desired
      ti_strategy = BuildImplicitStrategy(sdyn);
      // If there was no suitable implicit time integrator check for the
      // explicit case
      if (ti_strategy.is_null()) ti_strategy = BuildExplicitStrategy(sdyn);
      break;
    }
    default:
      dserror("Unknown integration strategy!");
      break;
  }

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildImplicitStrategy(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> ti_strategy = Teuchos::null;

  // get the dynamic type
  const enum INPAR::STR::DynamicType dyntype =
      INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");

  const bool is_prestress = Teuchos::getIntegralValue<INPAR::STR::PreStress>(
                                GLOBAL::Problem::Instance()->StructuralDynamicParams(),
                                "PRESTRESS") != INPAR::STR::PreStress::none;
  if (is_prestress or dyntype == INPAR::STR::dyna_statics or  // dynamic type
      dyntype == INPAR::STR::dyna_genalpha or dyntype == INPAR::STR::dyna_genalpha_liegroup or
      dyntype == INPAR::STR::dyna_onesteptheta or dyntype == INPAR::STR::dyna_gemm)
    ti_strategy = Teuchos::rcp(new STR::TIMINT::Implicit());

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildExplicitStrategy(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> ti_strategy = Teuchos::null;

  // what's the current problem type?
  GLOBAL::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();

  if (probtype == GLOBAL::ProblemType::fsi or probtype == GLOBAL::ProblemType::fsi_redmodels or
      probtype == GLOBAL::ProblemType::fsi_lung or probtype == GLOBAL::ProblemType::gas_fsi or
      probtype == GLOBAL::ProblemType::ac_fsi or probtype == GLOBAL::ProblemType::biofilm_fsi or
      probtype == GLOBAL::ProblemType::thermo_fsi)
    dserror("No explicit time integration with fsi");

  const enum INPAR::STR::DynamicType dyntype =
      INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");

  if (dyntype == INPAR::STR::dyna_expleuler or dyntype == INPAR::STR::dyna_centrdiff or
      dyntype == INPAR::STR::dyna_ab2 or dyntype == INPAR::STR::dyna_ab4)
    ti_strategy = Teuchos::rcp(new STR::TIMINT::Explicit());

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataSDyn> STR::TIMINT::Factory::BuildDataSDyn(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyndata_ptr = Teuchos::null;

  const enum INPAR::STR::DynamicType dyntype =
      INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");

  switch (dyntype)
  {
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_genalpha_liegroup:
      sdyndata_ptr = Teuchos::rcp(new STR::TIMINT::GenAlphaDataSDyn());
      break;
    case INPAR::STR::dyna_onesteptheta:
      sdyndata_ptr = Teuchos::rcp(new STR::TIMINT::OneStepThetaDataSDyn());
      break;
    case INPAR::STR::dyna_expleuler:
      sdyndata_ptr = Teuchos::rcp(new STR::TIMINT::ExplEulerDataSDyn());
      break;
    default:
      sdyndata_ptr = Teuchos::rcp(new STR::TIMINT::BaseDataSDyn());
      break;
  }

  return sdyndata_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> STR::TIMINT::Factory::BuildDataGlobalState() const
{
  return Teuchos::rcp(new STR::TIMINT::BaseDataGlobalState());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::BuildStrategy(const Teuchos::ParameterList& sdyn)
{
  Factory factory;
  return factory.BuildStrategy(sdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataSDyn> STR::TIMINT::BuildDataSDyn(
    const Teuchos::ParameterList& sdyn)
{
  Factory factory;
  return factory.BuildDataSDyn(sdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> STR::TIMINT::BuildDataGlobalState()
{
  Factory factory;
  return factory.BuildDataGlobalState();
}

BACI_NAMESPACE_CLOSE
