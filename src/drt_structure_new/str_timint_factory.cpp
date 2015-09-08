/*
 * str_timint_factory.cpp
 *
 *  Created on: Sep 4, 2015
 *      Author: hiermeier
 */

#include "str_timint_factory.H"

#include "../drt_inpar/inpar_structure.H"

#include <Teuchos_ParameterList.hpp>

// supported time integrator
#include "str_timint_implicit.H"
#include "str_timint_explicit.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildTimeIntegrator(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> timeintegrator = Teuchos::null;

  // Check first if a implicit time integrator is desired
  timeintegrator = BuildImplicitTimeIntegrator(sdyn);

  // If there was no suitable implicit time integrator check for the
  // explicit case
  if (timeintegrator.is_null())
    timeintegrator = BuildExplicitTimeIntegrator(sdyn);

  return timeintegrator;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildImplicitTimeIntegrator(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> timeintegrator = Teuchos::null;

  // get the prestress type
  const enum INPAR::STR::PreStress pstype =
      DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  // get the dynamic type
  const enum INPAR::STR::DynamicType dyntype =
      DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");

  if (pstype == INPAR::STR::prestress_mulf     or   // prestress type
      pstype == INPAR::STR::prestress_id       or
      dyntype == INPAR::STR::dyna_statics      or   // dynamic type
      dyntype == INPAR::STR::dyna_genalpha     or
      dyntype == INPAR::STR::dyna_onesteptheta or
      dyntype == INPAR::STR::dyna_gemm         or
      dyntype == INPAR::STR::dyna_statmech)
    timeintegrator = Teuchos::rcp(new STR::TIMINT::Implicit());

  return timeintegrator;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::Factory::BuildExplicitTimeIntegrator(
    const Teuchos::ParameterList& sdyn) const
{
  Teuchos::RCP<STR::TIMINT::Base> timeintegrator = Teuchos::null;

  // what's the current problem type?
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  if (probtype == prb_fsi           or
      probtype == prb_fsi_redmodels or
      probtype == prb_fsi_lung      or
      probtype == prb_gas_fsi       or
      probtype == prb_ac_fsi        or
      probtype == prb_biofilm_fsi   or
      probtype == prb_thermo_fsi)
    dserror("No explicit time integration with fsi");

  const enum INPAR::STR::DynamicType dyntype =
      DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");

  if (dyntype == INPAR::STR::dyna_expleuler or
      dyntype == INPAR::STR::dyna_centrdiff or
      dyntype == INPAR::STR::dyna_ab2)
    timeintegrator = Teuchos::rcp(new STR::TIMINT::Explicit());

  return timeintegrator;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::Base> STR::TIMINT::BuildTimeIntegrator(
    const Teuchos::ParameterList& sdyn)
{
  Factory factory;
  return factory.BuildTimeIntegrator(sdyn);
}
