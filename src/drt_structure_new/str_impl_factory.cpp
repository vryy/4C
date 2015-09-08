/*
 * str_impl_factory.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: hiermeier
 */


#include "str_impl_factory.H"

// supported implicit time integrators
#include "str_impl_prestress.H" // derived from statics
#include "str_impl_genalpha.H"
#include "str_impl_gemm.H"
#include "str_impl_statmech.H"  // derived from ost

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::IMPLICIT::Generic> STR::IMPLICIT::Factory::BuildImplicitIntegrator(
    const enum INPAR::STR::DynamicType& dynType,
    const enum INPAR::STR::PreStress& preStressType) const
{
  Teuchos::RCP<STR::IMPLICIT::Generic> implIntegrator = Teuchos::null;

  // check if we have a problem that needs to be prestressed
  if (preStressType==INPAR::STR::prestress_mulf || preStressType==INPAR::STR::prestress_id)
  {
    implIntegrator = Teuchos::rcp(new STR::IMPLICIT::PreStress());
    return implIntegrator;
  }

  switch (dynType)
  {
    // Static analysis
    case INPAR::STR::dyna_statics :
    {
      implIntegrator = Teuchos::rcp(new STR::IMPLICIT::Statics());
      break;
    }

    // Generalised-alpha time integration
    case INPAR::STR::dyna_genalpha :
    {
      implIntegrator = Teuchos::rcp(new STR::IMPLICIT::GenAlpha());
      break;
    }

    // One-step-theta (OST) time integration
    case INPAR::STR::dyna_onesteptheta :
    {
      implIntegrator = Teuchos::rcp(new STR::IMPLICIT::OneStepTheta());
      break;
    }

    // Generalised energy-momentum method (GEMM)
    case INPAR::STR::dyna_gemm :
    {
      implIntegrator = Teuchos::rcp(new STR::IMPLICIT::GEMM());
      break;
    }

    // Statistical Mechanics Time Integration
    case INPAR::STR::dyna_statmech :
    {
      implIntegrator = Teuchos::rcp(new STR::IMPLICIT::StatMech());
      break;
    }

    // Everything else
    default :
    {
      dserror("Unknow implicit time integrator!");
      break;
    }
  } // end of switch(dynType)

  return implIntegrator;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::IMPLICIT::Generic> STR::IMPLICIT::BuildImplicitIntegrator(
    const enum INPAR::STR::DynamicType& dynType,
    const enum INPAR::STR::PreStress& preStressType)
{
  Factory factory;
  return factory.BuildImplicitIntegrator(dynType,preStressType);
}
