/*-----------------------------------------------------------*/
/*!
\file str_factory.cpp

\brief factory for time integrator

\maintainer Michael Hiermeier

\date Nov 20, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_factory.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"

// supported implicit time integrators
#include "str_impl_statics.H"
#include "str_impl_prestress.H" // derived from statics
#include "str_impl_genalpha.H"
#include "str_impl_genalpha_liegroup.H"
#include "str_impl_gemm.H"
#include "str_impl_statmech.H"  // derived from ost

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator> STR::Factory::BuildIntegrator(
    const STR::TIMINT::BaseDataSDyn& datasdyn) const
{
  Teuchos::RCP<STR::Integrator> int_ptr = Teuchos::null;
  int_ptr = BuildImplicitIntegrator(datasdyn);
  if (int_ptr.is_null())
    int_ptr = BuildExplicitIntegrator(datasdyn);
  if (int_ptr.is_null())
    dserror("We could not find a suitable dynamic integrator (Dynamic Type).");

  return int_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator> STR::Factory::BuildImplicitIntegrator(
    const STR::TIMINT::BaseDataSDyn& datasdyn) const
{
  Teuchos::RCP<STR::IMPLICIT::Generic> impl_int_ptr = Teuchos::null;

  const enum INPAR::STR::DynamicType& dyntype = datasdyn.GetDynamicType();
  const enum INPAR::STR::PreStress& prestresstype = datasdyn.GetPreStressType();

  // check if we have a problem that needs to be prestressed
  if (prestresstype==INPAR::STR::prestress_mulf or prestresstype==INPAR::STR::prestress_id)
  {
    impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::PreStress());
    return impl_int_ptr;
  }

  switch (dyntype)
  {
    // Static analysis
    case INPAR::STR::dyna_statics :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::Statics());
      break;
    }

    // Generalised-alpha time integration
    case INPAR::STR::dyna_genalpha :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::GenAlpha());
      break;
    }

    // Generalised-alpha time integration for Lie groups (e.g. SO3 group of rotation matrices)
    case INPAR::STR::dyna_genalpha_liegroup :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::GenAlphaLieGroup());
      break;
    }

    // One-step-theta (OST) time integration
    case INPAR::STR::dyna_onesteptheta :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::OneStepTheta());
      break;
    }

    // Generalised energy-momentum method (GEMM)
    case INPAR::STR::dyna_gemm :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::Gemm());
      break;
    }

    // Statistical Mechanics Time Integration
    case INPAR::STR::dyna_statmech :
    {
      impl_int_ptr = Teuchos::rcp(new STR::IMPLICIT::StatMech());
      break;
    }

    // Everything else
    default :
    {
      /* Do nothing and return Techos::null. */
      break;
    }
  } // end of switch(dynType)

  return impl_int_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator> STR::Factory::BuildExplicitIntegrator(
    const STR::TIMINT::BaseDataSDyn& datasdyn) const
{
//  Teuchos::RCP<STR::EXPLICIT::Generic> expl_int_ptr = Teuchos::null;
  dserror("Not yet implemented!");

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator> STR::BuildIntegrator(
    const STR::TIMINT::BaseDataSDyn& datasdyn)
{
  STR::Factory factory;

  return factory.BuildIntegrator(datasdyn);
}
