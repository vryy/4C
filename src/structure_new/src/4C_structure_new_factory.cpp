/*-----------------------------------------------------------*/
/*! \file

\brief factory for time integrator


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_expl_ab2.hpp"
#include "4C_structure_new_expl_abx.hpp"
#include "4C_structure_new_expl_centrdiff.hpp"
#include "4C_structure_new_expl_forwardeuler.hpp"
#include "4C_structure_new_impl_genalpha.hpp"
#include "4C_structure_new_impl_genalpha_liegroup.hpp"
#include "4C_structure_new_impl_ost.hpp"        // derived from ost
#include "4C_structure_new_impl_prestress.hpp"  // derived from statics
#include "4C_structure_new_impl_statics.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Integrator> Solid::Factory::build_integrator(
    const Solid::TimeInt::BaseDataSDyn& datasdyn) const
{
  Teuchos::RCP<Solid::Integrator> int_ptr = Teuchos::null;
  int_ptr = build_implicit_integrator(datasdyn);
  if (int_ptr.is_null()) int_ptr = build_explicit_integrator(datasdyn);
  FOUR_C_ASSERT(
      !int_ptr.is_null(), "We could not find a suitable dynamic integrator (Dynamic Type).");

  return int_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Integrator> Solid::Factory::build_implicit_integrator(
    const Solid::TimeInt::BaseDataSDyn& datasdyn) const
{
  Teuchos::RCP<Solid::IMPLICIT::Generic> impl_int_ptr = Teuchos::null;

  const enum Inpar::Solid::DynamicType& dyntype = datasdyn.get_dynamic_type();
  const enum Inpar::Solid::PreStress& prestresstype = datasdyn.get_pre_stress_type();

  // check if we have a problem that needs to be prestressed
  const bool is_prestress = prestresstype != Inpar::Solid::PreStress::none;
  if (is_prestress)
  {
    impl_int_ptr = Teuchos::rcp(new Solid::IMPLICIT::PreStress());
    return impl_int_ptr;
  }

  switch (dyntype)
  {
    // Static analysis
    case Inpar::Solid::dyna_statics:
    {
      impl_int_ptr = Teuchos::rcp(new Solid::IMPLICIT::Statics());
      break;
    }

    // Generalised-alpha time integration
    case Inpar::Solid::dyna_genalpha:
    {
      impl_int_ptr = Teuchos::rcp(new Solid::IMPLICIT::GenAlpha());
      break;
    }

    // Generalised-alpha time integration for Lie groups (e.g. SO3 group of rotation matrices)
    case Inpar::Solid::dyna_genalpha_liegroup:
    {
      impl_int_ptr = Teuchos::rcp(new Solid::IMPLICIT::GenAlphaLieGroup());
      break;
    }

    // One-step-theta (OST) time integration
    case Inpar::Solid::dyna_onesteptheta:
    {
      impl_int_ptr = Teuchos::rcp(new Solid::IMPLICIT::OneStepTheta());
      break;
    }

    // Everything else
    default:
    {
      /* Do nothing and return Techos::null. */
      break;
    }
  }  // end of switch(dynType)

  return impl_int_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Integrator> Solid::Factory::build_explicit_integrator(
    const Solid::TimeInt::BaseDataSDyn& datasdyn) const
{
  Teuchos::RCP<Solid::EXPLICIT::Generic> expl_int_ptr = Teuchos::null;

  switch (datasdyn.get_dynamic_type())
  {
    // Forward Euler Scheme
    case Inpar::Solid::dyna_expleuler:
    {
      expl_int_ptr = Teuchos::rcp(new Solid::EXPLICIT::ForwardEuler());
      break;
    }

    // Central Difference Scheme
    case Inpar::Solid::dyna_centrdiff:
    {
      expl_int_ptr = Teuchos::rcp(new Solid::EXPLICIT::CentrDiff());
      break;
    }

    // Adams-Bashforth-2 Scheme
    case Inpar::Solid::dyna_ab2:
    {
      expl_int_ptr = Teuchos::rcp(new Solid::EXPLICIT::AdamsBashforth2());
      break;
    }

    // Adams-Bashforth-4 Scheme
    case Inpar::Solid::dyna_ab4:
    {
      expl_int_ptr = Teuchos::rcp(new Solid::EXPLICIT::AdamsBashforthX<4>());
      break;
    }

    // Everything else
    default:
    {
      /* Do nothing and return Techos::null. */
      break;
    }
  }  // end of switch(dynType)

  return expl_int_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Integrator> Solid::build_integrator(
    const Solid::TimeInt::BaseDataSDyn& datasdyn)
{
  Solid::Factory factory;

  return factory.build_integrator(datasdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Dbc> Solid::Factory::build_dbc(
    const Solid::TimeInt::BaseDataSDyn& datasdyn) const
{
  // if you want your model specific dbc object, check here if your model type is
  // active ( datasdyn.get_model_types() )and build your own dbc object
  Teuchos::RCP<Solid::Dbc> dbc = Teuchos::null;
  dbc = Teuchos::rcp(new Solid::Dbc());

  return dbc;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Solid::Dbc> Solid::build_dbc(const Solid::TimeInt::BaseDataSDyn& datasdyn)
{
  Solid::Factory factory;

  return factory.build_dbc(datasdyn);
}

FOUR_C_NAMESPACE_CLOSE
