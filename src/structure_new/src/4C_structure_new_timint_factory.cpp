// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_timint_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"
#include "4C_structure_new_timint_explicit.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::TimeInt::Factory::Factory()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::Base> Solid::TimeInt::Factory::build_strategy(
    const Teuchos::ParameterList& sdyn) const
{
  std::shared_ptr<Solid::TimeInt::Base> ti_strategy = nullptr;

  const auto intstrat =
      Teuchos::getIntegralValue<Inpar::Solid::IntegrationStrategy>(sdyn, "INT_STRATEGY");

  switch (intstrat)
  {
    case Inpar::Solid::int_standard:
    {
      // Check first if a implicit integration strategy is desired
      ti_strategy = build_implicit_strategy(sdyn);
      // If there was no suitable implicit time integrator check for the
      // explicit case
      if (!ti_strategy) ti_strategy = build_explicit_strategy(sdyn);
      break;
    }
    default:
      FOUR_C_THROW("Unknown integration strategy!");
      break;
  }

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::Base> Solid::TimeInt::Factory::build_implicit_strategy(
    const Teuchos::ParameterList& sdyn) const
{
  std::shared_ptr<Solid::TimeInt::Base> ti_strategy = nullptr;

  // get the dynamic type
  const auto dyntype = Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYPE");

  const bool is_prestress = Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
                                Global::Problem::instance()->structural_dynamic_params(),
                                "PRESTRESS") != Inpar::Solid::PreStress::none;
  if (is_prestress or dyntype == Inpar::Solid::dyna_statics or  // dynamic type
      dyntype == Inpar::Solid::dyna_genalpha or dyntype == Inpar::Solid::dyna_genalpha_liegroup or
      dyntype == Inpar::Solid::dyna_onesteptheta)
    ti_strategy = std::make_shared<Solid::TimeInt::Implicit>();

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::Base> Solid::TimeInt::Factory::build_explicit_strategy(
    const Teuchos::ParameterList& sdyn) const
{
  std::shared_ptr<Solid::TimeInt::Base> ti_strategy = nullptr;

  // what's the current problem type?
  Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();

  if (probtype == Core::ProblemType::fsi or probtype == Core::ProblemType::fsi_redmodels or
      probtype == Core::ProblemType::gas_fsi or probtype == Core::ProblemType::biofilm_fsi or
      probtype == Core::ProblemType::thermo_fsi)
    FOUR_C_THROW("No explicit time integration with fsi");

  const auto dyntype = Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYPE");

  if (dyntype == Inpar::Solid::dyna_expleuler or dyntype == Inpar::Solid::dyna_centrdiff or
      dyntype == Inpar::Solid::dyna_ab2 or dyntype == Inpar::Solid::dyna_ab4)
    ti_strategy = std::make_shared<Solid::TimeInt::Explicit>();

  return ti_strategy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::BaseDataSDyn> Solid::TimeInt::Factory::build_data_sdyn(
    const Teuchos::ParameterList& sdyn) const
{
  std::shared_ptr<Solid::TimeInt::BaseDataSDyn> sdyndata_ptr = nullptr;

  const auto dyntype = Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYPE");

  switch (dyntype)
  {
    case Inpar::Solid::dyna_genalpha:
    case Inpar::Solid::dyna_genalpha_liegroup:
      sdyndata_ptr = std::make_shared<Solid::TimeInt::GenAlphaDataSDyn>();
      break;
    case Inpar::Solid::dyna_onesteptheta:
      sdyndata_ptr = std::make_shared<Solid::TimeInt::OneStepThetaDataSDyn>();
      break;
    case Inpar::Solid::dyna_expleuler:
      sdyndata_ptr = std::make_shared<Solid::TimeInt::ExplEulerDataSDyn>();
      break;
    default:
      sdyndata_ptr = std::make_shared<Solid::TimeInt::BaseDataSDyn>();
      break;
  }

  return sdyndata_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>
Solid::TimeInt::Factory::build_data_global_state() const
{
  return std::make_shared<Solid::TimeInt::BaseDataGlobalState>();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::Base> Solid::TimeInt::build_strategy(
    const Teuchos::ParameterList& sdyn)
{
  Factory factory;
  return factory.build_strategy(sdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::BaseDataSDyn> Solid::TimeInt::build_data_sdyn(
    const Teuchos::ParameterList& sdyn)
{
  Factory factory;
  return factory.build_data_sdyn(sdyn);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::TimeInt::BaseDataGlobalState> Solid::TimeInt::build_data_global_state()
{
  Factory factory;
  return factory.build_data_global_state();
}

FOUR_C_NAMESPACE_CLOSE
