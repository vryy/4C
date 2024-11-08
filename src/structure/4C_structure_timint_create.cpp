// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_timint_create.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_timint_ab2.hpp"
#include "4C_structure_timint_centrdiff.hpp"
#include "4C_structure_timint_expleuler.hpp"
#include "4C_structure_timint_genalpha.hpp"
#include "4C_structure_timint_ost.hpp"
#include "4C_structure_timint_prestress.hpp"
#include "4C_structure_timint_statics.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* create marching time integrator */
std::shared_ptr<Solid::TimInt> Solid::tim_int_create(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, std::shared_ptr<Core::FE::Discretization>& actdis,
    std::shared_ptr<Core::LinAlg::Solver>& solver,
    std::shared_ptr<Core::LinAlg::Solver>& contactsolver,
    std::shared_ptr<Core::IO::DiscretizationWriter>& output)
{
  // set default output
  std::shared_ptr<Solid::TimInt> sti = nullptr;
  // try implicit integrators
  sti = tim_int_impl_create(
      timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
  // if nothing found try explicit integrators
  if (sti == nullptr)
  {
    sti = tim_int_expl_create(
        timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
  }

  // deliver
  return sti;
}

/*======================================================================*/
/* create implicit marching time integrator */
std::shared_ptr<Solid::TimIntImpl> Solid::tim_int_impl_create(
    const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& ioflags,
    const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
    std::shared_ptr<Core::FE::Discretization>& actdis,
    std::shared_ptr<Core::LinAlg::Solver>& solver,
    std::shared_ptr<Core::LinAlg::Solver>& contactsolver,
    std::shared_ptr<Core::IO::DiscretizationWriter>& output)
{
  std::shared_ptr<Solid::TimIntImpl> sti = nullptr;

  // TODO: add contact solver...

  // check if we have a problem that needs to be prestressed
  if (Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
          Global::Problem::instance()->structural_dynamic_params(), "PRESTRESS") !=
      Inpar::Solid::PreStress::none)
  {
    sti = std::make_shared<Solid::TimIntPrestress>(
        timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
    return sti;
  }

  // create specific time integrator
  switch (Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYPE"))
  {
    // Static analysis
    case Inpar::Solid::dyna_statics:
    {
      sti = std::make_shared<Solid::TimIntStatics>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }

    // Generalised-alpha time integration
    case Inpar::Solid::dyna_genalpha:
    {
      sti = std::make_shared<Solid::TimIntGenAlpha>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }

    // One-step-theta (OST) time integration
    case Inpar::Solid::dyna_onesteptheta:
    {
      sti = std::make_shared<Solid::TimIntOneStepTheta>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }

    // Everything else
    default:
    {
      // do nothing
      break;
    }
  }  // end of switch(sdyn->Type)

  // return the integrator
  return sti;
}

/*======================================================================*/
/* create explicit marching time integrator */
std::shared_ptr<Solid::TimIntExpl> Solid::tim_int_expl_create(
    const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& ioflags,
    const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
    std::shared_ptr<Core::FE::Discretization>& actdis,
    std::shared_ptr<Core::LinAlg::Solver>& solver,
    std::shared_ptr<Core::LinAlg::Solver>& contactsolver,
    std::shared_ptr<Core::IO::DiscretizationWriter>& output)
{
  std::shared_ptr<Solid::TimIntExpl> sti = nullptr;

  // what's the current problem type?
  Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();

  if (probtype == Core::ProblemType::fsi or probtype == Core::ProblemType::fsi_redmodels or
      probtype == Core::ProblemType::gas_fsi or probtype == Core::ProblemType::biofilm_fsi or
      probtype == Core::ProblemType::thermo_fsi)
  {
    FOUR_C_THROW("no explicit time integration with fsi");
  }

  // create specific time integrator
  switch (Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(sdyn, "DYNAMICTYPE"))
  {
    // forward Euler time integration
    case Inpar::Solid::dyna_expleuler:
    {
      sti = std::make_shared<Solid::TimIntExplEuler>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }
    // central differences time integration
    case Inpar::Solid::dyna_centrdiff:
    {
      sti = std::make_shared<Solid::TimIntCentrDiff>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }
    // Adams-Bashforth 2nd order (AB2) time integration
    case Inpar::Solid::dyna_ab2:
    {
      sti = std::make_shared<Solid::TimIntAB2>(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
      break;
    }

    // Everything else
    default:
    {
      // do nothing
      break;
    }
  }  // end of switch(sdyn->Type)

  // return the integrator
  return sti;
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
