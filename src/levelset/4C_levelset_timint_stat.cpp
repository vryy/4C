// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_levelset_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
ScaTra::LevelSetTimIntStationary::LevelSetTimIntStationary(
    std::shared_ptr<Core::FE::Discretization> actdis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Teuchos::ParameterList> params,
    std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
    std::shared_ptr<Teuchos::ParameterList> extraparams,
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      LevelSetAlgorithm(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}



/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetTimIntStationary::init()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntStationary::init();
  LevelSetAlgorithm::init();

  if (myrank_ == 0)
  {
    std::cout << "\n*------------------------------------------------------------------------*"
              << std::endl;
    std::cout << "|            stationary level set for coupled problems only              |"
              << std::endl;
    std::cout << "*------------------------------------------------------------------------*\n"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  setup time integration                                  rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetTimIntStationary::setup()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntStationary::setup();
  LevelSetAlgorithm::setup();

  if (myrank_ == 0)
  {
    std::cout << "\n*------------------------------------------------------------------------*"
              << std::endl;
    std::cout << "|            stationary level set for coupled problems only              |"
              << std::endl;
    std::cout << "*------------------------------------------------------------------------*\n"
              << std::endl;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
