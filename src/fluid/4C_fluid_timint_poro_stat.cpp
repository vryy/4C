// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_timint_poro_stat.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN


FLD::TimIntPoroStat::TimIntPoroStat(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntStationary(actdis, solver, params, output, alefluid),
      TimIntPoro(actdis, solver, params, output, alefluid)
{
}

void FLD::TimIntPoroStat::init()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntStationary::init();
  TimIntPoro::init();
}

void FLD::TimIntPoroStat::read_restart(int step)
{
  // call of base classes
  TimIntStationary::read_restart(step);
  TimIntPoro::read_restart(step);
}

FOUR_C_NAMESPACE_CLOSE
