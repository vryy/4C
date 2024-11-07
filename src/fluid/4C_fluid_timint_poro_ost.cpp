// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_timint_poro_ost.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN

FLD::TimIntPoroOst::TimIntPoroOst(const std::shared_ptr<Core::FE::Discretization>& actdis,
    const std::shared_ptr<Core::LinAlg::Solver>& solver,
    const std::shared_ptr<Teuchos::ParameterList>& params,
    const std::shared_ptr<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntPoro(actdis, solver, params, output, alefluid)
{
}

void FLD::TimIntPoroOst::init()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::init();
  TimIntPoro::init();
}

void FLD::TimIntPoroOst::read_restart(int step)
{
  // call of base classes
  TimIntOneStepTheta::read_restart(step);
  TimIntPoro::read_restart(step);
}
FOUR_C_NAMESPACE_CLOSE
