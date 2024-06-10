/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem driver for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

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

void FLD::TimIntPoroStat::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  TimIntPoro::Init();
}

void FLD::TimIntPoroStat::read_restart(int step)
{
  // call of base classes
  TimIntStationary::read_restart(step);
  TimIntPoro::read_restart(step);
}

FOUR_C_NAMESPACE_CLOSE
