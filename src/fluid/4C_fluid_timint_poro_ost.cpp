/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_poro_ost.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN

FLD::TimIntPoroOst::TimIntPoroOst(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
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
