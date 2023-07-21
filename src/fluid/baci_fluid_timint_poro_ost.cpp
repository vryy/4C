/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_timint_poro_ost.H"
#include "baci_io.H"

FLD::TimIntPoroOst::TimIntPoroOst(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntPoro(actdis, solver, params, output, alefluid)
{
}

void FLD::TimIntPoroOst::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  TimIntPoro::Init();
}

void FLD::TimIntPoroOst::ReadRestart(int step)
{
  // call of base classes
  TimIntOneStepTheta::ReadRestart(step);
  TimIntPoro::ReadRestart(step);
}