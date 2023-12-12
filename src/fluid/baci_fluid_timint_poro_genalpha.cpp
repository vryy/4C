/*-----------------------------------------------------------*/
/*! \file

\brief generalized alpha time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_timint_poro_genalpha.H"

#include "baci_io.H"

BACI_NAMESPACE_OPEN


FLD::TimIntPoroGenAlpha::TimIntPoroGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntGenAlpha(actdis, solver, params, output, alefluid),
      TimIntPoro(actdis, solver, params, output, alefluid)
{
}


void FLD::TimIntPoroGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  TimIntPoro::Init();
}

void FLD::TimIntPoroGenAlpha::GenAlphaIntermediateValues()
{
  // set intermediate values for acceleration and potential temporal
  // derivatives
  //
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)

  // consider both velocity and pressure degrees of freedom
  accam_->Update((alphaM_), *accnp_, (1.0 - alphaM_), *accn_, 0.0);

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  velaf_->Update((alphaF_), *velnp_, (1.0 - alphaF_), *veln_, 0.0);
}

void FLD::TimIntPoroGenAlpha::ReadRestart(int step)
{
  // call of base classes
  TimIntGenAlpha::ReadRestart(step);
  TimIntPoro::ReadRestart(step);
}

BACI_NAMESPACE_CLOSE
