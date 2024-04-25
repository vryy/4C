/*-----------------------------------------------------------*/
/*! \file

\brief OST time integrator for FS3I-AC problems. Is in diamond inhertance with TimIntOneStepTheta
       and TimIntAC


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_ac_ost.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntACOst::TimIntACOst(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntAC(actdis, solver, params, output, alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
 |  read restart data                                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntACOst::ReadRestart(int step)
{
  // call of base classes
  TimIntOneStepTheta::ReadRestart(step);
  TimIntAC::ReadRestart(step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
