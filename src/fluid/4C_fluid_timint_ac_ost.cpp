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
FLD::TimIntACOst::TimIntACOst(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntAC(actdis, solver, params, output, alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
 |  read restart data                                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntACOst::read_restart(int step)
{
  // call of base classes
  TimIntOneStepTheta::read_restart(step);
  TimIntAC::read_restart(step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
