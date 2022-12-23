/*-----------------------------------------------------------*/
/*! \file

\brief OST time integrator for FS3I-AC problems. Is in diamond inhertance with TimIntOneStepTheta
       and TimIntAC


\level 3

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_ac_ost.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntACOst::TimIntACOst(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntAC(actdis, solver, params, output, alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   Thon 12/14 |
*----------------------------------------------------------------------*/
FLD::TimIntACOst::~TimIntACOst() { return; }

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
