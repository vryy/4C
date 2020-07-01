/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_loma_genalpha.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLomaGenAlpha::TimIntLomaGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntGenAlpha(actdis, solver, params, output, alefluid),
      TimIntLoma(actdis, solver, params, output, alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLomaGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  TimIntLoma::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntLomaGenAlpha::~TimIntLomaGenAlpha() { return; }
