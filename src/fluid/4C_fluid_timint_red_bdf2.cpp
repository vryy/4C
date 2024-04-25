/*-----------------------------------------------------------*/
/*! \file

\brief BDF-2 time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_red_bdf2.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedModelsBDF2::TimIntRedModelsBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntBDF2(actdis, solver, params, output, alefluid),
      TimIntRedModels(actdis, solver, params, output, alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntBDF2::Init();
  TimIntRedModels::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  read restart data                                   rasthofer 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsBDF2::ReadRestart(int step)
{
  // call of base classes
  TimIntBDF2::ReadRestart(step);
  TimIntRedModels::ReadRestart(step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
