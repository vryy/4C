/*-----------------------------------------------------------*/
/*! \file

\brief BDF2 time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_loma_bdf2.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLomaBDF2::TimIntLomaBDF2(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntBDF2(actdis, solver, params, output, alefluid),
      TimIntLoma(actdis, solver, params, output, alefluid)
{
  std::cout << "\nWARNING: Loma has never been tested with BDF2 time integration!!\n" << std::endl;
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLomaBDF2::init()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntBDF2::init();
  TimIntLoma::init();

  return;
}

FOUR_C_NAMESPACE_CLOSE
