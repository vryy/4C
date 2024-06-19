/*-----------------------------------------------------------*/
/*! \file

\brief Stationary driver for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_red_stat.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedModelsStat::TimIntRedModelsStat(const Teuchos::RCP<Core::FE::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntStationary(actdis, solver, params, output, alefluid),
      TimIntRedModels(actdis, solver, params, output, alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsStat::init()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntStationary::init();
  TimIntRedModels::init();

  return;
}


/*----------------------------------------------------------------------*
 |  read restart data                                   rasthofer 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsStat::read_restart(int step)
{
  // call of base classes
  TimIntStationary::read_restart(step);
  TimIntRedModels::read_restart(step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
