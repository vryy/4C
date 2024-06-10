/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for thermal dynamics

\level 3

*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 01/12 |
 *----------------------------------------------------------------------*/
#include "4C_thermo_timint_expl.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_thermo_timint.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                               dano 01/12 |
 *----------------------------------------------------------------------*/
THR::TimIntExpl::TimIntExpl(const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& tdynparams,                        //!< input parameters
    const Teuchos::ParameterList& xparams,                           //!< extra flags
    Teuchos::RCP<Core::FE::Discretization> actdis,                   //!< current discretisation
    Teuchos::RCP<Core::LinAlg::Solver> solver,                       //!< the solver
    Teuchos::RCP<Core::IO::DiscretizationWriter> output              //!< the output
    )
    : TimInt(ioparams, tdynparams, xparams, actdis, solver, output)
{
  // get away
  return;
}  // TimIntExplEuler()


/*----------------------------------------------------------------------*
 | update time step                                          dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExpl::Update()
{
  // update temperature and temperature rate
  // after this call we will have tempn_ == temp_ (temp_{n+1} == temp_n), etc.
  UpdateStepState();
  // update time and step
  UpdateStepTime();
  // currently nothing, can include history dependency of materials
  UpdateStepElement();
  return;

}  // Update()


/*----------------------------------------------------------------------*
 | print step summary                                        dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExpl::PrintStep()
{
  // print out
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    print_step_text(stdout);
  }
}  // print_step()


/*----------------------------------------------------------------------*
 | print step summary                                        dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExpl::print_step_text(FILE* ofile)
{
  fprintf(ofile,
      "Finalised: step %6d"
      " | nstep %6d"
      " | time %-14.8E"
      " | dt %-14.8E"
      " | numiter %3d\n",
      //     " | wct %-14.8E\n",
      step_, stepmax_, (*time_)[0], (*dt_)[0], 0
      //       timer_->totalElapsedTime(true)
  );

  // print a beautiful line made exactly of 80 dashes
  fprintf(ofile,
      "--------------------------------------------------------------"
      "------------------\n");
  // do it, print now!
  fflush(ofile);

  // fall asleep
  return;

}  // print_step_text()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
