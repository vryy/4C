/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for thermal dynamics

\level 3

\maintainer Christoph Meier
*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 01/12 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_expl.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 | constructor                                               dano 01/12 |
 *----------------------------------------------------------------------*/
THR::TimIntExpl::TimIntExpl(const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& tdynparams,                        //!< input parameters
    const Teuchos::ParameterList& xparams,                           //!< extra flags
    Teuchos::RCP<DRT::Discretization> actdis,                        //!< current discretisation
    Teuchos::RCP<LINALG::Solver> solver,                             //!< the solver
    Teuchos::RCP<IO::DiscretizationWriter> output                    //!< the output
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
    PrintStepText(stdout);
  }

  if (printerrfile_)
  {
    PrintStepText(errfile_);
  }

  // fall asleep
  return;

}  // PrintStep()


/*----------------------------------------------------------------------*
 | print step summary                                        dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExpl::PrintStepText(FILE* ofile)
{
  fprintf(ofile,
      "Finalised: step %6d"
      " | nstep %6d"
      " | time %-14.8E"
      " | dt %-14.8E"
      " | numiter %3d\n",
      //     " | wct %-14.8E\n",
      step_, stepmax_, (*time_)[0], (*dt_)[0], 0
      //       timer_->ElapsedTime()
  );

  // print a beautiful line made exactly of 80 dashes
  fprintf(ofile,
      "--------------------------------------------------------------"
      "------------------\n");
  // do it, print now!
  fflush(ofile);

  // fall asleep
  return;

}  // PrintStepText()


/*----------------------------------------------------------------------*/
