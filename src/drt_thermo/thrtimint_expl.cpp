/*----------------------------------------------------------------------*/
/*!
\file thrtimint_expl.cpp
\brief Explicit time integration for thermal dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_expl.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/* constructor */
THR::TimIntExpl:: TimIntExpl
(
  const Teuchos::ParameterList& ioparams,  //!< ioflags
  const Teuchos::ParameterList& tdynparams,  //!< input parameters
  const Teuchos::ParameterList& xparams,  //!< extra flags
  Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
  Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
  Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
)
: TimInt
  (
    ioparams,
    tdynparams,
    xparams,
    actdis,
    solver,
    output
  )
{
  // get away
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void THR::TimIntExpl::PrintStep()
{
  // print out
  if ( (myrank_ == 0) and printscreen_ )
  {
    PrintStepText(stdout);
  }

  if (printerrfile_)
  {
    PrintStepText(errfile_);
  }

  // fall asleep
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void THR::TimIntExpl::PrintStepText
(
  FILE* ofile
)
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
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
