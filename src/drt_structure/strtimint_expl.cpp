/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.cpp
\brief Explicit time integration for spatial discretised 
       structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>

#include "strtimint.H"
#include "strtimint_expl.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl:: TimIntExpl
(
  const Teuchos::ParameterList& ioparams,  //!< ioflags
  const Teuchos::ParameterList& sdynparams,  //!< input parameters
  const Teuchos::ParameterList& xparams,  //!< extra flags
  Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
  Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
  Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
)
: TimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  )
{
  // explicit time integrators cannot handle constraints
  if (conman_->HaveConstraint())
    dserror("Explicit TIS cannot handle constraints");

  // cannot handle rotated DOFs
  if (locsysman_ != Teuchos::null)
    dserror("Explicit TIS cannot handle local co-ordinate systems");

  // get away
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntExpl::PrintStep()
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
void STR::TimIntExpl::PrintStepText
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
          step_, stepmax_, (*time_)[0], (*dt_)[0], 0);
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
