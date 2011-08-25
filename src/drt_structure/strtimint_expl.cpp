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
#include "../drt_constraint/constraint_manager.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_contact.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl:: TimIntExpl
(
  const Teuchos::ParameterList& ioparams,  //!< ioflags
  const Teuchos::ParameterList& sdynparams,  //!< input parameters
  const Teuchos::ParameterList& xparams,  //!< extra flags
  Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
  Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
  Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
  Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
)
: TimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  )
{
  // explicit time integrators cannot handle constraints
  if (conman_->HaveConstraint())
    dserror("Explicit TIS cannot handle constraints");

  // explicit time integrators can only handle penalty contact / meshtying
  if (cmtman_ != Teuchos::null)
  {
    INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
    if (soltype != INPAR::CONTACT::solution_penalty)
      dserror("Explicit TIS can only handle penalty contact / meshtying");
  }

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
/*----------------------------------------------------------------------*/
bool STR::TimIntExpl::UseContactSolver()
{
  // no contact possible -> return false
  if (!HaveContactMeshtying())
    return false;
  // contact possible -> check current status
  else
  {
    // currently not in contact -> return false
    if (!cmtman_->GetStrategy().IsInContact() &&
        !cmtman_->GetStrategy().WasInContact() &&
        !cmtman_->GetStrategy().WasInContactLastTimeStep())
      return false;
    // currently in contact -> return true
    else
      return true;
  }
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
          " | numiter %3d"
          " | wct %-14.8E\n",
          step_, stepmax_, (*time_)[0], (*dt_)[0], 0, timer_->ElapsedTime());
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
