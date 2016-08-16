/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.cpp

\brief Explicit time integration for structural dynamics

\level 1

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include <sstream>

#include "strtimint.H"
#include "strtimint_expl.H"

#include "../drt_cardiovascular0d/cardiovascular0d_manager.H"
#include "stru_aux.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/springdashpot_manager.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_contact/meshtying_contact_bridge.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl::TimIntExpl
(
  const Teuchos::ParameterList& timeparams,   //! time parameters
  const Teuchos::ParameterList& ioparams,     //!< ioflags
  const Teuchos::ParameterList& sdynparams,   //!< input parameters
  const Teuchos::ParameterList& xparams,      //!< extra flags
  Teuchos::RCP<DRT::Discretization> actdis,   //!< current discretisation
  Teuchos::RCP<LINALG::Solver> solver,        //!< the solver
  Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
  Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
)
: TimInt
  (
    timeparams,
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
  if (HaveContactMeshtying())
  {
    INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtbridge_->GetStrategy().Params(),"STRATEGY");
    if (soltype != INPAR::CONTACT::solution_penalty)
      dserror("Explicit TIS can only handle penalty contact / meshtying");
  }

  // cannot handle rotated DOFs
  if (locsysman_ != Teuchos::null)
    dserror("Explicit TIS cannot handle local co-ordinate systems");

  // explicit time integrators cannot handle nonlinear inertia forces
  if (HaveNonlinearMass())
    dserror("Explicit TIS cannot handle nonlinear inertia forces (flag: MASSLIN)");

  // get away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void STR::TimIntExpl::ApplyForceExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
  const Teuchos::RCP<Epetra_Vector> vel,  //!< velocity state
  Teuchos::RCP<Epetra_Vector>& fext  //!< external force
)
{
  Teuchos::ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0,"displacement", dis);
  discret_->SetState(0,"displacement new", dis);

  if (damping_ == INPAR::STR::damp_material)
    discret_->SetState(0,"velocity", vel);
  // get load vector
  discret_->EvaluateNeumann(p, *fext);

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntExpl::PrintStep()
{
  // print out
  if ( (myrank_ == 0) and printscreen_ and (StepOld()%printscreen_==0))
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
          " | numiter %3d"
          " | wct %-14.8E\n",
          step_, stepmax_, (*time_)[0], (*dt_)[0], 0, timer_->ElapsedTime());
  fprintf(ofile,
          "                      "
          " ( ts %-14.8E"
          " | te %-14.8E"
          " | tc %-14.8E)\n",
          dtsolve_,dtele_,dtcmt_);
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
