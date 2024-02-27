/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for structural dynamics

\level 1

*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_structure_timint_expl.hpp"

#include "baci_cardiovascular0d_manager.hpp"
#include "baci_constraint_manager.hpp"
#include "baci_constraint_springdashpot_manager.hpp"
#include "baci_contact_meshtying_contact_bridge.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_mortar_manager_base.hpp"
#include "baci_mortar_strategy_base.hpp"
#include "baci_structure_aux.hpp"
#include "baci_structure_timint.hpp"

#include <sstream>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl::TimIntExpl(const Teuchos::ParameterList& timeparams,  //! time parameters
    const Teuchos::ParameterList& ioparams,                            //!< ioflags
    const Teuchos::ParameterList& sdynparams,                          //!< input parameters
    const Teuchos::ParameterList& xparams,                             //!< extra flags
    Teuchos::RCP<DRT::Discretization> actdis,                          //!< current discretisation
    Teuchos::RCP<CORE::LINALG::Solver> solver,                         //!< the solver
    Teuchos::RCP<CORE::LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
    Teuchos::RCP<IO::DiscretizationWriter> output      //!< the output
    )
    : TimInt(timeparams, ioparams, sdynparams, xparams, actdis, solver, contactsolver, output)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntExpl::Init(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver)
{
  // call Init() in base class
  STR::TimInt::Init(timeparams, sdynparams, xparams, actdis, solver);

  // get away
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntExpl::Setup()
{
  // call Setup() in base class
  STR::TimInt::Setup();

  // explicit time integrators cannot handle constraints
  if (conman_->HaveConstraint())
    dserror("Currently, constraints cannot be done with explicit time integration.");

  // explicit time integrators can only handle penalty contact / meshtying
  if (HaveContactMeshtying())
  {
    INPAR::CONTACT::SolvingStrategy soltype =
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
            cmtbridge_->GetStrategy().Params(), "STRATEGY");
    if (soltype != INPAR::CONTACT::solution_penalty &&
        (soltype != INPAR::CONTACT::solution_multiscale))
      dserror(
          "Currently, only penalty or multi-scale contact / meshtying can be done with explicit "
          "time integration schemes.");
  }

  // cannot handle rotated DOFs
  if (locsysman_ != Teuchos::null)
    dserror("Explicit time integration schemes cannot handle local co-ordinate systems");

  // explicit time integrators cannot handle nonlinear inertia forces
  if (HaveNonlinearMass())
    dserror(
        "Explicit time integration schemes cannot handle nonlinear inertia forces (flag: MASSLIN)");

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void STR::TimIntExpl::ApplyForceExternal(const double time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector> dis,                   //!< displacement state
    const Teuchos::RCP<Epetra_Vector> vel,                   //!< velocity state
    Teuchos::RCP<Epetra_Vector>& fext                        //!< external force
)
{
  Teuchos::ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0, "displacement", dis);
  discret_->SetState(0, "displacement new", dis);

  if (damping_ == INPAR::STR::damp_material) discret_->SetState(0, "velocity", vel);
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
  if ((myrank_ == 0) and printscreen_ and (StepOld() % printscreen_ == 0))
  {
    PrintStepText(stdout);
  }
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntExpl::PrintStepText(FILE* ofile)
{
  fprintf(ofile,
      "Finalised: step %6d"
      " | nstep %6d"
      " | time %-14.8E"
      " | dt %-14.8E"
      " | numiter %3d"
      " | wct %-14.8E\n",
      step_, stepmax_, (*time_)[0], (*dt_)[0], 0, timer_->totalElapsedTime(true));
  fprintf(ofile,
      "                      "
      " ( ts %-14.8E"
      " | te %-14.8E"
      " | tc %-14.8E)\n",
      dtsolve_, dtele_, dtcmt_);
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

BACI_NAMESPACE_CLOSE
