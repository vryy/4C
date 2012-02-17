/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.cpp
\brief Explicit time integration for structural dynamics

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

#include "strtimint.H"
#include "strtimint_expl.H"
#include "stru_aux.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl::TimIntExpl
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
  ),
  fsisurface_(NULL),
  fifc_(Teuchos::null)
{
  // explicit time integrators cannot handle constraints
  if (conman_->HaveConstraint())
    dserror("Explicit TIS cannot handle constraints");

  // explicit time integrators can only handle penalty contact / meshtying
  if (HaveContactMeshtying())
  {
    INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cmtman_->GetStrategy().Params(),"STRATEGY");
    if (soltype != INPAR::CONTACT::solution_penalty)
      dserror("Explicit TIS can only handle penalty contact / meshtying");
  }

  // create empty interface force vector
  fifc_ = LINALG::CreateVector(*dofrowmap_, true);

  // cannot handle rotated DOFs
  if (locsysman_ != Teuchos::null)
    dserror("Explicit TIS cannot handle local co-ordinate systems");

  // get away
  return;
}


/*----------------------------------------------------------------------*/
/* introduce (robin) fsi surface extractor object */
void STR::TimIntExpl::SetSurfaceFSI
(
  const STR::AUX::MapExtractor* fsisurface  //!< the FSI surface
)
{
  fsisurface_ = fsisurface;
}


/*----------------------------------------------------------------------*/
/* Set forces due to interface with fluid */
void STR::TimIntExpl::SetForceInterface
(
  const STR::AUX::MapExtractor& extractor,
  Teuchos::RCP<Epetra_Vector> iforce  ///< the force on interface
)
{
  fifc_->PutScalar(0.0);
  extractor.AddFSICondVector(iforce, fifc_);
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
