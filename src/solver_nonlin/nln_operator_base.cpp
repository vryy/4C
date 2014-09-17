/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_base.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorBase::NlnOperatorBase()
: isinit_(false),
  issetup_(false),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  nlnproblem_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Initialization */
void NLNSOL::NlnOperatorBase::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem
    )
{
  // We need to call Setup() after Init()
  issetup_ = false;

  // fill member variables
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  nlnproblem_ = nlnproblem;

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/* Access to communicator */
const Epetra_Comm& NLNSOL::NlnOperatorBase::Comm() const
{
  if (comm_.is_null())
    dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/* Access to parameter list */
const Teuchos::ParameterList& NLNSOL::NlnOperatorBase::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}

/*----------------------------------------------------------------------------*/
/* Access to nonlinear problem */
Teuchos::RCP<NLNSOL::NlnProblem> NLNSOL::NlnOperatorBase::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror("The nonlinear problem 'nlnproblem_' has not been initialized, yet.");

  return nlnproblem_;
}

/*----------------------------------------------------------------------------*/
/* Print summary of current iteration */
void NLNSOL::NlnOperatorBase::PrintIterSummary(const int iter,
    const double fnorm2
    ) const
{
  // print only on one processor
  if (Comm().MyPID() == 0 && Params().get<bool>("Nonlinear Operator: Print Iterations"))
  {
    IO::cout << std::setprecision(6)
             << Label() << " iteration " << iter
             << ": |f| = " << fnorm2
             << IO::endl;
  }
}

/*----------------------------------------------------------------------------*/
/* Check if we still continue to iterate */
const bool NLNSOL::NlnOperatorBase::ContinueIterations(const int iter,
    const bool converged) const
{
  // initialize return value
  bool retval = false;

  if (iter < GetMaxIter() && not converged)
    retval = true;

  return retval;
}

/*----------------------------------------------------------------------------*/
/* Check if iteration loop converged successfully */
const bool NLNSOL::NlnOperatorBase::CheckSuccessfulConvergence(
    const int iter,
    const bool converged
    ) const
{
  // initialize return value
  bool successful = false;

  // make the decision
  if (IsSolver())
  {
    if (converged && iter < GetMaxIter()) // successful convergence
      return true;
    else // convergence failed
    {
      dserror("%s did not converge in %d iterations.", Label(), iter);
      return false;
    }
  }
  else // no solver, so we don't care for convergence
  {
    return true;
  }

  return successful;
}
