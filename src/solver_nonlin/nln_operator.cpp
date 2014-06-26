/*----------------------------------------------------------------------*/
/*!
\file nln_operator.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

// NOX

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperator::NlnOperator()
: isinit_(false),
  issetup_(false),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  nlnproblem_(Teuchos::null),
  issolver_(false)
{
  return;
}

/*----------------------------------------------------------------------*/
/* Initialization */
void NLNSOL::NlnOperator::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem
    )
{
  // fill member variables
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  nlnproblem_ = nlnproblem;

  if (Params().isParameter("Nonlinear Operator: Is Solver"))
    issolver_ = Params().get<bool>("Nonlinear Operator: Is Solver");

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------*/
/* Access to communicator */
const Epetra_Comm& NLNSOL::NlnOperator::Comm() const
{
  if (comm_.is_null())
    dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------*/
/* Access to parameter list */
const Teuchos::ParameterList& NLNSOL::NlnOperator::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}

/*----------------------------------------------------------------------*/
/* Access to nonlinear problem */
Teuchos::RCP<NLNSOL::NlnProblem> NLNSOL::NlnOperator::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror("The nonlinear problem 'nlnproblem_' has not been initialized, yet.");

  return nlnproblem_;
}
