/*----------------------------------------------------------------------------*/
/*!
\file linesearch_backtracking.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <iostream>

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

// NOX
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_MultiVector.H>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_backtracking.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::LineSearchBacktracking::LineSearchBacktracking()
 : NLNSOL::LineSearchBase(),
   trialstepsize_(1.0)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of line search object */
void NLNSOL::LineSearchBacktracking::Setup()
{
  // make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  trialstepsize_ = Params().sublist("Backtracking").get<double>("trial step size");

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Compute the line search parameter */
const double NLNSOL::LineSearchBacktracking::ComputeLSParam() const
{
  // make sure that Init() and Setup() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // the line search parameter
  double lsparam = 1.0;

  // take a full step
  Teuchos::RCP<Epetra_MultiVector> xnew = Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);

  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(xnew->Map(), true));
  Evaluate(*xnew, *fnew);

  double fnorm2 = 1.0e+12;
  ConvergenceCheck(*fnew, fnorm2);

  while (not IsSufficientDecrease(fnorm2, lsparam) && lsparam > 0.0)
  {
    // reduce trial line search parameter
    lsparam -= trialstepsize_;

    std::cout << "lsparam = " << lsparam;// << std::endl;

    // take a reduced step
    xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
    Evaluate(*xnew, *fnew);

    // evaluate and compute L2-norm of residual
    ConvergenceCheck(*fnew, fnorm2);

    std::cout << "\tfnorm2 = " << fnorm2 << "\tinitnorm = " << GetFNormOld() << std::endl;
  }

  // check for successful line search
  if (not IsSufficientDecrease(fnorm2, lsparam))
  {
    dserror("Sufficient decrease condition could not be satisfied.");

    lsparam = 0.0;
  }

  return lsparam;
}



