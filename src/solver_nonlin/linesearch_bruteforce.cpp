/*----------------------------------------------------------------------------*/
/*!
\file linesearch_bruteforce.cpp

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
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "linesearch_bruteforce.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchBruteForce::LineSearchBruteForce() : NLNSOL::LineSearchBase(), trialstepsize_(1.0)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBruteForce::Setup()
{
  // make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  trialstepsize_ = MyGetParameter<double>("line search: trial step size");

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBruteForce::ComputeLSParam(double& lsparam, bool& suffdecr) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::LineSearchBruteForce::ComputeLSParam");
  Teuchos::TimeMonitor monitor(*time);

  // make sure that Init() and Setup() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // ---------------------------------------------------------------------------
  // start with a full step
  // ---------------------------------------------------------------------------
  lsparam = 1.0;
  Teuchos::RCP<Epetra_MultiVector> xnew =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);

  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(xnew->Map(), true));
  ComputeF(*xnew, *fnew);

  double fnorm2 = 1.0e+12;
  ConvergenceCheck(*fnew, fnorm2);
  suffdecr = IsSufficientDecrease(fnorm2, lsparam);

  // ---------------------------------------------------------------------------
  // reduce step length in increments of #trialstepsize_
  // ---------------------------------------------------------------------------
  while (not suffdecr && lsparam > 0.0)
  {
    // reduce trial line search parameter
    lsparam -= trialstepsize_;

    if (getVerbLevel() > Teuchos::VERB_MEDIUM)
      *getOStream() << "Trying lsparam = " << lsparam << std::endl;

    // take a reduced step
    xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
    ComputeF(*xnew, *fnew);

    // check for sufficient decrease
    ConvergenceCheck(*fnew, fnorm2);
    suffdecr = IsSufficientDecrease(fnorm2, lsparam);

    if (getVerbLevel() > Teuchos::VERB_HIGH)
      *getOStream() << LabelShort() << ": lsparam = " << lsparam << std::endl;
  }

  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << LabelShort() << ": lsparam = " << lsparam << std::endl;

  // check for successful line search
  if (not suffdecr)
  {
    dserror("Sufficient decrease condition could not be satisfied.");

    lsparam = 0.0;
  }

  return;
}
