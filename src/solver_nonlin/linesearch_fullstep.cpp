/*----------------------------------------------------------------------------*/
/*!
\file linesearch_fullstep.cpp

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
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "linesearch_fullstep.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchFullStep::LineSearchFullStep() : NLNSOL::LineSearchBase() { return; }

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchFullStep::Setup()
{
  // make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchFullStep::ComputeLSParam(double& lsparam, bool& suffdecr) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::LineSearchFullStep::ComputeLSParam");
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

  // full step without caring for sufficient decrease
  ComputeLSParam(lsparam);

  // take the full step
  Teuchos::RCP<Epetra_MultiVector> xnew =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);

  // check for sufficient decrease
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(xnew->Map(), true));
  ComputeF(*xnew, *fnew);
  double fnorm2 = 0.0;
  ConvergenceCheck(*fnew, fnorm2);
  suffdecr = IsSufficientDecrease(fnorm2, lsparam);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchFullStep::ComputeLSParam(double& lsparam) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::LineSearchFullStep::ComputeLSParam");
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

  // full step
  lsparam = 1.0;

  return;
}
