/*----------------------------------------------------------------------------*/
/*!
\file linesearch_linear.cpp

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

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "linesearch_linear.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchLinear::LineSearchLinear()
 : NLNSOL::LineSearchBase()
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchLinear::Setup()
{
  // make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::LineSearchLinear::ComputeLSParam() const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::LineSearchLinear::ComputeLSParam");
  Teuchos::TimeMonitor monitor(*time);

  // make sure that Init() and Setup() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  double nominator = 0.0;
  double denominator = 0.0;

  Teuchos::RCP<Epetra_MultiVector> xnew = // new trial solution
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  Teuchos::RCP<Epetra_MultiVector> fnew = // residual at trial solution
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));

  // compute nominator
  GetFOld().Dot(GetXInc(), &nominator);

  // compute denominator
  xnew->Update(1.0, GetXOld(), 1.0, GetXInc(), 0.0);
  ComputeF(*xnew, *fnew);
  fnew->Dot(GetXInc(), &denominator);
  denominator -= nominator;

  /* return line search parameter with negative sign to account for ComputeF()
   * delivering a descending residual
   */
  return -nominator / denominator;
}
