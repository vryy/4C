/*----------------------------------------------------------------------------*/
/*!
\file linesearch_polynomial.cpp

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

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_polynomial.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::LineSearchPolynomial::LineSearchPolynomial()
 : NLNSOL::LineSearchBase(),
   itermax_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of line search object */
void NLNSOL::LineSearchPolynomial::Setup()
{
  // make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // max number of recursive polynomials
  itermax_ = Params().sublist("Polynomial2").get<int>("max number of recursive polynomials");

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Compute the line search parameter */
const double NLNSOL::LineSearchPolynomial::ComputeLSParam() const
{
  int err = 0;

  // make sure that Init() and Setup() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // the line search parameter
  double lsparam = 1.0;
  double lsparamold = 1.0;

  // try a full step first
  Teuchos::RCP<Epetra_MultiVector> xnew = Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
  if (err != 0) { dserror("Failed."); }

  Teuchos::RCP<Epetra_MultiVector> residual = Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  Evaluate(*xnew, *residual);

  double fnorm2fullstep = 1.0e+12;
  ConvergenceCheck(*residual, fnorm2fullstep);

  if (IsSufficientDecrease(fnorm2fullstep, lsparam))
    return lsparam;
  else
  {
    // try half step
    lsparamold = lsparam;
    lsparam = 0.5;

    err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
    if (err != 0) { dserror("Failed."); }

    Evaluate(*xnew, *residual);

    double fnorm2halfstep = 1.0e+12;
    ConvergenceCheck(*residual, fnorm2halfstep);

    if (IsSufficientDecrease(fnorm2halfstep, lsparam))
      return lsparam;
    else
    {
      // build polynomial model
      int iter = 0;

      double l1 = 0.0;
      double l2 = 1.0;
      double l3 = 0.5;

      double y1 = GetFNormOld();
      double y2 = fnorm2fullstep;
      double y3 = fnorm2halfstep;

      double fnorm2 = y3;
      double a = 0.0;
      double b = 0.0;

      while (not IsSufficientDecrease(fnorm2, lsparam) && iter < itermax_)
      {
        ++iter;

        // compute coefficients of 2nd order polynomial
        a = y3 - (y2-y1)/(l2-l1)*l3 - y1 + l1*(y2-y1)/(l2-l1);
        a = a / (l3*l3 - l3*(l2*l2-l1*l1)/(l2-l1) - l1*l1 + l1*(l2*l2-l1*l1)/(l2-l1));

        b = (y2-y1)/(l2-l1) - a*(l2*l2-l1*l1)/(l2-l1);

        lsparamold = lsparam;
        lsparam = - b / (2*a);

        // safeguard strategy
        Safeguard(lsparam, lsparamold);

        err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
        if (err != 0) { dserror("Failed."); }
        Evaluate(*xnew, *residual);

        ConvergenceCheck(*residual, fnorm2);

        l1 = l2;
        l2 = l3;
        l3 = lsparam;

        y1 = y2;
        y2 = y3;
        y3 = fnorm2;
      }

      if (not IsSufficientDecrease(fnorm2, lsparam) && iter > itermax_)
        dserror("Polynomial line search cannot satisfy sufficient decrease condition.");

      return lsparam;
    }
  }
}
