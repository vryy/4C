/*----------------------------------------------------------------------------*/
/*!

\brief Line search based on polynomials

\level 3

\maintainer Matthias Mayr
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
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "linesearch_polynomial.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchPolynomial::LineSearchPolynomial() : NLNSOL::LineSearchBase(), itermax_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchPolynomial::Setup()
{
  // make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // fill member variables
  itermax_ = MyGetParameter<int>("line search: max number of recursive polynomials");

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchPolynomial::ComputeLSParam(double& lsparam, bool& suffdecr) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::LineSearchPolynomial::ComputeLSParam");
  Teuchos::TimeMonitor monitor(*time);

  // create formatted output stream
  Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  Teuchos::OSTab tab(out);

  int err = 0;

  // make sure that Init() and Setup() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // the line search parameter
  lsparam = 1.0;
  double lsparamold = 1.0;

  // try a full step first
  Teuchos::RCP<Epetra_MultiVector> xnew =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
  if (err != 0)
  {
    dserror("Failed.");
  }

  // check for sufficient decrease
  Teuchos::RCP<Epetra_MultiVector> residual =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  ComputeF(*xnew, *residual);
  double fnorm2fullstep = 1.0e+12;
  bool converged = ConvergenceCheck(*residual, fnorm2fullstep);
  suffdecr = IsSufficientDecrease(fnorm2fullstep, lsparam);

  if (suffdecr)
  {
    if (getVerbLevel() > Teuchos::VERB_HIGH)
    {
      *getOStream() << LabelShort() << ": lsparam = " << lsparam << " after full step" << std::endl;
    }
    return;
  }
  else
  {
    // try half step
    lsparamold = lsparam;
    lsparam = 0.5;
    err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
    if (err != 0)
    {
      dserror("Failed.");
    }

    // check for sufficient decrease
    ComputeF(*xnew, *residual);
    double fnorm2halfstep = 1.0e+12;
    converged = ConvergenceCheck(*residual, fnorm2halfstep);
    suffdecr = IsSufficientDecrease(fnorm2halfstep, lsparam);

    if (converged or suffdecr)
    {
      if (getVerbLevel() > Teuchos::VERB_LOW)
      {
        *getOStream() << LabelShort() << ": lsparam = " << lsparam << " after half step"
                      << std::endl;
      }
      return;
    }
    else
    {
      // build polynomial model
      int iter = 0;

      // define three data points for a quadratic model
      double l1 = 0.0;
      double l2 = 1.0;
      double l3 = lsparam;

      double y1 = GetFNormOld();   // value at l1
      double y2 = fnorm2fullstep;  // value at l2
      double y3 = fnorm2halfstep;  // value at l3

      //      std::cout << "x_i\ty_i" << std::endl
      //                << l1 << "\t" << y1 << std::endl
      //                << l2 << "\t" << y2 << std::endl
      //                << l3 << "\t" << y3 << std::endl;

      double fnorm2 = y3;
      double a = 0.0;
      double b = 0.0;

      while (not converged and not suffdecr and iter < itermax_)
      {
        ++iter;

        // compute coefficients of 2nd order polynomial
        a = y3 - (y2 - y1) / (l2 - l1) * l3 - y1 + l1 * (y2 - y1) / (l2 - l1);
        a = a / (l3 * l3 - l3 * (l2 * l2 - l1 * l1) / (l2 - l1) - l1 * l1 +
                    l1 * (l2 * l2 - l1 * l1) / (l2 - l1));

        b = (y2 - y1) / (l2 - l1) - a * (l2 * l2 - l1 * l1) / (l2 - l1);

        /* Coefficient c is not needed for calculation of the minimizer of the
         * quadratic polynomial, but only for visualization of the polynomial.*/
        // const double c = y1 - a*l1*l1 - b*l1;

        // calculate new line search parameter as minimizer of polynomial model
        lsparamold = lsparam;
        if (a > 0.0)  // cf. [Kelley (1995), p. 143]
          lsparam = -b / (2 * a);
        else
          lsparam = 0.5 * lsparamold;

        // safeguard strategy
        Safeguard(lsparam, lsparamold);

        // take trial step
        err = xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);
        if (err != 0)
        {
          dserror("Failed.");
        }

        // check for sufficient decrease
        ComputeF(*xnew, *residual);
        converged = ConvergenceCheck(*residual, fnorm2);
        suffdecr = IsSufficientDecrease(fnorm2, lsparam);

        /* update interpolation points
         *
         * Note: Following [Kelley (2003), p. 12], we keep the value at
         * lsparam = 0.0 and just update those with lsparam > 0.0. */
        // l1 = l2;
        l2 = l3;
        l3 = lsparam;
        // y1 = y2;
        y2 = y3;
        y3 = fnorm2;
      }

      if (getVerbLevel() > Teuchos::VERB_LOW)
      {
        *getOStream() << LabelShort() << ": lsparam = " << lsparam << " after " << iter
                      << " iterations" << std::endl;
      }
    }
  }

  return;
}
