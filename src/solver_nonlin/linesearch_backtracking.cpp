/*----------------------------------------------------------------------------*/
/*!

\brief Linesearch algorithm based on backtracking

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
#include "linesearch_backtracking.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchBacktracking::LineSearchBacktracking()
    : NLNSOL::LineSearchBase(), itermax_(4), backtrackfac_(0.5)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBacktracking::Setup()
{
  // make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // fill member variables
  itermax_ = MyGetParameter<int>("line search: max number of backtracking steps");
  backtrackfac_ = MyGetParameter<double>("line search: factor for step size reduction");

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBacktracking::ComputeLSParam(double& lsparam, bool& suffdecr) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::LineSearchBacktracking::ComputeLSParam");
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

  // start backtracking with full step
  lsparam = 1.0;

  // take a full step
  Teuchos::RCP<Epetra_MultiVector> xnew =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);

  // check for sufficient decrease
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(xnew->Map(), true));
  ComputeF(*xnew, *fnew);
  double fnorm2 = 1.0e+12;
  bool converged = ConvergenceCheck(*fnew, fnorm2);
  suffdecr = IsSufficientDecrease(fnorm2, lsparam);

  int iter = 0;  // iteration index for multiple backtracking steps

  while (not converged and not suffdecr and iter < itermax_)
  {
    ++iter;

    // halve trial line search parameter
    lsparam *= backtrackfac_;

    //    *out << "lsparam = " << lsparam;// << std::endl;

    // take a reduced step
    xnew->Update(1.0, GetXOld(), lsparam, GetXInc(), 0.0);

    // check for sufficient decrease
    ComputeF(*xnew, *fnew);
    converged = ConvergenceCheck(*fnew, fnorm2);
    suffdecr = IsSufficientDecrease(fnorm2, lsparam);

    //    *out << "\tfnorm2 = " << fnorm2
    //         << "\tinitnorm = " << GetFNormOld()
    //         << std::endl;
  }

  if (getVerbLevel() > Teuchos::VERB_LOW)
  {
    *getOStream() << LabelShort() << ": lsparam = " << lsparam << " after " << iter
                  << " iterations." << std::endl;
  }

  return;
}
