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
#include "linesearch_linear.H"

#include "../drt_io/io_pstream.H"
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
  int err = 0;

  // make sure that Init() and Setup() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  double nominator = 0.0;
  double denominator = 0.0;

  // compute nominator
  Teuchos::RCP<Epetra_MultiVector> tmp =
      Teuchos::rcp(new Epetra_MultiVector(GetXOld().Map(), true));
  ComputeF(GetXOld(), *tmp);
  err = GetXInc().Dot(*tmp, &nominator);
  if (err != 0) { dserror("Dot product failed."); }

  // compute denominator
//  GetJacobianOperator()->Apply(*tmp, *tmp);
  err = GetXInc().Dot(*tmp, &denominator);
  if (err != 0) { dserror("Dot product failed."); }

  // return line search parameter
  return nominator / denominator;
}



