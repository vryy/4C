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

// baci
#include "linesearch_fullstep.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchFullStep::LineSearchFullStep()
 : NLNSOL::LineSearchBase()
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchFullStep::Setup()
{
  // make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // SetupLineSearch() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::LineSearchFullStep::ComputeLSParam() const
{
  // make sure that Init() and Setup() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  return 1.0;
}
