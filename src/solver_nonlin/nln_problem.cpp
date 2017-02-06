/*----------------------------------------------------------------------------*/
/*!
\file nln_problem.cpp

<pre>
\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>

\level 3
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard

// Epetra
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>

// NOX
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Vector.H>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "nln_problem.H"
#include "nln_utils.H"
#include "nln_utils_debugwriter.H"

#include "../drt_io/io.H"

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblem::NlnProblem()
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec,
    const std::string& description,
    const IO::VectorType vt) const
{
  if (HaveDebugWriter())
  {
    DebugWriter()->WriteVector(vec, description, vt);
  }
  else
  {
    if (getVerbLevel() > Teuchos::VERB_NONE)
    {
      *getOStream() << Label()
          << ": WARNING: Cant't write debug output of vector '"
          << description << "', since debug writer 'dbgwriter_' has not been "
              "set properly, yet."
          << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::WriteVector(const Epetra_MultiVector& vec,
    const std::string& description,
    const IO::VectorType vt) const
{
  WriteVector(Teuchos::rcp(&vec, false), description, vt);

  return;
}
