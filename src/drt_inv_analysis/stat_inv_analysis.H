/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.H
 *
<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANALYSIS_H
#define STAT_INV_ANALYSIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

// forward declarations

namespace COMM_UTILS
{
  class NestedParGroup;
}

namespace STR
{
  class TimInt;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
}

namespace DRT
{
  class Condition;
  class Discretization;
}

namespace STR
{


class StatInvAnalysis {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;




public:

  /// standard constructor
  StatInvAnalysis(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  //
  void Sample();
  //
  void SolveForwardProblem();

  //
  void CreateNewSample();

  //
  void SetParameters();


};  // class StatInvAnalysis

}  // namespace STR

#endif /*INV_ANALYSIS_H_*/
