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

  //! \brief Quick and dirty hack to compute an error based on a monitor file, copied from gen_inv_analysis
  double CalcErrorMonitorFile();

  //! \brief Solve nonlinear forward problem without using adapters
  Epetra_SerialDenseVector SolveForwardProblemOldStyle(bool outputtofile,int nmp_,std::vector<double>  timesteps_, int ndofs_,int nnodes_,std::vector<std::vector<int> > dofs_,std::vector<int>   nodes_);

  //! \brief Get displacement at dofs specified in monitor file
  Epetra_SerialDenseVector GetCalculatedCurve(Epetra_Vector& disp, int ndofs_, int nnodes_,std::vector<std::vector<int> > dofs_,std::vector<int>   nodes_ );

  //! \brief vector to store displacements of discretization
  Teuchos::RCP<const Epetra_Vector> displacement_;

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
