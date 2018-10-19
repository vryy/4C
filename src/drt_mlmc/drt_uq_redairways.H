/*!----------------------------------------------------------------------
\file drt_uq_redairways.H

\brief Class for performing UQ analysis with redairways

\level 2

\maintainer Jonas Biehler
*/
/*----------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#ifndef DRT_UQ_REDAIRWAYS_H_
#define DRT_UQ_REDAIRWAYS_H_

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include "Epetra_FECrsMatrix.h"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace ADAPTER
{
  class Structure;
}

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
  class Node;
  class Element;
}  // namespace DRT

// namespace LINALG
//{
//  class SparseMatrix;
//  class SparseOperator;
//}

namespace IO
{
  class DiscretizationWriter;
  class OutputControl;
  class ErrorFileControl;
}  // namespace IO

namespace UQ
{
  class RandomField;
  class MCMatParManager;
  class MCVarThicknessManager;

  /// UQ with redairway models with
  /// some sort of stochastic property
  class UQ_REDAIRWAYS
  {
   public:
    //! \brief standard constructor
    UQ_REDAIRWAYS(Teuchos::RCP<DRT::Discretization> dis);  ///< the discretisation

    //! \brief standard monte carlo analysis
    void Integrate();


   protected:
   private:
    //! \brief my structure discretisation to work at
    Teuchos::RCP<DRT::Discretization> discret_;

    //! \brief filename of standard control file
    std::string filename_;

    //! brief starting random seed
    int start_random_seed_;

    //! brief number of monte carlo runs
    int numb_run_;

    //! brief run to start with
    int start_run_;

    //! brief number of mc runs per npar group
    int numruns_pergroup_;

    //! brief flag to reduce output to files as well as stdout
    bool reduced_output_;

    Teuchos::RCP<MCMatParManager> my_matpar_manager_;
  };  // EOF class

}  // namespace UQ

#endif  /* DRT_UQ_REDAIRWAYS_H_ */
#endif  // HAVE_FFTw
