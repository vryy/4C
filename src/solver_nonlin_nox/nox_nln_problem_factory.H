/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem_factory.H

\maintainer Michael Hiermeier

\date Jul 23, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PROBLEM_FACTORY_H_
#define NOX_NLN_PROBLEM_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include <NOX_Common.H>

// Forward declarations
class Epetra_Vector;

namespace LINALG{
  class SparseOperator;
}

namespace NOX {
namespace NLN {
  class NoxProblem;
  class GlobalData;
namespace Problem {
  class Factory
  {
  public:
    //! Constructor
    Factory();

    //! Destructor
    ~Factory() {};

    //! Factory to build a NoxProblem object
    Teuchos::RCP<NOX::NLN::NoxProblem> BuildNoxProblem(
        Teuchos::RCP<Epetra_Vector>& x,
        Teuchos::RCP<LINALG::SparseOperator> A,
        Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData);

  };

  /*! Nonmember function to build a NoxProblem object.

  \relates NOX::NLN::Problem::Factory

  */
  Teuchos::RCP<NOX::NLN::NoxProblem> BuildNoxProblem(
          Teuchos::RCP<Epetra_Vector>& x,
          Teuchos::RCP<LINALG::SparseOperator> A,
          Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData);

} // namespace Problem
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_PROBLEM_FACTORY_H_ */
