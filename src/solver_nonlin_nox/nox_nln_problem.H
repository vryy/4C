/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem.H

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PROBLEM_H_
#define NOX_NLN_PROBLEM_H_

#include <Teuchos_RCP.hpp>

#include "nox_nln_globaldata.H"

#include <NOX_StatusTest_Generic.H>

// forward declaration
class Epetra_Vector;

namespace Teuchos {
  class ParameterList;
}

namespace LINALG {
  class Solver;
  class SparseOperator;
}

namespace NOX {
  class Utils;
namespace Abstract {
  class Group;
} // namespace Abstract
namespace Epetra {
  class LinearSystem;
  class Vector;
namespace Interface {
    class Required;
    class Jacobian;
    class Preconditioner;
} // namespace Interface
} // namespace Epetra
namespace StatusTest {
  class Generic;
} // namespace  StatusTest
namespace NLN {
namespace INNER {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
} // namespace INNER

  class NoxProblem
  {
  public:
    //! constructor
    NoxProblem(Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData,
        Teuchos::RCP<Epetra_Vector>& x,
        Teuchos::RCP<LINALG::SparseOperator> A);

    //! destructor
    virtual ~NoxProblem(){};

    //! initialize stuff (can be overloaded in derived classes)
    virtual void Initialize(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<LINALG::SparseOperator> A);

    //! create the linear system for the NOX framework
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem() const;

    //! create a nox group
    virtual Teuchos::RCP<NOX::Abstract::Group> CreateNoxGroup(
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys
        ) const;

    virtual void CreateStatusTests(
        Teuchos::RCP<NOX::StatusTest::Generic>& outerTest,
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTest) const;

    //! check final status of the non-linear solving procedure
    virtual int CheckFinalStatus(const NOX::StatusTest::StatusType& finalStatus,
        const bool& isAdaptiveTimeStep=false) const;

  private:
    //! throw error
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  protected:
    Teuchos::RCP<NOX::NLN::GlobalData> nlnGlobalData_;

    Teuchos::RCP<NOX::Epetra::Vector> xVector_;

    Teuchos::RCP<LINALG::SparseOperator> jac_;
  };
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_PROBLEM_H_ */
