/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem.H

\brief This class manages some of the necessary factory calls
       if a %NOX::NLN solver is supposed to be used. Therefore a
       lean function call becomes possible.

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PROBLEM_H_
#define NOX_NLN_PROBLEM_H_

#include <Teuchos_RCP.hpp>

#include "nox_nln_enum_lists.H"
#include "../drt_lib/drt_dserror.H"

#include <NOX_StatusTest_Generic.H>

// forward declaration
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
  class Scaling;
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
  class GlobalData;
namespace INNER {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
} // namespace INNER

  class Problem
  {
  public:

    //! minimal constructor
    Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData);

    //! standard constructor
    Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData,
        const Teuchos::RCP<NOX::Epetra::Vector>& x,
        const Teuchos::RCP<LINALG::SparseOperator>& A);

    //! destructor
    virtual ~Problem(){};

    //! initialize stuff (can be overloaded in derived classes)
    virtual void Initialize(const Teuchos::RCP<NOX::Epetra::Vector>& x,
        const Teuchos::RCP<LINALG::SparseOperator>& A);

    //! create the linear system for the NOX framework
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem() const;

    //! create a nox group
    virtual Teuchos::RCP<NOX::Abstract::Group> CreateGroup(
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys
        ) const;

    void CreateOuterStatusTest(
        Teuchos::RCP<NOX::StatusTest::Generic>& outerTests) const;

    virtual void CreateStatusTests(
        Teuchos::RCP<NOX::StatusTest::Generic>& outerTest,
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTest) const;

    //! check final status of the non-linear solving procedure
    virtual void CheckFinalStatus(const NOX::StatusTest::StatusType& finalStatus) const;

    /// access the global data object
    NOX::NLN::GlobalData & NlnGlobalData()
    {
      return * noxNlnGlobalData_;
    }

    /// access the global data object ptr
    Teuchos::RCP<NOX::NLN::GlobalData> NlnGlobalDataPtr()
    {
      return noxNlnGlobalData_;
    }

  protected:
    inline void CheckInit() const
    {
      if ( not isinit_ )
        dserror("You have to call Initialize() first, before you can use this"
            " function!");
    }

    inline const bool & IsJac() const { return isjac_; };

  protected:
    bool isinit_;

    bool isjac_;

    Teuchos::RCP<NOX::NLN::GlobalData> noxNlnGlobalData_;

    /** ptr to the state vector RCP. In this way the strong_count is neither lost
     *  nor increased. */
    const Teuchos::RCP<NOX::Epetra::Vector> * xVector_;

    /** ptr to the state matrix RCP. In this way the strong_count is neither lost
     *  nor increased. */
    const Teuchos::RCP<LINALG::SparseOperator> * jac_;

    Teuchos::RCP<LINALG::SparseOperator> precMat_;

    Teuchos::RCP<NOX::Epetra::Scaling> scalingObject_;

  };
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_PROBLEM_H_ */
