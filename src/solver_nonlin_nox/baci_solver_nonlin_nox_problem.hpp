/*-----------------------------------------------------------*/
/*! \file

\brief This class manages some of the necessary factory calls
       if a %NOX::NLN solver is supposed to be used. Therefore a
       lean function call becomes possible.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_PROBLEM_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_PROBLEM_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_enum_lists.hpp"
#include "baci_solver_nonlin_nox_forward_decl.hpp"
#include "baci_utils_exceptions.hpp"

#include <NOX_StatusTest_Generic.H>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declaration
namespace CORE::LINALG
{
  class Solver;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace NOX
{
  namespace NLN
  {
    class GlobalData;
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER

    class Problem
    {
     public:
      //! minimal constructor
      Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData);

      //! standard constructor
      Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData,
          const Teuchos::RCP<::NOX::Epetra::Vector>& x,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& A);

      //! destructor
      virtual ~Problem() = default;

      //! initialize stuff (can be overloaded in derived classes)
      virtual void Initialize(const Teuchos::RCP<::NOX::Epetra::Vector>& x,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& A);

      //! create the linear system for the NOX framework
      virtual Teuchos::RCP<::NOX::Epetra::LinearSystem> CreateLinearSystem() const;

      //! create a nox group
      virtual Teuchos::RCP<::NOX::Abstract::Group> CreateGroup(
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys) const;

      void CreateOuterStatusTest(Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests) const;

      virtual void CreateStatusTests(Teuchos::RCP<::NOX::StatusTest::Generic>& outerTest,
          Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTest) const;

      //! check final status of the non-linear solving procedure
      virtual void CheckFinalStatus(const ::NOX::StatusTest::StatusType& finalStatus) const;

      /// access the global data object
      NOX::NLN::GlobalData& NlnGlobalData() { return *noxNlnGlobalData_; }

      /// access the global data object ptr
      Teuchos::RCP<NOX::NLN::GlobalData> NlnGlobalDataPtr() { return noxNlnGlobalData_; }

     protected:
      inline void CheckInit() const
      {
        if (not isinit_)
          dserror(
              "You have to call Initialize() first, before you can use this"
              " function!");
      }

      inline const bool& IsJac() const { return isjac_; };

     protected:
      bool isinit_;

      bool isjac_;

      Teuchos::RCP<NOX::NLN::GlobalData> noxNlnGlobalData_;

      /** ptr to the state vector RCP. In this way the strong_count is neither lost
       *  nor increased. */
      const Teuchos::RCP<::NOX::Epetra::Vector>* xVector_;

      /** ptr to the state matrix RCP. In this way the strong_count is neither lost
       *  nor increased. */
      const Teuchos::RCP<CORE::LINALG::SparseOperator>* jac_;

      Teuchos::RCP<CORE::LINALG::SparseOperator> precMat_;
    };
  }  // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_PROBLEM_H
