/*-----------------------------------------------------------*/
/*!
\file loca_nln_problem.H

\maintainer Michael Hiermeier

\date Nov 20, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_LOCA_CONTINUATION_LOCA_NLN_PROBLEM_H_
#define SRC_LOCA_CONTINUATION_LOCA_NLN_PROBLEM_H_

#include "../solver_nonlin_nox/nox_nln_problem.H"

namespace LOCA {
  class GlobalData;
  class ParameterVector;
namespace StatusTest {
  class Abstract;
}
namespace NLN {
  class Problem : public NOX::NLN::Problem
  {
  public:
    //! constructor
    Problem(
        const Teuchos::RCP<NOX::NLN::GlobalData>& nox_nln_global_data_ptr,
        const Teuchos::RCP<LOCA::GlobalData>& loca_global_data_ptr,
        const Teuchos::RCP<NOX::Epetra::Vector>& x_ptr,
        const Teuchos::RCP<LINALG::SparseOperator>& jac_ptr,
        const Teuchos::RCP<LOCA::ParameterVector>& loca_param_vec_ptr);

    //! destructor
    virtual ~Problem() {};

    //! Create a LOCA group
    virtual Teuchos::RCP<NOX::Abstract::Group> CreateGroup(
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys) const;

    //! Create the loca-/inner-/outer-status test objects
    virtual void CreateStatusTests(Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
        Teuchos::RCP<LOCA::StatusTest::Abstract>& locaTests) const;

    virtual void CreateStatusTests(
        Teuchos::RCP<NOX::StatusTest::Generic>& outerTest,
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTest) const;

  protected:
    Teuchos::RCP<LOCA::GlobalData> loca_global_data_ptr_;

    Teuchos::RCP<LOCA::ParameterVector> loca_param_vec_ptr_;

    bool isLocaStatusTest_;
  };  // class Problem
} // namespace NLN
} // namespace LOCA



#endif /* SRC_LOCA_CONTINUATION_LOCA_NLN_PROBLEM_H_ */
