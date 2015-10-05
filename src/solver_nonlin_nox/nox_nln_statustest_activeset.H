/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_activeset.H

\maintainer Michael Hiermeier

\date Oct 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_

#include <NOX_StatusTest_Generic.H> // base class
#include "nox_nln_statustest_factory.H"

namespace NOX {
namespace NLN {
namespace StatusTest {
  class ActiveSet : public NOX::StatusTest::Generic
  {
  public:
    //! constructor
    ActiveSet(const NOX::NLN::StatusTest::QuantityType& qtype);

    //! destructor
    virtual ~ActiveSet() {};

    //! derived
    virtual NOX::StatusTest::StatusType checkStatus(
        const NOX::Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType);

    //! derived
    //! NOTE: returns the global status of all normF tests
    virtual NOX::StatusTest::StatusType getStatus() const;

    //! derived
    virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  private:
    //! current quantity type
    NOX::NLN::StatusTest::QuantityType qtype_;

    //! status of the active set test
    NOX::StatusTest::StatusType status_;

    //! cycle size if zigzagging is checked
    int cyclesize_;

    //! size of the active set
    int activesetsize_;
  };  // class ActiveSet
} // namespace StatusTest
} // namespace NLN
} // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_ */
