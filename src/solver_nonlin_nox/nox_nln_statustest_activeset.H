/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_activeset.H

\brief Check the active set for convergence. Only meaningful for
       inequality constrained problems.

\maintainer Michael Hiermeier

\date Oct 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_

#include <NOX_StatusTest_Generic.H> // base class

#include <Teuchos_RCP.hpp>
#include "nox_nln_enum_lists.H"

// forward declarations
class Epetra_Map;

namespace NOX {
namespace NLN {
namespace StatusTest {
  class ActiveSet : public NOX::StatusTest::Generic
  {
  public:
    //! constructor
    ActiveSet(const NOX::NLN::StatusTest::QuantityType& qtype,
        const int& max_cycle_size);

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

    //! maximal cycle size, which is checked
    int max_cycle_size_;

    //! cycle size if zigzagging is checked
    int cycle_size_;

    //! size of the active set
    int activesetsize_;

    std::deque<Teuchos::RCP<const Epetra_Map> > cycling_maps_;
  };  // class ActiveSet
} // namespace StatusTest
} // namespace NLN
} // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_ACTIVESET_H_ */
