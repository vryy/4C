/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_generic.H

\maintainer Michael Hiermeier

\date Aug 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INNER_STATUSTEST_GENERIC_H_
#define NOX_NLN_INNER_STATUSTEST_GENERIC_H_

#include <NOX_StatusTest_Generic.H>

// forward declaration ...
namespace NOX {
namespace Solver {
  class Generic;
} // namespace Solver
namespace Abstract {
  class Group;
} // namespace Abstract
namespace NLN {
namespace INNER {
namespace StatusTest {
namespace Interface {
  class Required;
} // namespace Interface

  enum StatusType
  {
    //! Unevaluated
    status_unevaluated = -4,
    //! Failed
    status_failed = -3,
    // No descent direction
    status_no_descent_direction=-2,
    //! Step too short
    status_step_too_short = -1,
    //! Step too long
    status_step_too_long = 0,
    //! Converged
    status_converged = 1
  };


  class Generic
  {
  public:
    //! constructor
    Generic() {};

    //! destructor
    virtual ~Generic() {};

    /** \brief %Test the inner stopping criterion
     *
     *  The test can (and should, if possible) be skipped if
     *  checkType is NOX::StatusType::None.  If the test is skipped, then
     *  the status should be set to NOX::StatusTest::Unevaluated. */
    virtual StatusType CheckStatus(
        const Interface::Required& interface,
        const NOX::Solver::Generic& solver,
        const NOX::Abstract::Group& grp,
        NOX::StatusTest::CheckType checkType) = 0;

    //! Return the result of the most recent inner checkStatus call
    virtual StatusType GetStatus() const = 0;

    //! Output formatted description of inner stopping test to output stream.
    virtual std::ostream& Print(std::ostream& stream, int indent = 0) const = 0;
  };

  // non-member function
  std::ostream& operator<<(std::ostream& os, StatusType type);
} // StatusTest
} // INNER
} // NLN
} // NOX


#endif /* NOX_NLN_INNER_STATUSTEST_GENERIC_H_ */
