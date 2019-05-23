/*-----------------------------------------------------------*/
/*!

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_generic.H"  // class definition
#include <iomanip>

std::ostream& NOX::NLN::INNER::StatusTest::operator<<(std::ostream& os, StatusType type)
{
  os << std::setiosflags(std::ios::left) << std::setw(13) << std::setfill('.');
  switch (type)
  {
    case status_failed:
      os << "Failed";
      break;
    case status_converged:
      os << "Converged";
      break;
    case status_unevaluated:
      os << "??";
      break;
    case status_no_descent_direction:
    case status_step_too_short:
    case status_step_too_long:
    default:
      os << "**";
      break;
  }
  os << std::resetiosflags(std::ios::adjustfield) << std::setfill(' ');
  return os;
}
