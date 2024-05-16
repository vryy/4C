/*----------------------------------------------------------------------------*/
/*! \file

\brief Bundle the floating point exception handling in one class.



\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_solver_nonlin_nox_floating_point_exception.hpp"

#include <NOX_Utils.H>
#ifdef FOUR_C_ENABLE_FE_TRAPPING
#include <fenv.h>
#endif

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::disable() const
{
  if (not shall_be_caught_) return;
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::clear() const
{
  if (not shall_be_caught_) return;
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  feclearexcept(FE_ALL_EXCEPT);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::FloatingPointException::checkAndPrint(std::ostream& os)
{
  int exception_occurred = 0;
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  if (fetestexcept(FE_INVALID))
  {
    os << "Floating Point Exception: INVALID\n";
    ++exception_occurred;
  }
  if (fetestexcept(FE_OVERFLOW))
  {
    os << "Floating Point Exception: OVERFLOW\n";
    ++exception_occurred;
  }
  if (fetestexcept(FE_DIVBYZERO))
  {
    os << "Floating Point Exception: DIVISON BY ZERO\n";
    ++exception_occurred;
  }
#endif
  return exception_occurred;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::enable() const
{
  if (not shall_be_caught_) return;
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::precompute() const { disable(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::FloatingPointException::postcompute(std::ostream& os) const
{
  const int err = checkAndPrint(os);
  if (err)
    os << ::NOX::Utils::fill(40, '-') << "\n"
       << "Caught floating point exceptions = " << err << std::endl;

  clear();
  enable();

  return err;
}

FOUR_C_NAMESPACE_CLOSE
