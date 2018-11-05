/*----------------------------------------------------------------------------*/
/*!
\file nox_nln_floating_point_exception.cpp

\brief Bundle the floating point exception handling in one class.

\maintainer Michael Hiermeier

\date Aug 27, 2018

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "nox_nln_floating_point_exception.H"

#include <NOX_Utils.H>
#include <fenv.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::disable() const
{
  if (not shall_be_caught_) return;

  fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::clear() const
{
  if (not shall_be_caught_) return;

  feclearexcept(FE_ALL_EXCEPT);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::FloatingPointException::checkAndPrint(std::ostream& os)
{
  int exception_occurred = 0;
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

  return exception_occurred;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::FloatingPointException::enable() const
{
  if (not shall_be_caught_) return;

  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
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
    os << NOX::Utils::fill(40, '-') << "\n"
       << "Caught floating point exceptions = " << err << std::endl;

  clear();
  enable();

  return err;
}
