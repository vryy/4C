/*!--------------------------------------------------------------------*
\file dserror_stub.cpp

\brief a dserror replacement, that does not crash, but throws an exception

\level 1

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "src/drt_lib/drt_dserror.H"

#include <stdexcept>
#include <string.h>
#include <sstream>

extern "C" void cpp_dsassert_func(
    const char* file, const int line, const int test, const char* text)
{
}

void cpp_dslatest(const std::string file, const int line) {}

extern "C" void cpp_dslatest(const char* file, const int line) {}

extern "C" void cpp_dserror_func(const char* text, ...) { throw std::runtime_error(text); }

void cpp_dserror_func(const std::string text, ...) { cpp_dserror_func(text.c_str()); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void run_time_error_latest(const std::string file, const std::string func, const int line) {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void run_time_error_func(const std::string& errorMsg, const bool& is_catch)
{
  throw std::runtime_error(errorMsg.c_str());
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void run_time_error_func(const std::string& errorMsg, const std::runtime_error& e)
{
  std::ostringstream msg;
  msg << errorMsg << "\n" << e.what();

  run_time_error_func(msg.str(), true);
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void run_time_error_func(const std::runtime_error& e)
{
  run_time_error_func("Caught runtime_error:", e);
}
