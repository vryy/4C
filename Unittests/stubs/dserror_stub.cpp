/*!--------------------------------------------------------------------*
\file dserror_stub.cpp

\brief a dserror replacement, that does not crash, but throws an exception

\level 1

\maintainer Karl-Robert Wichmann
*----------------------------------------------------------------------*/

#include "src/drt_lib/drt_dserror.H"

#include <stdexcept>
#include <string.h>

extern "C"
void cpp_dsassert_func(const char* file, const int line, const int test, const char* text)
{}

void cpp_dslatest(const std::string file, const int line)
{}

extern "C"
void cpp_dslatest(const char* file, const int line)
{}

extern "C"
void cpp_dserror_func(const char* text, ...)
{
  throw std::runtime_error(text);
}

void cpp_dserror_func(const std::string text, ...)
{
  cpp_dserror_func(text.c_str());
}
