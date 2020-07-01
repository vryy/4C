/*----------------------------------------------------------------------*/
/*! \file

\brief printing error messages

\level 2
*----------------------------------------------------------------------*/
#ifndef DSERROR_H
#define DSERROR_H

#ifdef __cplusplus

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <string>

namespace std
{
  class runtime_error;
}

void cpp_dsassert_func(
    const std::string file, const int line, const bool test, const std::string string);
void cpp_dserror_func(const std::string string, ...);
void cpp_dslatest(const std::string file, const int line);

void run_time_error_latest(const std::string file, const std::string func, const int line);
void run_time_error_func(const std::string& errorMsg, const bool& is_catch = false);
void run_time_error_func(const std::string& errorMsg, const std::runtime_error& e);
void run_time_error_func(const std::runtime_error& e);

#define run_time_error \
  run_time_error_latest(__FILE__, __PRETTY_FUNCTION__, __LINE__), run_time_error_func
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#ifdef THROWELEMENTERRORS

/// notify error subsystem of beginning of element loop
void EnterElementLoop();

/// notify error subsystem of end of element loop
void ExitElementLoop();

/// notify error subsystem of element error
void ElementError(int ele, const std::string& err);

#endif  // THROWELEMENTERRORS
#else

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

void cpp_dsassert_func(const char* file, const int line, const int test, const char* string);
void cpp_dserror_func(const char* string, ...);
void cpp_dslatest(const char* file, const int line);
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#endif
#endif
