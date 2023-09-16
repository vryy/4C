/*---------------------------------------------------------------------*/
/*! \file

\brief printing error messages

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_config.H"

#include "baci_utils_exceptions.H"

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <sstream>

BACI_NAMESPACE_OPEN

static int latest_line = -1;
static std::string latest_file = "{dserror_func call without prototype}";
static std::string latest_func = "{dummy latest function name}";
static int err_count = 0;

/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C" void cpp_dsassert_func(
    const char* file, const int line, const int test, const char* text)
{
#ifdef BACI_DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(const std::string file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
extern "C" void cpp_dslatest(const char* file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C" [[noreturn]] void cpp_dserror_func(const char* text, ...)
{
  int initialized;
  MPI_Initialized(&initialized);
  int myrank = 0;
  if (initialized > 0)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  }

  const std::size_t BUFLEN = 8192;
  char errbuf[BUFLEN];

  va_list ap;
  va_start(ap, text);

  snprintf(
      errbuf, BUFLEN, "PROC %d ERROR in %s, line %i:\n", myrank, latest_file.c_str(), latest_line);
  vsnprintf(&errbuf[strlen(errbuf)], BUFLEN - strlen(errbuf), text, ap);

  snprintf(&errbuf[strlen(errbuf)], BUFLEN - strlen(errbuf), "\n------------------\n");

  va_end(ap);

  throw CORE::Exception(errbuf);
} /* end of dserror_func */

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
[[noreturn]] void cpp_dserror_func(const std::string text, ...)
{
  cpp_dserror_func(text.c_str());
} /* end of dserror_func */

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void run_time_error_latest(const std::string file, const std::string func, const int line)
{
  latest_file = file;
  unsigned l = func.find("(", 0);
  latest_func = func.substr(0, l);
  latest_line = line;
}

namespace
{
  [[noreturn]] void run_time_error_func_internal(
      const std::string& errorMsg, const bool& is_catch = false)
  {
    if (not is_catch) err_count = 0;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    std::ostringstream msg_0;
    msg_0 << "[" << err_count << "] "
          << "PROC " << myrank << " ERROR in " << latest_file << ", line " << latest_line << ":\n";

    std::ostringstream msg;
    msg << msg_0.str() << "    " << latest_func << " - " << errorMsg;

    // increase level counter
    ++err_count;

    throw std::runtime_error(msg.str().c_str());
  }
}  // namespace

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
[[noreturn]] void run_time_error_func(const std::string& errorMsg, const std::runtime_error& e)
{
  std::ostringstream msg;
  msg << errorMsg << "\n" << e.what();

  run_time_error_func_internal(msg.str(), true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
[[noreturn]] void run_time_error_func(const std::string& errorMsg, const CORE::Exception& e)
{
  std::ostringstream msg;
  msg << errorMsg << "\n" << e.what();

  run_time_error_func_internal(msg.str(), true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
[[noreturn]] void run_time_error_func(const std::runtime_error& e)
{
  run_time_error_func("Caught runtime_error:", e);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
[[noreturn]] void run_time_error_func(const CORE::Exception& e)
{
  run_time_error_func("Caught runtime_error:", e);
}


const char* CORE::Exception::what() const noexcept
{
  what_message_ = message + to_string(*stacktrace);
  return what_message_.c_str();
}


CORE::Exception::Exception(std::string message)
    : message(std::move(message)), stacktrace(std::make_unique<boost::stacktrace::stacktrace>())
{
}

BACI_NAMESPACE_CLOSE
