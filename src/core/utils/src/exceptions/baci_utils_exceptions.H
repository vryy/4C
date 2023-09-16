/*---------------------------------------------------------------------*/
/*! \file

\brief central error printing functionality

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_UTILS_EXCEPTIONS_H
#define BACI_UTILS_EXCEPTIONS_H

#ifdef __cplusplus

#include "baci_config.H"

#include <boost/stacktrace.hpp>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>

BACI_NAMESPACE_OPEN

namespace CORE
{
  /**
   * @brief Base class for all BACI exceptions.
   *
   * Any exceptions generated directly by BACI will have this or a derived type. This allows to
   * catch BACI exceptions specifically by using this type in the `catch` clause.
   */
  class Exception : public std::exception
  {
   public:
    /**
     * Generate an Exception with the given message. A stacktrace is automatically attached to this
     * Exception.
     */
    explicit Exception(std::string message);

    /**
     * Return a message that describes what happened and includes a stack trace.
     */
    [[nodiscard]] const char* what() const noexcept override;

   private:
    /**
     * A user-defined message that explains what happened.
     */
    std::string message;

    /**
     * The generated stacktrace which is used to construct a nice error message when calling what().
     */
    std::unique_ptr<boost::stacktrace::stacktrace> stacktrace;

    /**
     * The full message that is returned by what. This message is composed of all the other
     * information stored in this class.
     *
     * @note This needs to be stored here since we only return a `const char*` from what().
     */
    mutable std::string what_message_;
  };

}  // namespace CORE

BACI_NAMESPACE_CLOSE


extern "C" [[noreturn]] void cpp_dserror_func(
    const char* text, ...);  // avoid constructing string from char*
[[noreturn]] void cpp_dserror_func(const std::string text, ...);
void cpp_dslatest(const std::string file, const int line);

inline void cpp_dsassert_func(const char* file, const int line, const bool test, const char* string)
{
  // inline to avoid function call in release mode
#ifdef BACI_DEBUG
  if (!test)
  {
    cpp_dslatest(file, line);
    cpp_dserror_func(string);
  }
#endif
} /* end of dsassert_func */

/** @name simplified version of the dserror implementation
 *
 *  \author hiermeier, \date 11/16 */
/// @{

/// gets the file name, the pretty function name and the line of the calling function
void run_time_error_latest(const std::string file, const std::string func, const int line);


/** throws a run time error during a catch block and forwards the existing
 *  extended error message */
[[noreturn]] void run_time_error_func(const std::string& errorMsg, const std::runtime_error& e);
[[noreturn]] void run_time_error_func(const std::string& errorMsg, const BACI::CORE::Exception& e);

/// throws a run time error during a catch block and forwards the existing error message
[[noreturn]] void run_time_error_func(const std::runtime_error& e);
[[noreturn]] void run_time_error_func(const BACI::CORE::Exception& e);

/// actual called macro definition
#define run_time_error \
  run_time_error_latest(__FILE__, __PRETTY_FUNCTION__, __LINE__), run_time_error_func

/// @}

// do not emit any code in release mode for dsassert (i.e., 'test' is not evaluated and optimized
// away)
#ifdef BACI_DEBUG
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#else
#define dsassert(test, string)
#endif

#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#else

#include "baci_config.H"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void cpp_dsassert_func(const char* file, const int line, const int test, const char* string);
void cpp_dserror_func(const char* string, ...);
void cpp_dslatest(const char* file, const int line);
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#endif


#endif
