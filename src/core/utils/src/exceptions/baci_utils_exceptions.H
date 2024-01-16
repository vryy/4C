/*---------------------------------------------------------------------*/
/*! \file

\brief central error printing functionality

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_UTILS_EXCEPTIONS_H
#define BACI_UTILS_EXCEPTIONS_H

#include "baci_config.H"

#include <boost/stacktrace/stacktrace_fwd.hpp>

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
     * Destructor.
     */
    ~Exception() override;

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

  namespace INTERNAL
  {
    /**
     * A helper struct taking the file name and line number from the error macro.
     */
    struct ErrorHelper
    {
      const char* file;
      int line_number;

      [[noreturn]] void operator()(const char* format, ...) const;
      [[noreturn]] void operator()(const std::string& format, ...) const;
    };

  }  // namespace INTERNAL
}  // namespace CORE

#ifdef BACI_DEBUG

/**
 * Assert that @p test is `true`. If not issue an error in the form of a CORE::Exception.
 * This macro is only active if BACI_DEBUG is set. In release mode, the @p test is not even
 * evaluated.
 *
 * This macro takes the test to evaluate and an error message.
 * For example:
 *
 * @code
 *   dsassert(vector.size() == dim, "Vector size does not equal dimension.");
 * @endcode
 */
#define dsassert(test, msg) \
  if (!(test))              \
  {                         \
    dserror(msg);           \
  }                         \
  static_assert(true, "Terminate dsassert with a comma.")

#else

/**
 * This macro would asserts that @p test is true, but only if BACI_DEBUG is set.
 */
#define dsassert(test, msg) static_assert(true, "Terminate dsassert with a comma.")

#endif

/**
 * Issue an error in the form of a CORE::Exception.
 *
 * This macro takes an error message, which may contain C-style formatting.
 * All format arguments are passed as additional arguments. For example:
 *
 * @code
 *   dserror("An error occured in iteration %d.", iter);
 * @endcode
 */
#define dserror \
  BACI::CORE::INTERNAL::ErrorHelper { __FILE__, __LINE__ }

BACI_NAMESPACE_CLOSE

#endif
