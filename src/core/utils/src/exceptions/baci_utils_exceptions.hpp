/*---------------------------------------------------------------------*/
/*! \file

\brief central error printing functionality

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_EXCEPTIONS_HPP
#define FOUR_C_UTILS_EXCEPTIONS_HPP

#include "baci_config.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace CORE
{
  namespace INTERNAL
  {
    class ExceptionImplementation;
  }

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
     * Return a message that describes what happened.
     */
    [[nodiscard]] const char* what() const noexcept override;

    /**
     * Return a message that describes what happened and includes a stack trace.
     *
     * @note Calling this function can be a lot more expensive than the what() function because the
     * stacktrace needs to be symbolyzed.
     */
    [[nodiscard]] std::string what_with_stacktrace() const noexcept;

   private:
    /**
     * Pointer to implementation. This technique is used to minimize the footprint of the exception
     * class that is put on the stack.
     */
    std::unique_ptr<INTERNAL::ExceptionImplementation> pimpl_;
  };

  namespace INTERNAL
  {

    [[noreturn]] void throw_error(const char* file, int line, const char* format, ...);
    [[noreturn]] void throw_error(const char* file, int line, const std::string& format, ...);

    /**
     * A helper struct taking the file name and line number from the error macro.
     */
    struct ErrorHelper
    {
      const char* file_name;
      int line_number;

      template <typename StringType, typename... Args>
      [[noreturn]] void operator()(const StringType& format, Args&&... args) const
      {
        static_assert(
            (... && std::is_trivial_v<std::decay_t<Args>>), "Can only format trivial types.");
        throw_error(file_name, line_number, format, std::forward<Args>(args)...);
      }
    };

  }  // namespace INTERNAL
}  // namespace CORE

#ifdef FOUR_C_ENABLE_ASSERTIONS

/**
 * Assert that @p test is `true`. If not issue an error in the form of a CORE::Exception.
 * This macro is only active if FOUR_C_ENABLE_ASSERTIONS is set. Otherwise, the @p test is not
 * even evaluated.
 *
 * This macro takes the test to evaluate and an error message.
 * For example:
 *
 * @code
 *   FOUR_C_ASSERT(vector.size() == dim, "Vector size does not equal dimension.");
 * @endcode
 */
#define FOUR_C_ASSERT(test, args...) \
  if (!(test))                       \
  {                                  \
    FOUR_C_THROW(args);              \
  }                                  \
  static_assert(true, "Terminate FOUR_C_ASSERT with a comma.")

#else

/**
 * This macro would asserts that @p test is true, but only if FOUR_C_ENABLE_ASSERTIONS is set.
 */
#define FOUR_C_ASSERT(test, args...) static_assert(true, "Terminate FOUR_C_ASSERT with a comma.")

#endif

/**
 * Throw an error in the form of a CORE::Exception.
 *
 * This macro takes an error message, which may contain C-style formatting.
 * All format arguments are passed as additional arguments. For example:
 *
 * @code
 *   FOUR_C_THROW("An error occured in iteration %d.", iter);
 * @endcode
 */
#define FOUR_C_THROW \
  FourC::CORE::INTERNAL::ErrorHelper { __FILE__, __LINE__ }

/**
 * Throw an error in the form of a CORE::Exception, unless the @p test is true.
 *
 * This macro takes an error message, which may contain C-style formatting.
 * All format arguments are passed as additional arguments. For example:
 *
 * @code
 *   FOUR_C_THROW_UNLESS(vector.size() == dim, "Vector size does not equal dimension d=%d.", dim);
 * @endcode
 */
#define FOUR_C_THROW_UNLESS(test, args...) \
  if (!(test))                             \
  {                                        \
    FOUR_C_THROW(args);                    \
  }                                        \
  static_assert(true, "Terminate FOUR_C_THROW_UNLESS with a comma.")

FOUR_C_NAMESPACE_CLOSE

#endif
