/*----------------------------------------------------------------------*/
/*! \file
\brief Internal classes to read lines from files
\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_IO_LINE_PARSER_HPP
#define FOUR_C_IO_LINE_PARSER_HPP

#include "baci_config.hpp"

#include "baci_utils_demangle.hpp"
#include "baci_utils_exceptions.hpp"

#include <istream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace IO
{
  /**
   * A helper to parse lines as defined in the dat file format into C++ data. This
   * is a low-level class intended for use inside more user-friendly input mechanisms.
   * The main reason why this is a class instead of a collection of functions, is to attach some
   * context in the constructor. For instance, you can pass the section name for better error
   * messages. This parser does not store any internal state related to the parsing because the
   * lines we want to parse follow a very simple grammar without recursion or nesting.
   */
  class LineParser
  {
   public:
    /**
     * Set up the Parser and give an optional additional scope message. This information is
     * prepended to all error messages. Example:
     *
     * @code
     *   LineParser parser("While reading section MY PARAMETERS: ");
     * @endcode
     */
    LineParser(std::string user_scope_message) : user_scope_(std::move(user_scope_message)) {}

    //! Read the next string in @p in and ensure it matches the expectation.
    void Consume(std::istream& in, const std::string& expected)
    {
      std::string read_string;
      in >> read_string;
      if (read_string != expected)
        FOUR_C_THROW(
            "%sCould not read expected string '%s'.", user_scope_.c_str(), expected.c_str());
    }

    //! Read a single value of given type.
    template <typename T>
    T Read(std::istream& in)
    {
      T read_object;
      in >> read_object;
      if (in.fail())
        FOUR_C_THROW("%sCould not read expected value of type '%s'.", user_scope_.c_str(),
            CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      return read_object;
    }

   private:
    //! Prepend a user message for better error messages.
    std::string user_scope_{};
  };
}  // namespace IO

FOUR_C_NAMESPACE_CLOSE

#endif
