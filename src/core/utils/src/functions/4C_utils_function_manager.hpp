/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_FUNCTION_MANAGER_HPP
#define FOUR_C_UTILS_FUNCTION_MANAGER_HPP


#include "4C_config.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"

#include <any>
#include <functional>
#include <typeindex>
#include <vector>

namespace Teuchos
{
  template <class T>
  class RCP;
}

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class DatFileReader;
}  // namespace Core::IO


namespace Core::UTILS
{
  /**
   * A class that collects various (mathematical) 4C Functions specified by the user.
   *
   * An instance of this class is available after reading the input file. The specific ways how to
   * parse Functions from input data are attached within the various modules.
   */
  class FunctionManager
  {
   public:
    /**
     * Type used to pass functions that create 4C Functions from a number of parsed lines.
     */
    using FunctionFactory = std::function<std::any(const std::vector<Input::LineDefinition>&)>;

    /// Return all known input lines that define a Function.
    std::vector<Input::LineDefinition> valid_function_lines();

    /// Read the 4C input file and set up all Functions.
    void ReadInput(Core::IO::DatFileReader& reader);

    /**
     * Tell the FunctionManager how to parse a set of @p possible_lines into a Function object.
     * When the ReadInput() function is called, this class will try to read lines encountered in
     * the input file according to the @p possible_lines that were previously passed in this
     * function. If one or more of these lines match, they are parsed and the @p function_factory
     * that was passed alongside these lines is called with the parsed data. @p function_factory
     * needs to return the actual function object wrapped in a Teuchos::RCP.
     */
    void add_function_definition(
        std::vector<Input::LineDefinition> possible_lines, FunctionFactory function_factory);

    /**
     * Get a Function by its @p id in the input file. In addition, you need to specify the type of
     * function interface that this function belongs to as the template argument. Note that, for
     * reasons lost to history, the @p id is a zero-based index and offset by -1 compared to the
     * number in the input file.
     *
     * @note This call performs potentially throwing, expensive casts. It is therefore recommended
     * to query required Function objects as early as possible and store them if they are required
     * repeatedly.
     */
    template <typename T>
    const T& FunctionById(int num) const;

    template <typename T>
    void SetFunctions(const std::vector<T>& functions)
    {
      functions_ = functions;
    };

   private:
    /// Internal storage for all functions. We use type erasure via std::any to store various
    /// Function objects belonging to distinct interfaces.
    std::vector<std::any> functions_;


    /**
     * Store the lines we can read and how to convert them into a 4C Function.
     */
    std::vector<std::pair<std::vector<Input::LineDefinition>, FunctionFactory>>
        attached_function_data_;
  };

  /**
   * Add valid built-in functions that are always available.
   * This includes space- and time-dependent functions and a generic function
   * depending on arbitrary variables.
   */
  void AddValidBuiltinFunctions(FunctionManager& function_manager);
}  // namespace Core::UTILS


// --- template and inline functions --- //

template <typename T>
const T& Core::UTILS::FunctionManager::FunctionById(int num) const
{
  const int input_id = num + 1;
  if (functions_.size() < (unsigned int)(input_id) || input_id < 1)
    FOUR_C_THROW("Function with index %d (i.e. input FUNCT%d) not available.", num, input_id);

  const auto& function_any = functions_[num];
  FOUR_C_ASSERT(function_any.has_value(), "Implementation error.");

  using StoredType = Teuchos::RCP<T>;
  if (typeid(StoredType) == function_any.type())
  {
    const auto function_rcp = std::any_cast<StoredType>(function_any);
    FOUR_C_ASSERT(function_rcp != Teuchos::null, "Implementation error.");
    return *function_rcp;
  }
  else
  {
    const std::string actual_type_name = std::invoke(
        [&function_any]()
        {
          const std::string actual_type_name_with_rcp_prefix =
              Core::UTILS::TryDemangle(function_any.type().name());

          // find the outermost pair of angle brackets which should enclose the type inside an RCP
          const std::size_t start = actual_type_name_with_rcp_prefix.find_first_of('<');
          const std::size_t end = actual_type_name_with_rcp_prefix.find_last_of('>');

          // if we actually have an RCP, return the contained type name, otherwise fall back to the
          // full type name
          if (std::string_view(actual_type_name_with_rcp_prefix.c_str(), start) == "Teuchos::RCP")
            return std::string(actual_type_name_with_rcp_prefix, start + 1, end - start - 1);
          else
            return actual_type_name_with_rcp_prefix;
        });

    FOUR_C_THROW(
        "You tried to query function %d as a function of type '%s'.\n"
        "Actually, it has type '%s'.",
        input_id, Core::UTILS::TryDemangle(typeid(T).name()).c_str(), actual_type_name.c_str());
  }
}

FOUR_C_NAMESPACE_CLOSE

#endif
