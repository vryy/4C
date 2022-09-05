/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FUNCTION_MANAGER_H
#define DRT_FUNCTION_MANAGER_H


#include <iostream>
#include <vector>
#include "drt_function.H"
#include "drt_dserror.H"
#include "demangle.H"
#include <any>


namespace DRT
{

  namespace INPUT
  {
    class Lines;
    class DatFileReader;

  }  // namespace INPUT

  namespace UTILS
  {
    /// class that manages all defined spatial functions
    class FunctionManager
    {
     public:
      /// Return all known input lines that define a DRT::Function.
      static Teuchos::RCP<DRT::INPUT::Lines> ValidFunctionLines();

      /// read the baci input file and set up all functions
      void ReadInput(DRT::INPUT::DatFileReader& reader);

      /**
       * Get a Function by its @p id in the input file. In addition, you need to specify the type of
       * function interface that this function belongs to as the template argument. Note that, for
       * reasons lost to history, the @p id is a zero-based index and offset by -1 compared to the
       * number in the input file.
       *
       * @note This call performs potentially throwing, expensive casts. It is therefore recommended
       * to query required Function objects as early as possible.
       */
      template <typename T>
      T& FunctionById(int num);

     private:
      /// Internal storage for all functions. We use type erasure via std::any to store various
      /// Function objects belonging to distinct interfaces.
      std::vector<std::any> functions_;
    };

    /// add valid function lines
    void AddValidFunctionFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);
  }  // namespace UTILS

}  // namespace DRT

/// Print the possible function definition lines to `std::cout`
void PrintFunctionDatHeader();


// --- template and inline functions --- //

template <typename T>
T& DRT::UTILS::FunctionManager::FunctionById(int num)
{
  const int input_id = num + 1;
  if (functions_.size() < (unsigned int)(input_id) || input_id < 1)
    dserror("Function with index %d (i.e. input FUNCT%d) not available.", num, input_id);

  auto& function_any = functions_[num];
  dsassert(function_any.has_value(), "Implementation error.");

  T* const queried_function = std::invoke(
      [&]() -> T*
      {
        try
        {
          const auto function_rcp = std::any_cast<Teuchos::RCP<T>>(function_any);
          dsassert(function_rcp != Teuchos::null, "Implementation error.");
          return function_rcp.get();
        }
        catch (const std::bad_any_cast& exc)
        {
          return nullptr;
        }
      });

  if (!queried_function)
  {
    const std::string actual_type_name = std::invoke(
        [&function_any]()
        {
          const std::string actual_type_name_with_rcp_prefix =
              TryDemangle(function_any.type().name());

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

    dserror(
        "You tried to query function %d as a function of type '%s'.\n"
        "Actually, it has type '%s'.",
        input_id, TryDemangle(typeid(T).name()).c_str(), actual_type_name.c_str());
  }

  return *queried_function;
}

#endif