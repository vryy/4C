/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of cell types traits.

In this file, the helper for defining type trait functionality is defined

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DETAILS_H
#define BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DETAILS_H

#include "baci_config.H"

#include "baci_discretization_fem_general_cell_type.H"
#include "baci_discretization_fem_general_cell_type_traits_def.H"
#include "baci_utils_exceptions.H"

#include <tuple>

BACI_NAMESPACE_OPEN

namespace CORE::FE::DETAILS
{
  //! @name helper functions for cell-type related type lists
  /// @{
  template <typename... t>
  struct Join
  {
    using type = decltype(std::tuple_cat(std::declval<t>()...));
  };


  template <template <typename...> typename TypeList, template <CellType> typename Base, typename T>
  struct apply_celltype_sequence;

  template <template <CellType> typename Base, template <typename...> typename TypeList,
      CellType... celltypes>
  struct apply_celltype_sequence<TypeList, Base, celltype_sequence<celltypes...>>
  {
    using type = TypeList<Base<celltypes>...>;
  };

  template <typename T>
  struct celltype_sequence_to_string_array;

  template <CellType... celltypes>
  struct celltype_sequence_to_string_array<celltype_sequence<celltypes...>>
  {
    static constexpr std::array value = {DETAILS::CellTypeInformation<celltypes>::name...};
  };

  template <typename T>
  struct celltype_sequence_to_array;

  template <CellType... celltypes>
  struct celltype_sequence_to_array<celltype_sequence<celltypes...>>
  {
    static constexpr std::array value = {celltypes...};
  };

  template <std::size_t first, typename T>
  struct celltype_sequence_generator;

  template <std::size_t first, std::size_t... ints>
  struct celltype_sequence_generator<first, std::integer_sequence<std::size_t, ints...>>
  {
    using type = celltype_sequence<static_cast<CellType>(ints + first)...>;
  };

  template <typename T>
  struct first_celltype_in_sequence;

  template <CellType celltype, CellType... celltypes>
  struct first_celltype_in_sequence<celltype_sequence<celltype, celltypes...>>
  {
    static constexpr CellType value = celltype;
  };

  template <CellType... celltypes>
  constexpr bool IsCellTypeInSequence(CellType celltype, celltype_sequence<celltypes...>)
  {
    return ((celltype == celltypes) || ...);
  }

  template <CellType first, CellType last>
  using make_celltype_sequence =
      typename celltype_sequence_generator<static_cast<std::size_t>(first),
          std::make_integer_sequence<std::size_t,
              static_cast<std::size_t>(last) + 1 - static_cast<std::size_t>(first)>>::type;
  /// @}

  //! @name Helper functions for the celltype switch
  /// @{
  template <CellType celltype, typename celltype_sequence, typename Function,
      typename UnsupportedCellTypeCallable>
  auto CellTypeSwitchItem(Function fct, UnsupportedCellTypeCallable unsupported_celltype_callable)
      -> std::invoke_result_t<Function,
          std::integral_constant<CellType, first_celltype_in_sequence<celltype_sequence>::value>>
  {
    if constexpr (IsCellTypeInSequence(celltype, celltype_sequence{}))
    {
      return fct(std::integral_constant<CellType, celltype>{});
    }

    constexpr bool is_convertible = std::is_convertible_v<
        std::invoke_result_t<UnsupportedCellTypeCallable,
            std::integral_constant<CellType, celltype>>,
        std::invoke_result_t<Function, std::integral_constant<CellType,
                                           first_celltype_in_sequence<celltype_sequence>::value>>>;
    constexpr bool is_void = std::is_void_v<std::invoke_result_t<UnsupportedCellTypeCallable,
        std::integral_constant<CellType, celltype>>>;
    static_assert(is_convertible || is_void,
        "The unsupported celltype callable must either return a convertible type to the return "
        "type of the regular callable or must return void (and throw an error)!");

    if constexpr (is_convertible)
    {
      return unsupported_celltype_callable(std::integral_constant<CellType, celltype>{});
    }

    unsupported_celltype_callable(std::integral_constant<CellType, celltype>{});
    dserror("Your unsupported celltype callable does not throw or return a compatible type!");
  }

  template <typename celltype_sequence>
  struct ThrowUnsupportedCellTypeError
  {
    template <typename CellTypeConstant>
    [[noreturn]] void operator()(CellTypeConstant celltype_t) const
    {
      constexpr std::array celltypes_str =
          celltype_sequence_to_string_array<celltype_sequence>::value;

      std::string celltypes_str_acc;

      for (const auto& item : celltypes_str)
      {
        if (!celltypes_str_acc.empty())
        {
          celltypes_str_acc += ", ";
        }
        celltypes_str_acc += item;
      }

      dserror(
          "The function you are calling is not implemented for the cell type %s. Supported "
          "celltypes are %s",
          DETAILS::CellTypeInformation<celltype_t()>::name, celltypes_str_acc.c_str());
    }
  };
  /// @}
}  // namespace CORE::FE::DETAILS

BACI_NAMESPACE_CLOSE

#endif