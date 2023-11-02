/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of cell types traits
\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_H
#define BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_H

#include "baci_discretization_fem_general_cell_type.H"

#include <tuple>
#include <type_traits>

namespace CORE::FE
{

  namespace DETAILS
  {
    template <CORE::FE::CellType celltype>
    using is_tet = std::bool_constant<celltype == CORE::FE::CellType::tet4 ||
                                      celltype == CORE::FE::CellType::tet10>;

    template <CORE::FE::CellType celltype>
    using is_hex = std::bool_constant<
        celltype == CORE::FE::CellType::hex8 || celltype == CORE::FE::CellType::hex16 ||
        celltype == CORE::FE::CellType::hex18 || celltype == CORE::FE::CellType::hex20 ||
        celltype == CORE::FE::CellType::hex27>;

    template <CORE::FE::CellType celltype>
    using is_wedge = std::bool_constant<celltype == CORE::FE::CellType::wedge6 ||
                                        celltype == CORE::FE::CellType::wedge15>;

    template <CORE::FE::CellType celltype>
    using is_pyramid = std::bool_constant<celltype == CORE::FE::CellType::pyramid5>;

    template <CORE::FE::CellType celltype>
    using is_quad = std::bool_constant<
        celltype == CORE::FE::CellType::quad4 || celltype == CORE::FE::CellType::quad6 ||
        celltype == CORE::FE::CellType::quad8 || celltype == CORE::FE::CellType::quad9>;

    template <CORE::FE::CellType celltype>
    using is_tri = std::bool_constant<celltype == CORE::FE::CellType::tri3 ||
                                      celltype == CORE::FE::CellType::tri6>;

    template <CORE::FE::CellType celltype>
    using is_line = std::bool_constant<
        celltype == CORE::FE::CellType::line2 || celltype == CORE::FE::CellType::line3 ||
        celltype == CORE::FE::CellType::line4 || celltype == CORE::FE::CellType::line5 ||
        celltype == CORE::FE::CellType::line6>;

    template <typename... t>
    struct Join
    {
      using type = decltype(std::tuple_cat(std::declval<t>()...));
    };
  }  // namespace DETAILS

  //! @name Type traits for cell shape based on cell types
  /// @{
  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_tet = DETAILS::is_tet<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_hex = DETAILS::is_hex<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_wedge = DETAILS::is_wedge<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_pyramid = DETAILS::is_pyramid<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_quad = DETAILS::is_quad<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_tri = DETAILS::is_tri<celltype>::value;

  template <CORE::FE::CellType celltype>
  inline static constexpr bool is_line = DETAILS::is_line<celltype>::value;
  /// @}


  //! @name Compile time lists of types with cell type as template argument
  ///@{
  template <typename... Ts>
  using BaseTypeList = std::tuple<Ts...>;

  /*!
   * @brief Join two compile time base type lists
   *
   * Example:
   * @code{.cpp}
   * using BaseHex = BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>>
   * using BaseTet = BaseTypeList<Base<CellType::tet4>, Base<CellType::tet10>>
   *
   * // result in BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>, Base<CellType::tet4>,
   * //   Base<CellType::tet10>>
   * using BaseHexAndTet = Join<BaseHex, BaseTet>;
   * @endcode
   *
   * @tparam t
   */
  template <typename... t>
  using Join = typename DETAILS::Join<t...>::type;
  ///@}

  //! @name Compile time lists of cell types
  //
  // You can create a compile time list of cell types via
  //
  // \code{.cpp}
  // using CellTypes = CellTypeList<CORE::FE::CellType::hex8, CORE::FE::CellType::hex27>
  // \endcode
  //
  ///@{
  template <CORE::FE::CellType... celltypes>
  struct CellTypeList
  {
    template <template <CORE::FE::CellType> typename Base>
    using Apply = BaseTypeList<Base<celltypes>...>;
  };

  /*!
   * @brief Apply cell type list to a class that takes the cell type as template parameter and
   * return a compile time list with the base type for each cell type in the cell type list.
   *
   * Suppose you have a template class @c Base<CellType> and you want to generate a compile time
   * list for different cell types. You can do it like so:
   *
   * @code{.cpp}
   * using CellTypes = CellTypeList<CellType::hex8, CellType::hex27>
   *
   * // results in BaseTypeList<Base<CellType::hex8>, Base<CellType::hex27>>
   * using BaseTypeList = ApplyToBase<Base, CellTypes>;
   *
   * @endcode
   *
   * You could use this BaseTypeList to create a variant of a the template class @c Base for each
   * cell type
   *
   * @code{.cpp}
   * namespace DETAILS
   * {
   *   template <typename Tuple>
   *   struct CreateVariant;
   *
   *   template <typename... Ts>
   *   struct CreateVariant<BaseTypeList<Ts...>>
   *   {
   *     using type = std::variant<Ts...>;
   *   };
   * }
   *
   * template <typename... Ts>
   * using CreateVariantType = typename DETAILS::CreateVariant<Ts...>::type;
   *
   * // results in std::variant<Base<CellType::hex8>, Base<CellType::hex27>>
   * using Variant = CreateVariantType<BaseTypeList>;
   * @endcode
   *
   * @tparam Base
   * @tparam CellTypeList
   */
  template <template <CORE::FE::CellType> typename Base, typename CellTypeList>
  using ApplyToBase = typename CellTypeList::template Apply<Base>;
  ///@}

}  // namespace CORE::FE

#endif