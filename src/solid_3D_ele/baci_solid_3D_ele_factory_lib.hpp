/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_FACTORY_LIB_HPP
#define FOUR_C_SOLID_3D_ELE_FACTORY_LIB_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"

#include <memory>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    /*!
     * @brief A simple container type that holds a pointer to the value of type T providing value
     * semantics.
     *
     * It is used to store types of different size in a std::variant.
     *
     * @tparam T
     */
    template <typename T>
    class VariantItem
    {
     public:
      VariantItem() { value_ = std::make_unique<T>(); }
      VariantItem(T&& obj) : value_(new T(std::move(obj))) {}
      VariantItem(const T& obj) : value_(new T(obj)) {}

      ~VariantItem() = default;

      VariantItem(const VariantItem& other) : VariantItem(*other.value_) {}
      VariantItem& operator=(const VariantItem& other)
      {
        *value_ = *other.value_;
        return *this;
      }

      // allow access via dereferencing operator propagating constness
      T& operator*() { return *value_; }
      const T& operator*() const { return *value_; }

      // allow class member access propagating constness
      T* operator->() { return value_.get(); }
      const T* operator->() const { return value_.get(); }

     private:
      std::unique_ptr<T> value_;
    };

    template <typename Tuple>
    struct CreateVariant;

    template <typename... Ts>
    struct CreateVariant<CORE::FE::BaseTypeList<Ts...>>
    {
      using type = std::variant<VariantItem<Ts>...>;
    };

  }  // namespace DETAILS

  /*!
   * @brief Meta function to create a std::variant of all provides types wrapped in a @p
   * DETAILS::VariantItem.
   */
  template <typename... Ts>
  using CreateVariantType = typename DETAILS::CreateVariant<Ts...>::type;

}  // namespace DRT::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
