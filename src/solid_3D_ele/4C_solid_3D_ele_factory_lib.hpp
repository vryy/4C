/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_FACTORY_LIB_HPP
#define FOUR_C_SOLID_3D_ELE_FACTORY_LIB_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_inpar_structure.hpp"

#include <memory>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Solid::ELEMENTS
{
  enum class EasType;
}
namespace Discret::ELEMENTS
{
  namespace Details
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
    struct CreateVariant<Core::FE::BaseTypeList<Ts...>>
    {
      using type = std::variant<VariantItem<Ts>...>;
    };

  }  // namespace Details

  /*!
   * @brief Meta function to create a std::variant of all provides types wrapped in a @p
   * Details::VariantItem.
   */
  template <typename... Ts>
  using CreateVariantType = typename Details::CreateVariant<Ts...>::type;

  namespace Details
  {
    /*!
     * @brief A struct that determines whether @p T is a valid template type
     *
     * The member variable value is true if the first template parameter is a valid type
     *
     * @tparam typename : Template parameter that may be a valid type or not
     */
    template <typename, typename = void>
    struct IsValidTypeTrait : std::false_type
    {
    };

    template <typename T>
    struct IsValidTypeTrait<T, std::void_t<typename T::type>> : std::true_type
    {
    };
  }  // namespace Details

  /*!
   * @brief Determines whether we have implemented a solid calculation formulation
   *
   * @tparam T typename: Template parameter that may be a valid type or not
   */
  template <typename T>
  constexpr bool is_valid_type = Details::IsValidTypeTrait<T>::value;

  /*!
   * @brief An automatic switch from runtime kinematic type to constexpr kinematic type.
   *
   * @tparam Function
   * @param kinem_type
   * @param fct : A callable function that has operators for each kinematic-type integral
   * constants
   * @return auto
   */
  template <typename Function>
  auto switch_kinematic_type(Inpar::Solid::KinemType kinem_type, Function fct)
  {
    switch (kinem_type)
    {
      case Inpar::Solid::KinemType::linear:
        return fct(
            std::integral_constant<Inpar::Solid::KinemType, Inpar::Solid::KinemType::linear>{});
      case Inpar::Solid::KinemType::nonlinearTotLag:
        return fct(std::integral_constant<Inpar::Solid::KinemType,
            Inpar::Solid::KinemType::nonlinearTotLag>{});
      case Inpar::Solid::KinemType::vague:
        return fct(
            std::integral_constant<Inpar::Solid::KinemType, Inpar::Solid::KinemType::vague>{});
    }

    FOUR_C_THROW("Your kinematic type is unknown: %d", kinem_type);
  }
}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
