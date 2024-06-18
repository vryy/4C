
/*! \file

\brief An interface defining an serializable class

\level 1
*/

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"

#include <variant>
#include <vector>


#ifndef FOUR_C_SOLID_3D_ELE_INTERFACE_SERIALIZABLE_HPP
#define FOUR_C_SOLID_3D_ELE_INTERFACE_SERIALIZABLE_HPP

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  /*!
   * @brief A type trait to check whether the type T supports the
   * pack(Core::Communication::PackBuffer&) method.
   *
   * @note T does not necessarily need to derive from some interface class
   * @{
   */
  template <typename T, typename AlwaysVoid = void>
  constexpr bool IsPackable = false;

  template <typename T>
  constexpr bool IsPackable<T, std::void_t<decltype(std::declval<const T>()->pack(
                                   std::declval<Core::Communication::PackBuffer&>()))>> = true;
  ///@}

  /*!
   * @brief A type trait to check whether the type T supports the
   * unpack(std::vector<char>::size_type&, const std::vector<char>&) method.
   *
   * @note T does not necessarily need to derive from some interface class
   * @{
   */
  template <typename T, typename AlwaysVoid = void>
  constexpr bool IsUnpackable = false;

  template <typename T>
  constexpr bool IsUnpackable<T,
      std::void_t<decltype(std::declval<T>()->unpack(std::declval<std::vector<char>::size_type&>(),
          std::declval<const std::vector<char>&>()))>> = true;
  ///@}


  namespace Details
  {
    /*!
     * @brief This struct should be used to serialize an item within a sum type (e.g.
     *std::variant).
     *
     * The pack(Core::Communication::PackBuffer&) is called for those types that provide the method.
     *Nothing is done for types that don't provide a Pack-method.
     *
     * @note If you have a type that needs to be serialized during execution, you could verify that
     * the Pack-member function is correct by asserting it with static_assert(IsPackable<T>); in
     *your cpp file
     */
    struct PackAction
    {
      PackAction(Core::Communication::PackBuffer& buffer) : data(buffer) {}

      template <typename T, std::enable_if_t<IsPackable<T>, bool> = true>
      void operator()(const T& packable)
      {
        packable->pack(data);
      }

      template <typename T, std::enable_if_t<!IsPackable<T>, bool> = true>
      void operator()(const T& other)
      {
        // do nothing if it is not packable
      }

      Core::Communication::PackBuffer& data;
    };

    /*!
     * @brief This struct should be used to deserialize an item within a sum type (e.g.
     * std::variant).
     *
     * unpack(std::size_t&, const std::vector<char>&) is called for those types that provide the
     * method. Nothing is done for types that don't provide a Unpack-method.
     *
     * @note If you have a type that needs to be serialized during execution, you could verify that
     * the Pack-member function is correct by asserting it with static_assert(IsUnpackable<T>); in
     * your cpp file
     */
    struct UnpackAction
    {
      UnpackAction(std::size_t& p, const std::vector<char>& d) : position(p), data(d) {}

      template <typename T, std::enable_if_t<IsUnpackable<T&>, bool> = true>
      void operator()(T& unpackable)
      {
        unpackable->unpack(position, data);
      }

      template <typename T, std::enable_if_t<!IsUnpackable<T&>, bool> = true>
      void operator()(T& other)
      {
        // do nothing if it is not unpackable
      }

      std::size_t& position;
      const std::vector<char>& data;
    };
  }  // namespace Details

  /*!
   * @brief Pack the item within the variant if it is packable
   *
   * @tparam VariantType A template argument to match all variant types
   * @param variant (in) : variant to pack
   * @param data (in/out) : the data to pack the variant to
   */
  template <typename VariantType>
  void pack(const VariantType& variant, Core::Communication::PackBuffer& data)
  {
    std::visit(Details::PackAction(data), variant);
  }

  /*!
   * @brief Unpack the data into the variant type
   *
   * @tparam VariantType A template argument to match all variant types
   * @param variant (in/out) : variant to unpack
   * @param position (in/out) : position where to start unpacking. Will be incremented within the
   * function
   * @param data (in) : data to unpack from
   */
  template <typename VariantType>
  void unpack(
      VariantType& variant, std::vector<char>::size_type& position, const std::vector<char>& data)
  {
    std::visit(Details::UnpackAction(position, data), variant);
  }
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif