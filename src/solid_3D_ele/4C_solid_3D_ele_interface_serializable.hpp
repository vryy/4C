// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_pack_helpers.hpp"

#include <type_traits>
#include <variant>


#ifndef FOUR_C_SOLID_3D_ELE_INTERFACE_SERIALIZABLE_HPP
#define FOUR_C_SOLID_3D_ELE_INTERFACE_SERIALIZABLE_HPP

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{

  namespace Internal
  {
    template <typename T>
    using VariantItemInternalType = decltype(*std::declval<T>());

    /*!
     * @brief This struct should be used to serialize an item within a sum type (e.g.
     * std::variant).
     *
     * The pack(Core::Communication::PackBuffer&) is called for those types that provide the method.
     * Nothing is done for types that don't provide a Pack-method.
     *
     * @note If you have a type that needs to be serialized during execution, you could verify that
     * the Pack-member function is correct by asserting it with
     * static_assert(Core::Communication::Packable<T>); in your cpp file
     */
    struct PackAction
    {
      PackAction(Core::Communication::PackBuffer& buffer) : data(buffer) {}

      template <typename T>
        requires(Core::Communication::Packable<const VariantItemInternalType<T>>)
      void operator()(const T& packable)
      {
        packable->pack(data);
      }

      template <typename T>
        requires(!Core::Communication::Packable<const VariantItemInternalType<T>>)
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
     * the Pack-member function is correct by asserting it with
     * static_assert(Core::Communication::is_unpackable<T>); in your cpp file
     */
    struct UnpackAction
    {
      UnpackAction(Core::Communication::UnpackBuffer& buffer) : buffer_(buffer) {}

      template <typename T>
        requires(Core::Communication::Unpackable<VariantItemInternalType<T>>)
      void operator()(T& unpackable)
      {
        unpackable->unpack(buffer_);
      }

      template <typename T>
        requires(!Core::Communication::Unpackable<VariantItemInternalType<T>>)
      void operator()(T& other)
      {
        // do nothing if it is not unpackable
      }

      Core::Communication::UnpackBuffer& buffer_;
    };
  }  // namespace Internal

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
    std::visit(Internal::PackAction(data), variant);
  }

  /*!
   * @brief Unpack the data into the variant type
   *
   * @tparam VariantType A template argument to match all variant types
   * @param variant (in/out) : variant to unpack
   * @param buffer (in/out) : the data to unpack the variant from
   */
  template <typename VariantType>
  void unpack(VariantType& variant, Core::Communication::UnpackBuffer& buffer)
  {
    std::visit(Internal::UnpackAction(buffer), variant);
  }
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif