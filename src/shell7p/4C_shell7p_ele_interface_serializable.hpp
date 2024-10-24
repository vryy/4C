// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SHELL7P_ELE_INTERFACE_SERIALIZABLE_HPP
#define FOUR_C_SHELL7P_ELE_INTERFACE_SERIALIZABLE_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements::Shell
{
  /*!
   * @brief An interface providing methods to serialize an object into a packbuffer and deserialize
   * from a char array
   *
   */
  class Serializable
  {
   public:
    virtual ~Serializable() = default;
    /*!
     * @brief Pack the whole state of the object into into the pack buffer
     *
     * @param data (out) : the buffer to pack into.
     */
    virtual void pack(Core::Communication::PackBuffer& data) const = 0;

    /*!
     * @brief Unpack the state of an object from a char vector starting from the position
     *
     * @param position (in/out) : position where to unpack the data
     * @param data     (in)     : data to be unpacked.
     */
    virtual void unpack(Core::Communication::UnpackBuffer& buffer) = 0;
  };
}  // namespace Discret::Elements::Shell

FOUR_C_NAMESPACE_CLOSE

#endif
