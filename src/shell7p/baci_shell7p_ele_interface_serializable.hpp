/*! \file

\brief An interface defining an serializable class

\level 1
*/
#ifndef FOUR_C_SHELL7P_ELE_INTERFACE_SERIALIZABLE_HPP
#define FOUR_C_SHELL7P_ELE_INTERFACE_SERIALIZABLE_HPP

#include "baci_config.hpp"

#include "baci_comm_pack_buffer.hpp"

#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS::SHELL
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
    virtual void Pack(CORE::COMM::PackBuffer& data) const = 0;

    /*!
     * @brief Unpack the state of an object from a char vector starting from the position
     *
     * @param position (in/out) : position where to unpack the data
     * @param data     (in)     : data to be unpacked.
     */
    virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data) = 0;
  };
}  // namespace DRT::ELEMENTS::SHELL

BACI_NAMESPACE_CLOSE

#endif  // BACI_SHELL7P_ELE_INTERFACE_SERIALIZABLE_H
