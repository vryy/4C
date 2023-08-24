
/*! \file

\brief An interface defining an serializable class

\level 1
*/

#include "baci_lib_pack_buffer.H"

#include <vector>


#ifndef BACI_SOLID_ELE_INTERFACE_SERIALIZABLE_H
#define BACI_SOLID_ELE_INTERFACE_SERIALIZABLE_H

namespace DRT::ELEMENTS
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
    virtual void Pack(DRT::PackBuffer& data) const = 0;

    /*!
     * @brief Unpack the state of an object from a char vector starting from the position
     *
     * @param position (in/out) : position where to unpack the data
     * @param data     (in)     : data to be unpacked.
     */
    virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data) = 0;
  };
}  // namespace DRT::ELEMENTS

#endif  // SOLID_ELE_INTERFACE_SERIALIZABLE_H