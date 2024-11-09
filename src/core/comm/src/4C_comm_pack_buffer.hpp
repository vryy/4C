// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_PACK_BUFFER_HPP
#define FOUR_C_COMM_PACK_BUFFER_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <cstring>
#include <typeinfo>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  /**
   * @brief A class to pack data into a buffer.
   *
   * This class allows to add various types of data to a buffer, which is essentially a vector of
   * characters. The buffer can be used to send data over MPI or serialize the data to disk.
   *
   * @note It is important that data added to the PackBuffer is extracted in the same order from
   * a corresponding UnpackBuffer.
   */
  class PackBuffer
  {
   public:
    PackBuffer() = default;

    std::vector<char>& operator()() { return buf_; }

    const std::vector<char>& operator()() const { return buf_; }

    /// Add a trivially copyable object, i.e., an object of a type that can be copied with memcpy.
    template <typename T, typename = std::enable_if_t<std::is_trivially_copyable_v<T>>>
    void add_to_pack(const T& stuff)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Write the type into the buffer
      const std::size_t hash = typeid(T).hash_code();
      scratch_buffer_.resize(sizeof(hash));
      std::memcpy(scratch_buffer_.data(), &hash, sizeof(hash));
      buf_.insert(buf_.end(), scratch_buffer_.begin(), scratch_buffer_.end());
#endif

      // Convert stuff into a vector of chars via a separate buffer.
      scratch_buffer_.resize(sizeof(T));
      std::memcpy(scratch_buffer_.data(), &stuff, sizeof(T));

      // Append the data to the actual buffer using insert() to get amortized constant complexity.
      buf_.insert(buf_.end(), scratch_buffer_.begin(), scratch_buffer_.end());
    }

    /// Add an array of trivially copyable objects, i.e., objects of a type that can be copied with
    /// memcpy.
    template <typename T, typename = std::enable_if_t<std::is_trivially_copyable_v<T>>>
    void add_to_pack(const T* stuff, std::size_t stuff_size)
    {
      FOUR_C_ASSERT(stuff_size % sizeof(T) == 0, "Size of stuff must be a multiple of sizeof(T).");

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Write the type into the buffer
      const std::size_t hash = typeid(T).hash_code();
      scratch_buffer_.resize(sizeof(hash));
      std::memcpy(scratch_buffer_.data(), &hash, sizeof(hash));
      buf_.insert(buf_.end(), scratch_buffer_.begin(), scratch_buffer_.end());
#endif

      // Convert stuff into a vector of chars via a separate buffer.
      scratch_buffer_.resize(stuff_size);
      std::memcpy(scratch_buffer_.data(), stuff, stuff_size);

      // Append the data to the actual buffer using insert() to get amortized constant complexity.
      buf_.insert(buf_.end(), scratch_buffer_.begin(), scratch_buffer_.end());
    }

   private:
    //! The actual buffer containing the packed data.
    std::vector<char> buf_;

    //! Scratch buffer used during packing to avoid reallocations.
    std::vector<char> scratch_buffer_;
  };


  /**
   * The counterpart of the PackBuffer class. This class is used to unpack data from a buffer that
   * was packed using the PackBuffer class. Internally tracks the position in the buffer every time
   * an object is extracted.
   */
  class UnpackBuffer
  {
   public:
    explicit UnpackBuffer(const std::vector<char>& data, std::size_t start_position = 0)
        : data_(data), position_(start_position)
    {
    }

    /*!
     * \brief Extract stuff from the UnpackBuffer.
     *
     * This method is a template for all POD types like int, char, double, ... . It extracts the
     * next object into @p stuff and assumes that the raw bytes actually came from an object of the
     * type T.
     */
    template <typename T>
    std::enable_if_t<std::is_trivially_copyable_v<T>, void> extract_from_pack(T& stuff)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, &data_[position_], sizeof(hash));
      position_ += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(),
          "Type mismatch during unpacking. Tried to extract type %s", typeid(T).name());
#endif

      memcpy(&stuff, &data_[position_], sizeof(T));
      position_ += sizeof(T);
    }

    /*!
     * \brief Extract multiple entries into @p stuff.
     *
     * Same as the other method but extracts @p stuff_size entries into the array @p stuff.
     */
    template <typename T>
    void extract_from_pack(T* stuff, const std::size_t stuff_size)
    {
      FOUR_C_ASSERT(stuff_size % sizeof(T) == 0, "Size of stuff must be a multiple of sizeof(T).");

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, &data_[position_], sizeof(hash));
      position_ += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(),
          "Type mismatch during unpacking. Tried to extract type %s", typeid(T).name());
#endif


      memcpy(stuff, &data_[position_], stuff_size);
      position_ += stuff_size;
    }

    /**
     * Get @p stuff but also leave it in the buffer for the next extraction.
     */
    template <typename T>
    std::enable_if_t<std::is_trivially_copyable_v<T>, void> peek(T& stuff) const
    {
      std::size_t position = position_;

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // Check that the type matches the type that was packed.
      std::size_t hash;
      std::memcpy(&hash, &data_[position_], sizeof(hash));
      position += sizeof(hash);
      FOUR_C_ASSERT(hash == typeid(T).hash_code(), "Type mismatch during unpacking.");
#endif

      memcpy(&stuff, &data_[position], sizeof(T));
    }

    /**
     * Returns true if this buffer contains no more data to be extracted.
     */
    [[nodiscard]] bool at_end() const { return position_ >= data_.size(); }

   private:
    const std::vector<char>& data_;
    std::vector<char>::size_type position_{0};
  };

}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
