#ifndef FOUR_C_COMM_PACK_BUFFER_HPP
#define FOUR_C_COMM_PACK_BUFFER_HPP

#include "4C_config.hpp"

#include <cstring>
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
    friend class SizeMarker;

   public:
    /**
     * This class is used to mark the size of an object in the buffer. The class uses RAII to
     * automatically prepend the size of the object in the buffer when an instance goes out of
     * scope. An instance of this class should be created when entering a `pack` method of a class.
     */
    class SizeMarker
    {
     public:
      SizeMarker(PackBuffer& data) : data_(data)
      {
        // add dummy object size, will be filled later
        data_.add_to_pack(-1);
        old_size_ = data_().size();
      }

      ~SizeMarker() { data_.set_object_size(old_size_); }

     private:
      PackBuffer& data_;
      std::size_t old_size_{};
    };

    PackBuffer() = default;

    std::vector<char>& operator()() { return buf_; }

    const std::vector<char>& operator()() const { return buf_; }

    /// Add a trivially copyable object, i.e., an object of a type that can be copied with memcpy.
    template <typename T, typename = std::enable_if_t<std::is_trivially_copyable_v<T>>>
    void add_to_pack(const T& stuff)
    {
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
      // Convert stuff into a vector of chars via a separate buffer.
      scratch_buffer_.resize(stuff_size);
      std::memcpy(scratch_buffer_.data(), stuff, stuff_size);

      // Append the data to the actual buffer using insert() to get amortized constant complexity.
      buf_.insert(buf_.end(), scratch_buffer_.begin(), scratch_buffer_.end());
    }

   private:
    /// set size of a ParObject after it has been inserted
    void set_object_size(std::size_t oldsize)
    {
      int osize = buf_.size() - oldsize;
      std::memcpy(&buf_[oldsize - sizeof(int)], &osize, sizeof(int));
    }

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
    std::enable_if_t<std::is_pod_v<T>, void> extract_from_pack(T& stuff)
    {
      peek(stuff);
      position_ += sizeof(T);
    }

    /*!
     * \brief Extract multiple entries into @p stuff.
     *
     * Same as the other method but extracts @p stuff_size entries into the array @p stuff.
     */
    template <typename Kind>
    void extract_from_pack(Kind* stuff, const std::size_t stuff_size)
    {
      memcpy(stuff, &data_[position_], stuff_size);
      position_ += stuff_size;
    }

    /**
     * Get @p stuff but also leave it in the buffer for the next extraction.
     */
    template <typename T>
    std::enable_if_t<std::is_pod_v<T>, void> peek(T& stuff) const
    {
      memcpy(&stuff, &data_[position_], sizeof(T));
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
