/*---------------------------------------------------------------------*/
/*! \file

\brief Data packing for sending over MPI

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_COMM_PACK_BUFFER_HPP
#define FOUR_C_COMM_PACK_BUFFER_HPP

#include "4C_config.hpp"

#include <cstring>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer
  {
    friend class SizeMarker;

   public:
    /**
     * This class is used to mark the size of an object in the buffer. The class uses RAII to
     * automatically prepend the size of the object in the buffer when an instance goes out of
     * scope.
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

    PackBuffer() : size_(0), grow_(true) {}

    void StartPacking()
    {
      grow_ = false;
      buf_.reserve(size_);
    }

    std::vector<char>& operator()() { return buf_; }

    const std::vector<char>& operator()() const { return buf_; }

    /// add POD object
    template <typename kind>
    void add_to_pack(const kind& stuff)
    {
      std::size_t osize = sizeof(kind);
      if (grow_)
      {
        size_ += osize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize(oldsize + osize);
        std::memcpy(&buf_[oldsize], &stuff, osize);
      }
    }

    /// add array of POD objects
    template <typename kind>
    void add_to_pack(const kind* stuff, std::size_t stuffsize)
    {
      if (grow_)
      {
        size_ += stuffsize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize(oldsize + stuffsize);
        std::memcpy(&buf_[oldsize], stuff, stuffsize);
      }
    }

   private:
    /// set size of a ParObject after it has been inserted
    void set_object_size(std::size_t oldsize)
    {
      if (not grow_)
      {
        int osize = buf_.size() - oldsize;
        std::memcpy(&buf_[oldsize - sizeof(int)], &osize, sizeof(int));
      }
    }

    std::vector<char> buf_;
    std::size_t size_;
    bool grow_;
  };
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
