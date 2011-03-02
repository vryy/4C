#ifndef DRT_PACK_BUFFER_H
#define DRT_PACK_BUFFER_H

#include <vector>

namespace DRT
{

  class PackBuffer
  {
    friend class ParObject;

  public:

    PackBuffer() : size_( 0 ), grow_( true ) {}

    void StartPacking()
    {
      grow_ = false;
      buf_.reserve( size_ );
    }

    std::vector<char> & operator()() { return buf_; }

    const std::vector<char> & operator()() const { return buf_; }

    /// add POD object
    template <typename kind>
    void AddtoPack( const kind& stuff )
    {
      std::size_t osize = sizeof(kind);
      if ( grow_ )
      {
        size_ += osize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize( oldsize+osize );
        memcpy( &buf_[oldsize], &stuff, osize );
      }
    }

    /// add array of POD objects
    template <typename kind>
    void AddtoPack( const kind* stuff, std::size_t stuffsize )
    {
      if ( grow_ )
      {
        size_ += stuffsize;
      }
      else
      {
        std::size_t oldsize = buf_.size();
        buf_.resize( oldsize+stuffsize );
        memcpy( &buf_[oldsize], stuff, stuffsize );
      }
    }

  private:

    /// set size of a ParObject after it has been inserted
    void SetObjectSize( std::size_t oldsize )
    {
      if ( not grow_ )
      {
        std::size_t osize = buf_.size() - oldsize;
        memcpy( &buf_[oldsize], &osize, sizeof( std::size_t ) );
      }
    }

    std::vector<char> buf_;
    std::size_t size_;
    bool grow_;
  };
}

#endif
