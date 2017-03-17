/*----------------------------------------------------------------------*/
/**
\file pairedvector.H

\brief This class is meant as a drop-in replacement for std::maps, when other
       storage and access characteristics are needed.

<pre>
\maintainer Karl-Robert Wichmann
            wichmann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef PAIREDVECTOR_H_
#define PAIREDVECTOR_H_

#include <vector>
#include <utility>
#include <cstddef>
#include <stdexcept>

namespace GEN
{

/**
 * @brief A substitute for std::maps, that has different storage and access
 * characteristics.
 *
 * @tparam Key Type of key
 * @tparam T   Type of element
 *
 * This class is meant as a drop-in replacement for std::maps, when other
 * storage and access characteristics are needed.
 * The memory is allocated beforehand to eliminate the overhead of repeated
 * memory allocation. This requires the knowledge of an upper bound on the
 * number of entries and the size is not meant to be changed after
 * initialization. There are however some instances where that is inevitable.
 * The access characteristics are equivalent to those of a vector, which is
 * the container it is based on. Note especially that the elements are not
 * sorted.
 */
template<typename Key, typename T>
class pairedvector
{
 private:
  typedef std::vector<std::pair<Key,T> > pairedvector_type;

 public:
  typedef typename pairedvector_type::iterator       iterator;
  typedef typename pairedvector_type::const_iterator const_iterator;

  /**
   *  @brief  constructor creates no elements, but reserves the maximum
   *          number of entries.
   *  @param reserve The number of elements that are preallocated
   */
  pairedvector(size_t reserve):
    _m(reserve,std::pair<Key,T>()),
    _entries(0)
  {}

  /**
   *  @brief  constructor creates no elements, but reserves the maximum
   *          number of entries.
   *  @param reserve The number of elements that are preallocated
   *  @param default_key default value for the key within the pair
   *  @param default_T   default value for the data within the pair
   */
  pairedvector(size_t reserve, Key default_key, T default_T):
    _m(reserve,std::pair<Key,T>(default_key,default_T)),
    _entries(0)
  {}

  /**
   *  The dtor only erases the elements, and note that if the
   *  elements themselves are pointers, the pointed-to memory is
   *  not touched in any way.  Managing the pointer is the user's
   *  responsibility.
   */
  ~pairedvector()
  {}

  /**
   *  Returns a read/write iterator that points to the first
   *  element in the %pairedvector.  Iteration is done in ordinary
   *  element order.
   */
  iterator begin()
  {
    return _m.begin();
  }

  /**
   *  Returns a read-only (constant) iterator that points to the first pair
   *  in the %pairedvector.  Iteration is done in ordinary
   *  element order.
   */
  const_iterator begin() const
  {
    return _m.begin();
  }

  /**
   *  Returns a read/write iterator that points one past the last
   *  pair in the %pairedvector.  Iteration is done in ordinary
   *  element order.
   */
  iterator end()
  {
    return _m.begin()+_entries;
  }

  /**
   *  Returns a read-only (constant) iterator that points one past the last
   *  pair in the %pairedvector.  Iteration is done in ordinary
   *  element order.
   */
  const_iterator end() const
  {
    return _m.begin()+_entries;
  }

  /**
   *  @brief  Tries to locate an element in a %pairedvector.
   *  @param  k  Key of (key, value) %pair to be located.
   *  @return Iterator pointing to sought-after element, or end() if not
   *          found.
   *
   *  This function takes a key and tries to locate the element with which
   *  the key matches.  If successful the function returns an iterator
   *  pointing to the sought after %pair.  If unsuccessful it returns the
   *  past-the-end ( @c end() ) iterator.
   */
  iterator find(const Key k)
  {
    for (iterator it = this->begin(); it != this->end(); ++it)
    {
      if (it->first == k)
        return it;
    }
    return this->end();
  }

  /**
   *  @brief  Tries to locate an element in a %pairedvector.
   *  @param  k  Key of (key, value) %pair to be located.
   *  @return Read-only (constant) iterator pointing to sought-after
   *          element, or end() if not found.
   *
   *  This function takes a key and tries to locate the element with which
   *  the key matches.  If successful the function returns a constant
   *  iterator pointing to the sought after %pair. If unsuccessful it
   *  returns the past-the-end ( @c end() ) iterator.
   */
  const_iterator find(const Key k) const
  {
    for (const_iterator it = this->begin(); it != this->end(); ++it)
    {
      if (it->first == k)
        return it;
    }
    return this->end();
  }

  /**
   *  Erases all elements in a %pairedvector.  Note that this function only
   *  erases the elements, and that if the elements themselves are
   *  pointers, the pointed-to memory is not touched in any way.
   *  Managing the pointer is the user's responsibility.
   */
  void clear()
  {
    _entries = 0;

    // The vector must be overwritten explicitly, to avoid holding pointers
    // to smart- or reference counting pointers.
    for (typename pairedvector_type::iterator it = _m.begin(); it != _m.end(); ++it)
      *it = std::pair<Key,T>();
  }

  /**
   *  @brief  Resizes the %pairedvector to the specified number of elements.
   *  @param  new_size  Number of elements the %vector should contain.
   *  @param  x  Data with which new elements should be populated.
   *
   *  This function will %resize the %pairedvector to the specified
   *  number of elements.  If the number is smaller than the
   *  %pairedvector's current size the %pairedvector is truncated, otherwise
   *  the %pairedvector is extended and new elements are populated with
   *  given data.
   */
  void resize(size_t new_size, std::pair<Key,T> x = std::pair<Key,T>())
  {
    _m.resize(new_size,x);

    // If vector is truncated to new_size, adapt number of entries.
    if(new_size<_entries)
      _entries=new_size;
  }

  /**
   *  @brief  Subscript ( @c [] ) access to %pairedvector data.
   *  @param  k  The key for which data should be retrieved.
   *  @return A reference to the data of the (key,data) %pair.
   *
   *  Allows for easy lookup with the subscript ( @c [] )
   *  operator.  Returns data associated with the key specified in
   *  subscript.  If the key does not exist, a pair with that key
   *  is created using default values, which is then returned.
   */
  T& operator[](const Key k)
  {
    iterator it = this->find(k);
    if (it != this->end())
      return it->second;

    if (_entries >= _m.size())
      throw std::length_error("pairedvector::operator[]");

    _entries++;
    _m[_entries-1].first = k;
    return _m[_entries-1].second;
  }

  /**
   *  @brief  Access to %pairedvector data.
   *  @param  k  The key for which data should be retrieved.
   *  @return A reference to the data whose key is equivalent to @a k, if
   *          such a data is present in the %pairedvector.
   *  @throw  std::runtime_error("invalid key")  If no such data is present.
   */
  T& at(const Key k)
  {
    iterator it = this->find(k);
    if (it == this->end())
      throw std::runtime_error("pairedvector::at(): invalid key");

    return it->second;
  }

  const T& at(const Key k) const
  {
    const_iterator it = this->find(k);
    if (it == this->end())
      throw std::runtime_error("pairedvector::at(): invalid key");

    return it->second;
  }

  /** Returns the current capacity of the %pairedvector.  */
  size_t capacity() const
  {
    return _m.size();
  }

  /**  Returns the number of elements in the %pairedvector.  */
  size_t size() const
  {
    return _entries;
  }

  /**
   *  @brief  Erasing elements is not allowed for this data type!
   *  @param  k  Key of element to be erased.
   *
   */
  void erase(const Key k)
  {
    throw std::runtime_error("entries cannot be removed");
    return;
  }

  /**
   *  @brief  %pairedvector assignment operator.
   *  @param  x  A %vector of identical element and allocator types.
   *
   *  All the elements of @a x are copied, but unlike the copy constructor,
   *  the allocator object is not copied.
   */
  pairedvector_type& operator=(const pairedvector_type& x)
  {
    if (x.size() > _m.size())
      throw std::length_error("pairedvector::operator=");

    size_t i = 0;
    for (const_iterator src = x.begin(); src != x.end(); ++src)
    {
      _m[i].first  = src->first;
      _m[i++].second = src->second;
    }
    _entries = x.size();
    return *this;
  }

  /** Returns true if the %pairedvector is empty.  (Thus begin() would equal
   *  end().)
  */
  bool empty() const
  {
    return (_entries <= 0);
  }

 private:
  pairedvector() {};

  pairedvector_type _m;
  size_t _entries;

};

}


#endif /* PAIREDVECTOR_H_ */
