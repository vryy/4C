/*---------------------------------------------------------------------*/
/*! \file

\brief This class is meant as a replacement for std::maps, when other
       storage and access characteristics are needed.

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef PAIREDVECTOR_H_
#define PAIREDVECTOR_H_

#include <vector>
#include <utility>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include <map>

#include "pairedobj_insert_policy.H"

namespace GEN
{
  /// copy types
  enum CopyType : char
  {
    DeepCopy,  ///< copy the whole paired vector
    ShapeCopy  ///< copy the keys of the paired vector but clear all values.
  };

  /**
   * @brief A substitute for std::maps, that has different storage and access
   * characteristics.
   *
   * @tparam Key Type of key
   * @tparam T   Type of element
   *
   * \note This class is no longer a pure drop-in solution for std::map. Actually
   * it can still be used as one but with restricted functionality.
   *
   * The memory is allocated beforehand to eliminate the overhead of repeated
   * memory allocation. This requires the knowledge of an upper bound on the
   * number of entries and the size is not meant to be changed after
   * initialization. There are however some instances where that is inevitable.
   * The access characteristics are equivalent to those of a vector, which is
   * the container it is based on. Note especially that the elements are not
   * sorted.
   */
  template <typename Key, typename T, typename insert_policy = default_insert_policy<Key, T>>
  class pairedvector : protected insert_policy
  {
   private:
    typedef pairedvector<Key, T, insert_policy> class_type;
    typedef insert_policy base_type;
    typedef typename base_type::pairedvector_type pairedvector_type;
    typedef typename base_type::pair_type pair_type;

   public:
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;

    /**
     *  @brief  constructor creates no elements, but reserves the maximum
     *          number of entries.
     *  @param reserve The number of elements that are preallocated
     */
    pairedvector(size_t reserve)
        : _m(reserve + insert_policy::capacity_offset(), pair_type()), _entries(0)
    {
    }

    /**
     *  @brief  empty constructor creates no elements and does not reserve any
     *          number of entries. Use resize as soon as you know the necessary
     *          number of elements.
     */
    pairedvector() : _m(insert_policy::capacity_offset(), pair_type()), _entries(0) {}

    /**
     *  @brief  constructor creates no elements, but reserves the maximum
     *          number of entries.
     *  @param reserve The number of elements that are preallocated
     *  @param default_key default value for the key within the pair
     *  @param default_T   default value for the data within the pair
     */
    pairedvector(size_t reserve, Key default_key, T default_T)
        : _m(reserve + insert_policy::capacity_offset(), pair_type(default_key, default_T)),
          _entries(0)
    {
    }

    /**
     *  @brief  copy constructor
     *
     *  @param[in] source %pairedmatrix object we want to copy.
     *  @param[in] type   Apply this copy type.
     *
     *  \author hiermeier \date 05/17 */
    pairedvector(const pairedvector<Key, T>& source, enum GEN::CopyType type = DeepCopy)
        : _m(0, pair_type()), _entries(0)

    {
      clone(source);

      switch (type)
      {
        case ShapeCopy:
        {
          for (std::pair<Key, T>& entry : *this) entry.second = T();

          break;
        }
        default:
          break;
      }
    }

    /**
     *  The dtor only erases the elements, and note that if the
     *  elements themselves are pointers, the pointed-to memory is
     *  not touched in any way.  Managing the pointer is the user's
     *  responsibility.
     */
    ~pairedvector() {}

    /**
     *  Returns a read/write iterator that points to the first
     *  element in the %pairedvector.  Iteration is done in ordinary
     *  element order.
     */
    iterator begin() { return base_type::begin(_m.begin()); }

    /**
     *  Returns a read-only (constant) iterator that points to the first pair
     *  in the %pairedvector.  Iteration is done in ordinary
     *  element order.
     */
    const_iterator begin() const { return base_type::begin(_m.begin()); }

    /**
     *  Returns a read/write iterator that points one past the last
     *  pair in the %pairedvector.  Iteration is done in ordinary
     *  element order.
     */
    iterator end() { return base_type::end(_m.begin() + _entries); }

    /**
     *  Returns a read-only (constant) iterator that points one past the last
     *  pair in the %pairedvector.  Iteration is done in ordinary
     *  element order.
     */
    const_iterator end() const { return base_type::end(_m.begin() + _entries); }

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
    iterator find(const Key k) { return base_type::find(k, _m, _entries); }

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
    const_iterator find(const Key k) const { return base_type::find(k, _m, _entries); }

    /**
     * @param  x  Data with which old elements are overwritten.
     *
     *  Erases all elements in a %pairedvector.  Note that this function only
     *  erases the elements, and that if the elements themselves are
     *  pointers, the pointed-to memory is not touched in any way.
     *  Managing the pointer is the user's responsibility.
     */
    void clear(const pair_type& x = pair_type())
    {
      if (empty()) return;

      _entries = 0;

      // The vector must be overwritten explicitly, to avoid holding pointers
      // to smart- or reference counting pointers.
      for (typename pairedvector_type::iterator it = _m.begin(); it != _m.end(); ++it) *it = x;
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
    void resize(size_t new_size, pair_type x = pair_type())
    {
      // adapt sentinel value thresholds
      if (_m.size() >= insert_policy::capacity_offset())
        std::fill(_m.end() - insert_policy::capacity_offset(), _m.end(), x);

      _m.resize(new_size + insert_policy::capacity_offset(), x);

      // If vector is truncated to new_size, adapt number of entries.
      if (new_size < _entries) _entries = new_size;
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
    T& operator[](const Key k) { return base_type::get(k, _m, _entries); }

    /** @brief In the default case same behavior as %operator[] */
    T& operator()(const Key k) { return base_type::operator()(k, _m, _entries); }

    /** @brief In the default case same behavior as %operator[] with minor overhead
     *
     *  @note If the quick_insert_policy is chosen, speed-up of look-up for
     *  repetitive accesses with constant access pattern. See the policy for more
     *  information.
     *
     *  @param[in] k          The key for which data should be retrieved.
     *  @param[in] rep_count  Repetition counter. First repetition (rep_count==0)
     *                        is used to initialize the access pattern.
     *  @return A reference to the data of the (key,data) %pair.*/
    T& repetitive_access(const Key k, const int rep_count)
    {
      return base_type::repetitive_access(k, rep_count, _m, _entries);
    }

    /** @brief complete special internally used data structures.
     *
     *  This routine is without functionality for the default_insert_policy, but
     *  plays an important roll for the quick_insert_policy as well as for the
     *  insert_and_sort_policy. See the repetitive_access and operator() methods
     *  for more information.
     *
     *  @author hiermeier @date 05/17 */
    void complete() { _entries = insert_policy::complete(_m, _entries); }

    /**
     *  @brief  Access to %pairedvector data.
     *  @param  k  The key for which data should be retrieved.
     *  @return A reference to the data whose key is equivalent to @a k, if
     *          such a data is present in the %pairedvector.
     *  @throw  std::runtime_error("invalid key")  If no such data is present.
     */
    T& at(const Key k) { return base_type::at(k, _m, _entries); }

    /**
     *  @brief  Access to %pairedvector data.
     *  @param  k  The key for which data should be retrieved.
     *  @return A reference to the data whose key is equivalent to @a k, if
     *          such a data is present in the %pairedvector.
     *  @throw  std::runtime_error("invalid key")  If no such data is present.
     */
    const T& at(const Key k) const { return base_type::at(k, _m, _entries); }

    /** Returns the current capacity of the %pairedvector.  */
    size_t capacity() const
    {
      return (_m.size() > 0 ? _m.size() - insert_policy::capacity_offset() : 0);
    }

    /**  Returns the number of elements in the %pairedvector.  */
    size_t size() const { return _entries; }

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

    /** Returns true if the %pairedvector is empty.  (Thus begin() would equal
     *  end().)
     */
    bool empty() const { return (_entries <= 0); }

    /**
     *  @brief  Swaps data with another %pairedvector.
     *  @param  x  A %pairedvector of identical allocator type.
     *
     *  This exchanges the elements between two vectors in constant time.
     *  (Three pointers and the entries information, so it should be quite fast.)
     *
     *  \author hiermeier \date 06/17 */
    void swap(class_type& x)
    {
      // swap internal data structure
      _m.swap(x._m);
      base_type::swap(x);

      // swap entry number
      const size_t my_entries = _entries;
      _entries = x._entries;
      x._entries = my_entries;
    }

    /**
     *  @brief  assign a source %pairedvector to this %pairedvector
     *  @param  a %pairedvector of identical allocator type.
     *
     *  If necessary the capacity of this %pairedvector will be modified. All
     *  previous elements are lost and will be overwritten by the source values.
     */
    class_type& operator=(const class_type& source)
    {
      clone(source, DeepCopy);
      return *this;
    }

    /**
     *  @brief  assign a source %std::map to this %pairedvector
     *  @param  a %std_::map of identical allocator type.
     *
     *  If necessary the capacity of this %pairedvector will be modified. All
     *  previous elements are lost and will be overwritten by the source values.
     */
    class_type& operator=(const std::map<Key, T>& source)
    {
      clear();
      clone(source);

      return *this;
    }

    /** @brief clone the source object
     *
     *  If necessary the capacity of this %pairedvector will be modified. All
     *  previous elements are lost and will be overwritten by the source values.
     *
     *  @param[in] source  Clone the given object.
     *  @param[in] type    type for the clone procedure. If ShapeCopy is chosen,
     *                     only the key but not the values are copied.
     *
     *  @author hiermeier @date 07/17 */
    void clone(const class_type& source, const enum CopyType type)
    {
      clear();
      clone(source);

      switch (type)
      {
        case ShapeCopy:
        {
          for (pair_type& entry : *this) entry.second = T();

          break;
        }
        default:
          break;
      }
    }

    /**
     *  @brief  print the %pairedvector
     *  @param  output stream.
     *  @param  sort the entries before the actual print is performed.
     *
     *  Print the %pairedvector information in column format. By default the
     *  entries are sorted with respect to their KEY entries.
     *
     *  @author hiermeier @date 07/17 */
    void print(std::ostream& os, bool sort = true) const
    {
      pairedvector_type sorted_m(_m.begin(), _m.begin() + _entries);
      if (sort) std::sort(sorted_m.begin(), sorted_m.end(), pair_comp<pair_type>);

      os << "GEN::pairedvector [size= " << size() << ", capacity=" << capacity() << "]\n";
      if (sort) os << "sorted ";
      os << "entries {KEY, T}:\n";
      for (auto& p : sorted_m) os << "{" << p.first << ", " << p.second << "}\n";
    }

    /** @brief access the raw internally stored data (read-only)
     *
     *  @note Do NOT modify this object, otherwise the paired vector may be
     *  compromised.
     */
    const pairedvector_type& data() const { return _m; }

    /** @brief Activate and deactivate the debug functionality
     *
     *  @note This is only working, if the underlying data structures are
     *  compiled with included DEBUG output.
     *
     *  @param[in] isdebug  bool value for activation/deactivation
     *
     *  @author hiermeier @date 07/17 */
    void setDebugMode(bool isdebug) { base_type::setDebugMode(isdebug); }

   protected:
    /// access the internally stored data
    pairedvector_type& data() { return _m; }

    /** @brief internal clone method
     *
     *  @param[in] source  copy the source object into this.
     *
     *  @author hiermeier @date 05/17 */
    void clone(const class_type& source)
    {
      base_type::clone(source);

      const size_t src_capacity = source.capacity();
      if (capacity() < src_capacity) resize(src_capacity);

      _entries = source.size();

      iterator it = _m.begin();
      for (const auto& src_i : source) *(it++) = src_i;
    }

   private:
    /// raw data vector
    pairedvector_type _m;

    /// number of entries (without the pre-allocated entries, _entries != _m.size())
    size_t _entries;

  };  // class pairedvector

  /** @brief print sorted %pairedvector (row format)
   *
   *  @param[in] os  output stream
   *  @param[in] vec paired vector object which is going to be printed
   *  @return The modified output stream.
   *
   *  @author hiermeier @date 05/17 */
  template <typename Key, typename T0, typename... Ts>
  std::ostream& operator<<(std::ostream& os, const pairedvector<Key, T0, Ts...>& vec)
  {
    std::vector<std::pair<Key, T0>> sorted_vec(vec.begin(), vec.end());
    std::sort(sorted_vec.begin(), sorted_vec.end(), pair_comp<std::pair<Key, T0>>);

    os << "[size= " << vec.size() << ", capacity=" << vec.capacity() << "]";
    for (auto& p : sorted_vec)
      os << " (" << p.first << ", " << std::setprecision(3) << std::scientific << p.second << ")";

    return os;
  }

  /// template alias for the pairedvector using the default_insert_policy
  template <typename Key, typename T>
  using default_pairedvector = pairedvector<Key, T, default_insert_policy<Key, T>>;

  /// template alias for the pairedvector using the quick_insert_policy
  template <typename Key, typename T>
  using quick_pairedvector = pairedvector<Key, T, quick_insert_policy<Key, T>>;

}  // namespace GEN

#endif /* PAIREDVECTOR_H_ */
