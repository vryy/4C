/*---------------------------------------------------------------------*/
/*! \file

\brief These classes are meant as policies for the paired objects and
       are NOT supposed to be accessed directly. Use the methods of the derived
       classes, instead.

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_HEADERS_PAIREDOBJ_INSERT_POLICY_H_
#define SRC_HEADERS_PAIREDOBJ_INSERT_POLICY_H_

//#define DEBUG_INSERT_POLICY

#include <vector>
#include <utility>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <iostream>
#include <limits>

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
#include "../drt_lib/drt_dserror.H"
#endif

namespace GEN
{
  /// binary function
  template <typename T>
  inline bool pair_comp(const T& i, const T& j)
  {
    return i.first < j.first;
  }

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class default_insert_policy
  {
    typedef default_insert_policy<Key, T> class_type;

   protected:
    typedef std::pair<Key, T> pair_type;
    typedef std::vector<pair_type> pairedvector_type;
    typedef typename pairedvector_type::iterator iterator;
    typedef typename pairedvector_type::const_iterator const_iterator;

   public:
    default_insert_policy(){/* empty */};

    default_insert_policy(const class_type& source){/*empty*/};

    ~default_insert_policy() = default;

    iterator begin(const iterator& first) { return first; }

    const_iterator begin(const const_iterator& first) const { return first; }

    iterator end(const iterator& last) { return last; }

    const_iterator end(const const_iterator& last) const { return last; }

    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      for (iterator it = data.begin(); it != last; ++it)
      {
        if (it->first == k) return it;
      }
      return last;
    }

    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      for (const_iterator it = data.begin(); it != last; ++it)
      {
        if (it->first == k) return it;
      }
      return last;
    }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it != last) return it->second;

      if (entries >= data.size()) throw std::length_error("pairedvector::operator[]");

      ++entries;
      last->first = k;
      return last->second;
    }

    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it == last) throw std::runtime_error("pairedvector::at(): invalid key");

      return it->second;
    }

    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      const_iterator it = find(k, data, entries);
      if (it == last) throw std::runtime_error("pairedvector::at(): invalid key");

      return it->second;
    }

    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      return unique_entries;
    }

    void swap(class_type& x)
    { /* empty */
    }

    void clone(const class_type& source)
    { /* empty */
    }

    void setDebugMode(bool isdebug)
    {
#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      _isdebug = isdebug;
#else
      std::stringstream msg;
      msg << __LINE__ << " - " << __FUNCTION__ << " in " << __FILE__ << ":\n";
      msg << "Set the define \"DEBUG_INSERT_POLICY\" first to"
             " use the debug functionality!";
      throw std::runtime_error(msg.str());
#endif
    }

    /** @brief return a constant number of offset values
     *
     *   The default strategy does not use this feature. See quick_insert_policy
     *   for an example.
     *
     *  @author hiermeier @date 11/17 */
    static size_t capacity_offset() { return 0; }

   protected:
    bool _isdebug = false;
  };

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class quick_insert_policy : public default_insert_policy<Key, T>
  {
    typedef quick_insert_policy<Key, T> class_type;

   protected:
    typedef default_insert_policy<Key, T> base_type;
    typedef typename base_type::pair_type pair_type;
    typedef typename base_type::pairedvector_type pairedvector_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;

   public:
    quick_insert_policy() = default;

    quick_insert_policy(const class_type& source) { clone(source); }

    ~quick_insert_policy() = default;

    /** @brief return a constant number of offset values
     *
     *  This number will be added to a new capacity and can be used for
     *  sentinel values for example.
     *
     *  @author hiermeier @date 11/17 */
    static size_t capacity_offset() { return 1; }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it != last) return it->second;

      if (entries + capacity_offset() >= data.size())
        throw std::length_error("pairedvector::get()");

      ++entries;
      return last->second;
    }

    /// copy from the base class
    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);
      if (it == last) throw std::runtime_error("pairedvector::at(): invalid key");

      return it->second;
    }

    /// copy from the base class
    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      const_iterator last = data.begin() + entries;
      const_iterator it = find(k, data, entries);
      if (it == last) throw std::runtime_error("pairedvector::at(): invalid key");

      return it->second;
    }

    /// copy from the base class
    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      return get(k, data, entries);
    }

    /** @brief Quick linear search routine
     *
     *  The unnecessary bound check of the for loop is avoided by introduction
     *  of a sentinel value at the one past the last position.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Internally stored data.
     *  @param[in] entries Number of current entries.
     *  @return A iterator to the pair with the given key is returned. If the
     *          key is not a part of the stored data, an iterator to one past
     *          the last position will be returned.
     *
     *  @author hiermeier @date 11/17 */
    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      // set sentinel values
      (data.begin() + entries)->first = k;

      iterator it = data.begin();
      for (;;)
      {
        if (it->first == k) return it;
        ++it;
      }

      exit(EXIT_FAILURE);
    }

    /** @brief Quick linear search routine (const version)
     *
     *  The unnecessary bound check of the for loop is avoided by introduction
     *  of a sentinel value at the one past the last position. Therefore a
     *  technically const_cast becomes necessary. Nevertheless, the real data
     *  set is not touched.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Internally stored data.
     *  @param[in] entries Number of current entries.
     *  @return A const_iterator to the pair with the given key is returned. If
     *          the key is not a part of the stored data, an iterator to one past
     *          the last position will be returned.
     *
     *  @author hiermeier @date 11/17 */
    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      // set sentinel values
      pairedvector_type& mutable_data = const_cast<pairedvector_type&>(data);
      (mutable_data.begin() + entries)->first = k;

      const_iterator cit = data.begin();
      for (;;)
      {
        if (cit->first == k) return cit;
        ++cit;
      }

      exit(EXIT_FAILURE);
    }

    /** @brief repetitive access of the paired_vector entries in the same order
     *
     *  This method is very helpful, if a big paired_vector structure must be
     *  accessed repetitively in the same order. In this case only the first
     *  repetition (i.e. rep_count == 0) uses the internal linear search
     *  routine, while all following access repetitions (i.e. rep_count > 0)
     *  use the internal id-map under the assumption of an constant access
     *  pattern.
     *
     *  \note This method can not be used, if the access pattern changes during
     *  successive repetitions!
     *
     *  @param[in] k         Unique identification key of the stored data.
     *  @param[in] rep_count Repetition counter.
     *  @param[in] data      Internal stored data set.
     *  @param[in] entries   Updated number of entries.
     *  @return The value corresponding to the given key.
     *
     *  @author hiermeier @date 11/17 */
    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
#if defined(DEBUG) || defined(DEBUG_INSERT_POLICY)
      if (rep_count < -1)
        dserror(
            "The repetition counter is not allowed to "
            "be smaller than -1!");
#endif

      switch (rep_count)
      {
        case -1:
        {
          return get(k, data, entries);
        }
        case 0:
        {
#if defined(DEBUG) || defined(DEBUG_INSERT_POLICY)
          if (_prev_repcount != -1 and _prev_repcount != rep_count and !_isfilled)
            dserror(
                "At the end of each repetition Complete "
                "must be called!");
#endif
          init_id_map();
          _isfilled = false;
          return first_access(k, data, entries);
        }
        default:
        {
#if defined(DEBUG) || defined(DEBUG_INSERT_POLICY)
          if (_prev_repcount != rep_count and !_isfilled)
            dserror(
                "At the end of each repetition Complete "
                "must be called!");
#endif
          _isfilled = false;
          _prev_repcount = rep_count;
          return subsequent_access(k, data, entries);
        }
      }
    }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (_isfilled) return unique_entries;

      // These things must be done only once after the first complete call,
      // i.e. repetition counter equal to zero.
      if (_prev_repcount == 0)
      {
        _id_map_size = _id_map.size();
        if (_id_index != _id_map_size)
        {
          std::stringstream msg;
          msg << "Size (== " << _id_map.size() << " ) <--> count (== " << _id_index
              << " ) mismatch!\n";
          throw std::runtime_error(msg.str());
        }
      }

      _isfilled = true;
      // reset id-index after each single complete call
      _id_index = 0;

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "\nCOMPLETE\n";
      }
#endif

      return unique_entries;
    }

    void clone(const class_type& source)
    {
      base_type::clone(source);

      _id_map = source._id_map;
      _id_map_size = source._id_map_size;
      _id_index = source._id_index;

      _prev_repcount = source._prev_repcount;
      _isfilled = source._isfilled;
    }

   private:
    /** @brief Initiate the internally stored ID map
     *
     *  The old map is cleared, while the old size is used as an estimate for
     *  the necessary storage space.
     *
     *  @author hiermeier @date 11/17 */
    void init_id_map()
    {
      // nothing to do, if it is a successive access in the first repetition
      // Note: The _isfilled check guarantees that the maps are cleared, even if
      // only one GP had been considered in the previous attempt.
      if (_prev_repcount == 0 and not _isfilled) return;

      _id_map.clear();
      _id_map.reserve(_id_map_size);
      _prev_repcount = 0;

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << __FUNCTION__ << std::endl;
      }
#endif
    }

    T& first_access(const Key k, pairedvector_type& data, size_t& entries)
    {
      iterator last = data.begin() + entries;
      iterator it = find(k, data, entries);

      // By using push_back we follow the std::vector policy if the set
      // capacity is exceeded. See previous reserve call in init_id_map().
      _id_map.push_back(it - data.begin());
      ++_id_index;

      if (it == last)
      {
        if (entries + capacity_offset() >= data.size())
          throw std::length_error("quick_insert_policy::get()");

        ++entries;
      }

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "repetition count  = " << _prev_repcount << "\n";
        std::cout << "_id_map index     = " << _id_index << "\n";
        std::cout << "data id           = " << _id_map.back() << "\n";
        std::cout << "size of _id_map   = " << _id_map.size() << "\n";
        std::cout << "key               = " << k << "\n\n";
      }
#endif

      return it->second;
    }

    T& subsequent_access(const Key k, pairedvector_type& data, size_t& entries)
    {
      pair_type& curr_pair = data[get_id()];
#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "index #" << _id_index << std::endl;
        std::cout << "repetition count = " << _prev_repcount << "\n";
        std::cout << "stored key       = " << curr_pair.first << "\n";
        std::cout << "received key     = " << k << "\n\n";
      }
#endif
      if (curr_pair.first != k)
      {
        std::stringstream msg;
        msg << "quick_insert_policy::subsequent_access(): Mismatch of stored "
               "key and given key: "
            << curr_pair.first << " <--> " << k << "\n";
        throw std::runtime_error(msg.str());
      }

      return curr_pair.second;
    }

    size_t get_id() { return _id_map[_id_index++]; }

   private:
    std::vector<size_t> _id_map;
    size_t _id_map_size = 512;
    size_t _id_index = 0;

    int _prev_repcount = -1;
    bool _isfilled = false;
  };

  /*--------------------------------------------------------------------------*/
  template <typename Key, typename T>
  class insert_and_sort_policy : public default_insert_policy<Key, T>
  {
    typedef insert_and_sort_policy<Key, T> class_type;

   protected:
    typedef default_insert_policy<Key, T> base_type;
    typedef typename base_type::pair_type pair_type;
    typedef typename base_type::pairedvector_type pairedvector_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;

   public:
    insert_and_sort_policy()
        : _non_unique_vec(0),
          _non_unique_entries(0),
          _dyn_max_allowed_capacity(CONST_MAX_ALLOWED_CAPACITY),
          _isfilled(true),
          _max_value(0.0){};

    ~insert_and_sort_policy(){};

    iterator begin(const iterator& first)
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::begin(first);
    }

    const_iterator begin(const const_iterator& first) const
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::begin(first);
    }

    iterator end(const iterator& last)
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::end(last);
    }

    const_iterator end(const const_iterator& last) const
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::end(last);
    }

    iterator find(const Key k, pairedvector_type& data, size_t entries)
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::find(k, data, entries);
    }

    const_iterator find(const Key k, const pairedvector_type& data, size_t entries) const
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::find(k, data, entries);
    }

    T& get(const Key k, pairedvector_type& data, size_t& entries)
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::get(k, data, entries);
    }

    T& at(const Key k, pairedvector_type& data, size_t entries)
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::at(k, data, entries);
    }

    const T& at(const Key k, const pairedvector_type& data, size_t entries) const
    {
      throwIfNotFilled(__LINE__, __FUNCTION__);
      return base_type::at(k, data, entries);
    }

    T& repetitive_access(const Key k, const int rep_count, pairedvector_type& data, size_t& entries)
    {
      return this->operator()(k, data, entries);
    }

    /** @brief Insert values into the non-unique vector
     *
     *  @note This routine makes it necessary that you call complete before you can
     *  access the inserted values.
     *  @note Values with an absolute value smaller than the relative machine
     *  precision are ignored.
     *
     *  @param[in] k       Unique key value.
     *  @param[in] data    Unique paired data vector (unused at this point)
     *  @param[in] entries Number of unique data in the data vector.
     *  @return Reference to a new slot to insert data.
     *
     *  @author hiermeier @date 07/17 */
    T& operator()(const Key k, pairedvector_type& data, size_t& entries)
    {
      if (_non_unique_entries == _dyn_max_allowed_capacity and
          (not setDynamicMaxAllowedCapacity(2 * entries)))
      {
        entries = midComplete(data, entries);
      }

      if (_isfilled)
      {
        _isfilled = false;
        initCapacity();
        ++_non_unique_entries;
      }

      pair_type* pair_ptr = _non_unique_vec.data() + currentId();

      setMaxValue(pair_ptr->second);

      /* If reasonable values have been inserted, we check the capacity and go
       * to the next entry. Otherwise, we stick to the current entry. */
      if (std::abs(pair_ptr->second) > getRelativeMachinePrecision())
      {
        increaseCapacity();
        ++_non_unique_entries;
        pair_ptr = _non_unique_vec.data() + currentId();
      }

      pair_ptr->first = k;
      return pair_ptr->second;
    }

    inline bool isFilled() const { return _isfilled; }

    size_t complete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (_isfilled) return unique_entries;

      const size_t new_unique_entries = groupAndMerge(unique_vec, unique_entries);

      finalPostComplete();

      return new_unique_entries;
    }

    void swap(class_type& x)
    {
      base_type::swap(x);

      const bool my_isfilled = _isfilled;
      _isfilled = x._isfilled;
      x._isfilled = my_isfilled;

      const size_t my_non_unique_entries = _non_unique_entries;
      _non_unique_entries = x._non_unique_entries;
      x._non_unique_entries = my_non_unique_entries;

      const size_t my_dyn_max_allowed_capacity = _dyn_max_allowed_capacity;
      _dyn_max_allowed_capacity = x._dyn_max_allowed_capacity;
      x._dyn_max_allowed_capacity = my_dyn_max_allowed_capacity;

      const double my_max_value = _max_value;
      _max_value = x._max_value;
      x._max_value = my_max_value;

      _non_unique_vec.swap(x._non_unique_vec);
    }

    void clone(const class_type& source) { throwIfNotFilled(__LINE__, __FUNCTION__); }

   private:
    /** @brief Internal complete call
     *
     *  This is a complete call done during the a standard insert call into the
     *  non-unique vector. Its purpose is to group and merge the non-unique data
     *  as soon as some bound for the size of the non-unique vector is reached.
     *  In contrast to a default complete call the already allocated memory of
     *  the non-unique vector is not freed afterwards, since the insertion of
     *  new values is not yet finished.
     *
     *  @param[in/out] unique_vec     This object will contain the new
     *                                completed unique data.
     *  @param[in]     unique_entries Previous/old number of unique data.
     *  @return The number of the unique data in the end of this routine.
     *
     *  @author hiermeier @data 07/17 */
    size_t midComplete(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      if (_isfilled) return unique_entries;

      const size_t new_unique_entries = groupAndMerge(unique_vec, unique_entries);

      midPostComplete();

      return new_unique_entries;
    }

    /** @brief Throw a run time error if any method is called which
     *  presumes a filled state
     *
     *  @param[in] linenumber  Line number in the calling function.
     *  @param[in] functname   Function name of the calling function.
     *
     *  @author hiermeier @date 07/17 */
    inline void throwIfNotFilled(int linenumber, const std::string& functname) const
    {
      if (not _isfilled)
      {
        std::stringstream msg;
        msg << "LINE " << linenumber << " in " << functname
            << ": access denied! "
               "call complete first.";
        throw std::runtime_error(msg.str());
      }
    }

    /** @brief Group all entries in the %_non_unique_vec and insert the
     *  merged values of each group into the %unique_vec
     *
     *  @param[out] unique_vec     Vector containing the current unique entries.
     *                             This vector will be modified by inserting
     *                             the entries collected in the non-unique vector.
     *  @param[in]  unique_entries Number of old unique entries.
     *  @return Number of the new unique entries contained in unique_vec.
     *
     *  @author hiermeier @date 07/17 */
    size_t groupAndMerge(pairedvector_type& unique_vec, const size_t unique_entries)
    {
      // get total number of entries and merge all entries temporal in the
      // _non_unique_vec member
      const size_t merged_entries = unique_entries + _non_unique_entries;
      if (_non_unique_vec.size() < merged_entries + 1)
        _non_unique_vec.resize(merged_entries + 1, pair_type());

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "unique_vec:\n";
        print(unique_vec);
        std::cout << "_non_unique_vec:\n";
        print(_non_unique_vec);
      }
#endif

      std::copy(unique_vec.begin(), unique_vec.begin() + unique_entries,
          _non_unique_vec.begin() + _non_unique_entries);

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "merged _non_unique_vec + unique_vec:\n";
        print(_non_unique_vec);
      }
#endif

      // group the entries
      const size_t num_grps = groupData(_non_unique_vec.begin(), merged_entries);

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "sorted merged _non_unique_vec + unique_vec:\n";
        print(_non_unique_vec);
      }
#endif

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "number of groups = " << num_grps << "\n";
      }
#endif

      if (unique_vec.size() < num_grps) unique_vec.resize(num_grps);

      // copy and merge all non-unique entries into one unique vector
      const typename pairedvector_type::iterator last_result = mergeGroupData(
          _non_unique_vec.begin(), _non_unique_vec.begin() + merged_entries, unique_vec.begin());

      const size_t corrected_length = static_cast<size_t>(last_result - unique_vec.begin());

      // reset the remaining part of the unique vector
      if (corrected_length < unique_entries)
        std::fill(last_result, unique_vec.begin() + unique_entries, pair_type());

#if defined(DEBUG_INSERT_POLICY) || defined(DEBUG)
      if (this->_isdebug)
      {
        std::cout << "corrected num_unique = " << corrected_length << "\n";
        std::cout << "new unique data:\n";
        print(unique_vec);
        std::cout << "\n\n";
      }
#endif

      return corrected_length;
    }

    /** @brief reset stuff after an internal complete call
     *
     *  In contrast to the %finalPostComplete method we do not free the
     *  allocated memory for the %_non_unique_vec, instead we just clear the
     *  content. The %_max_value and the current necessary dynamic
     *  maximal allowed capacity are also kept.
     *
     *  @author hiermeier @date 07/17 */
    void midPostComplete()
    {
      std::fill(
          _non_unique_vec.begin(), _non_unique_vec.begin() + _non_unique_entries, pair_type());

      if (_max_capacity < _non_unique_entries) _max_capacity = _non_unique_entries;

      _non_unique_entries = 0;
      _isfilled = true;
    }

    /** @brief reset all class members at the end of a final externally
     *  initialized complete call
     *
     *  The allocated memory for the %_non_unique_vec is freed. Furthermore,
     *  all class members which are only important during one call are reset.
     *
     *  @author hiermeier @date 07/17 */
    void finalPostComplete()
    {
      _non_unique_vec.clear();

      if (_max_capacity < _non_unique_entries) _max_capacity = _non_unique_entries;

      _dyn_max_allowed_capacity = CONST_MAX_ALLOWED_CAPACITY;
      _max_value = 0.0;
      _non_unique_entries = 0;
      _isfilled = true;
    }

    /** @brief Group data into clusters with same KEY values
     *
     *  Group data in the given interval by creating clusters of pairs with the
     *  same KEY value.
     *
     *  @param[in] sbegin  iterator pointing at the first element of the interval
     *  @return Number of groups.
     *
     *  @author hiermeier @date 07/17 */
    size_t groupData(typename pairedvector_type::iterator sbegin, const size_t num_entries) const
    {
      switch (num_entries)
      {
        case 0:
          return 0;
        case 1:
          return 1;
        case 2:
        {
          if (sbegin->first == (sbegin + 1)->first) return 1;
          return 2;
        }
        default:
          return groupBigData(sbegin, sbegin + num_entries);
      }
    }

    /** @brief Group data into clusters with same KEY values
     *
     *  Group data in the given interval by creating clusters of pairs with the
     *  same KEY value. No extra memory is allocated since the grouping works
     *  by a simple swap routine.
     *
     *  @param[in] sbegin  iterator pointing at the first element of the interval
     *  @param[in] slast   iterator pointing one past the last element of the interval
     *  @return Number of groups.
     *
     *  @author hiermeier @date 07/17 */
    size_t groupBigData(typename pairedvector_type::iterator sbegin,
        const typename pairedvector_type::iterator slast) const
    {
      typename pairedvector_type::iterator result = sbegin;
      const typename pairedvector_type::iterator prev_last = slast - 1;

      size_t num_grps = 1;

      for (;;)
      {
        // set sentinel accordingly
        *slast = *result;
        sbegin = result + 1;

        for (;;)
        {
          pair_type& val = *sbegin;
          if (val.first == result->first)
          {
            if (sbegin == slast) break;

            // swap values
            (++result)->swap(val);
          }
          ++sbegin;
        }

        if (result == prev_last) break;

        ++result;
        ++num_grps;
      }

      return num_grps;
    }

    /** @brief Merge the values in each group
     *
     *  @pre The input values in between %sbegin and %slast must be grouped.
     *
     *  @param[in]  sbegin  constant iterator pointing to the first entry of
     *                      already grouped interval
     *  @param[in]  slast   iterator pointing to one past the last
     *                      entry of the already grouped interval (used as
     *                      sentinel value)
     *  @param[out] result  iterator pointing to the merged result interval
     *  @return Return the iterator pointing to one past the last entry of the
     *          result interval
     *
     *  @author hiermeier @date 07/17 */
    typename pairedvector_type::iterator mergeGroupData(
        typename pairedvector_type::const_iterator sbegin,
        typename pairedvector_type::iterator slast,
        typename pairedvector_type::iterator result) const
    {
      if (sbegin == slast) return result;

      // add sentinel
      slast->first = ((slast - 1)->first + 1);

      *result = *sbegin;
      for (;;)
      {
        const std::pair<Key, T>& val = *(++sbegin);
        if (val.first == result->first)
          result->second += val.second;
        else
        {
          if (sbegin == slast) break;

          // erase/overwrite zero entries
          if (std::abs(result->second) > getRelativeMachinePrecision()) ++result;

          *result = val;
        }
      }

      // Kick the last element out if the absolute value of the last result is
      // less than the relative machine precision.
      return (std::abs(result->second) > getRelativeMachinePrecision() ? ++result : result);
    }

    /// Return the local id of the last element
    inline size_t currentId() const { return (_non_unique_entries - 1); }

    /** @brief Set a new maximal allowed capacity of the non-unique vector
     *
     *  The maximal allowed capacity is only modified if the new dynamic bound
     *  value is larger as the current one.
     *
     *  @param[in] dyn_bound New dynamic bound value.
     *  @return If the threshold is modified, TRUE is returned, otherwise FALSE.
     *
     *  @author hiermeier @date 07/17 */
    bool setDynamicMaxAllowedCapacity(const size_t dyn_bound)
    {
      if (dyn_bound > _dyn_max_allowed_capacity)
      {
        _dyn_max_allowed_capacity = dyn_bound;
        return true;
      }
      return false;
    }

    /** @brief Set maximal absolute value occurring during the current insertion
     *         interval, i.e. between the last and next complete call.
     *
     *  @param[in] curr_val This is the current value which is going to be inserted.
     *
     *  @author hiermeier @date 07/17 */
    void setMaxValue(const double curr_val)
    {
      const double abs_curr_val = std::abs(curr_val);
      if (abs_curr_val > _max_value) _max_value = abs_curr_val;
    }

    /** @brief Get the relative machine precision
     *
     *  Relative to the max value occurring during the current insertion
     *  interval. See setMaxValue for more information.
     *
     *  @author hiermeier @date 07/17 */
    inline double getRelativeMachinePrecision() const { return _max_value * MACHINE_PRECISION; }

    inline void initCapacity()
    {
      _non_unique_vec.resize(std::max(1, static_cast<int>(_max_capacity)), pair_type());
    }

    void increaseCapacity()
    {
      const size_t curr_length = _non_unique_vec.size();
      if (_non_unique_entries < curr_length) return;

      _non_unique_vec.resize(2 * curr_length, pair_type());
    }

    void print(const pairedvector_type& pvec, std::ostream& os = std::cout) const
    {
      os << "[Size = " << pvec.size() << "]\n"
         << "{ ";
      size_t c = 0;
      typename pairedvector_type::const_iterator first = pvec.begin();
      while (first != pvec.end())
      {
        const pair_type& pair = *first;
        os << "[" << c++ << "](" << pair.first << ", " << pair.second << ")";
        if (++first != pvec.end()) os << "; ";
      }
      os << " };\n" << std::flush;
    }

   private:
    /** temporal allocated vector containing the unsorted and therefore
     *  non-unique entries */
    pairedvector_type _non_unique_vec;

    /** @brief number of the inserted non-unique entries
     *
     *  Note that this number is in general not equal to the size of the
     *  _non_unique_vec, since memory is pre-allocated. */
    size_t _non_unique_entries;

    /** maximal necessary capacity for all paired-vector objects of the same
     *  types */
    static size_t _max_capacity;

    /// absolute machine precision
    static constexpr double MACHINE_PRECISION = 10 * std::numeric_limits<double>::epsilon();

    /// constant maximal allowed capacity
    static constexpr size_t CONST_MAX_ALLOWED_CAPACITY = 512;

    /** @brief Dynamically corrected maximal allowed capacity.
     *
     *  If the unique vector exceeds the constant capacity, the upper capacity
     *  bound must be increased to stay meaningful.
     *
     *  See the method %setDynamicMaxAllowedCapacity() for more information. */
    size_t _dyn_max_allowed_capacity;

    bool _isfilled;
    double _max_value;
  };
}  // namespace GEN

template <typename Key, typename T>
size_t GEN::insert_and_sort_policy<Key, T>::_max_capacity = 0;



#endif /* SRC_HEADERS_PAIREDOBJ_INSERT_POLICY_H_ */
