/*---------------------------------------------------------------------*/
/*! \file

\brief This class is meant as a replacement for std::maps of std::maps, when
       other storage and access characteristics are needed.

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_HEADERS_PAIREDMATRIX_H_
#define SRC_HEADERS_PAIREDMATRIX_H_

#include "pairedvector.H"

namespace GEN
{
  /// struct containing the base type of the pairedmatrix class
  template <typename Key, typename T, typename inner_insert_policy, typename outer_insert_policy>
  struct pairedmatrix_base
  {
    typedef pairedvector<Key, pairedvector<Key, T, inner_insert_policy>, outer_insert_policy> type;
  };

  template <typename Key, typename T, typename inner_insert_policy = default_insert_policy<Key, T>,
      typename outer_insert_policy =
          default_insert_policy<Key, pairedvector<Key, T, inner_insert_policy>>>
  class pairedmatrix
      : public pairedmatrix_base<Key, T, inner_insert_policy, outer_insert_policy>::type
  {
    typedef pairedvector<Key, T, inner_insert_policy> inner_pairedvector_type;
    typedef pairedvector<Key, inner_pairedvector_type, outer_insert_policy> base_type;
    typedef std::pair<Key, inner_pairedvector_type> pair_type;
    typedef std::vector<pair_type> pairedmatrix_type;
    typedef pairedmatrix<Key, T, inner_insert_policy, outer_insert_policy> class_type;

   public:
    /**
     *  @brief  constructor creates no elements, but reserves the maximum
     *          number of entries.
     *  @param reserve The number of elements that are preallocated
     */
    pairedmatrix(size_t reserve)
        : base_type(reserve, Key(), inner_pairedvector_type(reserve)), _max_row_capacity(reserve)
    {
    }

    /**
     *  @brief  empty constructor creates no elements and does not reserve any
     *          number of entries. Use resize as soon as you know the necessary
     *          number of elements.
     */
    pairedmatrix() : base_type(), _max_row_capacity(0) {}

    /**
     *  @brief  constructor creates no elements, but reserves the maximum
     *          number of entries.
     *  @param reserve The number of elements that are preallocated
     *  @param default_key default value for the key within the pair
     *  @param default_T   default value for the data within the pair
     */
    pairedmatrix(size_t reserve, Key default_key, inner_pairedvector_type default_T)
        : base_type(reserve, default_key, default_T), _max_row_capacity(default_T.capacity())
    {
    }

    /**
     *  @brief  copy constructor
     *
     *  @param[in] source %pairedmatrix object we want to copy.
     *  @param[in] type   Apply this copy type.
     */
    pairedmatrix(const inner_pairedvector_type& source, enum GEN::CopyType type = DeepCopy)
        : base_type(), _max_row_capacity(0)
    {
      const size_t src_max_row_capacity = maxRowCapacity(source);
      clear(default_pair(src_max_row_capacity));
      resize(source.capacity(), src_max_row_capacity);

      base_type::clone(source);

      switch (type)
      {
        case ShapeCopy:
        {
          for (pair_type& row : *this)
            for (std::pair<Key, T>& entry : row.second) entry.second = T();

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
    ~pairedmatrix() {}

    /**
     * @brief  clear the %pairedmatrix content
     *
     *  Erases all elements in a %pairedmatrix. Note that this function only
     *  erases the elements, and that if the elements themselves are
     *  pointers, the pointed-to memory is not touched in any way.
     *  Managing the pointer is the user's responsibility.
     *
     *  \note This method keeps the current max row capacity. If you want to
     *  reset the capacity as well, call the second clear routine using a
     *  default constructor as input argument.
     *
     *  \author hiermeier \date 07/17
     */
    void clear() { clear(_max_row_capacity); }

    void clear(const pair_type& x) { base_type::clear(x); }

    void resize(size_t new_size)
    {
      if (new_size > _max_row_capacity) _max_row_capacity = new_size;

      base_type::resize(new_size, default_pair(_max_row_capacity));
    }

    void resize(size_t new_size, const pair_type& x) { base_type::resize(new_size, x); }

    /** @brief assign operator
     *
     *  Perform a deep copy of the input object. */
    class_type& operator=(const base_type& source)
    {
      clone(source, DeepCopy);
      return *this;
    }

    void clone(const base_type& source, const enum CopyType type)
    {
      const size_t src_row_capacity = maxRowCapacity(source);
      clear(src_row_capacity);
      resize(source.capacity(), src_row_capacity);

      base_type::clone(source);

      switch (type)
      {
        case ShapeCopy:
        {
          for (pair_type& row : *this)
            for (std::pair<Key, T>& entry : row.second) entry.second = T();

          break;
        }
        default:
          break;
      }
    }

    void complete()
    {
      base_type::complete();
      for (pair_type& vec : *this) vec.second.complete();
    }

    /**
     *  @brief  print the %pairedmatrix
     *
     *  @param[in] os    Output stream.
     *  @param[in] sort  Sort the entries before the actual print is performed.
     *
     *  Print the %pairedmatrix information in column format. By default the
     *  entries are sorted with respect to their KEY entries.
     */
    void print(std::ostream& os = std::cout, bool sort = true) const
    {
      pairedmatrix_type sorted_m(base_type::begin(), base_type::end());
      if (sort) std::sort(sorted_m.begin(), sorted_m.end(), pair_comp<pair_type>);

      os << "GEN::pairedmatrix [size= " << base_type::size()
         << ", capacity=" << base_type::capacity()
         << ", max. capacity per row=" << maxRowCapacity(*this) << "]\n";
      if (sort) os << "sorted ";
      os << "entries {KEY, T}:\n";
      for (auto& p : sorted_m) os << "{" << p.first << ", " << p.second << "}\n";
    }

    /** @brief Activate and deactivate the debug functionality
     *
     *  @note This is only working, if the underlying data structures are
     *  compiled with included DEBUG output.
     *
     *  @param[in] isdebug
     */
    void setDebugMode(bool isdebug)
    {
      base_type::setDebugMode(isdebug);
      for (pair_type& pair : this->data()) pair.second.setDebugMode(isdebug);
    }

    /// return a default pair of correct types with the specified capacity
    static pair_type default_pair(const size_t row_capacity)
    {
      const Key empty_first = Key();
      const inner_pairedvector_type empty_second(row_capacity);

      return pair_type(empty_first, empty_second);
    }

   private:
    inline void clear(const size_t row_capacity)
    {
      base_type::clear(pair_type(Key(), inner_pairedvector_type(row_capacity)));
    }

    inline void resize(const size_t new_size, const size_t row_capacity)
    {
      _max_row_capacity = row_capacity;
      base_type::resize(new_size, default_pair(row_capacity));
    }

    /** @brief Detect the maximal row capacity of a given %pairedmatrix
     *
     *  @param mat %pairedmatrix object we want to investigate.
     *
     *  @author hiermeier @date 05/17 */
    inline size_t maxRowCapacity(const base_type& mat) const
    {
      size_t max_row_capacity = 0;
      for (const pair_type& row : mat.data())
      {
        const size_t row_capacity = row.second.capacity();
        if (max_row_capacity < row_capacity) max_row_capacity = row_capacity;
      }

      return max_row_capacity;
    }

   private:
    size_t _max_row_capacity;
  };  // class pairedmatrix

  /*--------------------------------------------------------------------------*/

  /** @brief (Re)set a %pairedvector
   *
   *  @param[in]  reserve_size  Maximal necessary capacity of the %pairedvector.
   *  @param[out] paired_vec    Reset this %pairedvector.
   *
   *  \author hiermeier \date 03/17 */
  template <typename Key, typename... Ts>
  inline void reset(const unsigned reserve_size, pairedvector<Key, Ts...>& paired_vec)
  {
    if (not paired_vec.empty()) paired_vec.clear();

    if (paired_vec.capacity() < reserve_size) paired_vec.resize(reserve_size);
  }

  /** @brief (Re)set a %pairedmatrix
   *
   *  @param[in]  reserve_size        New capacity of the %pairedmatrix.
   *  @param[in]  inner_reserve_size  New capacity of each row.
   *  @param[out] paired_mat          Reset this %pairedmatrix.
   *
   *  \author hiermeier \date 03/17 */
  template <typename Key, typename... Ts>
  inline void reset(const unsigned reserve_size, const unsigned row_reserve_size,
      pairedmatrix<Key, Ts...>& paired_mat)
  {
    if (not paired_mat.empty())
    {
      paired_mat.clear(pairedmatrix<Key, Ts...>::default_pair(row_reserve_size));
    }

    if (paired_mat.capacity() < reserve_size)
    {
      paired_mat.resize(reserve_size, pairedmatrix<Key, Ts...>::default_pair(row_reserve_size));
    }
  }

  /** @brief (Re)set a %pairedmatrix
   *
   *  @param[in]  reserve_size  New capacity of the %pairedmatrix,
   *                            used for rows and columns.
   *  @param[out] paired_mat    Reset this %pairedmatrix.
   *
   *  @author hiermeier @date 03/17 */
  template <typename Key, typename... Ts>
  inline void reset(const unsigned reserve_size, pairedmatrix<Key, Ts...>& paired_mat)
  {
    reset<Key, Ts...>(reserve_size, reserve_size, paired_mat);
  }

  /** @brief (Re)set a std::vector of paired-matrices
   *
   *  @param[in]  num_vec         New size of the %std::vector.
   *  @param[in]  reserve_size    New capacity of the paired object,
   *                              used for rows and columns.
   *  @param[out] vec_paired_obj  Reset this vector of paired objects.
   *
   *  @author hiermeier @date 03/17 */
  template <typename paired_type>
  inline void reset(
      const unsigned num_vec, const unsigned reserve_size, std::vector<paired_type>& vec_paired_obj)
  {
    if (vec_paired_obj.size() != num_vec)
    {
      vec_paired_obj.clear();
      const paired_type empty_paired_obj;
      vec_paired_obj.resize(num_vec, empty_paired_obj);
    }

    for (auto& paired_obj : vec_paired_obj) reset(reserve_size, paired_obj);
  }

  /** @brief (Re)set a std::vector of paired-objects
   *
   *  The vector size is kept constant.
   *
   *  @param[in]  reserve_size    New capacity of the paired object.
   *  @param[out] vec_paired_obj  Reset this vector of paired objects.
   *
   *  @author hiermeier @date 03/17 */
  template <typename paired_type>
  inline void reset(const unsigned reserve_size, std::vector<paired_type>& vec_paired_obj)
  {
    reset(vec_paired_obj.size(), reserve_size, vec_paired_obj);
  }

  /** @brief Weak reset a paired-vector
   *
   *  The capacity is kept constant. Furthermore, the already set key values
   *  are kept and the number of entries stays untouched. Only the corresponding
   *  set values will be reseted.
   *
   *  @param[out] vec_paired_obj  Weak reset this vector of paired objects.
   *
   *  @author hiermeier @date 11/17 */
  template <typename Key, typename T0, typename... Ts>
  inline void weak_reset(pairedvector<Key, T0, Ts...>& paired_vec)
  {
    for (auto& pair : paired_vec) pair.second = T0();
  }

  /** @brief Weak reset a paired-matrix
   *
   *  The capacities are kept constant. Furthermore, the already set key values
   *  are kept and the number of entries stays untouched. Only the corresponding
   *  set values will be reseted.
   *
   *  @param[out] vec_paired_obj  Weak reset this vector of paired objects.
   *
   *  @author hiermeier @date 11/17 */
  template <typename Key, typename... Ts>
  inline void weak_reset(pairedmatrix<Key, Ts...>& paired_mat)
  {
    for (auto& pair : paired_mat) weak_reset(pair.second);
  }

  /** @brief weak reset a std::vector of paired objects
   *
   *  The vector size as well as the capacities are kept constant. Furthermore,
   *  the already set key values are kept and the number of entries stays
   *  untouched. Only the correpsonding set values will be reseted.
   *
   *  @param[out] vec_paired_obj  Weak reset this vector of paired objects.
   *
   *  @author hiermeier @date 11/17 */
  template <typename paired_type>
  inline void weak_reset(std::vector<paired_type>& vec_paired_obj)
  {
    for (auto& paired_obj : vec_paired_obj) weak_reset(paired_obj);
  }

  /** @brief Increase the capacity of a paired_type object  if necessary
   *
   *  @param paired_obj Check capacity of this object and increase it if necessary
   *
   *  @return Capacity of the possibly modified paired object. If no modification
   *          took place, the old capacity is returned.
   *
   *  @author hiermeier @date 03/17 */
  template <typename paired_type>
  inline size_t increaseCapacity(paired_type& paired_obj)
  {
    size_t new_capacity = paired_obj.capacity();

    // if the capacity is still sufficient, do nothing
    if (new_capacity > paired_obj.size()) return new_capacity;

    switch (new_capacity)
    {
      // special case: capacity is zero
      case 0:
        ++new_capacity;
        break;
      // default case: increase old capacity by a factor of 2
      default:
        new_capacity *= 2;
        break;
    }

    paired_obj.resize(new_capacity);

    return new_capacity;
  }

  /** @brief Copy Source %std::vector of paired_type objects into a new %std::vector
   *
   *  @param[in]  source Source vector which is going to be copied.
   *  @param[out] target Target of the copy operation.
   *  @param[in]  type   Optional copy type. Per default shape and values are copied.
   *
   *  @author hiermeier @date 05/17 */
  template <typename paired_type>
  inline void copy(const std::vector<paired_type>& source, std::vector<paired_type>& target,
      const enum CopyType type = DeepCopy)
  {
    const unsigned vec_dim = source.size();

    target.resize(vec_dim);

    auto src_it = source.begin();
    for (auto it = target.begin(); it != target.end(); ++it, ++src_it) it->clone(*src_it, type);
  }

  /** @brief Copy Source of paired_type into a new object of paired_type
   *
   *  @param[in]  source Source object which is going to be copied.
   *  @param[out] target Target of the copy operation.
   *  @param[in]  type   Optional copy type. Per default shape and values are copied.
   *
   *  @author hiermeier @date 05/17 */
  template <typename paired_type>
  inline void copy(
      const paired_type& source, paired_type& target, const enum CopyType type = DeepCopy)
  {
    target.clone(source, type);
  }

  /** @brief Print a set of paired objects
   *
   *  @param[in]  vec_paired_obj  Print this object.
   *  @param[out] os              Use this output stream.
   *
   *  @author hiermeier @date 07/17 */
  template <typename paired_type>
  inline void print(const std::vector<paired_type>& vec_paired_obj, std::ostream& os = std::cout)
  {
    size_t i = 0;
    const size_t vec_size = vec_paired_obj.size();

    for (const paired_type& paired_obj : vec_paired_obj)
    {
      os << "component #" << ++i << " of " << vec_size << ":\n";
      paired_obj.print(os);
    }
  }

  /** @brief Complete a set of paired objects
   *
   *  @param[in]  vec_paired_obj  Complete all paired objects contained in this
   *                              vector.
   *
   *  @author hiermeier @date 07/17 */
  template <typename paired_type>
  inline void complete(std::vector<paired_type>& vec_paired_obj)
  {
    for (paired_type& paired_obj : vec_paired_obj) paired_obj.complete();
  }

  /** @brief Set debug mode status in each of the contained paired objects
   *
   *  @pre You have to activate the compiler define flag first!
   *
   *  @param[in] vec_paired_obj  Set of paired objects.
   *  @param[in] isdebug         New debug status.
   *
   *  @author hiermeier @date 07/17 */
  template <typename paired_type>
  inline void setDebugMode(std::vector<paired_type>& vec_paired_obj, bool isdebug)
  {
    for (paired_type& paired_obj : vec_paired_obj)
    {
      paired_obj.setDebugMode(isdebug);
    }
  }

  /// template alias for the pairedmatrix using the default_insert_policy
  template <typename Key, typename T>
  using default_pairedmatrix = pairedmatrix<Key, T, default_insert_policy<Key, T>,
      default_insert_policy<Key, default_pairedvector<Key, T>>>;

  /// template alias for the pairedmatrix using the quick_insert_policy
  template <typename Key, typename T>
  using quick_pairedmatrix = pairedmatrix<Key, T, quick_insert_policy<Key, T>,
      quick_insert_policy<Key, quick_pairedvector<Key, T>>>;
}  // namespace GEN


#endif /* SRC_HEADERS_PAIREDMATRIX_H_ */
