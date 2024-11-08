// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_PACK_HELPERS_HPP
#define FOUR_C_COMM_PACK_HELPERS_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_pairedvector.hpp"

#include <array>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{

  namespace Internal
  {
    template <class T>
    using is_enum_class =
        std::integral_constant<bool, !std::is_convertible<T, int>::value && std::is_enum<T>::value>;


  }
  /*!
   * @brief A type trait to check whether the type T supports the
   * pack(Core::Communication::PackBuffer&) method.
   *
   * @note T does not necessarily need to derive from some interface class
   * @{
   */
  template <typename T, typename AlwaysVoid = void>
  constexpr bool is_packable = false;

  template <typename T>
  constexpr bool is_packable<T, std::void_t<decltype(std::declval<const std::decay_t<T>>().pack(
                                    std::declval<Core::Communication::PackBuffer&>()))>> = true;
  ///@}

  /*!
   * @brief A type trait to check whether the type T supports the
   * unpack(std::vector<char>::size_type&, const std::vector<char>&) method.
   *
   * @note T does not necessarily need to derive from some interface class
   * @{
   */
  template <typename T, typename AlwaysVoid = void>
  constexpr bool is_unpackable = false;

  template <typename T>
  constexpr bool is_unpackable<T, std::void_t<decltype(std::declval<std::decay_t<T>>().unpack(
                                      std::declval<Core::Communication::UnpackBuffer&>()))>> = true;
  ///@}

  //! @name Routines to help pack stuff into a char vector

  /*!
   * \brief Add stuff to the PackBuffer.
   *
   * This function works for all trivially copyable types, i.e. types that can be copied with
   * memcpy(). This includes all POD types, but also some user-defined types.
   */
  template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>, int> = 0>
  inline void add_to_pack(PackBuffer& data, const T& stuff)
  {
    data.add_to_pack(stuff);
  }

  /*!
   * \brief Add stuff to a char vector data
   *
   * This method adds an array to data
   * \param[in,out] data  char vector stuff shall be added to
   * \param[in] stuff  ptr to stuff that has length stuffsize (in byte)
   * \param[in] stuffsize length of stuff in byte
   */
  template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>, int> = 0>
  void add_to_pack(PackBuffer& data, const T* stuff, const int stuffsize)
  {
    data.add_to_pack(stuff, stuffsize);
  }

  /**
   * Add an object that implements a `pack()` method to the buffer.
   */
  template <typename T, std::enable_if_t<is_packable<T>, int> = 0>
  void add_to_pack(PackBuffer& data, const T& obj)
  {
    obj.pack(data);
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::vector<T>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::vector<T> that get's added to stuff
   */
  template <typename T>
  void add_to_pack(PackBuffer& data, const std::vector<T>& stuff)
  {
    int numele = stuff.size();
    add_to_pack(data, numele);

    // If T is trivially copyable, we can just copy the bytes. Otherwise, recursively call the
    // pack function for every element. Note that vector<bool> is a special case, as it does not
    // provide the data() method.
    if constexpr (std::is_trivially_copyable_v<T> && !std::is_same_v<T, bool>)
    {
      add_to_pack(data, stuff.data(), numele * sizeof(T));
    }
    else
    {
      for (const auto& elem : stuff) add_to_pack(data, elem);
    }
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::array<T, n>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::array<T, n> that get's added to stuff
   */
  template <typename T, std::size_t numentries>
  void add_to_pack(PackBuffer& data, const std::array<T, numentries>& stuff)
  {
    add_to_pack(data, stuff.data(), numentries * sizeof(T));
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::map<T,U>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::map<T,U> that get's added to stuff
   */
  template <typename T, typename U>
  void add_to_pack(PackBuffer& data, const std::map<T, U>& stuff)
  {
    int numentries = (int)stuff.size();
    add_to_pack(data, numentries);

    for (const auto& entry : stuff)
    {
      add_to_pack(data, entry.first);
      add_to_pack(data, entry.second);
    }
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::unordered_map<T,U>
   *
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::unordered_map<T,U> that get's added to stuff
   */
  template <typename T, typename U>
  void add_to_pack(PackBuffer& data, const std::unordered_map<T, U>& stuff)
  {
    int numentries = (int)stuff.size();
    add_to_pack(data, numentries);

    for (const auto& entry : stuff)
    {
      add_to_pack(data, entry.first);
      add_to_pack(data, entry.second);
    }
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::pair<T,U>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::pair<T,U> that get's added to stuff
   */
  template <typename T, typename U>
  void add_to_pack(PackBuffer& data, const std::pair<T, U>& stuff)
  {
    add_to_pack(data, stuff.first);
    add_to_pack(data, stuff.second);
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is a template for Pairedvector<Ts...>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff Pairedvector<Ts...> that get's added to stuff
   */
  template <typename... Ts>
  void add_to_pack(PackBuffer& data, const Core::Gen::Pairedvector<Ts...>& stuff)
  {
    int numentries = (int)stuff.size();
    add_to_pack(data, numentries);

    int i = 0;
    for (const auto& colcurr : stuff)
    {
      add_to_pack(data, colcurr.first);
      add_to_pack(data, colcurr.second);
      ++i;
    }

    if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
  }

  /*!
   * \brief Add stuff to the end of a char vector data first
   *
   * This method is a template for std::vector< Pairedvector<Ts...> >
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::vector<Pairedvector<Ts...> > that get's added to stuff
   */
  template <typename... Ts>
  void add_to_pack(PackBuffer& data, const std::vector<Core::Gen::Pairedvector<Ts...>>& stuff)
  {
    int numentries = (int)stuff.size();
    add_to_pack(data, numentries);

    int i = 0;
    for (auto& paired_vec : stuff)
    {
      add_to_pack(data, paired_vec);
      ++i;
    }

    if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is templated for std::set<T, U>
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff std::set<T> that get's added to stuff
   */
  template <typename T, typename U>
  void add_to_pack(PackBuffer& data, const std::set<T, U>& stuff)
  {
    int numentries = (int)stuff.size();
    add_to_pack(data, numentries);

    // iterator
    typename std::set<T, U>::const_iterator colcurr;

    int i = 0;
    for (colcurr = stuff.begin(); colcurr != stuff.end(); ++colcurr)
    {
      add_to_pack(data, *colcurr);
      ++i;
    }

    if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is an overload of the above template for Core::LinAlg::SerialDenseMatrix
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff Core::LinAlg::SerialDenseMatrix that get's added to stuff
   */
  void add_to_pack(PackBuffer& data, const Core::LinAlg::SerialDenseMatrix& stuff);

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is an overload of the above template for Core::LinAlg::SerialDenseMatrix
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff Core::LinAlg::SerialDenseVector that get's added to stuff
   */
  void add_to_pack(PackBuffer& data, const Core::LinAlg::SerialDenseVector& stuff);

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is an overload of the above template for Matrix
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff Matrix that get's added to data
   */
  template <unsigned i, unsigned j>
  void add_to_pack(PackBuffer& data, const Core::LinAlg::Matrix<i, j>& stuff)
  {
    add_to_pack(data, i);
    add_to_pack(data, j);
    add_to_pack(data, stuff.data(), stuff.m() * stuff.n() * sizeof(double));
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is an overload of the above template for an STL vector containing matrices
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff vector of matrices that get's added to data
   */
  template <unsigned i, unsigned j>
  void add_to_pack(PackBuffer& data, const std::vector<Core::LinAlg::Matrix<i, j>>& stuff)
  {
    // add length of vector to be packed so that later the vector can be restored with correct
    // length when unpacked
    int vectorlength = stuff.size();
    add_to_pack(data, vectorlength);

    for (int p = 0; p < vectorlength; ++p)
    {
      const double* A = stuff[p].data();

      // add all data in vector to pack
      add_to_pack(data, A, i * j * sizeof(double));
    }
  }

  /*!
   * \brief Add stuff to the end of a char vector data
   *
   * This method is an overload of the above template for string
   * \param[in,out] data char string stuff shall be added to
   * \param[in] stuff string that get's added to stuff
   */
  void add_to_pack(PackBuffer& data, const std::string& stuff);

  /**
   * Template to forward to the implementation on UnpackBuffer.
   */
  template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>, int> = 0>
  void extract_from_pack(UnpackBuffer& buffer, T& stuff)
  {
    buffer.extract_from_pack(stuff);
  }

  /**
   * Template to forward to the implementation on UnpackBuffer.
   */
  template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>, int> = 0>
  void extract_from_pack(UnpackBuffer& buffer, T* stuff, std::size_t stuff_size)
  {
    buffer.extract_from_pack(stuff, stuff_size);
  }

  double extract_double(UnpackBuffer& buffer);

  /**
   * Extract an object that implements an `unpack()` method from the buffer.
   */
  template <typename T, std::enable_if_t<is_unpackable<T>, int> = 0>
  void extract_from_pack(UnpackBuffer& data, T& obj)
  {
    obj.unpack(data);
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::vector<T>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::vector<T> to extract from data
   */
  template <typename T>
  void extract_from_pack(UnpackBuffer& buffer, std::vector<T>& stuff)
  {
    int dim = 0;
    buffer.extract_from_pack(dim);
    stuff.resize(dim);

    if constexpr (std::is_same_v<T, bool>)
    {
      // Note: this loop cannot be range-based due to the quirks of std::vector<bool>.
      for (int i = 0; i < dim; ++i)
      {
        bool value;
        extract_from_pack(buffer, value);
        stuff[i] = value;
      }
    }
    else if constexpr (std::is_trivially_copyable_v<T>)
    {
      extract_from_pack(buffer, stuff.data(), dim * sizeof(T));
    }
    else
    {
      for (auto& elem : stuff) extract_from_pack(buffer, elem);
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::array<T, n>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::map<T,n> to extract from data
   */
  template <typename T, std::size_t numentries>
  void extract_from_pack(UnpackBuffer& buffer, std::array<T, numentries>& stuff)
  {
    int size = numentries * sizeof(T);
    extract_from_pack(buffer, stuff.data(), size);
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::map<T,U>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented by
   * this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::map<T,U> to extract from data
   */
  template <typename T, typename U>
  void extract_from_pack(UnpackBuffer& buffer, std::map<T, U>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();

    for (int i = 0; i < numentries; i++)
    {
      T first;
      U second;
      extract_from_pack(buffer, first);
      extract_from_pack(buffer, second);

      // add to map
      stuff.insert(std::pair<T, U>(first, second));
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::unordered_map<T,U>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data  char vector where stuff is extracted from
   * \param[out] stuff std::unordered_map<T,U> to extract from data
   */
  template <typename T, typename U>
  void extract_from_pack(UnpackBuffer& buffer, std::unordered_map<T, U>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();

    for (int i = 0; i < numentries; i++)
    {
      T first;
      U second;
      extract_from_pack(buffer, first);
      extract_from_pack(buffer, second);

      // add to map
      stuff.insert({first, second});
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::pair<T,U>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::pair<T,U> to extract from data
   */
  template <typename T, typename U>
  void extract_from_pack(UnpackBuffer& buffer, std::pair<T, U>& stuff)
  {
    extract_from_pack(buffer, stuff.first);
    extract_from_pack(buffer, stuff.second);
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::vector<std::pair<T,U>>*
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::vector<std::pair<T,U>> to extract from data
   */
  template <typename T, typename U>
  void extract_from_pack(UnpackBuffer& buffer, std::vector<std::pair<T, U>>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();

    for (int i = 0; i < numentries; i++)
    {
      T first;
      U second;
      extract_from_pack(buffer, first);
      extract_from_pack(buffer, second);

      // add to map
      stuff.push_back(std::pair<T, U>(first, second));
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is a template for stuff of type Pairedvector<Key,T0,Ts...>*
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff Pairedvector<Key,T0,Ts...> to extract from data
   */
  template <typename Key, typename T0, typename... Ts>
  void extract_from_pack(UnpackBuffer& buffer, Core::Gen::Pairedvector<Key, T0, Ts...>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();
    stuff.resize(numentries);

    for (int i = 0; i < numentries; i++)
    {
      Key first;
      T0 second;
      extract_from_pack(buffer, first);
      extract_from_pack(buffer, second);

      // add to map
      stuff[first] = second;
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is a template for stuff of type std::vector< Pairedvector<Ts...> >*
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::vector< Pairedvector<Ts...> > to extract from data
   */
  template <typename... Ts>
  void extract_from_pack(UnpackBuffer& buffer, std::vector<Core::Gen::Pairedvector<Ts...>>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();
    stuff.resize(numentries);

    Core::Gen::Pairedvector<Ts...> paired_vec;
    for (int i = 0; i < numentries; i++)
    {
      extract_from_pack(buffer, paired_vec);


      // add to map
      stuff[i] = paired_vec;
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is templated for stuff of type std::set<T>
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff std::set<T> to extract from data
   */
  template <typename T, typename U>
  void extract_from_pack(UnpackBuffer& buffer, std::set<T, U>& stuff)
  {
    int numentries = 0;
    extract_from_pack(buffer, numentries);

    stuff.clear();

    for (int i = 0; i < numentries; i++)
    {
      T value;
      extract_from_pack(buffer, value);

      // add to set
      stuff.insert(value);
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is an overload of the above template for stuff of type
   * Core::LinAlg::SerialDenseMatrix
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff Core::LinAlg::SerialDenseMatrix to extract from data
   */
  void extract_from_pack(UnpackBuffer& buffer, Core::LinAlg::SerialDenseMatrix& stuff);

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is an overload of the above template for stuff of type
   * Core::LinAlg::SerialDenseVector
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff Core::LinAlg::SerialDenseVector to extract from data
   */
  void extract_from_pack(UnpackBuffer& buffer, Core::LinAlg::SerialDenseVector& stuff);

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is an overload of the above template for stuff of type Matrix
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff Matrix to extract from data
   */
  template <unsigned int i, unsigned int j>
  void extract_from_pack(UnpackBuffer& buffer, Core::LinAlg::Matrix<i, j>& stuff)
  {
    int m = 0;
    extract_from_pack(buffer, m);
    if (m != i) FOUR_C_THROW("first dimension mismatch");
    int n = 0;
    extract_from_pack(buffer, n);
    if (n != j) FOUR_C_THROW("second dimension mismatch");
    extract_from_pack(buffer, stuff.data(), stuff.m() * stuff.n() * sizeof(double));
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is an overload of the above template for stuff of type STL vector of matrices
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff STL vector of matrices to extract from data
   */
  template <unsigned int i, unsigned int j>
  void extract_from_pack(UnpackBuffer& buffer, std::vector<Core::LinAlg::Matrix<i, j>>& stuff)
  {
    // get length of vector to be extracted and allocate according amount of memory for all
    // extracted data
    int vectorlength;
    extract_from_pack(buffer, vectorlength);

    // resize vector stuff appropriately
    stuff.resize(vectorlength);

    for (int p = 0; p < vectorlength; ++p)
    {
      double* A = stuff[p].data();

      // actual extraction of data
      buffer.extract_from_pack(A, i * j * sizeof(double));
    }
  }

  /*!
   * \brief Extract stuff from a char vector data and increment position
   *
   * This method is an overload of the above template for stuff of type string
   *
   * \param[in,out] position place in data where to extract stuff. Position will be incremented
   * by this method
   * \param[in] data char vector where stuff is extracted from
   * \param[out] stuff string to extract from data
   */
  void extract_from_pack(UnpackBuffer& buffer, std::string& stuff);
  //@}

  /*!
   *  \brief Extract the type id at a given @p position from a @p data vector, assert if the
   * extracted type id matches the @p desired_type_id, and return the extracted type id.
   *
   * During extraction the position is incremented.
   *
   * \param[in,out] position Position in data vector where type id is extracted
   * \param[in] data Data vector where type id is extracted from
   * \param[in] desired_type_id Id of the desired type
   */
  int extract_and_assert_id(UnpackBuffer& buffer, const int desired_type_id);
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
