// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_MPI_UTILS_HPP
#define FOUR_C_COMM_MPI_UTILS_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <mpi.h>

#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  //! Merge map @p map_in (key of type @p T and value of type @p U) from all procs to a merged
  //! map (key of type @p T and value of type @p U). It is distributed to all procs.
  template <typename T, typename U>
  std::map<T, U> all_reduce(const std::map<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge map @p unordered map_in (key of type @p T and value of type @p U) from all procs to a
  //! merged unordered map (key of type @p T and value of type @p U). It is distributed to all
  //! procs.
  template <typename T, typename U>
  std::unordered_map<T, U> all_reduce(
      const std::unordered_map<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge unordered multimap @p map_in (key of type @p T and value of type @p U) from all procs
  //! to a merged unordered multimap (key of type @p T and value of type @p U). It is distributed
  //! to all procs.
  template <typename T, typename U>
  std::unordered_multimap<T, U> all_reduce(
      const std::unordered_multimap<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge vector of pairs @p pairs_in (items of type @p T and @p U) from all procs to a merged
  //! vector (items of type @p T and @p U). The merged items are in an unspecified order. It is
  //! distributed to all procs.
  template <typename T, typename U>
  std::vector<std::pair<T, U>> all_reduce(
      const std::vector<std::pair<T, U>>& pairs_in, const Epetra_Comm& comm);

  //! Merge @p set_in (items of type @p T) from all procs to a merged set (items of type @p T). It
  //! is distributed to all procs.
  template <typename T>
  std::set<T> all_reduce(const std::set<T>& set_in, const Epetra_Comm& comm);

  //! Merge vector @p vec_in (items of type @p T) from all procs to a merged vector (items of type
  //! @p T). The items of are in an unspecified order. The result is distributed to all procs.
  template <typename T>
  std::vector<T> all_reduce(const std::vector<T>& vec_in, const Epetra_Comm& comm);

  /**
   * Perform an all-reduce operation on a value @p value using the reduction operation @p
   * reduction_op. The result is distributed to all procs.
   */
  template <typename T>
  T all_reduce(const T& value, const std::function<T(const T&, const T&)>& reduction_op,
      const Epetra_Comm& comm);

  /**
   * Gather a value @p value from all procs to a vector of values. The result is distributed to all
   * procs.
   */
  template <typename T>
  std::vector<T> all_gather(const T& value, const Epetra_Comm& comm);

}  // namespace Core::Communication

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace Core::Communication::Internal
{
  //! Broadcast a map or vector<pair>
  template <typename T, typename U, typename M>
  void all_reduce_map_like_to_vectors(
      const M& map_in, std::vector<T>& vec_out1, std::vector<U>& vec_out2, const Epetra_Comm& comm)
  {
    // split map or std::vector<std::pair> into two vectors
    std::vector<T> my_gid_vec1;
    std::vector<U> my_gid_vec2;
    for (const auto& pair : map_in)
    {
      my_gid_vec1.emplace_back(pair.first);
      my_gid_vec2.emplace_back(pair.second);
    }
    std::vector<T> vec1 = all_reduce(my_gid_vec1, comm);
    std::vector<U> vec2 = all_reduce(my_gid_vec2, comm);

    FOUR_C_ASSERT(vec1.size() == vec2.size(), "Vectors must have the same length.");

    // reconstruct map-like object
    for (unsigned i = 0; i < vec1.size(); ++i)
    {
      vec_out1.emplace_back(vec1[i]);
      vec_out2.emplace_back(vec2[i]);
    }
  }
}  // namespace Core::Communication::Internal

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::map<T, U> Core::Communication::all_reduce(
    const std::map<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  Internal::all_reduce_map_like_to_vectors<T, U>(map_in, vec1, vec2, comm);
  std::map<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::unordered_map<T, U> Core::Communication::all_reduce(
    const std::unordered_map<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  Internal::all_reduce_map_like_to_vectors<T, U>(map_in, vec1, vec2, comm);
  std::unordered_map<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::unordered_multimap<T, U> Core::Communication::all_reduce(
    const std::unordered_multimap<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  Internal::all_reduce_map_like_to_vectors<T, U>(map_in, vec1, vec2, comm);
  std::unordered_multimap<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::vector<std::pair<T, U>> Core::Communication::all_reduce(
    const std::vector<std::pair<T, U>>& pairs_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  Internal::all_reduce_map_like_to_vectors<T, U>(pairs_in, vec1, vec2, comm);
  std::vector<std::pair<T, U>> pairs_out;
  for (unsigned i = 0; i < vec1.size(); ++i)
    pairs_out.emplace_back(std::make_pair(vec1[i], vec2[i]));
  return pairs_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
std::set<T> Core::Communication::all_reduce(const std::set<T>& set_in, const Epetra_Comm& comm)
{
  std::vector<T> vec_in, vec_out;
  for (const auto& val : set_in) vec_in.emplace_back(val);
  vec_out = all_reduce(vec_in, comm);
  std::set<T> set_out;
  for (const auto& val : vec_out) set_out.insert(val);
  return set_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
std::vector<T> Core::Communication::all_reduce(
    const std::vector<T>& vec_in, const Epetra_Comm& comm)
{
  std::vector<std::vector<T>> gathered = all_gather(vec_in, comm);

  // Merge the data from all processors into a single vector
  std::vector<T> all_vectors;
  for (auto& received : gathered)
  {
    all_vectors.insert(all_vectors.end(), received.begin(), received.end());
  }

  return all_vectors;
}


template <typename T>
T Core::Communication::all_reduce(const T& value,
    const std::function<T(const T&, const T&)>& reduction_op, const Epetra_Comm& comm)
{
  std::vector<T> all_values = all_gather(value, comm);
  return std::accumulate(all_values.cbegin() + 1, all_values.cend(), all_values[0], reduction_op);
}


template <typename T>
std::vector<T> Core::Communication::all_gather(const T& value, const Epetra_Comm& comm)
{
  MPI_Comm mpi_comm = dynamic_cast<const Epetra_MpiComm&>(comm).Comm();
  const int n_procs = comm.NumProc();

  // Use PackBuffer to serialize the data into a vector<char>
  PackBuffer buffer;
  add_to_pack(buffer, value);

  int local_size = buffer().size();
  std::vector<int> size_all_data(n_procs, 0);
  [[maybe_unused]] int err =
      MPI_Allgather(&local_size, 1, MPI_INT, size_all_data.data(), 1, MPI_INT, mpi_comm);
  FOUR_C_ASSERT(err == MPI_SUCCESS, "MPI_Allgather failed.");

  // Compute offset of the data from the i-th rank
  std::vector<int> accumulated_offset(n_procs);
  std::partial_sum(size_all_data.begin(), size_all_data.end() - 1, accumulated_offset.begin() + 1);

  std::vector<char> combined_buffer(accumulated_offset.back() + size_all_data.back());

  err = MPI_Allgatherv(buffer().data(), local_size, MPI_CHAR,
      /*full received buffer*/ combined_buffer.data(),
      /*sizes of processor contributions*/ size_all_data.data(),
      /*displacement of processor contributions in full buffer*/ accumulated_offset.data(),
      MPI_CHAR, mpi_comm);
  FOUR_C_ASSERT(err == MPI_SUCCESS, "MPI_Allgatherv failed.");

  // Merge the data from all processors into a single vector
  std::vector<T> all_data(n_procs);
  for (int i = 0; i < n_procs; ++i)
  {
    // Since UnpackBuffer cannot view data, we need to make a copy of the data.
    std::vector<char> local_buffer(combined_buffer.begin() + accumulated_offset[i],
        combined_buffer.begin() + accumulated_offset[i] + size_all_data[i]);

    UnpackBuffer unpack_buffer(local_buffer);
    extract_from_pack(unpack_buffer, all_data[i]);
  }

  return all_data;
}

FOUR_C_NAMESPACE_CLOSE

#endif