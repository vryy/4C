/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods to broadcast data across procs

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_COMM_BROADCAST_UTILS_HPP
#define FOUR_C_COMM_BROADCAST_UTILS_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>

#include <map>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  //! Merge map @p map_in (key of type @p T and value of type @p U) from all procs to a merged
  //! map (key of type @p T and value of type @p U). It is distributed to to all procs.
  template <typename T, typename U>
  std::map<T, U> broadcast(const std::map<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge map @p unordered map_in (key of type @p T and value of type @p U) from all procs to a
  //! merged unordered map (key of type @p T and value of type @p U). It is distributed to to all
  //! procs.
  template <typename T, typename U>
  std::unordered_map<T, U> broadcast(
      const std::unordered_map<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge unordered multimap @p map_in (key of type @p T and value of type @p U) from all procs
  //! to a merged unordered multimap (key of type @p T and value of type @p U). It is distributed
  //! to to all procs.
  template <typename T, typename U>
  std::unordered_multimap<T, U> broadcast(
      const std::unordered_multimap<T, U>& map_in, const Epetra_Comm& comm);

  //! Merge vector of pairs @p pairs_in (items of type @p T and @p U) from all procs to a merged
  //! vector (items of type @p T and @p U). The merged items are in an unspecified order. It is
  //! distributed to to all procs.
  template <typename T, typename U>
  std::vector<std::pair<T, U>> broadcast(
      const std::vector<std::pair<T, U>>& pairs_in, const Epetra_Comm& comm);

  //! Merge @p set_in (items of type @p T) from all procs to a merged set (items of type @p T). It
  //! is distributed to to all procs.
  template <typename T>
  std::set<T> broadcast(const std::set<T>& set_in, const Epetra_Comm& comm);

  //! Merge vector @p vec_in (items of type @p T) from all procs to a merged vector (items of type
  //! @p T). The items of are in an unspecified order. It is distributed to to all procs.
  template <typename T>
  std::vector<T> broadcast(const std::vector<T>& vec_in, const Epetra_Comm& comm);
}  // namespace Core::Communication

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace Core::Communication::DETAIL
{
  //! Broadcast a map or vector<pair>
  template <typename T, typename U, typename M>
  void BroadcastMapLikeToVetors(
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
    std::vector<T> vec1;
    std::vector<U> vec2;
    vec1 = broadcast(my_gid_vec1, comm);
    vec2 = broadcast(my_gid_vec2, comm);

    FOUR_C_ASSERT(vec1.size() == vec2.size(), "Vectors must have the same length.");

    // reconstruct map-like object
    for (unsigned i = 0; i < vec1.size(); ++i)
    {
      vec_out1.emplace_back(vec1[i]);
      vec_out2.emplace_back(vec2[i]);
    }
  }
}  // namespace Core::Communication::DETAIL

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::map<T, U> Core::Communication::broadcast(const std::map<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(map_in, vec1, vec2, comm);
  std::map<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::unordered_map<T, U> Core::Communication::broadcast(
    const std::unordered_map<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(map_in, vec1, vec2, comm);
  std::unordered_map<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::unordered_multimap<T, U> Core::Communication::broadcast(
    const std::unordered_multimap<T, U>& map_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(map_in, vec1, vec2, comm);
  std::unordered_multimap<T, U> map_out;
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
  return map_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::vector<std::pair<T, U>> Core::Communication::broadcast(
    const std::vector<std::pair<T, U>>& pairs_in, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(pairs_in, vec1, vec2, comm);
  std::vector<std::pair<T, U>> pairs_out;
  for (unsigned i = 0; i < vec1.size(); ++i)
    pairs_out.emplace_back(std::make_pair(vec1[i], vec2[i]));
  return pairs_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
std::set<T> Core::Communication::broadcast(const std::set<T>& set_in, const Epetra_Comm& comm)
{
  std::vector<T> vec_in, vec_out;
  for (const auto& val : set_in) vec_in.emplace_back(val);
  vec_out = broadcast(vec_in, comm);
  std::set<T> set_out;
  for (const auto& val : vec_out) set_out.insert(val);
  return set_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
std::vector<T> Core::Communication::broadcast(const std::vector<T>& vec_in, const Epetra_Comm& comm)
{
  std::vector<T> vec_out;
  for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
  {
    // communicate size of vector
    int size = static_cast<int>(vec_in.size());
    comm.Broadcast(&size, 1, iproc);

    // new vectors to be filled (by this proc, if MyPID == iproc or other procs by
    // communication)
    std::vector<T> vec_broadcast;
    if (iproc == comm.MyPID()) vec_broadcast = vec_in;

    // communicate vector
    vec_broadcast.resize(size);
    comm.Broadcast(vec_broadcast.data(), size, iproc);

    // append communicated vector to vec_out
    for (const auto& item : vec_broadcast) vec_out.emplace_back(item);
  }
  return vec_out;
}

FOUR_C_NAMESPACE_CLOSE

#endif