/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_GID_VECTOR_H
#define BACI_LIB_UTILS_GID_VECTOR_H

#include "baci_config.H"

#include "baci_utils_exceptions.H"

#include <Epetra_Comm.h>
#include <Teuchos_RCPDecl.hpp>

#include <map>
#include <set>
#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    /*!
     * \brief Add nodal GID on this processor to existing list of GIDs
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] nodegid               nodal GID
     * @param[out] my_gid_vec           vector/set with my node GIDs
     */
    void AddOwnedNodeGID(const Discretization& dis, int nodegid, std::vector<int>& my_gid_vec);
    void AddOwnedNodeGID(const Discretization& dis, int nodegid, std::set<int>& my_gid_set);

    /*!
     * \brief Add nodal GIDs on this processor to existing list from list with global GIDs
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] global_node_gid_vec   vector/set with all node GIDs
     * @param[out] my_gid_vec           vector/set with my node GIDs
     */
    template <typename T, typename U>
    void AddOwnedNodeGIDFromList(
        const Discretization& dis, const T& global_node_gid_vec, U& my_gid_list)
    {
      for (const int nodegid : global_node_gid_vec) AddOwnedNodeGID(dis, nodegid, my_gid_list);
    }

    /*!
     * \brief check, whether node with GID is owned by this processor
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] node_gid              GID of node to be checked
     * @return                          indicates, whether node is owned by this processor
     */
    bool IsNodeGIDOnThisProc(const DRT::Discretization& dis, int node_gid);

    //! Merge map @p map_in (key of type @p T and value of type @p U) from all procs to a merged
    //! map (key of type @p T and value of type @p U). It is distributed to to all procs.
    template <typename T, typename U>
    std::map<T, U> BroadcastMap(const std::map<T, U>& map_in, const Epetra_Comm& comm);

    //! Merge vector of pairs @p pairs_in (items of type @p T and @p U) from all procs to a merged
    //! vector (items of type @p T and @p U). The merged items are in an unspecified order. It is
    //! distributed to to all procs.
    template <typename T, typename U>
    std::vector<std::pair<T, U>> BroadcastPairVector(
        const std::vector<std::pair<T, U>>& pairs_in, const Epetra_Comm& comm);

    //! Merge @p set_in (items of type @p T) from all procs to a merged set (items of type @p T). It
    //! is distributed to to all procs.
    template <typename T>
    std::set<T> BroadcastSet(const std::set<T>& set_in, const Epetra_Comm& comm);

    //! Merge vector @p vec_in (items of type @p T) from all procs to a merged vector (items of type
    //! @p T). The items of are in an unspecified order. It is distributed to to all procs.
    template <typename T>
    std::vector<T> BroadcastVector(const std::vector<T>& vec_in, const Epetra_Comm& comm);
  }  // namespace UTILS
}  // namespace DRT


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace DRT::UTILS::DETAIL
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
    vec1 = BroadcastVector(my_gid_vec1, comm);
    vec2 = BroadcastVector(my_gid_vec2, comm);

    dsassert(vec1.size() == vec2.size(), "Vectors must have the same length.");

    // reconstruct map-like object
    for (unsigned i = 0; i < vec1.size(); ++i)
    {
      vec_out1.emplace_back(vec1[i]);
      vec_out2.emplace_back(vec2[i]);
    }
  }
}  // namespace DRT::UTILS::DETAIL

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
std::map<T, U> DRT::UTILS::BroadcastMap(const std::map<T, U>& map_in, const Epetra_Comm& comm)
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
std::vector<std::pair<T, U>> DRT::UTILS::BroadcastPairVector(
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
std::set<T> DRT::UTILS::BroadcastSet(const std::set<T>& set_in, const Epetra_Comm& comm)
{
  std::vector<T> vec_in, vec_out;
  for (const auto& val : set_in) vec_in.emplace_back(val);
  vec_out = BroadcastVector(vec_in, comm);
  std::set<T> set_out;
  for (const auto& val : vec_out) set_out.insert(val);
  return set_out;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
std::vector<T> DRT::UTILS::BroadcastVector(const std::vector<T>& vec_in, const Epetra_Comm& comm)
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

BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_GID_VECTOR_H
