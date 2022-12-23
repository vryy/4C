/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef UTILS_H_GID_VECTOR
#define UTILS_H_GID_VECTOR

#include "dserror.H"
#include <Epetra_Comm.h>
#include <map>
#include <Teuchos_RCPDecl.hpp>
#include <vector>

namespace DRT
{
  class DiscretizationInterface;

  namespace UTILS
  {
    /*!
     * \brief Add nodal GID on this processor to existing vector of GIDs
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] nodegid               nodal GID
     * @param[out] my_gid_vec           vector with my node GIDs
     */
    void AddOwnedNodeGID(
        const DiscretizationInterface& dis, int nodegid, std::vector<int>& my_gid_vec);

    /*!
     * \brief Add nodal GIDs on this processor to existing vector from vector with global GIDs
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] gloabl_node_gid_vec   vector with all node GIDs
     * @param[out] my_gid_vec           vector with my node GIDs
     */
    void AddOwnedNodeGIDVector(const DiscretizationInterface& dis,
        const std::vector<int>& global_node_gid_vec, std::vector<int>& my_gid_vec);

    /*!
     * \brief check, whether node with GID is owned by this processor
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] node_gid              GID of node to be checked
     * @return                          indicates, whether node is owned by this processor
     */
    bool IsNodeGIDOnThisProc(const DRT::DiscretizationInterface& dis, int node_gid);

    /*!
     * \brief remove nodal GIDs from vector when they are owned by this processor
     *
     * @param[in] dis                   discretization, that holds nodes with gids
     * @param[in] node_gids_to_remove   vector with node GIDs to remove
     * @param[in/out] node_gid_vec      vector to remove GIDs from
     */
    void RemoveNodeGIDsFromVector(const DiscretizationInterface& dis,
        const std::vector<int>& node_gids_to_remove, std::vector<int>& node_gid_vec);

    //! Merge map @p map_in (key of type @p T and value of type @p U) from all procs to a merged
    //! map @p map_out (key of type @p T and value of type @p U). It is distributed to to all procs.
    template <typename T, typename U>
    void BroadcastMap(
        const std::map<T, U>& map_in, std::map<T, U>& map_out, const Epetra_Comm& comm);

    //! Merge vector of pairs @p pairs_in (items of type @p T and @p U) from all procs to a merged
    //! vector @p pairs_out (items of type @p T and @p U). The items of @p pairs_out are in an
    //! unspecified order. It is distributed to to all procs.
    template <typename T, typename U>
    void BroadcastPairVector(const std::vector<std::pair<T, U>>& pairs_in,
        std::vector<std::pair<T, U>>& pairs_out, const Epetra_Comm& comm);

    //! Merge vector @p vec_in (items of type @p T) from all procs to a merged vector @p vec_out
    //! (items of type @p T). The items of @p vec_out are in an unspecified order. It is distributed
    //! to to all procs.
    template <typename T>
    void BroadcastVector(
        const std::vector<T>& vec_in, std::vector<T>& vec_out, const Epetra_Comm& comm);
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
    BroadcastVector(my_gid_vec1, vec1, comm);
    BroadcastVector(my_gid_vec2, vec2, comm);

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
void DRT::UTILS::BroadcastMap(
    const std::map<T, U>& map_in, std::map<T, U>& map_out, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(map_in, vec1, vec2, comm);
  for (unsigned i = 0; i < vec1.size(); ++i) map_out.insert(std::make_pair(vec1[i], vec2[i]));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T, typename U>
void DRT::UTILS::BroadcastPairVector(const std::vector<std::pair<T, U>>& pairs_in,
    std::vector<std::pair<T, U>>& pairs_out, const Epetra_Comm& comm)
{
  std::vector<T> vec1;
  std::vector<U> vec2;
  DETAIL::BroadcastMapLikeToVetors<T, U>(pairs_in, vec1, vec2, comm);
  for (unsigned i = 0; i < vec1.size(); ++i)
    pairs_out.emplace_back(std::make_pair(vec1[i], vec2[i]));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void DRT::UTILS::BroadcastVector(
    const std::vector<T>& vec_in, std::vector<T>& vec_out, const Epetra_Comm& comm)
{
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
    comm.Broadcast(&vec_broadcast[0], size, iproc);

    // append communicated vector to vec_out
    for (const auto& item : vec_broadcast) vec_out.emplace_back(item);
  }
}

#endif  // #ifndef UTILS_H_VECTOR
