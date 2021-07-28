/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef UTILS_H_GID_VECTOR
#define UTILS_H_GID_VECTOR

namespace DRT
{
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
        Teuchos::RCP<DRT::Discretization> dis, int nodegid, std::vector<int>& my_gid_vec);

    /*!
     * \brief Add nodal GIDs on this processor to existing vector from vector with global GIDs
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] gloabl_node_gid_vec   vector with all node GIDs
     * @param[out] my_gid_vec           vector with my node GIDs
     */
    void AddOwnedNodeGIDVector(Teuchos::RCP<DRT::Discretization> dis,
        const std::vector<int>& global_node_gid_vec, std::vector<int>& my_gid_vec);

    /*!
     * \brief check, whether node with GID is owned by this processor
     *
     * @param[in] dis                   discretization, that holds nodes with GIDs
     * @param[in] node_gid              GID of node to be checked
     * @return                          indicates, whether node is owned by this processor
     */
    bool IsNodeGIDOnThisProc(Teuchos::RCP<DRT::Discretization> dis, int node_gid);

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
    void RemoveNodeGIDsFromVector(Teuchos::RCP<DRT::Discretization> dis,
        const std::vector<int>& node_gids_to_remove, std::vector<int>& node_gid_vec);
  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_H_VECTOR
