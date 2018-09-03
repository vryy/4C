/*----------------------------------------------------------------------------*/
/**
\file xfield_field_coupling_dofset.H

\brief DoF set for coupling a xfield and a field discretization at a common
       interface

\maintainer Michael Hiermeier

\date Sep 26, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_DOFSET_H_
#define SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_DOFSET_H_

#include "../drt_lib/drt_dofset_fixed_size.H"


namespace XFEM
{
  namespace XFieldField
  {
    /** \brief DoF set for coupling a xfield and a field discretization at a common interface
     *
     *  This coupling dof set helps to couple two discretizations (e.g. XFEM and standard)
     *  by returning the maximum number of DoF's at one coupling node.
     *
     *  \author hiermeier \date 09/16 */
    class CouplingDofSet : public DRT::FixedSizeDofSet
    {
     public:
      /** \brief constructor of the coupling DoF set
       *
       *  \param my_num_reserve_dof_per_node (in): reserve this number of DoF's for all nodes
       *                                           (e.g. standard DoFs + enrichment DoFs)
       *  \param g_node_index_range          (in): index range of nodes
       *  \param g_num_std_dof_per_node      (in): number of standard DoF's per node (w/o enriched
       * DoFs) \param my_num_dofs_per_node        (in): map containing the actual number of DoFs per
       * node
       *
       *  \author hiermeier
       *  \date 09/16 */
      CouplingDofSet(const int& numreservedofpernode, const int& g_node_index_range,
          const int& g_num_std_dof_per_node, const std::map<int, int>& my_num_dofs_per_node);

      /// destructor
      virtual ~CouplingDofSet(){};

      /** \brief Get the GIDs of all DoFs of a node.
       *
       *  Ask the current DofSet for the gids of the dofs of this node. The
       *  required vector is created and filled on the fly. So better keep it
       *  if you need more than one dof gid.
       *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
       *
       *  Additional input nodal_dof_set_id: If the node contains more than one set of dofs, which
       * can be evaluated, the ID of the set needs to be given. Currently only the case for XFEM.
       *
       *  \param dof             (out): vector of dof gids (to be filled)
       *  \param node            (in) : the node
       *  \param nodal_dofset_id (in) : id of the nodal dofset
       *
       *  \author hiermeier \date 10/16 */
      virtual void Dof(
          std::vector<int>& dofs, const DRT::Node* node, unsigned nodal_dofset_id) const;

      /** \brief Get the number of standard DoFs per coupling node
       *
       *  This value is supposed to be constant over all nodes!
       *
       *  \author hiermeier
       *  \date 10/16 */
      int NumStandardDofPerNode() const;

     protected:
      virtual int NumDofPerNode(const DRT::Node& node) const;

     private:
      /** \brief Get the number of DoFs of the node with the given nodal global ID
       *
       *  \param node_gid (in): nodal global id
       *
       *  \author hiermeier
       *  \date 09/16 */
      int MyNumDofPerNode(const int& node_gid) const;

     private:
      /// map containing the number of DoF's per node
      std::map<int, int> my_num_dof_per_node_;

      int g_num_std_dof_per_node_;
    };
  }  // namespace XFieldField
}  // namespace XFEM



#endif /* SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_DOFSET_H_ */
