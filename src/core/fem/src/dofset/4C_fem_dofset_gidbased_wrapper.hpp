/*----------------------------------------------------------------------*/
/*! \file

 \brief subproxy functionality to dofsets

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_FEM_DOFSET_GIDBASED_WRAPPER_HPP
#define FOUR_C_FEM_DOFSET_GIDBASED_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_base.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::DOFSets
{
  /*! \brief This class extends a standard DofSet to a DofSet whose mapping uses global IDs

    This class is used like a proxy. It is meant to be used as secondary DofSet in a
    discretization if there are two only partially coupled Discretizations with matching
    nodes. If you have only partial coupling, the standard Lid matching based on the
    column map lids does not work anymore, since the secondary discretization may have
    a different amount of elements and entirely different maps. Therefore, we have to
    retrieve the dofs by searching for the matching GIDs, instead of the LIDs.

    This class holds a source discretization sourcedis_ and a source DofSet sourcedofset_.
    Within all Dof(..) and NumDof(..) methods, first, the source node/element is searched within the
    source discretization (based on the GID of the input node/element) and then used as a new input
    for the source DofSet.

    \note This class inherits from the dofset base class and has to provide the
          whole dofset functionality.

    \warning No GID based mapping is performed for face elements, but the
          standard versions are called instead. The reason for this is, that boundary elements
          (which are face elements) built from a condition do not have unique GIDs, so
          a unique GID based mapping is not possible.

    \author Andreas Rauch
    \author Anh-Tu Vuong
    \date 10/16    */
  class DofSetGIDBasedWrapper : public DofSetBase
  {
   public:
    //! Standard Constructor
    DofSetGIDBasedWrapper(Teuchos::RCP<Core::FE::Discretization> sourcedis,
        Teuchos::RCP<DofSetInterface> sourcedofset);

    //! Destructor
    ~DofSetGIDBasedWrapper() override;

    //! original DofSet has new dofs
    int assign_degrees_of_freedom(
        const Core::FE::Discretization& dis, const unsigned dspos, const int start) override;

    //! original DofSet has been reset
    void reset() override;

    /// Notify original dofset of new proxies
    void notify_assigned() override;

    /// our original DofSet dies
    void disconnect(DofSetInterface* dofset) override;

    /// Returns true if filled
    bool filled() const override { return sourcedofset_->filled(); };

    //! @name Access methods

    /// get number of nodal dofs
    int num_dof_per_node(
        const Core::Nodes::Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      if (not sourcedis_->have_global_node(node.id())) return 0;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node.id());
      check_is_assigned();
      return sourcedofset_->num_dof_per_node(*sourcenode);
    };

    /// Get number of dofs for given node
    int num_dof(const Core::Nodes::Node* node) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return 0;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->num_dof(sourcenode);
    }

    /// Get number of dofs for given element
    int num_dof(const Core::Elements::Element* element) const override
    {
      check_is_assigned();
      if (element->is_face_element()) return sourcedofset_->num_dof(element);
      if (not sourcedis_->have_global_element(element->id())) return 0;
      Core::Elements::Element* sourceele = sourcedis_->g_element(element->id());
      return sourcedofset_->num_dof(sourceele);
    }

    /// Get the gid of a dof for given node
    int dof(const Core::Nodes::Node* node, int dof) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return -1;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->dof(sourcenode, dof);
    }

    /// Get the gid of all dofs of a node
    std::vector<int> dof(const Core::Nodes::Node* node) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return std::vector<int>();
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->dof(sourcenode);
    }

    /// Get the gid of all dofs of a node
    void dof(std::vector<int>& dof,     ///< vector of dof gids (to be filled)
        const Core::Nodes::Node* node,  ///< the node
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->dof(dof, sourcenode, nodaldofset);
    }

    /// Get the gid of all dofs of a element
    std::vector<int> dof(const Core::Elements::Element* element) const override
    {
      check_is_assigned();
      if (element->is_face_element()) return sourcedofset_->dof(element);
      if (not sourcedis_->have_global_element(element->id())) return std::vector<int>();
      Core::Elements::Element* sourceele = sourcedis_->g_element(element->id());
      return sourcedofset_->dof(sourceele);
    }

    /// Get the gid of all dofs of a node
    void dof(const Core::Nodes::Node* node, std::vector<int>& lm) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->dof(sourcenode, lm);
    }

    /// Get the gid of all dofs of a node
    void dof(const Core::Nodes::Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,  ///< first index of vector at which will be written to end
        std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      return sourcedofset_->dof(sourcenode, startindex, lm);
    }

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void dof(const Core::Elements::Element*
                 element,  ///< element which provides its expected number of DOFs per node
        const Core::Nodes::Node* node,  ///< node, for which you want the DOF positions
        std::vector<int>& lm  ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->have_global_node(node->id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->g_node(node->id());
      Core::Elements::Element* sourceele = sourcedis_->g_element(element->id());
      return sourcedofset_->dof(sourceele, sourcenode, lm);
    }

    /// Get the gid of a dof for given element
    int dof(const Core::Elements::Element* element, int dof) const override
    {
      check_is_assigned();
      if (element->is_face_element()) return sourcedofset_->dof(element, dof);
      if (not sourcedis_->have_global_element(element->id())) return -1;
      Core::Elements::Element* sourceele = sourcedis_->g_element(element->id());
      return sourcedofset_->dof(sourceele, dof);
    }


    /// Get the gid of all dofs of a element
    void dof(const Core::Elements::Element* element, std::vector<int>& lm) const override
    {
      check_is_assigned();
      if (element->is_face_element()) return sourcedofset_->dof(element, lm);
      if (not sourcedis_->have_global_element(element->id())) return;
      Core::Elements::Element* sourceele = sourcedis_->g_element(element->id());
      return sourcedofset_->dof(sourceele, lm);
    }

    /// Print this class
    void print(std::ostream& os) const override { sourcedofset_->print(os); };

    /// Print the dofsets in the static_dofsets_ list
    void print_all_dofsets(const Epetra_Comm& comm) const override
    {
      sourcedofset_->print_all_dofsets(comm);
    };

    /// Get Number of Global Elements of degree of freedom row map
    int num_global_elements() const override
    {
      check_is_assigned();
      return sourcedofset_->num_global_elements();
    };

    /// Get degree of freedom row map
    const Epetra_Map* dof_row_map() const override
    {
      check_is_assigned();
      return sourcedofset_->dof_row_map();
    };

    /// Get degree of freedom column map
    const Epetra_Map* dof_col_map() const override
    {
      check_is_assigned();
      return sourcedofset_->dof_col_map();
    };

    /// Get maximum GID of degree of freedom row map
    int max_all_gid() const override
    {
      check_is_assigned();
      return sourcedofset_->max_all_gid();
    };

    /// Get minimum GID of degree of freedom row map
    int min_all_gid() const override
    {
      check_is_assigned();
      return sourcedofset_->min_all_gid();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int max_gi_din_list(const Epetra_Comm& comm) const override
    {
      check_is_assigned();
      return sourcedofset_->max_gi_din_list(comm);
    };

    /// are the dof maps already initialized?
    bool initialized() const override
    {
      check_is_assigned();
      return sourcedofset_->initialized();
    };

    //@}

   private:
    /// check if \ref assign_degrees_of_freedom was called on parent dofset
    void check_is_assigned() const;

    //! source discretization
    Teuchos::RCP<Core::FE::Discretization> sourcedis_;

    //! source dofset wrapped in this class
    Teuchos::RCP<DofSetInterface> sourcedofset_;

    bool isassigned_;
  };
}  // namespace Core::DOFSets


FOUR_C_NAMESPACE_CLOSE

#endif
