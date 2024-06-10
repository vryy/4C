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
    void Reset() override;

    /// Notify original dofset of new proxies
    void NotifyAssigned() override;

    /// our original DofSet dies
    void Disconnect(DofSetInterface* dofset) override;

    /// Returns true if filled
    bool Filled() const override { return sourcedofset_->Filled(); };

    //! @name Access methods

    /// get number of nodal dofs
    int NumDofPerNode(
        const Core::Nodes::Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      if (not sourcedis_->HaveGlobalNode(node.Id())) return 0;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node.Id());
      check_is_assigned();
      return sourcedofset_->NumDofPerNode(*sourcenode);
    };

    /// Get number of dofs for given node
    int NumDof(const Core::Nodes::Node* node) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return 0;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->NumDof(sourcenode);
    }

    /// Get number of dofs for given element
    int NumDof(const Core::Elements::Element* element) const override
    {
      check_is_assigned();
      if (element->IsFaceElement()) return sourcedofset_->NumDof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return 0;
      Core::Elements::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->NumDof(sourceele);
    }

    /// Get the gid of a dof for given node
    int Dof(const Core::Nodes::Node* node, int dof) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return -1;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, dof);
    }

    /// Get the gid of all dofs of a node
    std::vector<int> Dof(const Core::Nodes::Node* node) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return std::vector<int>();
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode);
    }

    /// Get the gid of all dofs of a node
    void Dof(std::vector<int>& dof,     ///< vector of dof gids (to be filled)
        const Core::Nodes::Node* node,  ///< the node
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(dof, sourcenode, nodaldofset);
    }

    /// Get the gid of all dofs of a element
    std::vector<int> Dof(const Core::Elements::Element* element) const override
    {
      check_is_assigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return std::vector<int>();
      Core::Elements::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele);
    }

    /// Get the gid of all dofs of a node
    void Dof(const Core::Nodes::Node* node, std::vector<int>& lm) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, lm);
    }

    /// Get the gid of all dofs of a node
    void Dof(const Core::Nodes::Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,  ///< first index of vector at which will be written to end
        std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, startindex, lm);
    }

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void Dof(const Core::Elements::Element*
                 element,  ///< element which provides its expected number of DOFs per node
        const Core::Nodes::Node* node,  ///< node, for which you want the DOF positions
        std::vector<int>& lm  ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      check_is_assigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      Core::Nodes::Node* sourcenode = sourcedis_->gNode(node->Id());
      Core::Elements::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, sourcenode, lm);
    }

    /// Get the gid of a dof for given element
    int Dof(const Core::Elements::Element* element, int dof) const override
    {
      check_is_assigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, dof);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return -1;
      Core::Elements::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, dof);
    }


    /// Get the gid of all dofs of a element
    void Dof(const Core::Elements::Element* element, std::vector<int>& lm) const override
    {
      check_is_assigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, lm);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return;
      Core::Elements::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, lm);
    }

    /// Print this class
    void Print(std::ostream& os) const override { sourcedofset_->Print(os); };

    /// Print the dofsets in the static_dofsets_ list
    void PrintAllDofsets(const Epetra_Comm& comm) const override
    {
      sourcedofset_->PrintAllDofsets(comm);
    };

    /// Get Number of Global Elements of degree of freedom row map
    int NumGlobalElements() const override
    {
      check_is_assigned();
      return sourcedofset_->NumGlobalElements();
    };

    /// Get degree of freedom row map
    const Epetra_Map* dof_row_map() const override
    {
      check_is_assigned();
      return sourcedofset_->dof_row_map();
    };

    /// Get degree of freedom column map
    const Epetra_Map* DofColMap() const override
    {
      check_is_assigned();
      return sourcedofset_->DofColMap();
    };

    /// Get maximum GID of degree of freedom row map
    int MaxAllGID() const override
    {
      check_is_assigned();
      return sourcedofset_->MaxAllGID();
    };

    /// Get minimum GID of degree of freedom row map
    int MinAllGID() const override
    {
      check_is_assigned();
      return sourcedofset_->MinAllGID();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int MaxGIDinList(const Epetra_Comm& comm) const override
    {
      check_is_assigned();
      return sourcedofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    bool Initialized() const override
    {
      check_is_assigned();
      return sourcedofset_->Initialized();
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
