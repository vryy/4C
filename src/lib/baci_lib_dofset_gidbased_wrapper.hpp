/*----------------------------------------------------------------------*/
/*! \file

 \brief subproxy functionality to dofsets

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_LIB_DOFSET_GIDBASED_WRAPPER_HPP
#define FOUR_C_LIB_DOFSET_GIDBASED_WRAPPER_HPP

#include "baci_config.hpp"

#include "baci_lib_discret.hpp"
#include "baci_lib_dofset_base.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  /*! \brief This class extends a standard DofSet to a DofSet whose mapping uses global IDs

    This class is used like a proxy. It is meant to be used as secondary DofSet in a
    Discretization if there are two only partially coupled Discretizations with matching
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
    DofSetGIDBasedWrapper(Teuchos::RCP<DRT::Discretization> sourcedis,
        Teuchos::RCP<DRT::DofSetInterface> sourcedofset);

    //! Destructor
    ~DofSetGIDBasedWrapper() override;

    //! original DofSet has new dofs
    int AssignDegreesOfFreedom(
        const Discretization& dis, const unsigned dspos, const int start) override;

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
    int NumDofPerNode(const Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      if (not sourcedis_->HaveGlobalNode(node.Id())) return 0;
      DRT::Node* sourcenode = sourcedis_->gNode(node.Id());
      CheckIsAssigned();
      return sourcedofset_->NumDofPerNode(*sourcenode);
    };

    /// Get number of dofs for given node
    int NumDof(const Node* node) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return 0;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->NumDof(sourcenode);
    }

    /// Get number of dofs for given element
    int NumDof(const Element* element) const override
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->NumDof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return 0;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->NumDof(sourceele);
    }

    /// Get the gid of a dof for given node
    int Dof(const Node* node, int dof) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return -1;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, dof);
    }

    /// Get the gid of all dofs of a node
    std::vector<int> Dof(const Node* node) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return std::vector<int>();
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode);
    }

    /// Get the gid of all dofs of a node
    void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const Node* node,            ///< the node
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(dof, sourcenode, nodaldofset);
    }

    /// Get the gid of all dofs of a element
    std::vector<int> Dof(const Element* element) const override
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return std::vector<int>();
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele);
    }

    /// Get the gid of all dofs of a node
    void Dof(const Node* node, std::vector<int>& lm) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, lm);
    }

    /// Get the gid of all dofs of a node
    void Dof(const Node* node,      ///< node, for which you want the dof positions
        const unsigned startindex,  ///< first index of vector at which will be written to end
        std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
    ) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, startindex, lm);
    }

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void Dof(
        const Element* element,  ///< element which provides its expected number of DOFs per node
        const Node* node,        ///< node, for which you want the DOF positions
        std::vector<int>& lm     ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, sourcenode, lm);
    }

    /// Get the gid of a dof for given element
    int Dof(const Element* element, int dof) const override
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, dof);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return -1;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, dof);
    }


    /// Get the gid of all dofs of a element
    void Dof(const Element* element, std::vector<int>& lm) const override
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, lm);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
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
      CheckIsAssigned();
      return sourcedofset_->NumGlobalElements();
    };

    /// Get degree of freedom row map
    const Epetra_Map* DofRowMap() const override
    {
      CheckIsAssigned();
      return sourcedofset_->DofRowMap();
    };

    /// Get degree of freedom column map
    const Epetra_Map* DofColMap() const override
    {
      CheckIsAssigned();
      return sourcedofset_->DofColMap();
    };

    /// Get maximum GID of degree of freedom row map
    int MaxAllGID() const override
    {
      CheckIsAssigned();
      return sourcedofset_->MaxAllGID();
    };

    /// Get minimum GID of degree of freedom row map
    int MinAllGID() const override
    {
      CheckIsAssigned();
      return sourcedofset_->MinAllGID();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int MaxGIDinList(const Epetra_Comm& comm) const override
    {
      CheckIsAssigned();
      return sourcedofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    bool Initialized() const override
    {
      CheckIsAssigned();
      return sourcedofset_->Initialized();
    };

    //@}

   private:
    /// check if \ref AssignDegreesOfFreedom was called on parent dofset
    void CheckIsAssigned() const;

    //! source discretization
    Teuchos::RCP<DRT::Discretization> sourcedis_;

    //! source dofset wrapped in this class
    Teuchos::RCP<DRT::DofSetInterface> sourcedofset_;

    bool isassigned_;
  };
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif
