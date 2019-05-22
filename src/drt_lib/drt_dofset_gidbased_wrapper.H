/*----------------------------------------------------------------------*/
/*!

 \brief subproxy functionality to dofsets

\level 2

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/


#ifndef DRT_DOFSET_SUBSET_H_
#define DRT_DOFSET_SUBSET_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

#include "drt_dofset_base.H"

#include "drt_element.H"
#include "drt_node.H"
#include "drt_discret.H"

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
    virtual ~DofSetGIDBasedWrapper();

    //! original DofSet has new dofs
    virtual int AssignDegreesOfFreedom(
        const Discretization& dis, const unsigned dspos, const int start);

    //! original DofSet has been reset
    virtual void Reset();

    /// Notify original dofset of new proxies
    virtual void NotifyAssigned();

    /// our original DofSet dies
    virtual void Disconnect(DofSetInterface* dofset);

    /// Returns true if filled
    virtual bool Filled() const { return sourcedofset_->Filled(); };

    //! @name Access methods

    /// get number of nodal dofs
    virtual int NumDofPerNode(
        const Node& node  ///< node, for which you want to know the number of dofs
        ) const
    {
      if (not sourcedis_->HaveGlobalNode(node.Id())) return 0;
      DRT::Node* sourcenode = sourcedis_->gNode(node.Id());
      CheckIsAssigned();
      return sourcedofset_->NumDofPerNode(*sourcenode);
    };

    /// Get number of dofs for given node
    virtual int NumDof(const Node* node) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return 0;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->NumDof(sourcenode);
    }

    /// Get number of dofs for given element
    virtual int NumDof(const Element* element) const
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->NumDof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return 0;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->NumDof(sourceele);
    }

    /// Get the gid of a dof for given node
    virtual int Dof(const Node* node, int dof) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return -1;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, dof);
    }

    /// Get the gid of all dofs of a node
    virtual std::vector<int> Dof(const Node* node) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return std::vector<int>();
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode);
    }

    /// Get the gid of all dofs of a node
    virtual void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const Node* node,                    ///< the node
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
        ) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(dof, sourcenode, nodaldofset);
    }

    /// Get the gid of all dofs of a element
    virtual std::vector<int> Dof(const Element* element) const
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return std::vector<int>();
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele);
    }

    /// Get the gid of all dofs of a node
    virtual void Dof(const Node* node, std::vector<int>& lm) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, lm);
    }

    /// Get the gid of all dofs of a node
    virtual void Dof(const Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,      ///< first index of vector at which will be written to end
        std::vector<int>& lm  ///< already allocated vector to be filled with dof positions
        ) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      return sourcedofset_->Dof(sourcenode, startindex, lm);
    }

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    virtual void Dof(
        const Element* element,  ///< element which provides its expected number of DOFs per node
        const Node* node,        ///< node, for which you want the DOF positions
        std::vector<int>& lm     ///< already allocated vector to be filled with DOF positions
        ) const
    {
      CheckIsAssigned();
      if (not sourcedis_->HaveGlobalNode(node->Id())) return;
      DRT::Node* sourcenode = sourcedis_->gNode(node->Id());
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, sourcenode, lm);
    }

    /// Get the gid of a dof for given element
    virtual int Dof(const Element* element, int dof) const
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, dof);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return -1;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, dof);
    }


    /// Get the gid of all dofs of a element
    virtual void Dof(const Element* element, std::vector<int>& lm) const
    {
      CheckIsAssigned();
      if (element->IsFaceElement()) return sourcedofset_->Dof(element, lm);
      if (not sourcedis_->HaveGlobalElement(element->Id())) return;
      DRT::Element* sourceele = sourcedis_->gElement(element->Id());
      return sourcedofset_->Dof(sourceele, lm);
    }

    /// Print this class
    virtual void Print(std::ostream& os) const { sourcedofset_->Print(os); };

    /// Print the dofsets in the static_dofsets_ list
    virtual void PrintAllDofsets(const Epetra_Comm& comm) const
    {
      sourcedofset_->PrintAllDofsets(comm);
    };

    /// Get Number of Global Elements of degree of freedom row map
    virtual int NumGlobalElements() const
    {
      CheckIsAssigned();
      return sourcedofset_->NumGlobalElements();
    };

    /// Get degree of freedom row map
    virtual const Epetra_Map* DofRowMap() const
    {
      CheckIsAssigned();
      return sourcedofset_->DofRowMap();
    };

    /// Get degree of freedom column map
    virtual const Epetra_Map* DofColMap() const
    {
      CheckIsAssigned();
      return sourcedofset_->DofColMap();
    };

    /// Get maximum GID of degree of freedom row map
    virtual int MaxAllGID() const
    {
      CheckIsAssigned();
      return sourcedofset_->MaxAllGID();
    };

    /// Get minimum GID of degree of freedom row map
    virtual int MinAllGID() const
    {
      CheckIsAssigned();
      return sourcedofset_->MinAllGID();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    virtual int MaxGIDinList(const Epetra_Comm& comm) const
    {
      CheckIsAssigned();
      return sourcedofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    virtual bool Initialized() const
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


#endif /* DRT_DOFSET_SUBSET_H_ */
