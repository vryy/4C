/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_DOFSET_PROXY_HPP
#define FOUR_C_DISCRETIZATION_DOFSET_PROXY_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CORE::Dofsets
{
  /*! \brief Proxy to a DofSet that does not own dofs itself

    As the Discretization handles DofSets a little implicit, a DofSetProxy is
    needed to change the DofSet behaviour. We need a DofSet that

    - returns dof numbers owned by a different DofSet
    - does not reset those dof numbers
    - does not assign dof numbers itself

    The DofSetProxy is meant to be used as secondary DofSet in a Discretization
    if there are two fully volumetric coupled Discretizations with
    matching nodes. Think of Structure-Thermo coupling.

    \author u.kue
    \author Andreas Rauch
    \date 10/16    */
  class DofSetProxy : public DofSetBase
  {
   public:
    //! @name Construction

    /// Constructor
    explicit DofSetProxy(DofSetInterface* dofset);

    /// Destructor
    ~DofSetProxy() override;

    /// create a copy of this object
    virtual Teuchos::RCP<DofSetProxy> Clone() { return Teuchos::rcp(new DofSetProxy(*this)); }

    /// Add Dof Set to list #static_dofsets_
    void AddDofSettoList() override;

    /// Replace a Dof Set in list #static_dofsets_ with this
    void replace_in_static_dofsets(Teuchos::RCP<DofSetInterface> olddofset) override
    {
      dofset_->replace_in_static_dofsets(olddofset);
    };

    /// Assign dof numbers using all elements and nodes of the discretization.
    int assign_degrees_of_freedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start) override;

    /// returns true if \ref dofset_ is filled
    bool Filled() const override;

    /// reset all internal variables
    void Reset() override;

    //@}

    //! @name Communication
    /// The original DofSet sends notifications if it changes.

    /// original DofSet has new dofs
    void NotifyAssigned() override;

    /// our original DofSet dies
    void Disconnect(DofSetInterface* dofset) override;

    //@}


    //! @name Access methods

    /// Get number of dofs for given node
    int NumDof(const DRT::Node* node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      CheckIsAssigned();
      return dofset_->NumDof(node);
    };

    /// Get number of dofs for given element
    int NumDof(
        const DRT::Element* element  ///< element, for which you want to know the number of dofs
    ) const override
    {
      CheckIsAssigned();
      return dofset_->NumDof(element);
    };

    /// get number of nodal dofs
    int NumDofPerNode(
        const DRT::Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      CheckIsAssigned();
      return dofset_->NumDofPerNode(node);
    };

    /// Get the gid of a dof for given node
    int Dof(const DRT::Node* node, int dof) const override
    {
      CheckIsAssigned();
      return dofset_->Dof(node, dof);
    };

    /// Get the gid of a dof for given element
    int Dof(const DRT::Element* element,  ///< element, for which you want the dof positions
        int dof) const override
    {
      CheckIsAssigned();
      return dofset_->Dof(element, dof);
    };

    /// Get the gid of all dofs of a node
    std::vector<int> Dof(const DRT::Node* node  ///< node, for which you want the dof positions
    ) const override
    {
      CheckIsAssigned();
      return dofset_->Dof(node);
    };

    /// Get the gid of all dofs of a node
    void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const DRT::Node* node,       ///< node, for which you want the dof positions
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      CheckIsAssigned();
      dofset_->Dof(dof, node, nodaldofset);
    };

    /// Get the gid of all dofs of a element
    std::vector<int> Dof(const DRT::Element* element) const override
    {
      CheckIsAssigned();
      return dofset_->Dof(element);
    };

    /// Get the gid of all dofs of a node and the location matrix
    void Dof(const DRT::Node* node, std::vector<int>& lm) const override
    {
      CheckIsAssigned();
      dofset_->Dof(node, lm);
    };

    /// Get the gid of all dofs of a node
    void Dof(const DRT::Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,   ///< first index of vector at which will be written to end
        std::vector<int>& lm         ///< already allocated vector to be filled with dof positions
    ) const override
    {
      CheckIsAssigned();
      dofset_->Dof(node, startindex, lm);
    };

    /// Get the gid of all dofs of a element and the location matrix
    void Dof(const DRT::Element* element, std::vector<int>& lm) const override
    {
      CheckIsAssigned();
      dofset_->Dof(element, lm);
    };

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void Dof(const DRT::Element*
                 element,       ///< element which provides its expected number of DOFs per node
        const DRT::Node* node,  ///< node, for which you want the DOF positions
        std::vector<int>& lm    ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      CheckIsAssigned();
      dofset_->Dof(element, node, lm);
    };

    /// Print this class
    void Print(std::ostream& os) const override { dofset_->Print(os); };

    /// Print the dofsets in the static_dofsets_ list
    void PrintAllDofsets(const Epetra_Comm& comm) const override
    {
      dofset_->PrintAllDofsets(comm);
    };

    /// Get Number of Global Elements of degree of freedom row map
    int NumGlobalElements() const override
    {
      CheckIsAssigned();
      return dofset_->NumGlobalElements();
    };

    /// Get degree of freedom row map
    const Epetra_Map* DofRowMap() const override
    {
      CheckIsAssigned();
      return dofset_->DofRowMap();
    };

    /// Get degree of freedom column map
    const Epetra_Map* DofColMap() const override
    {
      CheckIsAssigned();
      return dofset_->DofColMap();
    };

    /// Get maximum GID of degree of freedom row map
    int MaxAllGID() const override
    {
      CheckIsAssigned();
      return dofset_->MaxAllGID();
    };

    /// Get minimum GID of degree of freedom row map
    int MinAllGID() const override
    {
      CheckIsAssigned();
      return dofset_->MinAllGID();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int MaxGIDinList(const Epetra_Comm& comm) const override
    {
      CheckIsAssigned();
      return dofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    bool Initialized() const override
    {
      CheckIsAssigned();
      return dofset_->Initialized();
    };

   protected:
    /// check if \ref assign_degrees_of_freedom was called on parent dofset
    void CheckIsAssigned() const;

    /// pointer to the parent dofset represented by this proxy
    DofSetInterface* dofset_;

   private:
    /// assigned flag
    bool isassigned_;

  };  // class DofSetProxy

}  // namespace CORE::Dofsets


FOUR_C_NAMESPACE_CLOSE

#endif
