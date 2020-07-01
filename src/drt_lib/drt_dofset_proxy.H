/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef DRT_DOFSET_PROXY_H
#define DRT_DOFSET_PROXY_H

#include "drt_dofset_base.H"


namespace DRT
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
    virtual ~DofSetProxy();

    /// create a copy of this object
    virtual Teuchos::RCP<DofSetProxy> Clone() { return Teuchos::rcp(new DofSetProxy(*this)); }

    /// Add Dof Set to list #static_dofsets_
    virtual void AddDofSettoList();

    /// Replace a Dof Set in list #static_dofsets_ with this
    virtual void ReplaceInStaticDofsets(Teuchos::RCP<DofSetInterface> olddofset)
    {
      dofset_->ReplaceInStaticDofsets(olddofset);
    };

    /// Assign dof numbers using all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
        const Discretization& dis, const unsigned dspos, const int start);

    /// returns true if \ref dofset_ is filled
    virtual bool Filled() const;

    /// reset all internal variables
    virtual void Reset();

    //@}

    //! @name Communication
    /// The original DofSet sends notifications if it changes.

    /// original DofSet has new dofs
    virtual void NotifyAssigned();

    /// our original DofSet dies
    virtual void Disconnect(DofSetInterface* dofset);

    //@}


    //! @name Access methods

    /// Get number of dofs for given node
    virtual int NumDof(const Node* node  ///< node, for which you want to know the number of dofs
        ) const
    {
      CheckIsAssigned();
      return dofset_->NumDof(node);
    };

    /// Get number of dofs for given element
    virtual int NumDof(
        const Element* element  ///< element, for which you want to know the number of dofs
        ) const
    {
      CheckIsAssigned();
      return dofset_->NumDof(element);
    };

    /// get number of nodal dofs
    virtual int NumDofPerNode(
        const Node& node  ///< node, for which you want to know the number of dofs
        ) const
    {
      CheckIsAssigned();
      return dofset_->NumDofPerNode(node);
    };

    /// Get the gid of a dof for given node
    virtual int Dof(const Node* node, int dof) const
    {
      CheckIsAssigned();
      return dofset_->Dof(node, dof);
    };

    /// Get the gid of a dof for given element
    virtual int Dof(const Element* element,  ///< element, for which you want the dof positions
        int dof) const
    {
      CheckIsAssigned();
      return dofset_->Dof(element, dof);
    };

    /// Get the gid of all dofs of a node
    virtual std::vector<int> Dof(const Node* node  ///< node, for which you want the dof positions
        ) const
    {
      CheckIsAssigned();
      return dofset_->Dof(node);
    };

    /// Get the gid of all dofs of a node
    virtual void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const Node* node,                    ///< node, for which you want the dof positions
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
        ) const
    {
      CheckIsAssigned();
      dofset_->Dof(dof, node, nodaldofset);
    };

    /// Get the gid of all dofs of a element
    virtual std::vector<int> Dof(const Element* element) const
    {
      CheckIsAssigned();
      return dofset_->Dof(element);
    };

    /// Get the gid of all dofs of a node and the location matrix
    virtual void Dof(const Node* node, std::vector<int>& lm) const
    {
      CheckIsAssigned();
      dofset_->Dof(node, lm);
    };

    /// Get the gid of all dofs of a node
    virtual void Dof(const Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,      ///< first index of vector at which will be written to end
        std::vector<int>& lm  ///< already allocated vector to be filled with dof positions
        ) const
    {
      CheckIsAssigned();
      dofset_->Dof(node, startindex, lm);
    };

    /// Get the gid of all dofs of a element and the location matrix
    virtual void Dof(const Element* element, std::vector<int>& lm) const
    {
      CheckIsAssigned();
      dofset_->Dof(element, lm);
    };

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    virtual void Dof(
        const Element* element,  ///< element which provides its expected number of DOFs per node
        const Node* node,        ///< node, for which you want the DOF positions
        std::vector<int>& lm     ///< already allocated vector to be filled with DOF positions
        ) const
    {
      CheckIsAssigned();
      dofset_->Dof(element, node, lm);
    };

    /// Print this class
    virtual void Print(std::ostream& os) const { dofset_->Print(os); };

    /// Print the dofsets in the static_dofsets_ list
    virtual void PrintAllDofsets(const Epetra_Comm& comm) const { dofset_->PrintAllDofsets(comm); };

    /// Get Number of Global Elements of degree of freedom row map
    virtual int NumGlobalElements() const
    {
      CheckIsAssigned();
      return dofset_->NumGlobalElements();
    };

    /// Get degree of freedom row map
    virtual const Epetra_Map* DofRowMap() const
    {
      CheckIsAssigned();
      return dofset_->DofRowMap();
    };

    /// Get degree of freedom column map
    virtual const Epetra_Map* DofColMap() const
    {
      CheckIsAssigned();
      return dofset_->DofColMap();
    };

    /// Get maximum GID of degree of freedom row map
    virtual int MaxAllGID() const
    {
      CheckIsAssigned();
      return dofset_->MaxAllGID();
    };

    /// Get minimum GID of degree of freedom row map
    virtual int MinAllGID() const
    {
      CheckIsAssigned();
      return dofset_->MinAllGID();
    };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    virtual int MaxGIDinList(const Epetra_Comm& comm) const
    {
      CheckIsAssigned();
      return dofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    virtual bool Initialized() const
    {
      CheckIsAssigned();
      return dofset_->Initialized();
    };

   protected:
    /// check if \ref AssignDegreesOfFreedom was called on parent dofset
    void CheckIsAssigned() const;

    /// pointer to the parent dofset represented by this proxy
    DofSetInterface* dofset_;

   private:
    /// assigned flag
    bool isassigned_;

  };  // class DofSetProxy

}  // namespace DRT


#endif
