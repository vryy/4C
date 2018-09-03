/*!----------------------------------------------------------------------
\file drt_dofset_interface.H

\brief Common interface class for all sets of degrees of freedom.

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSET_INTERFACE_H
#define DOFSET_INTERFACE_H

#include <vector>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Node;
  class Element;
  class Discretization;

  /*! \brief Common interface class for all sets of degrees of freedom.
   *
   * This is a pure virtual class all classes managing sets of degrees of freedom
   * should inherit from.
   *
   * \date 10/2016
   * \author Andreas Rauch    */
  class DofSetInterface
  {
   public:
    //! @name Construction

    /// Standard Constructor
    DofSetInterface(){};

    /// Destructor
    virtual ~DofSetInterface(){};

    //@}


    //! @name Public Access Methods

    /// Get number of dofs for given node
    virtual int NumDof(const Node* node  ///< node, for which you want to know the number of dofs
        ) const = 0;

    /// Get number of dofs for given element
    virtual int NumDof(
        const Element* element  ///< element, for which you want to know the number of dofs
        ) const = 0;

    /// get number of nodal dofs
    virtual int NumDofPerNode(
        const Node& node  ///< node, for which you want to know the number of dofs
        ) const = 0;

    /// Get the gid of a dof for given node
    virtual int Dof(const Node* node,  ///< node, for which you want the dof positions
        int dof                        ///< number of dof for which you want the dof position
        ) const = 0;

    /// Get the gid of a dof for given element
    virtual int Dof(const Element* element,  ///< element, for which you want the dof positions
        int dof                              ///< number dof for which you want the dof position
        ) const = 0;

    /// Get the gid of all dofs of a node
    virtual std::vector<int> Dof(const Node* node  ///< node, for which you want the dof positions
        ) const = 0;

    /// Get the gid of all dofs of a node
    virtual void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const Node* node,                    ///< node, for which you want the dof positions
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
        ) const = 0;

    /// Get the gid of all dofs of a element
    virtual std::vector<int> Dof(const Element* element) const = 0;

    /// Get the gid of all dofs of a node and the location matrix
    virtual void Dof(const Node* node, std::vector<int>& lm) const = 0;

    /// Get the gid of all dofs of a node
    virtual void Dof(const Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,      ///< first index of vector at which will be written to end
        std::vector<int>& lm  ///< already allocated vector to be filled with dof positions
        ) const = 0;

    /// Get the gid of all dofs of a element and the location matrix
    virtual void Dof(const Element* element, std::vector<int>& lm) const = 0;

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    virtual void Dof(
        const Element* element,  ///< element which provides its expected number of DOFs per node
        const Node* node,        ///< node, for which you want the DOF positions
        std::vector<int>& lm     ///< already allocated vector to be filled with DOF positions
        ) const = 0;

    //@}


    //! @name Utility Methods

    /// Print this class
    virtual void Print(std::ostream& os) const = 0;

    /// Print the dofsets in the static_dofsets_ list
    virtual void PrintAllDofsets(const Epetra_Comm& comm) const = 0;

    /// Returns true if filled
    virtual bool Filled() const = 0;

    /// Add Dof Set to list #static_dofsets_
    virtual void AddDofSettoList() = 0;

    /// Replace a Dof Set in list #static_dofsets_ with this
    virtual void ReplaceInStaticDofsets(Teuchos::RCP<DofSetInterface> olddofset) = 0;

    /// Get Number of Global Elements of degree of freedom row map
    virtual int NumGlobalElements() const = 0;

    /// Get degree of freedom row map
    virtual const Epetra_Map* DofRowMap() const = 0;

    /// Get degree of freedom column map
    virtual const Epetra_Map* DofColMap() const = 0;

    /// Get maximum GID of degree of freedom row map
    virtual int MaxAllGID() const = 0;

    /// Get minimum GID of degree of freedom row map
    virtual int MinAllGID() const = 0;

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    virtual int MaxGIDinList(const Epetra_Comm& comm) const = 0;

    /// are the dof maps already initialized?
    virtual bool Initialized() const = 0;

    //@}



    //! @name Setup and Initialization

    /// Assign dof numbers using all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
        const Discretization& dis, const unsigned dspos, const int start) = 0;

    /// reset all internal variables
    virtual void Reset() = 0;

    //@}


    //! @name Proxy management
    /// Proxies need to know about changes to the DofSet.

    /// Notify proxies of new dofs
    virtual void NotifyAssigned() = 0;

    /// Notify proxies of reset
    virtual void NotifyReset() = 0;

    /// Register new proxy to notify
    virtual void Register(DofSetInterface* dofset) = 0;

    /// Remove proxy from list
    virtual void Unregister(DofSetInterface* dofset) = 0;

    /// our original DofSet dies
    virtual void Disconnect(DofSetInterface* dofset) = 0;

    //@}

  };  // class DofSetInterface

}  // namespace DRT


#endif  // #ifndef DOFSET_INTERFACE_H
