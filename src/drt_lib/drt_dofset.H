/*!----------------------------------------------------------------------
\file drt_dofset.H

\brief A set of degrees of freedom

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSET_H
#define DOFSET_H

#include <vector>
#include <list>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "drt_dofset_base.H"
#include "drt_dserror.H"
#include "drt_element.H"
#include "drt_node.H"

namespace DRT
{

class Discretization;
class DofSetProxy;

/*!
\brief A set of degrees of freedom

\note This is an internal class of the discretization module that one
should not need to touch on an ordinary day. It is here to support the
Discretization class. Everyday use should happen via the
Discretization class only.

<h3>Purpose</h3>

This class represents one set of degrees of freedom for the
Discretization class in the usual parallel fashion. That is there is a
DofRowMap() and a DofColMap() that return the maps of the global FE
system of equation in row and column setting respectively. These maps
are used by the algorithm's Epetra_Vector classes among others.

There can be dofs in nodes, faces, and elements. And additionally to the
above maps this class needs to know the global dof ids of all nodes
and elements. In order to provide this information nodal and elemental
column vectors are provided that store the number of dofs and the
local column map id of the first dof for each node and element. Since
dof numbers to one object are always consecutive this is all that's
needed. So the methods NumDof() and Dof() can be provided for nodes,
faces, and elements.

One has to keep in mind, however, that the lookup of one dof gid
involves some table lookups. Therefore there is a special version of
Dof() that gathers and returns all dof gids to one node or element at
once. This is to be preferred if more that one lookup is needed.

The point of holding these maps and vectors in a class of its own is
to enable multiple sets of dofs on the same mesh. But judging from
past experience this feature will not be used that often. So effort
has been made to hide the possibility of multiple DofSets.

The setup is done by AssignDegreesOfFreedom(). This method uses two
redundant nodal and elemental vectors. It would be hard to avoid
those. Lets hope we can always afford them.

\note It is guaranteed that the same mesh (nodes and elements are all
the same with the same global ids) is assigned the same set of dofs
all the time independent of its parallel distribution. That's crucial
to be able to redistribute a mesh without losing the old vectors.

<h3>Invariants</h3>

There are two possible states in this class: Reset and setup. To
change back and forth use AssignDegreesOfFreedom() and Reset().

<h3>Dof number uniqueness</h3>

Each DofSet assigns unique dof numbers that do not occur in any other
DofSet. This is true as long as the number of dofs per DofSet does not
change. To achieve this we keep a list of dof sets internally.

<h3>Copying behaviour</h3>

Please note that even though Michael does not like it this class
contains neither copy constructor nor assignment operator. This is
intended. It is legal to copy this objects of class. The internal
variables (all Teuchos::RCPs) know how to copy themselves. So the
default versions will do just fine. (Far better than buggy hand
written versions.) And due to the two possible states there is no
reason to deep copy any of the local map and vector variables.

\author u.kue
*/
class DofSet: public DofSetBase
{
  friend class DofSetProxy;
public:

  /*!
  \brief Standard Constructor

  */
  DofSet();

  /*!
  \brief Destructor

  */
  virtual ~DofSet();

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new DofSet(*this)); }

  //! @name Access methods

  /// Get number of dofs for given node
  virtual int NumDof(const Node* node) const
    {
      int lid = node->LID();
      if (lid==-1)
        return 0;
      return (*numdfcolnodes_)[lid];
    }

  /// Get number of dofs for given element
  virtual int NumDof(const Element* element) const
    {
      // check if this is a face element
      int lid = element->LID();
      if (lid==-1)
        return 0;
      if (element->IsFaceElement())
        return (*numdfcolfaces_)[lid];
      else
        return (*numdfcolelements_)[lid];
    }

  /// Get the gid of a dof for given node
  virtual int Dof(const Node* node, int dof) const
    {
      int lid = node->LID();
      if (lid==-1)
        return -1;
      if (pccdofhandling_)
        return dofscolnodes_->GID((*shiftcolnodes_)[lid]+dof);
      else
        return (*idxcolnodes_)[lid]+dof;
    }

  /// Get the gid of a dof for given element
  virtual int Dof(const Element* element, int dof) const
    {
      int lid = element->LID();
      if (lid==-1)
        return -1;
      if (element->IsFaceElement())
        return (*numdfcolfaces_)[lid]+dof;
      else
        return (*idxcolelements_)[lid]+dof;
    }

  /// Get the gid of all dofs of a node
  virtual std::vector<int> Dof(const Node* node) const
    {
      const int lid = node->LID();
      if (lid==-1)
        return std::vector<int>();
      const int idx = (*idxcolnodes_)[lid];
      std::vector<int> dof((*numdfcolnodes_)[lid]);
      for (unsigned i=0; i<dof.size(); ++i)
      {
        if(pccdofhandling_)
          dof[i] = dofscolnodes_->GID((*shiftcolnodes_)[lid]+i);
        else
          dof[i] = idx+i;
      }
      return dof;
    }

  /// Get the gid of all dofs of a node
  virtual void Dof( std::vector<int>& dof, ///< vector of dof gids (to be filled)
                    const Node* node,      ///< the node
                    unsigned nodaldofset   ///< number of nodal dof set of the node (currently !=0 only for XFEM)
                    ) const
    {
      dsassert(nodaldofset==0,"only one nodal dofset supported!");
      dof = Dof(node);
    }

  /// Get the gid of all dofs of a element
  virtual std::vector<int> Dof(const Element* element) const
    {
      int lid = element->LID();
      if (lid==-1)
        return std::vector<int>();
      int idx = element->IsFaceElement() ? (*idxcolfaces_)[lid] : (*idxcolelements_)[lid];
      std::vector<int> dof(element->IsFaceElement() ?
          (*numdfcolfaces_)[lid] : (*numdfcolelements_)[lid]);
      for (unsigned i=0; i<dof.size(); ++i)
        dof[i] = idx+i;
      return dof;
    }

  /// Get the gid of all dofs of a node
  virtual void Dof(const Node* node, std::vector<int>& lm) const
    {
      int lid = node->LID();
      if (lid==-1)
        return;
      int idx = (*idxcolnodes_)[lid];
      int size = (*numdfcolnodes_)[lid];
      for (int i=0; i<size; ++i)
      {
        if (pccdofhandling_)
          lm.push_back(dofscolnodes_->GID((*shiftcolnodes_)[lid]+i));
        else
          lm.push_back(idx+i);
      }
    }

  /// Get the gid of all dofs of a node
  virtual void Dof(
      const Node* node,           ///< node, for which you want the dof positions
      const unsigned startindex,  ///< first index of vector at which will be written to end
      std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
      ) const
  {
    const int lid = node->LID();
    if (lid==-1)
      return;
    const int idx = (*idxcolnodes_)[lid];
    const int size = (*numdfcolnodes_)[lid];
    dsassert(lm.size() >= (startindex+size), "vector<int> lm too small");
    for (int i=0; i<size; ++i)
    {
      if (pccdofhandling_)
        lm[startindex+i] = dofscolnodes_->GID((*shiftcolnodes_)[lid]+i);
      else
        lm[startindex+i] = idx+i;
    }
  }

  /// Get the gid of all dofs of a element
  virtual void Dof(const Element* element, std::vector<int>& lm) const
    {
      int lid = element->LID();
      if (lid==-1)
        return;
      int idx = element->IsFaceElement() ? (*idxcolfaces_)[lid] : (*idxcolelements_)[lid];
      int size = element->IsFaceElement() ? (*numdfcolfaces_)[lid] : (*numdfcolelements_)[lid];
      for (int i=0; i<size; ++i)
        lm.push_back(idx+i);
    }

  /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
  virtual void Dof(
      const Element* element,     ///< element which provides its expected number of DOFs per node
      const Node* node,           ///< node, for which you want the DOF positions
      std::vector<int>& lm        ///< already allocated vector to be filled with DOF positions
      ) const
  {
    const int lid = node->LID();
    if (lid==-1)
      return;
    const int idx = (*idxcolnodes_)[lid];
    // this method is used to setup the vector of number of dofs, so we cannot ask numdfcolelements_ here as above.
    // Instead we have to ask the node itself
    const int size = NumDofPerNode(*element,*node);
    for (int i=0; i<size; ++i)
    {
      if (pccdofhandling_)
        lm.push_back(dofscolnodes_->GID((*shiftcolnodes_)[lid]+i));
      else
        lm.push_back(idx+i);
    }
  }

  /*!
  \brief Print this class
  */
  virtual void Print(std::ostream& os) const ;

  virtual bool Filled() const { return filled_; }

  //@}

  //! @name Construction

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

  /// reset all internal variables
  virtual void Reset();

  //@}

private:

  //! @name Proxy management
  /// Proxies need to know about changes to the DofSet.

  /// Register new proxy to notify
  void RegisterProxy(DofSetProxy* proxy);

  /// Remove proxy from list
  void UnregisterProxy(DofSetProxy* proxy);

  //@}

protected:

  /// get number of nodal dofs
  virtual int NumDofPerNode( const Node & node) const
  {
    const int numele = node.NumElement();
    const DRT::Element* const * myele = node.Elements();
    int numdf=0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,NumDofPerNode(*myele[j],node));
    return numdf;
  }

  /// get number of nodal dofs for this element at this node
  virtual int NumDofPerNode( const Element & element, const Node & node) const
  {
    return element.NumDofPerNode( node);
  }

  /// get number of element dofs for this element
  virtual int NumDofPerElement( const Element & element) const
  {
    return element.NumDofPerElement( );
  }

  virtual int NumDofPerFace( const Element & element, int face) const
  {
    return element.NumDofPerFace(face);
  }

  //! @name Proxy management

  /// Notify proxies of new dofs
  void NotifyAssigned();

  /// Notify proxies of reset
  void NotifyReset();

  ///Get Reserved Max Number Dofs per Node
  virtual void GetReservedMaxNumDofperNode(int& maxnodenumdf){return;};

  //@}

  /// filled flag
  bool filled_;

  /// position of dofset inside its discretization
  unsigned dspos_;

  /// number of dofs for each node
  Teuchos::RCP<Epetra_IntVector>                 numdfcolnodes_;

  /// number of dofs for each face
  Teuchos::RCP<Epetra_IntVector>                 numdfcolfaces_;

  /// number of dofs for each element
  Teuchos::RCP<Epetra_IntVector>                 numdfcolelements_;

  /// column map gid of first dof for each node
  Teuchos::RCP<Epetra_IntVector>                 idxcolnodes_;

  /// column map gid of first dof for each face
  Teuchos::RCP<Epetra_IntVector>                 idxcolfaces_;

  /// column map gid of first dof for each element
  Teuchos::RCP<Epetra_IntVector>                 idxcolelements_;

  //***************************************************************************
  // Remark:                                                       popp 02/2016
  // --------------------------------------------------------------------------
  // In the logic outlined above and implemented in AssignDegreesOfFreedom() we
  // assume that each nodal dof only appears once and therefore that we can
  // define a unique access to all dofs of one node by knowing (1) its number of
  // nodal dofs and (2) its first dof gid. The second, third,... dof gid are
  // then simply defined by incrementation +1 starting from the first dof gid.
  // However, in the case where we have point coupling conditions, e.g. due to
  // joints in structural mechanics (where certain dofs of two nodes at the same
  // geometrical position are to be coupled and certain others not), this basic
  // assumption does not hold anymore. Instead, we want to assign the same(!)
  // dof to both nodes in this case, i.e. the same dof is then shared by more
  // than one node and the numbering logic mentioned above is not sufficient
  // anymore. For such cases in particular (but maybe also in general?), we
  // should store not only the first dof for each node, but all(!) dofs for
  // each node. The idxcolnodes_ vector defined above is then replaced by the
  // two following data containers:
  //    dofscolnodes_ (Epetra_Map) and shiftcolnodes_ (Epetra_IntVector)
  // Nevertheless, for the time being, the old version with idxcolnodes_ stays
  // around, since the new version is only implemented for the standard DofSet
  // class declared here (but not for special TransparentDofSets, Proxy, etc.).
  //***************************************************************************
  /// special dof handling due to point coupling conditions
  bool pccdofhandling_;

  /// column map of all dofs for each node (possibly non-unique)
  Teuchos::RCP<Epetra_Map>                       dofscolnodes_;

  /// shift value for access to column map of all dofs for each node
  Teuchos::RCP<Epetra_IntVector>                 shiftcolnodes_;
  //***************************************************************************

private:

  /// list of my proxies
  std::list<DofSetProxy*> proxies_;
}; // class DofSet
} // namespace DRT


// << operator
std::ostream& operator << (std::ostream& os, const DRT::DofSet& dofset);


#endif  // #ifndef DOFSET_H
