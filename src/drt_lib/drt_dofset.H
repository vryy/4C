/*!----------------------------------------------------------------------
\file drt_dofset.H

\brief A set of degrees of freedom

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
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
are used by the algorithm's Epetra_Vector classes amoung others.

There can be dofs in nodes and in elements. And additionally to the
above maps this class needs to know the global dof ids of all nodes
and elements. In order to provide this information nodal and elemental
column vectors are provided that store the number of dofs and the
local column map id of the first dof for each node and element. Since
dof numbers to one object are always consecutive this is all that's
needed. So the methods NumDof() and Dof() can be provided for nodes
and elements.

One has to keep in mind, however, that the lookup of one dof gid
involves some table lookups. Therefore there is a special version of
Dof() that gathers and returns all dof gids to one node or element at
once. This is to be prefered if more that one lookup is needed.

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
contains neither copy constructor nor assignment operatior. This is
intended. It is legal to copy this objects of class. The internal
variables (all RCPs) know how to copy themselves. So the
default versions will do just fine. (Far better that buggy hand
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
      int lid = element->LID();
      if (lid==-1)
        return 0;
      return (*numdfcolelements_)[lid];
    }

  /// Get the gid of a dof for given node
  virtual int Dof(const Node* node, int dof) const
    {
      int lid = node->LID();
      if (lid==-1)
        return -1;
      return (*idxcolnodes_)[lid]+dof;
    }

  /// Get the gid of a dof for given element
  virtual int Dof(const Element* element, int dof) const
    {
      int lid = element->LID();
      if (lid==-1)
        return -1;
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
        dof[i] = idx+i;
      return dof;
    }

  /// Get the gid of all dofs of a element
  virtual std::vector<int> Dof(const Element* element) const
    {
      int lid = element->LID();
      if (lid==-1)
        return std::vector<int>();
      int idx = (*idxcolelements_)[lid];
      std::vector<int> dof((*numdfcolelements_)[lid]);
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
        lm.push_back(idx+i);
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
      lm[startindex+i] = idx+i;
  }

  /// Get the gid of all dofs of a element
  virtual void Dof(const Element* element, std::vector<int>& lm) const
    {
      int lid = element->LID();
      if (lid==-1)
        return;
      int idx = (*idxcolelements_)[lid];
      int size = (*numdfcolelements_)[lid];
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
    const int size = NumDofPerNode(*element,*node,0);
    for (int i=0; i<size; ++i)
      lm.push_back(idx+i);
  }

  /*!
  \brief Print this class
  */
  virtual void Print(ostream& os) const ;

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
  virtual int NumDofPerNode( const Node & node, unsigned dspos ) const
  {
    const int numele = node.NumElement();
    const DRT::Element* const * myele = node.Elements();
    int numdf=0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,NumDofPerNode(*myele[j],node,dspos));
    return numdf;
  }

  /// get number of nodal dofs for this element at this node
  virtual int NumDofPerNode( const Element & element, const Node & node, unsigned dspos ) const
  {
    return element.NumDofPerNode( dspos, node );
  }

  /// get number of element dofs for this element
  virtual int NumDofPerElement( const Element & element, unsigned dspos ) const
  {
    return element.NumDofPerElement( dspos );
  }

  //! @name Proxy management

  /// Notify proxies of new dofs
  void NotifyAssigned();

  /// Notify proxies of reset
  void NotifyReset();

  //@}

  /// filled flag
  bool filled_;

  /// position of dofset inside its discretization
  unsigned dspos_;

  /// number of dofs for each node
  Teuchos::RCP<Epetra_IntVector>                 numdfcolnodes_;

  /// number of dofs for each element
  Teuchos::RCP<Epetra_IntVector>                 numdfcolelements_;

  /// column map gid of first dof for each node
  Teuchos::RCP<Epetra_IntVector>                 idxcolnodes_;

  /// column map gid of first dof for each element
  Teuchos::RCP<Epetra_IntVector>                 idxcolelements_;

private:

  /// list of my proxies
  std::list<DofSetProxy*> proxies_;
}; // class DofSet
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::DofSet& dofset);


#endif  // #ifndef DOFSET_H
