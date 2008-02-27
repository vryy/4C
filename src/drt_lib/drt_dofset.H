/*!----------------------------------------------------------------------
\file drt_dofset.H

\brief A set of degrees of freedom

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DOFSET_H
#define DOFSET_H

#include <vector>
#include <list>

#include <Epetra_IntVector.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

using namespace std;

#include "drt_dserror.H"
#include "drt_element.H"
#include "drt_node.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

class Discretization;


/*!
\class DRT::DofSet

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
variables (all RefCountPtrs) know how to copy themselves. So the
default versions will do just fine. (Far better that buggy hand
written versions.) And due to the two possible states there is no
reason to deep copy any of the local map and vector variables.

\author u.kue
*/
class DofSet
{
public:

  /*!
  \brief Standard Constructor

  */
  DofSet();

  /*!
  \brief Destructor

  */
  virtual ~DofSet();


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
      return dofcolmap_->GID((*idxcolnodes_)[lid]+dof);
    }

  /// Get the gid of a dof for given element
  virtual int Dof(const Element* element, int dof) const
    {
      int lid = element->LID();
      if (lid==-1)
        return -1;
      return dofcolmap_->GID((*idxcolelements_)[lid]+dof);
    }

  /// Get the gid of all dofs of a node
  virtual vector<int> Dof(const Node* node) const
    {
      int lid = node->LID();
      if (lid==-1)
        return vector<int>();
      int idx = (*idxcolnodes_)[lid];
      vector<int> dof((*numdfcolnodes_)[lid]);
      for (unsigned i=0; i<dof.size(); ++i)
        dof[i] = dofcolmap_->GID(idx+i);
      return dof;
    }

  /// Get the gid of all dofs of a element
  virtual vector<int> Dof(const Element* element) const
    {
      int lid = element->LID();
      if (lid==-1)
        return vector<int>();
      int idx = (*idxcolelements_)[lid];
      vector<int> dof((*numdfcolelements_)[lid]);
      for (unsigned i=0; i<dof.size(); ++i)
        dof[i] = dofcolmap_->GID(idx+i);
      return dof;
    }

  /// Get degree of freedom row map
  virtual const Epetra_Map* DofRowMap() const { return dofrowmap_.get(); }

  /// Get degree of freedom column map
  virtual const Epetra_Map* DofColMap() const { return dofcolmap_.get(); }

  /*!
  \brief Print this class
  */
  virtual void Print(ostream& os) const ;

  //@}

  //! @name Construction

  /// Assign dof numbers to all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const int start);

  /// reset all internal variables
  virtual void Reset();

  //@}

protected:

  /// unique row map of degrees of freedom (node and element dofs))
  RefCountPtr<Epetra_Map>                       dofrowmap_;

  /// unique column map of degrees of freedom (node and element dofs))
  RefCountPtr<Epetra_Map>                       dofcolmap_;

  /// number of dofs for each node
  RefCountPtr<Epetra_IntVector>                 numdfcolnodes_;

  /// number of dofs for each element
  RefCountPtr<Epetra_IntVector>                 numdfcolelements_;

  /// column map lid of first dof for each node
  RefCountPtr<Epetra_IntVector>                 idxcolnodes_;

  /// column map lid of first dof for each element
  RefCountPtr<Epetra_IntVector>                 idxcolelements_;

  /// store all dofset ever created
  /*!
    This is hack to get unique dof numbers on all dof sets. With this in place
    we can combine the maps from any dof sets to form block systems and the
    like.
   */
  static std::list<DofSet*>                     dofsets_;

}; // class DofSet
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::DofSet& dofset);


#endif  // #ifndef DOFSET_H
#endif  // #ifdef CCADISCRET
