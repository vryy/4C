/*!----------------------------------------------------------------------
\file drt_xdofset.H

\brief dofset to manage the varying number and the varying meaning of 
       XFEM degrees of freedom

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef D_XFEM
#ifndef XDOFSET_H
#define XDOFSET_H

#include <vector>

#include <Epetra_IntVector.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

using namespace std;

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

class Discretization;


/*!
\class DRT::PBCDofSet

\brief A set of degrees of freedom

\author gammi
*/
class XDofSet: public DofSet
{
public:



  /*!
  \brief Standard Constructor


  <pre>

  create a dofset that allows coupled nodes for periodic boundary
  conditions                                         gammi 05/07
  
  </pre>

  \param couplednodes (i) list of coupled nodes


  \return void
  
  */
  XDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes);

  /*!
  \brief Destructor

  */
  virtual ~XDofSet();


  /// Get number of dofs for given node
  virtual int NumDof(Node* node) const
    {
      int lid = node->LID();
      if (lid==-1)
        return 0;
      return (*numdfcolnodes_)[lid];
    }
  
  /// Get number of dofs for given element
  virtual int NumDof(Element* element) const
    {
      int lid = element->LID();
      if (lid==-1)
        return 0;
      return (*numdfcolelements_)[lid];
    }

  /// Get the gid of a dof for given node
  virtual int Dof(Node* node, int dof) const
    {
      int lid = node->LID();
      if (lid==-1)
        return -1;
      return dofcolmap_->GID((*idxcolnodes_)[lid]+dof);
    }
  /// Get the gid of a dof for given element
  virtual int Dof(Element* element, int dof) const
    {
      int lid = element->LID();
      if (lid==-1)
        return -1;
      return dofcolmap_->GID((*idxcolelements_)[lid]+dof);
    }

  /// Get the gid of all dofs of a node
  virtual vector<int> Dof(Node* node) const
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
  virtual vector<int> Dof(Element* element) const
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
  
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const int start);

  /// reset all internal variables
  virtual void Reset();

protected:

  //!\brief master and slave node connectivity for periodic boundary conditions
  RefCountPtr<map<int,vector<int> > >           perbndcouples_;
}; // class XDofSet
} // namespace DRT



#endif  // #ifndef XDOFSET_H
#endif  // #ifdef D_XFEM
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
