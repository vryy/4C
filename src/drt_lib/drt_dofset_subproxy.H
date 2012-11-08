/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_subproxy.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef DRT_DOFSET_SUBSET_H_
#define DRT_DOFSET_SUBSET_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Map.h>

#include "drt_dofset_proxy.H"

namespace DRT
{
  class DofSetSubProxy : public DofSetProxy
  {
  public:

    /*!
    \brief Standard Constructor

    */
    explicit DofSetSubProxy(DofSet* dofset,
                          Teuchos::RCP<const Epetra_Map> subcolnodes,
                          Teuchos::RCP<const Epetra_Map> subcoleles );

    /*!
    \brief Destructor

    */
     ~DofSetSubProxy(){};

     /// Get the gid of all dofs of a node
     virtual vector<int> Dof(const Node* node) const
       {
         const int gid = node->Id();
         const int lid = subcolnodes_->LID(gid);
         if (lid==-1)
           return vector<int>();
         const int idx = (*idxcolnodes_)[lid];
         vector<int> dof((*numdfcolnodes_)[lid]);
         for (unsigned i=0; i<dof.size(); ++i)
           dof[i] = idx+i;
         return dof;
       }

     /// Get the gid of all dofs of a element
     virtual vector<int> Dof(const Element* element) const
       {
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
         if (lid==-1)
           return vector<int>();
         int idx = (*idxcolelements_)[lid];
         vector<int> dof((*numdfcolelements_)[lid]);
         for (unsigned i=0; i<dof.size(); ++i)
           dof[i] = idx+i;
         return dof;
       }

  private:

    Teuchos::RCP<const Epetra_Map>                 subcolnodes_;

    Teuchos::RCP<const Epetra_Map>                 subcoleles_;
  };
}


#endif /* DRT_DOFSET_SUBSET_H_ */
