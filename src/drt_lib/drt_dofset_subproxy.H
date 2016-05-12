/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_subproxy.H

 \brief subproxy functionality to dofsets

 <pre>
   \level 1

   \maintainer Anh-Tu Vuong
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef DRT_DOFSET_SUBSET_H_
#define DRT_DOFSET_SUBSET_H_

#include <Teuchos_RCP.hpp>
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
                            Teuchos::RCP<Epetra_Map>& subcolnodes,
                            Teuchos::RCP<Epetra_Map>& subcolele);

    /*!
    \brief Destructor

    */
     ~DofSetSubProxy(){};

     /// original DofSet has been reset
     void NotifyReset();

     //! @name Access methods

     /// Get the gid of all dofs of a node
     virtual std::vector<int> Dof(const Node* node) const
       {
         const int gid = node->Id();
         const int lid = subcolnodes_->LID(gid);
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
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
         if (lid==-1)
           return std::vector<int>();
         int idx = (*idxcolelements_)[lid];
         std::vector<int> dof((*numdfcolelements_)[lid]);
         for (unsigned i=0; i<dof.size(); ++i)
           dof[i] = idx+i;
         return dof;
       }

     /// Get number of dofs for given node
     virtual int NumDof(const Node* node) const
       {
         const int gid = node->Id();
         const int lid = subcolnodes_->LID(gid);
         if (lid==-1)
           return 0;
         return (*numdfcolnodes_)[lid];
       }

     /// Get number of dofs for given element
     virtual int NumDof(const Element* element) const
       {
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
         if (lid==-1)
           return 0;
         return (*numdfcolelements_)[lid];
       }

     /// Get the gid of a dof for given node
     virtual int Dof(const Node* node, int dof) const
       {
         const int gid = node->Id();
         const int lid = subcolnodes_->LID(gid);
         if (lid==-1)
           return -1;
         return (*idxcolnodes_)[lid]+dof;
       }

     /// Get the gid of a dof for given element
     virtual int Dof(const Element* element, int dof) const
       {
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
         if (lid==-1)
           return -1;
         return (*idxcolelements_)[lid]+dof;
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(const Node* node, std::vector<int>& lm) const
       {
         const int gid = node->Id();
         const int lid = subcolnodes_->LID(gid);
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
       const int gid = node->Id();
       const int lid = subcolnodes_->LID(gid);
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
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
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
         const int gid = element->Id();
         const int lid = subcoleles_->LID(gid);
       if (lid==-1)
         return;
       const int idx = (*idxcolnodes_)[lid];
       const int size = NumDofPerNode(*element,*node);
       for (int i=0; i<size; ++i)
         lm.push_back(idx+i);
     }

  private:

    //! map containing all node GIDs of this (sub)dofset
     Teuchos::RCP<Epetra_Map>                 subcolnodes_;

    //! map containing all element GIDs of this (sub)dofset
     Teuchos::RCP<Epetra_Map>                 subcoleles_;

  };
}


#endif /* DRT_DOFSET_SUBSET_H_ */
