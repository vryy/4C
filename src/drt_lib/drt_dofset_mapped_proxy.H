/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_mapped_proxy.H

 \brief A proxy of a dofset that does not rely on same GID/LID numbers but uses
         a defined node mapping instead (not implemented for element DOFs).

   \level 3

   \maintainer Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_DOFSET_MAPPED_PROXY_H_
#define SRC_DRT_LIB_DRT_DOFSET_MAPPED_PROXY_H_

#include "drt_dofset_proxy.H"

#include "drt_discret.H"

namespace DRT
{
  class DofSetMappedProxy : public DofSetProxy
  {
  public:

    /*!
    \brief Standard Constructor

    */
    explicit DofSetMappedProxy(Teuchos::RCP<DofSet>  dofset,
                               const Teuchos::RCP<const DRT::Discretization> sourcedis,
                               const std::string& couplingcond);

    /*!
    \brief Destructor

    */
     ~DofSetMappedProxy(){};

     /// Assign dof numbers using all elements and nodes of the discretization.
     virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

     /// original DofSet has been reset
     //virtual void NotifyReset();

     //! @name Access methods

     /// Get the gid of all dofs of a node
     virtual std::vector<int> Dof(const Node* node) const
       {
         const DRT::Node* sourcenode = GetSourceNode(node->LID());
         return DofSetProxy::Dof(sourcenode);
       }

     /// Get the gid of all dofs of a element
     virtual std::vector<int> Dof(const Element* element) const
       {
        return std::vector<int>();
         //return DofSetProxy::Dof(element);
       }

     /// Get number of dofs for given node
     virtual int NumDof(const Node* node) const
       {
       const DRT::Node* sourcenode = GetSourceNode(node->LID());

         return DofSetProxy::NumDof(sourcenode);
       }

     /// Get number of dofs for given element
     virtual int NumDof(const Element* element) const
       {
         // element dofs not yet supported
         return 0;
         //return DofSetProxy::NumDof(element);
       }

     /// Get the gid of a dof for given node
     virtual int Dof(const Node* node, int dof) const
       {
           const DRT::Node* sourcenode = GetSourceNode(node->LID());
          return DofSetProxy::Dof(sourcenode,dof);
       }

     /// Get the gid of a dof for given element
     virtual int Dof(const Element* element, int dof) const
       {
       return 0;
         return DofSetProxy::Dof(element,dof);
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(const Node* node, std::vector<int>& lm) const
       {
         const DRT::Node* sourcenode = GetSourceNode(node->LID());
         DofSetProxy::Dof(sourcenode,lm);

         return;
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(
         const Node* node,           ///< node, for which you want the dof positions
         const unsigned startindex,  ///< first index of vector at which will be written to end
         std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
         ) const
     {
       const DRT::Node* sourcenode = GetSourceNode(node->LID());
       DofSetProxy::Dof(sourcenode,startindex,lm);
     }

     /// Get the gid of all dofs of a element
     virtual void Dof(const Element* element, std::vector<int>& lm) const
     {
       // element dofs not yet supported
       return;
       //DofSetProxy::Dof(element,lm);
     }

     /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
     virtual void Dof(
         const Element* element,     ///< element which provides its expected number of DOFs per node
         const Node* node,           ///< node, for which you want the DOF positions
         std::vector<int>& lm        ///< already allocated vector to be filled with DOF positions
         ) const
     {
       const DRT::Node* sourcenode = GetSourceNode(node->LID());
       DofSetProxy::Dof(element,sourcenode,lm);
     }

  private:

     //! get corresponding source node from source discretization
     const DRT::Node* GetSourceNode(int targetLid) const
     {
       if(targetLid>sourcetotargetnodemapping_->MyLength())
         dserror("Target Lid %i not found in discretization!",targetLid);
       int sourcegid = (*sourcetotargetnodemapping_)[targetLid];

       if(not sourcedis_->HaveGlobalNode(sourcegid))
         dserror("Master Gid %i not found in source discretization!",sourcegid);
       return sourcedis_->gNode(sourcegid);
     }

    //! map containing the mapping of the target node GID to the corresponding source node GID (value)
     Teuchos::RCP<Epetra_IntVector>  sourcetotargetnodemapping_;

    //! underlying actual dofset
    Teuchos::RCP<DofSet> dofset_;

    //! source discretization
    const Teuchos::RCP<const DRT::Discretization> sourcedis_;

    //! condition string defining the coupling
    const std::string couplingcond_;
  };
}



#endif /* SRC_DRT_LIB_DRT_DOFSET_MAPPED_PROXY_H_ */
