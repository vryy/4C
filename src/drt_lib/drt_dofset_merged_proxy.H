/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_merged_proxy.H

 \brief A proxy of a dofset that adds additional, existing degrees of freedom from the same
        discretization to nodes (not yet to elements).

   \level 3

   \maintainer Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_DOFSET_MERGED_PROXY_H_
#define SRC_DRT_LIB_DRT_DOFSET_MERGED_PROXY_H_


#include "drt_dofset_proxy.H"

#include "drt_discret.H"

namespace DRT
{
  class DofSetMergedProxy : public DofSetProxy
  {
  public:

    /*!
    \brief Standard Constructor

    */
    explicit DofSetMergedProxy(Teuchos::RCP<DofSet>  dofset,
                               const Teuchos::RCP<const DRT::Discretization> sourcedis,
                               const std::string& couplingcond_master,
                               const std::string& couplingcond_slave);

    /*!
    \brief Destructor

    */
     ~DofSetMergedProxy(){};

     /// Assign dof numbers using all elements and nodes of the discretization.
     virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

     //! @name Access methods

     /*!
     \brief Get the gid of all dofs of a node
      \note convention: First the slave dofs and then the master dofs are inserted into full dof vector!
                        Thus all definitions in the input file concerning dof numbering have to be set accordingly
                        (e.g. for reactions in MAT_scatra_reaction and MAT_scatra_reaction,
                         see test case 'ssi_3D_tet4_tet4_tri3.dat')
     */
     virtual std::vector<int> Dof(const Node* node) const
       {
         std::vector<int> slavedof = DofSetProxy::Dof(node);
         const DRT::Node* masternode = GetMasterNode(node->LID());
         std::vector<int> masterdof = DofSetProxy::Dof(masternode);

         std::vector<int> dof;
         dof.reserve( slavedof.size() + masterdof.size() ); // preallocate memory
         dof.insert( dof.end(), slavedof.begin(), slavedof.end() );
         dof.insert( dof.end(), masterdof.begin(), masterdof.end() );

         return dof;
       }

     /// Get the gid of all dofs of a element
     virtual std::vector<int> Dof(const Element* element) const
       {
         return DofSetProxy::Dof(element);
       }

     /// Get number of dofs for given node
     virtual int NumDof(const Node* node) const
       {
       const DRT::Node* masternode = GetMasterNode(node->LID());

         return DofSetProxy::NumDof(node)+DofSetProxy::NumDof(masternode);
       }

     /// Get number of dofs for given element
     virtual int NumDof(const Element* element) const
       {
         return DofSetProxy::NumDof(element);
       }

     /// Get the gid of a dof for given node
     virtual int Dof(const Node* node, int dof) const
       {
         const int numslavedofs = DofSetProxy::NumDof(node);
         if(dof<numslavedofs)
           return DofSetProxy::Dof(node,dof);
         else
         {
           const DRT::Node* masternode = GetMasterNode(node->LID());
           return DofSetProxy::Dof(masternode,dof-numslavedofs);
         }
       }

     /// Get the gid of a dof for given element
     virtual int Dof(const Element* element, int dof) const
       {
         return DofSetProxy::Dof(element,dof);
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(const Node* node, std::vector<int>& lm) const
       {
         DofSetProxy::Dof(node,lm);
         const DRT::Node* masternode = GetMasterNode(node->LID());
         DofSetProxy::Dof(masternode,lm);

         return;
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(
         const Node* node,           ///< node, for which you want the dof positions
         const unsigned startindex,  ///< first index of vector at which will be written to end
         std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
         ) const
     {
       const int numslavedofs = DofSetProxy::NumDof(node);
       DofSetProxy::Dof(node,startindex,lm);

       const DRT::Node* masternode = GetMasterNode(node->LID());
       DofSetProxy::Dof(masternode,startindex+numslavedofs,lm);
     }

     /// Get the gid of all dofs of a element
     virtual void Dof(const Element* element, std::vector<int>& lm) const
       {
         DofSetProxy::Dof(element,lm);
       }

     /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
     virtual void Dof(
         const Element* element,     ///< element which provides its expected number of DOFs per node
         const Node* node,           ///< node, for which you want the DOF positions
         std::vector<int>& lm        ///< already allocated vector to be filled with DOF positions
         ) const
     {
       DofSetProxy::Dof(element,node,lm);

       const DRT::Node* masternode = GetMasterNode(node->LID());
       DofSetProxy::Dof(element,masternode,lm);
     }

  private:

     const DRT::Node* GetMasterNode(int slaveLid) const
     {
       if(slaveLid>slavetomasternodemapping_->MyLength())
         dserror("Target Lid %i not found in discretization!",slaveLid);
       int mastergid = (*slavetomasternodemapping_)[slaveLid];
       if(not sourcedis_->HaveGlobalNode(mastergid))
         dserror("Master Gid %i not found in source discretization!",mastergid);
       return sourcedis_->gNode(mastergid);
     }

    //! map containing the mapping of the slave node ID (key) to the corresponding master node (value)
    Teuchos::RCP<Epetra_IntVector>  slavetomasternodemapping_;

    //! underlying proxy of actual dofset
    Teuchos::RCP<DofSet> dofset_;

    //! source discretization
    const Teuchos::RCP<const DRT::Discretization> sourcedis_;

    //! condition strings defining the coupling
    const std::string couplingcond_master_;
    const std::string couplingcond_slave_;
  };
}


#endif /* SRC_DRT_LIB_DRT_DOFSET_MERGED_PROXY_H_ */
