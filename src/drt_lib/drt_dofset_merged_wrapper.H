/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_merged_wrapper.H

 \brief A dofset that adds additional, existing degrees of freedom from the same
        discretization to nodes (not yet to elements).

 \level 2

 \maintainer  Andreas Rauch
              rauch@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_DOFSET_MERGED_WRAPPER_H_
#define SRC_DRT_LIB_DRT_DOFSET_MERGED_WRAPPER_H_


#include "drt_dofset_proxy.H"

#include "drt_discret.H"

#include <Epetra_IntVector.h>

namespace DRT
{
/*!
\brief A dofset that adds additional, existing degrees of freedom from the same
       discretization to nodes.

\warning not implemented for element DOFs

  The Dofs of the nodes to be merged are defined by master and slave side conditions given as input.
  Overlapping nodes are identified using a search tree and this dofset will handle the dofs of one node
  as if there were one in the Dof(..) and NumDof(..) methods (see the implementation of theses methods).

\warning For the Dof(..) methods providing the full dof vector an ordering of the nodes is assumed. That
         is, first the Dofs from the slave node are filled into the dof vector followed by the Dofs of the
         master node.

*/

  class DofSetMergedWrapper : public DofSetBase
  {
  public:

    //! Standard Constructor
    DofSetMergedWrapper(
        Teuchos::RCP<DofSetInterface>  dofset,
        Teuchos::RCP<const DRT::Discretization> sourcedis,
        const std::string& couplingcond_master,
        const std::string& couplingcond_slave);

     //! Destructor
     ~DofSetMergedWrapper();

     /// Returns true if filled
     virtual bool Filled() const
     {return filled_ and sourcedofset_->Filled();};

     /// Assign dof numbers using all elements and nodes of the discretization.
     virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

     /// reset all internal variables
     virtual void Reset();

     //! @name Proxy management
     /// Proxies need to know about changes to the DofSet.

     /// our original DofSet dies
     virtual void Disconnect(DofSetInterface* dofset);

     //@}

     /// Get degree of freedom row map
     virtual const Epetra_Map* DofRowMap() const;

     /// Get degree of freedom column map
     virtual const Epetra_Map* DofColMap() const;

     //! @name Access methods

     /// Get number of dofs for given node
     virtual int NumDof(
         const Node* node
         ) const
     {
       const DRT::Node* masternode = GetMasterNode(node->LID());
       const DRT::Node* slavenode = GetSlaveNode(node->LID());
       return sourcedofset_->NumDof(slavenode)+sourcedofset_->NumDof(masternode);
     }

     /// Get number of dofs for given element
     virtual int NumDof(
         const Element* element
         ) const
      {return sourcedofset_->NumDof(element);}

     /// get number of nodal dofs
     virtual int NumDofPerNode(
         const Node & node ///< node, for which you want to know the number of dofs
         ) const
     {
       const DRT::Node* masternode = GetMasterNode(node.LID());
       const DRT::Node* slavenode  = GetSlaveNode(node.LID());
       return sourcedofset_->NumDofPerNode(*masternode)+sourcedofset_->NumDofPerNode(*slavenode);
     };

     /*! \brief Get the gid of all dofs of a node
      *
      \note convention: First the slave dofs and then the master dofs are inserted into full dof vector!
                        Thus all definitions in the input file concerning dof numbering have to be set accordingly
                        (e.g. for reactions in MAT_scatra_reaction and MAT_scatra_reaction,
                         see test case 'ssi_3D_tet4_tet4_tri3.dat')  */
     virtual std::vector<int> Dof(const Node* node) const
       {
         const DRT::Node* slavenode = GetSlaveNode(node->LID());
         std::vector<int> slavedof = sourcedofset_->Dof(slavenode);
         const DRT::Node* masternode = GetMasterNode(node->LID());
         std::vector<int> masterdof = sourcedofset_->Dof(masternode);

         std::vector<int> dof;
         dof.reserve( slavedof.size() + masterdof.size() ); // preallocate memory
         dof.insert( dof.end(), slavedof.begin(), slavedof.end() );
         dof.insert( dof.end(), masterdof.begin(), masterdof.end() );

         return dof;
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(
         std::vector<int>& dof, ///< vector of dof gids (to be filled)
         const Node* node,      ///< node, for which you want the dof positions
         unsigned nodaldofset   ///< number of nodal dof set of the node (currently !=0 only for XFEM)
         ) const
     {
       const DRT::Node* masternode = GetMasterNode(node->LID());
       const DRT::Node* slavenode  = GetSlaveNode(node->LID());
       std::vector<int> slavedof;
       sourcedofset_->Dof(slavedof,slavenode,nodaldofset);
       std::vector<int> masterdof;
       sourcedofset_->Dof(masterdof,masternode,nodaldofset);

       dof.reserve( slavedof.size() + masterdof.size() ); // preallocate memory
       dof.insert( dof.end(), slavedof.begin(), slavedof.end() );
       dof.insert( dof.end(), masterdof.begin(), masterdof.end() );
     }

     /// Get the gid of all dofs of a element
     virtual std::vector<int> Dof(const Element* element) const
       {
         return sourcedofset_->Dof(element);
       }

     /// Get the gid of a dof for given node
     virtual int Dof(const Node* node, int dof) const
       {
         const DRT::Node* slavenode = GetSlaveNode(node->LID());
         const int numslavedofs = sourcedofset_->NumDof(slavenode);
         if(dof<numslavedofs)
           return sourcedofset_->Dof(slavenode,dof);
         else
         {
           const DRT::Node* masternode = GetMasterNode(node->LID());
           return sourcedofset_->Dof(masternode,dof-numslavedofs);
         }
       }

     /// Get the gid of a dof for given element
     virtual int Dof(const Element* element, int dof) const
       {
         return sourcedofset_->Dof(element,dof);
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(const Node* node, std::vector<int>& lm) const
       {
         const DRT::Node* masternode = GetMasterNode(node->LID());
         const DRT::Node* slavenode = GetSlaveNode(node->LID());
         sourcedofset_->Dof(slavenode,lm);
         sourcedofset_->Dof(masternode,lm);

         return;
       }

     /// Get the gid of all dofs of a node
     virtual void Dof(
         const Node* node,           ///< node, for which you want the dof positions
         const unsigned startindex,  ///< first index of vector at which will be written to end
         std::vector<int>& lm        ///< already allocated vector to be filled with dof positions
         ) const
     {
       const DRT::Node* slavenode = GetSlaveNode(node->LID());
       const int numslavedofs = sourcedofset_->NumDof(slavenode);
       sourcedofset_->Dof(slavenode,startindex,lm);

       const DRT::Node* masternode = GetMasterNode(node->LID());
       sourcedofset_->Dof(masternode,startindex+numslavedofs,lm);
     }

     /// Get the gid of all dofs of a element
     virtual void Dof(const Element* element, std::vector<int>& lm) const
     {sourcedofset_->Dof(element,lm); }

     /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
     virtual void Dof(
         const Element* element,     ///< element which provides its expected number of DOFs per node
         const Node* node,           ///< node, for which you want the DOF positions
         std::vector<int>& lm        ///< already allocated vector to be filled with DOF positions
         ) const
     {
       const DRT::Node* slavenode = GetSlaveNode(node->LID());
       sourcedofset_->Dof(element,slavenode,lm);

       const DRT::Node* masternode = GetMasterNode(node->LID());
       sourcedofset_->Dof(element,masternode,lm);
     }

     /// Get maximum GID of degree of freedom row map
     virtual int MaxAllGID() const
     {return sourcedofset_->MaxAllGID();};

     /// Get minimum GID of degree of freedom row map
     virtual int MinAllGID() const
     {return sourcedofset_->MinAllGID();};

     /// Get Max of all GID assigned in the DofSets in front of current one in the list #static_dofsets_
     virtual int MaxGIDinList(const Epetra_Comm& comm) const
     {return sourcedofset_->MaxGIDinList(comm);};

     /// are the dof maps already initialized?
     virtual bool Initialized() const
     {return sourcedofset_->Initialized();};

     /// Get Number of Global Elements of degree of freedom row map
     virtual int NumGlobalElements() const
     {return sourcedofset_->NumGlobalElements();};

  private:

     //! get the master node to a corresponding slave node (given by LID)
     const DRT::Node* GetMasterNode(int slaveLid) const
     {
       dsassert(slaveLid<=master_nodegids_col_layout_->MyLength(),"Slave node Lid out of range!");
       int mastergid = (*master_nodegids_col_layout_)[slaveLid];
       return sourcedis_->gNode(mastergid);
     }

     //! get the slave node to a corresponding master node (given by LID)
     const DRT::Node* GetSlaveNode(int masterLid) const
     {
       dsassert(masterLid<=slave_nodegids_col_layout_->MyLength(),"Master node Lid out of range!");
       int slavegid = (*slave_nodegids_col_layout_)[masterLid];
       return sourcedis_->gNode(slavegid);
     }

    //! master node gids in col layout matching conditioned slave nodes
    Teuchos::RCP<Epetra_IntVector>  master_nodegids_col_layout_;

    //! slave node gids in col layout matching conditioned master nodes
    Teuchos::RCP<Epetra_IntVector>  slave_nodegids_col_layout_;

    //! underlying actual dofset
    Teuchos::RCP<DofSetInterface> sourcedofset_;

    //! source discretization
    Teuchos::RCP<const DRT::Discretization> sourcedis_;

    //! condition strings defining the coupling
    const std::string couplingcond_master_;
    const std::string couplingcond_slave_;

    /// filled flag
    bool filled_;
  };
}


#endif /* SRC_DRT_LIB_DRT_DOFSET_MERGED_WRAPPER_H_ */
