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

/*!
\brief A proxy of a dofset that does not rely on same GID/LID numbers but uses
       a defined node mapping instead (not implemented for element DOFs).

This dofset proxy was originally written to handle heterogeneous reactions :
\ref SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy .
This means reactions between volume-bound scalars and surface-bound scalars in the
context of scalar transport simulations.
It is also used in solid-scatra interaction \ref SSI::SSI_Base::SSI_Base .
One main application is the Cell Migration Model \ref CellMigrationControlAlgorithm() .
For further information see also the documentation of the class :
\ref SSI::SSI_Part2WC_PROTRUSIONFORMATION::SSI_Part2WC_PROTRUSIONFORMATION .


<h3>AssignDegreesOfFreedom</h3>

The important work is done in \ref AssignDegreesOfFreedom . If you want to use
or extend this class, please be aware of the documentation of this method. Also
have in mind the needs of those models and methods in BACI this class was originally
intended to support (see above).


<h3>Usage</h3>

After construction the newly created dofset has to be added to the target discretization.
The target discretization then has to be filled. During \ref FillComplete() the method
\ref AssignDegreesOfFreedom is called.

One example for the proper use of this class:

\code

  Teuchos::RCP<DRT::DofSet> newdofset_for_target_dis =
      Teuchos::rcp(new DRT::DofSetMappedProxy(sourcedis->GetDofSetProxy(),sourcedis,"CouplingSourceToTarget",setofcouplingids));

  target_dis->AddDofSet(newdofset_for_target_dis);

  target_dis->FillComplete();

\endcode

See also the example '.dat'-file 'immersed_cell_biochemo_mechano_pureProtrusion_h8' in
the BACI 'Input' folder.

\note Date of this documentation 08/16 by author rauch ( rauch@lnm.mw.tum.de )

\date 08/2016
\author vuong (vuong@lnm.mw.tum.de)
*/
  class DofSetMappedProxy : public DofSetProxy
  {
  public:

    /*!
    \brief Standard Constructor

    \param dofset       (in) : dofset proxy of source discretization \ref sourcedis
    \param sourcedis    (in) : source discretization for the new dofset
    \param couplingcond (in) : condition to be coupled
    \param condids      (in) : set of condition ids to merge into one dofset

    */
    DofSetMappedProxy(Teuchos::RCP<DofSet>  dofset,
                      const Teuchos::RCP<const DRT::Discretization> sourcedis,
                      const std::string& couplingcond,
                      const std::set<int> condids);

    /*!
    \brief Destructor

    */
     ~DofSetMappedProxy(){};


     /*!
     \brief Assign dof numbers using all elements and nodes of the discretization.

            This method is called during \ref FillComplete() on the target
            discretization.


           <h3>Usage</h3>

           An important feature of this dofset proxy is, that portions of source and
           target discretizations, whic hare supposed to be coupled by this dofset do
           not necessarily need to have the same number of nodes, as long as the mapping
           is still unique. We give an example to show what this means:

           Imagine, we have given the following conditions in our .dat-file:

           -----DESIGN SSI COUPLING SOLIDTOSCATRA VOL CONDITIONS
           DVOL   2
           // scatra volume matching to struct volume
           E 2 - 1
           // struct volume matching to scatra volume
           E 1 - 1
           ----DESIGN SSI COUPLING SOLIDTOSCATRA SURF CONDITIONS
           DSURF   2
          // separate surface discretization matching struct volume boundary
          E 1 - 2
          // struct volume boundary
          E 2 - 2

          We have two matching volume discretizations for structure and scatra, and
          we have a scatra surface discretization matching the boundary of the struct
          volume. Thus, we have more scatra nodes, than structure nodes, because we
          have NumNodeScatraVolume + NumNodeScatraSurface scatra nodes, but only
          NumNodeStructVolume (= NumNodeScatraVolume) structure nodes.
          However, since we are matching, each scatra node can be mapped uniquiely to
          a structure node. Note that we could not map each structure node uniquely
          to a scatra node because at the boundary two scatra nodes reside at the same
          location (volume + surface)!
          So in the first case, if we want to add the structure dofset to the scatradis,
          the structure is the source and the scatra is the target.
          We can now create the dofset proxy from the struct dis via

          \code
           Teuchos::RCP<DRT::DofSet> newdofset =
               Teuchos::rcp(new DRT::DofSetMappedProxy(structdis->GetDofSetProxy(),
                                                       structdis,
                                                       "SSICouplingSolidToScatra",
                                                       couplingids));
          \endcode

          The std::set 'couplingids' contains the ids from both, VOL and SURF conditions
          given above. This way, one unique struct dofset for all scatra nodes is created.

     */
     virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

     //! @name Access methods

     /// Get the gid of all dofs of a node
     virtual std::vector<int> Dof(const Node* node) const
       {
         const DRT::Node* sourcenode = GetSourceNode(node->LID());
         if(sourcenode==NULL)
           return std::vector<int>();
         else
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
         if(sourcenode==NULL)
           return 0;
         else
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
          if(sourcenode==NULL)
            return -1;
          else
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
         if(sourcenode==NULL)
             return;
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
       if(sourcenode==NULL)
           return;
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
       if(sourcenode==NULL)
           return;
       DofSetProxy::Dof(element,sourcenode,lm);
     }

  private:

     //! get corresponding source node from source discretization
     const DRT::Node* GetSourceNode(int targetLid) const
     {
       if(targetLid>targetlidtosourcegidmapping_->MyLength())
         dserror("Target Lid %i not found in discretization!",targetLid);
       int sourcegid = (*targetlidtosourcegidmapping_)[targetLid];

       if(sourcegid==-1)
         return NULL;
       if(not sourcedis_->HaveGlobalNode(sourcegid))
         dserror("Master Gid %i not found in source discretization!",sourcegid);
       return sourcedis_->gNode(sourcegid);
     }

    //! map containing the mapping of the target node GID to the corresponding source node GID (value)
     Teuchos::RCP<Epetra_IntVector>  targetlidtosourcegidmapping_;

    //! underlying actual dofset
    Teuchos::RCP<DofSet> dofset_;

    //! source discretization
    const Teuchos::RCP<const DRT::Discretization> sourcedis_;

    //! condition string defining the coupling
    const std::string couplingcond_;

    //! ID of condition the dofset is build from
    const std::set<int> condids_;
  };
}



#endif /* SRC_DRT_LIB_DRT_DOFSET_MAPPED_PROXY_H_ */
