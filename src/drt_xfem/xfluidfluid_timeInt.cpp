/*----------------------------------------------------------------------*/
/*!
\file xfluidfluid_timeInt.classes

\brief class provides functionalities for xfluid-fluid time integration
<pre>
Maintainer:
             Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <EpetraExt_MatrixMatrix.h>

#include "xfluidfluid_timeInt.H"

#include "../drt_xfem/xfem_dofset.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"

#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

#include "../drt_xfem/xfem_dofset_transparent_independent.H"

// -------------------------------------------------------------------
//  constructor
// -------------------------------------------------------------------
XFEM::XFluidFluidTimeIntegration::XFluidFluidTimeIntegration(
  const Teuchos::RCP<DRT::Discretization> &  bgdis,
  const Teuchos::RCP<DRT::Discretization> &  embdis,
  Teuchos::RCP<GEO::CutWizard>               wizard,
  int                                        step,
  enum INPAR::XFEM::XFluidFluidTimeInt       xfem_timeintapproach,
  const Teuchos::ParameterList&              params
  ) :
  bgdis_(bgdis),
  embdis_(embdis),
  step_(step),
  myrank_(bgdis->Comm().MyPID()),
  numproc_(bgdis->Comm().NumProc()),
  currentbgdofmap_(bgdis->DofRowMap()),
  timeintapproach_(xfem_timeintapproach),
  params_(params)

{
  // create node-to-dof-maps of the background fluid
  CreateBgNodeMaps(wizard);

  Teuchos::ParameterList& params_xfem  = params_.sublist("XFEM");
  gmsh_debug_out_ = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT");

  Teuchos::ParameterList& params_xf_gen = params_.sublist("XFLUID DYNAMIC/GENERAL");
  searchradius_fac_= params_xf_gen.get<double>("XFLUIDFLUID_SEARCHRADIUS");

  // in case the embedded discretization is empty on this proc
  if (!embdis_->NumMyRowElements())
  {
    minradius_ = 0.0;
    return;
  }

  // determine the radius of the search tree
  switch (embdis_->lRowElement(0)->Shape())
  {
  case DRT::Element::hex8:
    FindSearchRadius<DRT::Element::hex8>();
    break;
  case DRT::Element::hex20:
    FindSearchRadius<DRT::Element::hex20>();
    break;
  case DRT::Element::hex27:
    FindSearchRadius<DRT::Element::hex27>();
    break;
  default:
    dserror("Unsupported element shape %s!", DRT::DistypeToString(embdis_->lRowElement(0)->Shape()).c_str()); break;
  }
} // end constructor

// -------------------------------------------------------------------
// build maps of node ids to their dof-gids in this time step
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CreateBgNodeMaps(Teuchos::RCP<GEO::CutWizard> wizard)
{
  // map of standard nodes and their dof-ids

  // run through nodes on this proc
  for (int lid=0; lid<bgdis_->NumMyRowNodes(); lid++)
  {
    // get the node from the wizard and the discretization
    DRT::Node * node = bgdis_->lRowNode(lid);
    const int node_gid = node->Id();
    GEO::CUT::Node * n = wizard->GetNode(node_gid);

    if (n!=NULL) // xfem nodes
    {
      // set of volumecells which belong to the node;
      // the vector index corresponds to a dofset, with which the set of volume-cell sets are associated.
      // (for hex20 we have sets of volumecells, otherwise for linear
      // elements the size of the set is one, so 1 dofset --> 1 volume-cell set!)
      const std::vector<GEO::CUT::NodalDofSet* > & vcs_vec = n->NodalDofSets();


      // get the id from the parent node's parent element for every dofset
      for (size_t idofset = 0; idofset < vcs_vec.size(); ++ idofset)
      {
        // we just need the first vc to find out the parent element:
        // - get the volume-cell set
        // - get the first volume-cell (vc) of the set
        const GEO::CUT::VolumeCell * vc = *(vcs_vec.at(idofset)->VolumeCellComposite().begin()->begin());

        // get parent element id
        const int parent_id = vc->GetParentElementId();
        nodeToParentEle_[node_gid].insert(parent_id);
        DRT::Element * actele = bgdis_->gElement(parent_id);

        // get the element handle and the nds vector
        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector< int > > nds_sets;
        GEO::CUT::ElementHandle * e = wizard->GetElement( actele );
        e->GetVolumeCellsDofSets( cell_sets, nds_sets, false); //(include_inner=false)

        parentEleToNds_[parent_id] = nds_sets;
      }

      // Get the point from the cut, to determine node position w.r.t. interface position
      GEO::CUT::Point::PointPosition pos = n->point()->Position();

      switch (pos)
      {
        case GEO::CUT::Point::outside:
        {
          if (bgdis_->NumDof(node) > 0)
          {
            std::vector<int> gdofs = bgdis_->Dof(node);
            nodeToDof_std_np_[node_gid] = gdofs;
          }
          else
          {
            dserror("Node is located in physical fluid domain but has no dofs!");
          }
          break;
        }
        case GEO::CUT::Point::oncutsurface:
        {
          if (bgdis_->NumDof(node) > 0)
          {
            std::vector<int> gdofs = bgdis_->Dof(node);
            nodeToDof_std_np_[node_gid] = gdofs;
          }
          else
          {
            // Node is located on fluid-fluid interface, but has no dofs.
            // It happens if the embedded fluid is not fully encapsulated in
            // the background fluid, but shares common side surface - this can happen in pseudo-2D-examples.
          }
          break;
        }
        case GEO::CUT::Point::inside:
        {
          // enriched?
          if (bgdis_->NumDof(node) > 0)
          {
            std::vector<int> gdofs = bgdis_->Dof(node);
            nodeToDof_enriched_np_[node_gid] = gdofs;
          }
          // void otherwise!
          break;
        }
        default:
          dserror("Invalid point position %d.", pos); break;
      }
    }
    else if( bgdis_->NumDof(node) > 0) // no xfem node
    {
      std::vector<int> gdofs = bgdis_->Dof(node);
      nodeToDof_std_np_[node_gid] = gdofs;
    }
    else
      dserror("Non-xfem-node with ID %d has a invalid number of %d DOF!", node_gid, bgdis_->NumDof(node));
  }

#ifdef PARALLEL
  // build a reduced map of all noderowmaps of all processors
  const Epetra_Map * noderowmap = bgdis_->NodeRowMap();
  Teuchos::RCP<Epetra_Map> allnoderowmap = LINALG::AllreduceEMap(*noderowmap);
  // gather the information of all processors
  DRT::Exporter ex(*noderowmap,*allnoderowmap,bgdis_->Comm());
  ex.Export(nodeToDof_std_np_);
  ex.Export(nodeToDof_enriched_np_);
#endif

}

// -------------------------------------------------------------------
// save the old node-to-dofset maps of the background fluid
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SaveBgNodeMaps()
{
  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  nodeToDof_std_n_ = nodeToDof_std_np_;
  nodeToDof_enriched_n_ = nodeToDof_enriched_np_;
  nodeToDof_std_np_.clear();
  nodeToDof_enriched_np_.clear();
}

// -------------------------------------------------------------------
// - save the old node-to-dofset maps of the background fluid
// - create node- & dof-maps for the current step
// - give gmsh-output
// - return true, if the maps haven't changed
// -------------------------------------------------------------------
bool XFEM::XFluidFluidTimeIntegration::SaveBgNodeMapsAndCreateNew(Teuchos::RCP<GEO::CutWizard> wizard)
{

  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  SaveBgNodeMaps();

  // Create new node-to-dofset maps
  CreateBgNodeMaps(wizard);

  // Save the old dof row maps and get the new ones
  oldbgdofmap_ = currentbgdofmap_;
  currentbgdofmap_ =  bgdis_->DofRowMap();

  // determine a change in the background fluid state by comparing
  // maps from subsequent time steps
  if ((oldbgdofmap_->SameAs(*currentbgdofmap_) and (nodeToDof_std_n_ == nodeToDof_std_np_) and
       (nodeToDof_enriched_n_ == nodeToDof_enriched_np_)))
    samemaps_ = true;
  else samemaps_ = false;

  if (gmsh_debug_out_)
    GmshOutput();

  return samemaps_;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void  XFEM::XFluidFluidTimeIntegration::CreateBgNodeMapsForRestart(Teuchos::RCP<GEO::CutWizard>  wizard)
{
  // Create node-to-dof maps for the background fluid
  CreateBgNodeMaps(wizard);
  currentbgdofmap_ =  bgdis_->DofRowMap();
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void  XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectors(
  const Teuchos::RCP<Epetra_Vector>        bgstate_velnp,
  const Teuchos::RCP<Epetra_Vector>        bgstate_veln,
  const Teuchos::RCP<Epetra_Vector>        bgstate_velnm,
  const Teuchos::RCP<Epetra_Vector>        bgstate_accnp,
  const Teuchos::RCP<Epetra_Vector>        bgstate_accn,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnp,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_veln,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnm,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_accnp,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_accn,
  const Teuchos::RCP<const Epetra_Vector>  embstate_velnp,
  const Teuchos::RCP<const Epetra_Vector>  embstate_veln,
  const Teuchos::RCP<const Epetra_Vector>  embstate_velnm,
  const Teuchos::RCP<const Epetra_Vector>  embstate_accnp,
  const Teuchos::RCP<const Epetra_Vector>  embstate_accn,
  const Teuchos::RCP<const Epetra_Vector>  aledispn)
{
  // form state containers
  const TargetState bgfluid_state(
      bgstate_velnp,bgstate_veln,bgstate_velnm,bgstate_accnp,bgstate_accn
      );
  const SourceState bgfluid_staten(
      bgstaten_velnp,bgstaten_veln,bgstaten_velnm,bgstaten_accnp,bgstaten_accn
      );
  const SourceState emb_state(
      embstate_velnp,embstate_veln,embstate_velnm,embstate_accnp,embstate_accn
      );
  // REMARK: this method can write to the target Epetra_Vectors, but the pointers are const,
  // as the TargetState is const (const pointer to non-const pointee)
  SetNewBgStateVectors(
      bgfluid_state, bgfluid_staten, emb_state, aledispn);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectors(const TargetState &                       bgfluid_state,
                                                            const SourceState &                       bgfluid_state_n,
                                                            const SourceState &                       embfluid_state,
                                                            const Teuchos::RCP<const Epetra_Vector> & aledisp
)
{
  // do projection from embedded element to background fluid nodes
  if (timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
      timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved or
      timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
  {
    SetNewBgStateVectorFullProjection(bgfluid_state,
                                      bgfluid_state_n,
                                      embfluid_state,
                                      aledisp);
  } // no projection, use the ghost values instead
  else if (timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues)
  {
    SetNewBgStateVectorKeepGhostValues(bgfluid_state,
                                       bgfluid_state_n,
                                       embfluid_state,
                                       aledisp);
  }
  else
    dserror("Xfem time integration approach unknown!");
}

// -------------------------------------------------------------------
// If enriched values are available, keep them (no projection from
// embedded fluid in this case)
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectorKeepGhostValues(
  const TargetState &                       bgfluid_state,
  const SourceState &                       bgfluid_state_n,
  const SourceState &                       embfluid_state,
  const Teuchos::RCP<const Epetra_Vector> & aledispn
)
{
  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_xyz;
  // the vector containing interpolated values for each state vector from embedded dis. For each node one interpolated value..
  std::vector<LINALG::Matrix<20,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodes_nohistory;

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis_->lRowNode(lnid);
    std::map<int, std::vector<int> >::const_iterator iterstn  = nodeToDof_std_n_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterstnp = nodeToDof_std_np_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iteren   = nodeToDof_enriched_n_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterenp  = nodeToDof_enriched_np_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != nodeToDof_std_n_.end() and iterstnp != nodeToDof_std_np_.end()) or
        (iterstn != nodeToDof_std_n_.end() and iterenp != nodeToDof_enriched_np_.end()))
    {
      //int numsets = bgdis_->NumDof(bgnode)/4;
      std::vector<int> gdofsn = iterstn->second;

#ifdef DEBUG
      if (iterstn->second.size()>4)
        IO::cout << " INFO: more standard sets!!!! "<< "Node GID " << bgnode->Id() << IO::endl;
#endif

      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnp_,bgfluid_state_n.velnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.veln_, bgfluid_state_n.veln_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnm_,bgfluid_state_n.velnm_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accnp_,bgfluid_state_n.accnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accn_, bgfluid_state_n.accn_);

    }
    // Project dofs from embdis to bgdis:
    // n:void -> n+1:std, n:void ->  n+1:enriched
    else if (((iterstn == nodeToDof_std_n_.end() and iteren == nodeToDof_enriched_n_.end()) and iterstnp != nodeToDof_std_np_.end()) or
             ((iterstn == nodeToDof_std_n_.end() and iteren == nodeToDof_enriched_n_.end()) and iterenp  != nodeToDof_enriched_np_.end()))
    {
      // save the information needed from bg-nodes to call the method CommunicateNodes
      // set of are projected nodes
      projectedNodes_.insert(bgnode->Id());

      LINALG::Matrix<3,1> bgnode_xyz(bgnode->X());
      // collect the coordinates of bg-nodes without history values
      bgnodes_xyz.push_back(bgnode_xyz);

      // collect the bg-node ids which need a projection
      bgnodes_nohistory.push_back(bgnode->Id());

      // the vector of interpolated values
      LINALG::Matrix<20,1>    interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);
    }

    //keep the ghost dofs:
    //n: enriched -> n+1:enriched, n: enriched -> n+1: std
    else if ((iteren != nodeToDof_enriched_n_.end() and iterenp != nodeToDof_enriched_np_.end()) or
             (iteren != nodeToDof_enriched_n_.end() and iterstnp != nodeToDof_std_np_.end()))
    {
      const int numsets = bgdis_->NumDof(bgnode)/4;
      if (numsets > 1)
        IO::cout << "ghost-fluid-approach just available for one dofset!" << IO::endl;

      std::vector<int> gdofsn = iteren->second;
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnp_,bgfluid_state_n.velnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.veln_, bgfluid_state_n.veln_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnm_,bgfluid_state_n.velnm_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accnp_,bgfluid_state_n.accnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accn_, bgfluid_state_n.accn_);

    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn  == nodeToDof_std_n_.end()       and iteren    == nodeToDof_enriched_n_.end() and
               iterstnp == nodeToDof_std_np_.end()      and iterenp   == nodeToDof_enriched_np_.end() ) or
              (iterstn  != nodeToDof_std_n_.end()       and (iterstnp == nodeToDof_std_np_.end() and iterenp == nodeToDof_enriched_np_.end())) or
              (iteren   != nodeToDof_enriched_n_.end()  and (iterenp  == nodeToDof_enriched_np_.end() and iterstnp == nodeToDof_std_np_.end())) )
    {

      if ( bgdis_->NumDof(bgnode) > 0 )
        dserror("Node %d is not supposed to carry any dof!", bgnode->Id());
    }
    else
      dserror("Background node %d has reached an invalid state. Check the query!", bgnode->Id());
  }

  // build search tree for projection from embedded discretization
  SetupSearchTree(aledispn);

  // call the Round Robin Communicator
  CommunicateNodes(bgnodes_xyz,interpolated_vecs,bgnodes_nohistory,
                   aledispn,bgfluid_state,bgfluid_state_n,embfluid_state,bgdis_);

}//XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectorKeepGhostValues

// -------------------------------------------------------------------
// project values from embedded discretization
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectorFullProjection(
  const TargetState &                       bgfluid_state,
  const SourceState &                       bgfluid_state_n,
  const SourceState &                       embfluid_state,
  const Teuchos::RCP<const Epetra_Vector> & aledispn
)
{
  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_xyz;
  // the vector containing interpolated values for each state vector from embedded dis. For each node one interpolated value..
  // Entries 0 to 3: velnp_, Entries 4 to 7: veln_,  Entries 8 to 11: velnm_, Entries 12 to 15: accn_, Entries 16 to 19: accnp_
  std::vector<LINALG::Matrix<20,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodes_nohistory;

  projectedNodes_.clear();

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis_->lRowNode(lnid);
    std::map<int, std::vector<int> >::const_iterator iterstn  = nodeToDof_std_n_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterstnp = nodeToDof_std_np_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iteren   = nodeToDof_enriched_n_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterenp  = nodeToDof_enriched_np_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != nodeToDof_std_n_.end() and iterstnp != nodeToDof_std_np_.end()) or
        (iterstn != nodeToDof_std_n_.end() and iterenp  != nodeToDof_enriched_np_.end()))
    {
      std::vector<int> gdofsn = iterstn->second;

#ifdef DEBUG
      if (iterstn->second.size() > 4)
        IO::cout << "INFO: more standard sets!!!!"<< "Node GID " << bgnode->Id() << " size " << iterstn->second.size()  << IO::endl;
#endif
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnp_,bgfluid_state_n.velnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.veln_, bgfluid_state_n.veln_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.velnm_,bgfluid_state_n.velnm_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accnp_,bgfluid_state_n.accnp_);
      WriteValuesToBgStateVector(bgnode,gdofsn,bgfluid_state.accn_, bgfluid_state_n.accn_);

    }
    // Project dofs from embdis to bgdis:
    // n:void -> n+1:std, n:enriched -> n+1:enriched,
    // n:enriched -> n+1: std, n:void ->  n+1:enriched
    else if (((iterstn == nodeToDof_std_n_.end()      and iteren   == nodeToDof_enriched_n_.end()) and iterstnp != nodeToDof_std_np_.end()) or
             (iteren   != nodeToDof_enriched_n_.end() and iterenp  != nodeToDof_enriched_np_.end()) or
             (iteren   != nodeToDof_enriched_n_.end() and iterstnp != nodeToDof_std_np_.end()) or
             ((iterstn == nodeToDof_std_n_.end()      and iteren   == nodeToDof_enriched_n_.end()) and iterenp != nodeToDof_enriched_np_.end()))
    {

      // save the information needed from bg-nodes to call the method CommunicateNodes
      // set of are projected nodes
      projectedNodes_.insert(bgnode->Id());

      // collect the coordinates of bg-nodes without history values
      bgnodes_xyz.push_back(LINALG::Matrix<3,1>(bgnode->X()));

      // collect the bg-node ids which need a projection
      bgnodes_nohistory.push_back(bgnode->Id());

      // the vector of interpolated values of all state vector (max five state vectors) for every node
      LINALG::Matrix<20,1> interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);

    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn  == nodeToDof_std_n_.end()      and  iteren   == nodeToDof_enriched_n_.end()   and
               iterstnp == nodeToDof_std_np_.end()     and  iterenp  == nodeToDof_enriched_np_.end()) or
              (iterstn  != nodeToDof_std_n_.end()      and (iterstnp == nodeToDof_std_np_.end()       and iterenp == nodeToDof_enriched_np_.end())) or
              (iteren   != nodeToDof_enriched_n_.end() and (iterenp  == nodeToDof_enriched_np_.end()  and iterstnp == nodeToDof_std_np_.end())) )
    {

      if ( bgdis_->NumDof(bgnode) > 0 )
        dserror("No dofsets expected for node %d!", bgnode->Id());
    }
    else
      dserror("Unreasonable state change of background fluid node. Check the query.");
  }

  // build search tree for projection from embedded discretization
  SetupSearchTree(aledispn);

  // call the Round Robin Communicator
  CommunicateNodes(bgnodes_xyz,interpolated_vecs,
                   bgnodes_nohistory,aledispn,
                   bgfluid_state,bgfluid_state_n,
                   embfluid_state,
                   bgdis_);

}//XFEM::XFluidFluidTimeIntegration::SetNewBgStateVectorFullProjection

//-----------------------------------------------------------------------------
//------------------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CommunicateNodes(std::vector<LINALG::Matrix<3,1> >  &      nodes_xyz,
                                                        std::vector<LINALG::Matrix<20,1> > &      interpolated_vecs,
                                                        std::vector<int>  &                       nodes_nohistory,
                                                        const Teuchos::RCP<const Epetra_Vector> & aledispn,
                                                        const TargetState &                       target_state,
                                                        const SourceState &                       bgfluid_state_n,
                                                        const SourceState &                       embfluid_state,
                                                        const Teuchos::RCP<DRT::Discretization>&  discret
)
{

  // get number of processors and the current processors id
  const int numproc=embdis_->Comm().NumProc();

  //information how many processors work at all
  std::vector<int> allproc(numproc);

  // create an exporter for point to point comunication
  DRT::Exporter exporter(embdis_->Comm());

  // necessary variables
  MPI_Request request;

  // define send and receive blocks
  std::vector<char> sblock;
  std::vector<char> rblock;

  // vector which identifies if a node has already interpolated values (initialize to false)
  std::vector<int> have_values(nodes_nohistory.size(),0);

  //----------------------------------------------------------------------
  // communication is done in a round robin loop
  //----------------------------------------------------------------------
  for (int np=0; np<numproc+1; ++np)
  {
    // in the first step, we cannot receive anything
    if (np > 0)
    {
      ReceiveBlock(rblock,exporter,request);

      std::vector<char>::size_type position = 0;
      DRT::ParObject::ExtractfromPack(position,rblock,nodes_xyz);
      DRT::ParObject::ExtractfromPack(position,rblock,interpolated_vecs);
      DRT::ParObject::ExtractfromPack(position,rblock,nodes_nohistory);
      DRT::ParObject::ExtractfromPack(position,rblock,have_values);
    }

    // in the last step, we keep everything on this proc
    if (np < numproc)
    {
      // -----------------------
      // do what we wanted to do
      FindEmbEleAndInterpolateValues(nodes_xyz,interpolated_vecs,have_values,
                                     aledispn,embfluid_state);

      // Pack info into block to send it
      PackValues(nodes_xyz,interpolated_vecs,nodes_nohistory,have_values,sblock);

      // add size to sendblock
      SendBlock(sblock,exporter,request);
    }
  } // end of loop over processors

  //----------------------------------------------------------------------------------------------
  //set the interpolated values to statevec_np
  //---------------------------------------------------------------------------------------------
  for (size_t i=0; i<nodes_nohistory.size(); ++i)
  {
    const DRT::Node* node = discret->gNode(nodes_nohistory.at(i));
    // number of dof-sets
    const unsigned int numsets = discret->NumDof(node)/4;

    // if interpolated values are available
    if (have_values.at(i) == 1)
    {
      unsigned int offset = 0;
      for (unsigned int set=0; set<numsets; set++)
      {
        // offset for different state vectors
        for (unsigned int isd = 0; isd < 4; ++isd)
        {
          unsigned int countvecs = 0;
          (*target_state.velnp_)[target_state.velnp_->Map().LID(discret->Dof(node)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
          if (target_state.veln_ != Teuchos::null)
          {
            countvecs += 4;
            (*target_state.veln_) [target_state.veln_->Map().LID(discret->Dof(node)[offset+isd])]  = interpolated_vecs.at(i)(isd+countvecs);
          }
          if (target_state.velnm_ != Teuchos::null)
          {
            countvecs += 4;
            (*target_state.velnm_)[target_state.velnm_->Map().LID(discret->Dof(node)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
          }
          if (target_state.accnp_ != Teuchos::null)
          {
           countvecs += 4;
           (*target_state.accnp_)[target_state.accnp_->Map().LID(discret->Dof(node)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
          }
          if (target_state.accn_ != Teuchos::null)
          {
            countvecs += 4;
            (*target_state.accn_) [target_state.accn_->Map().LID(discret->Dof(node)[offset+isd])]  = interpolated_vecs.at(i)(isd+countvecs);
          }
        }
        offset += 4;
      }
    }
    // if no embedded element is found, try to find an enriched value (only when interpolating for background fluid)
    else if (discret->Name() == "xfluid")
    {
      std::map<int, std::vector<int> >::const_iterator iterstnp = nodeToDof_std_np_.find(node->Id());
      std::map<int, std::vector<int> >::const_iterator iteren   = nodeToDof_enriched_n_.find(node->Id());
      std::map<int, std::vector<int> >::const_iterator iterenp  = nodeToDof_enriched_np_.find(node->Id());

      if ((iteren != nodeToDof_enriched_n_.end() and iterenp != nodeToDof_enriched_np_.end())
          or ((iteren != nodeToDof_enriched_n_.end() and iterstnp != nodeToDof_std_np_.end())))
      {
        IO::cout << "CHECK: Took enriched values !!" << " Node GID " << node->Id() << IO::endl;
        IO::cout << " Warning: You may need to make your search radius bigger in the dat-file!" << IO::endl;
        std::vector<int> gdofsn = iteren->second;

        WriteValuesToBgStateVector(node,gdofsn,target_state.velnp_,bgfluid_state_n.velnp_);
        WriteValuesToBgStateVector(node,gdofsn,target_state.veln_, bgfluid_state_n.veln_);
        WriteValuesToBgStateVector(node,gdofsn,target_state.velnm_,bgfluid_state_n.velnm_);
        WriteValuesToBgStateVector(node,gdofsn,target_state.accnp_,bgfluid_state_n.accnp_);
        WriteValuesToBgStateVector(node,gdofsn,target_state.accn_, bgfluid_state_n.accn_);
      }
      else
      {
#ifdef DEBUG
        std::map<int, std::vector<int> >::const_iterator iterstn  = nodeToDof_std_n_.find(node->Id());

        IO::cout << " Warning: No patch element found for the node " << node->Id();
        if ((iterstn == nodeToDof_std_n_.end() and iteren == nodeToDof_enriched_n_.end()) and iterstnp != nodeToDof_std_np_.end())
          IO::cout << " n:void -> n+1:std  " << IO::endl;
        else if (iteren != nodeToDof_enriched_n_.end() and iterenp != nodeToDof_enriched_np_.end())
          IO::cout << " n:enriched -> n+1:enriched " << IO::endl;
        else if (iteren != nodeToDof_enriched_n_.end() and iterstnp != nodeToDof_std_np_.end())
          IO::cout << " n:enriched -> n+1: std " << IO::endl;
        else if ((iterstn == nodeToDof_std_n_.end() and iteren == nodeToDof_enriched_n_.end()) and iterenp != nodeToDof_enriched_np_.end())
          IO::cout << " n:void ->  n+1:enriched " << IO::endl;
#endif
        if (embdis_->Comm().MyPID() == 0)
          IO::cout << " Warning: You may need to make your search radius bigger in the dat-file!" << IO::endl;
      }
    }
  }// end of loop over nodes without history

 }//XFEM::XFluidFluidTimeIntegration::CommunicateNodes

//---------------------------------------------------------
// receive a block in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::ReceiveBlock(std::vector<char> &   rblock,
                                                    DRT::Exporter  &      exporter,
                                                    MPI_Request    &      request)
{
  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();
  int myrank =embdis_->Comm().MyPID();

  // necessary variables
  int length =-1;
  int frompid=(myrank+numproc-1)%numproc;
  int tag    =frompid;

  // receive from predecessor
  exporter.ReceiveAny(frompid,tag,rblock,length);

#ifdef DEBUG
  IO::cout << "----receiving " << rblock.size() <<  " bytes: to proc " << myrank << " from proc " << frompid << IO::endl;
#endif

  if (tag!=(myrank+numproc-1)%numproc)
  {
    dserror("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();

  return;
} // XFEM::XFluidFluidTimeIntegration::ReceiveBlock

//---------------------------------------------------------
// send a block in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SendBlock(std::vector<char>  & sblock  ,
                                                 DRT::Exporter & exporter,
                                                 MPI_Request   & request )
{
  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();
  int myrank =embdis_->Comm().MyPID();

  // Send block to next proc.
  int tag    =myrank;
  int frompid=myrank;
  int topid  =(myrank+1)%numproc;

#ifdef DEBUG
   IO::cout << "----sending " << sblock.size() <<  " bytes: from proc " << myrank << " to proc " << topid << IO::endl;
#endif

  exporter.ISend(frompid,topid,
                 &(sblock[0]),sblock.size(),
                 tag,request);

  // for safety
  exporter.Comm().Barrier();

  return;
} // XFluidFluidTimeIntegration::SendBlock

//---------------------------------------------------------
// pack values in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PackValues(std::vector<LINALG::Matrix<3,1> >  & nodes_xyz,
                                                  std::vector<LINALG::Matrix<20,1> > & interpolatedvec,
                                                  std::vector<int>                   & nodes_nohistory,
                                                  std::vector<int>                   & have_values,
                                                  std::vector<char>                  & sblock)
{
  // Pack info into block to send
  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data,nodes_xyz);
  DRT::ParObject::AddtoPack(data,interpolatedvec);
  DRT::ParObject::AddtoPack(data,nodes_nohistory);
  DRT::ParObject::AddtoPack(data,have_values);
  data.StartPacking();

  DRT::ParObject::AddtoPack(data,nodes_xyz);
  DRT::ParObject::AddtoPack(data,interpolatedvec);
  DRT::ParObject::AddtoPack(data,nodes_nohistory);
  DRT::ParObject::AddtoPack(data,have_values);
  swap( sblock, data() );

} // XFluidFluidTimeIntegration::PackValue

//--------------------------------------------------------
//
//--------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::FindEmbEleAndInterpolateValues(std::vector<LINALG::Matrix<3,1> >        & nodes_xyz,
                                                                      std::vector<LINALG::Matrix<20,1> >       & interpolated_vecs,
                                                                      std::vector<int>                         & have_values,
                                                                      const Teuchos::RCP<const Epetra_Vector>  & aledispn,
                                                                      const SourceState &                        embfluid_state)

{
  // loop over the nodes (coordinates)
  for (size_t i=0; i<nodes_xyz.size(); ++i)
  {
    // indicates that we found the embedded element, the background node is covered by
    bool insideelement = false;

    // node coordinate
    const LINALG::Matrix<3,1> & node_xyz = nodes_xyz.at(i);
    // interpolated vector which is zero at the beginning
    LINALG::Matrix<20,1> interpolatedvec = interpolated_vecs.at(i);

    //search for near elements
    std::map<int,std::set<int> > closeeles =
        searchTree_->searchElementsInRadius(*embdis_,emb_nodepositions_n_,node_xyz,minradius_,0);

    // Remark: it could be that closeles is empty on one processor but still has elements on other processors.

    if (closeeles.empty())
    {
      IO::cout << "The search radius is empty on one processor! You may need to change the XFLUIDFLUID_SEARCHRADIUS is your dat-file."<< IO::endl;
      continue;
    }

    // loop over the map of background node-IDs and elements within the search radius
    for (std::map<int, std::set<int> >::const_iterator closele = closeeles.begin(); closele != closeeles.end(); closele++)
    {
      if (insideelement) break;
      // loop over the set of embedded elements within the search radius
      for (std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
      {
        if (insideelement) break;
        DRT::Element* pele = embdis_->gElement(*eleIter); // eleIter is the gid of the pele
        // determine values for background fluid node (either projected from embedded element or ghost values)
        switch (pele->Shape())
        {
        case DRT::Element::hex8:
          insideelement = ComputeSpatialToElementCoordAndProject<DRT::Element::hex8>(
             pele,embfluid_state,aledispn,*embdis_,node_xyz,interpolatedvec);
          break;
        case DRT::Element::hex20:
          insideelement = ComputeSpatialToElementCoordAndProject<DRT::Element::hex20>(
             pele,embfluid_state,aledispn,*embdis_,node_xyz,interpolatedvec);
          break;
        case DRT::Element::hex27:
          insideelement = ComputeSpatialToElementCoordAndProject<DRT::Element::hex27>(
             pele,embfluid_state,aledispn,*embdis_,node_xyz,interpolatedvec);
          break;
        default:
          dserror("Unsupported element shape %s!", DRT::DistypeToString(pele->Shape()).c_str()); break;
        }

        if (insideelement and have_values.at(i)==0)
        {
          // mark the node id, if the embedded-element is found
          have_values.at(i) = 1;
          for (size_t j=0; j<20; ++j)
          {
            //set the interpolated values
            interpolated_vecs.at(i)(j) = interpolatedvec(j);
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------
//
//--------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetupSearchTree(const Teuchos::RCP<const Epetra_Vector>  & aledispn)
{
  // init of 3D search tree
  searchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // find previous node positions for embedded fluid discretization
  for (int lid = 0; lid < embdis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = embdis_->lColNode(lid);
    std::vector<int> src_dofs(4);
    std::vector<double> mydisp(4);

    // get the current displacement
    embdis_->Dof(node,0,src_dofs);
    DRT::UTILS::ExtractMyValues(*aledispn,mydisp,src_dofs);

    emb_nodepositions_n_[node->Id()](0) = node->X()[0]+mydisp.at(0);
    emb_nodepositions_n_[node->Id()](1) = node->X()[1]+mydisp.at(1);
    emb_nodepositions_n_[node->Id()](2) = node->X()[2]+mydisp.at(2);
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*embdis_,emb_nodepositions_n_);
  searchTree_->initializeTree(rootBox,*embdis_,GEO::TreeType(GEO::OCTTREE));
}

//-------------------------------------------------------------------
// Write the values of node from bgstatevec_n to bgstatevec_np
//--------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::WriteValuesToBgStateVector(
  const DRT::Node*                           node,
  std::vector<int>                           gdofs_n,
  const Teuchos::RCP<Epetra_Vector> &        statevec,
  const Teuchos::RCP<const Epetra_Vector> &  statevec_n)
{
  const unsigned int dofset_size = 4;
  const unsigned int numsets = bgdis_->NumDof(node)/dofset_size;

#ifdef DEBUG
  if (numsets > 1)
    IO::cout << "Info: more dofsets in transfer.. " <<  "Node GID " << node->Id() << IO::endl;
#endif

  unsigned int offset = 0;
  for (unsigned int set=0; set<numsets; set++)
  {
    for (unsigned int d=0; d<dofset_size; ++d)
    {
      (*statevec)[statevec->Map().LID(bgdis_->Dof(node)[offset+d])] =
        (*statevec_n)[statevec_n->Map().LID(gdofs_n[d])];
    }
    offset += dofset_size;
  }
}//XFEM::XFluidFluidTimeIntegration::WriteValuesToBgStateVector

// -------------------------------------------------------------------
//  Check whether the unknown node lies in the source element.
//  If yes fill the interpolatedvec and return true.
// -------------------------------------------------------------------
template<DRT::Element::DiscretizationType distype>
bool XFEM::XFluidFluidTimeIntegration::ComputeSpatialToElementCoordAndProject(
  const DRT::Element *                      src_ele,
  const SourceState &                       src_state,
  const Teuchos::RCP<const Epetra_Vector> & src_disp,
  const DRT::Discretization &               src_dis,
  const LINALG::Matrix<3,1> &               node_xyz,
  LINALG::Matrix<20,1> &                    interpolatedvec
)
{
  // problem dimension
  const unsigned int dim = DRT::UTILS::DisTypeToDim<distype>::dim;
  // number of dof per node
  const unsigned int numdofpernode = dim + 1;

  // embedded element's nodes
  const DRT::Node*const* src_elenodes = src_ele->Nodes();
  // number of embedded element's nodes
  const unsigned int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  // embedded element's coordinates
  LINALG::Matrix<dim,numnodes> src_xyze(true);

  // flag, whether the node is covered by src_ele or not
  bool inside = false;

  std::vector<double> myval(numdofpernode);
  std::vector<double> mydisp(dim);
  std::map<unsigned int, std::vector<int> > src_dofs;

  // loop over the nodes of the embedded element and save the node coordinates!
  // the nodal dof are saved afterwards, as we don't need them, unless the background fluid node was covered
  for (unsigned int inode = 0; inode < numnodes; ++inode)
  {
    src_dofs[inode].resize(numdofpernode,0);
    src_dis.Dof(src_elenodes[inode],0,src_dofs[inode]);
    DRT::UTILS::ExtractMyValues(*src_disp,mydisp,src_dofs[inode]);

    for (unsigned int d=0; d<dim; ++d)
    {
      // get the coordinates of patch element
      // add the current displacement to the coordinates
      src_xyze(d,inode) = src_elenodes[inode]->X()[d] + mydisp[d];
    }
  }

  // check whether the node with missing values is located in the projection source element
  GEO::CUT::Position<distype> pos(src_xyze,node_xyz);
  const double tol = 1e-10;
  inside = pos.ComputeTol(tol);

  // von ursula
  // bool in  = GEO::currentToVolumeElementCoordinates(DRT::Element::hex8, pxyze, x, xsi);
  // in  = GEO::checkPositionWithinElementParameterSpace(xsi, DRT::Element::hex8);

  if (inside)
  {
    // get the coordinates of x in element coordinates of associated embedded element (xsi)
    LINALG::Matrix<3,1> xsi(true);
    xsi = pos.LocalCoordinates();

    // evaluate shape function
    LINALG::SerialDenseVector shp(numnodes);
    DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), distype );
    // extract the embedded element state values and interpolate
    for (unsigned inode = 0; inode<numnodes; ++inode)
    {
      unsigned int offset = 0;
      if (src_state.velnp_ != Teuchos::null)
      {
        DRT::UTILS::ExtractMyValues(*src_state.velnp_,myval,src_dofs[inode]);
        for (unsigned int isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(inode);
        }
        myval.clear();
      }

      offset += numdofpernode;
      if (src_state.veln_ != Teuchos::null)
      {
        DRT::UTILS::ExtractMyValues(*src_state.veln_,myval,src_dofs[inode]);
        for (unsigned int isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(inode);
        }
        myval.clear();
      }

      offset += numdofpernode;
      if (src_state.velnm_ != Teuchos::null)
      {
        DRT::UTILS::ExtractMyValues(*src_state.velnm_,myval,src_dofs[inode]);
        for (unsigned int isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(inode);
        }
        myval.clear();
      }

      offset += numdofpernode;
      if (src_state.accnp_ != Teuchos::null)
      {
        DRT::UTILS::ExtractMyValues(*src_state.accnp_,myval,src_dofs[inode]);
        for (unsigned int isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(inode);
        }
        myval.clear();
      }

      offset += numdofpernode;
      if (src_state.accn_ != Teuchos::null)
      {
        DRT::UTILS::ExtractMyValues(*src_state.accn_,myval,src_dofs[inode]);
        for (unsigned int isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(inode);
        }
        myval.clear();
      }
    }
  }
  return inside;
}//ComputeSpatialToElementCoordAndProject

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewEmbStateVectors(
  const Teuchos::RCP<Epetra_Vector>        embstate_velnp,
  const Teuchos::RCP<Epetra_Vector>        embstate_veln,
  const Teuchos::RCP<Epetra_Vector>        embstate_velnm,
  const Teuchos::RCP<Epetra_Vector>        embstate_accnp,
  const Teuchos::RCP<Epetra_Vector>        embstate_accn,
  const Teuchos::RCP<const Epetra_Vector>  embstaten_velnp,
  const Teuchos::RCP<const Epetra_Vector>  embstaten_veln,
  const Teuchos::RCP<const Epetra_Vector>  embstaten_velnm,
  const Teuchos::RCP<const Epetra_Vector>  embstaten_accnp,
  const Teuchos::RCP<const Epetra_Vector>  embstaten_accn,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnp,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_veln,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnm,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_accnp,
  const Teuchos::RCP<const Epetra_Vector>  bgstaten_accn,
  const Teuchos::RCP<const Epetra_Vector>  aledispnp,
  const Teuchos::RCP<const Epetra_Vector>  aledispn)
{
  if (gmsh_debug_out_)
    GmshOutputForInterpolateFSI(aledispnp,aledispn);

  // our first option is to find another embedded element -
  // we then interpolate from it's values
  // form state containers
  const TargetState embfluid_state(
      embstate_velnp,embstate_veln,embstate_velnm,embstate_accnp,embstate_accn
      );
  const SourceState embfluid_staten(
      embstaten_velnp,embstaten_veln,embstaten_velnm,embstaten_accnp,embstaten_accn
      );
  const SourceState bgfluid_staten(
      bgstaten_velnp,bgstaten_veln,bgstaten_velnm,bgstaten_accnp,bgstaten_accn
      );

  // this method can write to the target Epetra_Vectors, but the pointers are const,
  // as the TargetState is const (const pointer to non-const pointee)
  SetNewEmbStateVectors(
      embfluid_state, embfluid_staten, bgfluid_staten, aledispnp, aledispn);
}//SetNewEmbStateVectors

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewEmbStateVectors(
  const TargetState &                       embfluid_state,
  const SourceState &                       embfluid_state_n,
  const SourceState &                       bgfluid_state_n,
  const Teuchos::RCP<const Epetra_Vector> & aledispnp,
  const Teuchos::RCP<const Epetra_Vector> & aledispn
)
{
  // new position of embedded nodes
  std::vector<LINALG::Matrix<3,1> > embnodes_xyz;
  // the vector containing interpolated values for each state vector from embedded dis. For each node one interpolated value..
  // Entries 0 to 3: velnp_, Entries 4 to 7: veln_,  Entries 8 to 11: velnm_, Entries 12 to 15: accn_, Entries 16 to 19: accnp_
  std::vector<LINALG::Matrix<20,1> > interpolated_vecs;
  LINALG::Matrix<20,1> interpolatedvec(true);

  // ids of embedded nodes, that need a projection (currently, we take all)
  std::vector<int> embnodes_nohistory;

  // loop over row nodes of embedded discretization of each processor
  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);

    // Project dofs from old ALE-mesh location to new node position

    LINALG::Matrix<3,1> pos;
    std::vector<int> src_dofs(4);
    std::vector<double> mydisp(4);

    // get the current displacement
    embdis_->Dof(embnode,0,src_dofs);
    DRT::UTILS::ExtractMyValues(*aledispnp,mydisp,src_dofs);

    pos(0) = embnode->X()[0]+mydisp.at(0);
    pos(1) = embnode->X()[1]+mydisp.at(1);
    pos(2) = embnode->X()[2]+mydisp.at(2);

    embnodes_xyz.push_back(pos);
    embnodes_nohistory.push_back(embnode->Id());

    // the vector of interpolated values of all state vector (max five state vectors) for every node
    interpolated_vecs.push_back(interpolatedvec);
  }

  // build search tree for projection from embedded discretization
  SetupSearchTree(aledispn);

  // call the Round Robin Communicator
  CommunicateNodes(embnodes_xyz,interpolated_vecs,
                   embnodes_nohistory,aledispn,
                   embfluid_state,bgfluid_state_n,
                   embfluid_state_n,embdis_);
}

// ------------------------------------------------------------------------
// Find an appropriate radius for the search tree. The minimum radius is the
// max. diameter of an embedded element, scaled with a user-defined safety
// factor.
// ------------------------------------------------------------------------
template<DRT::Element::DiscretizationType distype>
void XFEM::XFluidFluidTimeIntegration::FindSearchRadius()
{
  DRT::Element* actele = embdis_->lRowElement(0);
  const DRT::Node* const* nodes = actele->Nodes();

  // problem dimension
  const unsigned int dim = DRT::UTILS::DisTypeToDim<distype>::dim;

  // we are looking for the maximum diameter of the embedded element
  // as an estimate for the search radius
  //REMARK: the selection of the embedded element for this estimate is still
  //arbitrary --> choose a sufficiently large safety factor in the input file
  double max_diameter = 0.0;

  // build connectivity matrix for every surface of the embedded element
  std::vector< std::vector<int> > connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(distype);

  //-----------------------------------------------------------------------------------
  // We have hex elements & the faces are quads:
  // the first 4 nodes in the element node numbering vector for a given surface are the
  // corner nodes (equally numbered for hex8/20/hex27), followed by the central nodes
  // in case of hex20/27; in approx. diameter estimation, mid nodes are neglected
  //-----------------------------------------------------------------------------------

  // loop over element surfaces
  for (std::vector< std::vector<int> >::const_iterator ic = connectivity.begin();
       ic != connectivity.end(); ++ ic)
  {
    // get the set of nodes (connected in sequence) for the current surface
    const std::vector<int> & surf_nodeset = *ic;

    // compute the connections 0th->2nd, 1st->3rd corner node
    for (unsigned int icn = 0; icn < 2; ++ icn)
    {
      // next but one node position in vector
      const unsigned icnn = icn+2;

      // compute the distance
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[surf_nodeset[icnn]]->X()[isd] - nodes[icn]->X()[isd];
        dist_square += dx*dx;
      }

      double dist = sqrt(dist_square);

      // new maximum?
      if (dist > max_diameter)
        max_diameter = dist;
    }
  } // done with the surface elements

  // the spatial diagonals

  const unsigned ncn_face = 4;
  for (unsigned icn = 0; icn < 1; ++ icn)
  {
    // diagonally opposite (0-6, 1-7)
    {
      const unsigned icn_opp = icn + 2 + ncn_face;
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[icn_opp]->X()[isd] - nodes[icn]->X()[isd];
        dist_square += dx*dx;
      }
      double dist = sqrt(dist_square);
      if (dist > max_diameter)
        max_diameter = dist;
    }

    // diagonally opposite (2-4, 3-5)
    {
      const unsigned icn_opp = icn + ncn_face;
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[icn_opp]->X()[isd] - nodes[icn+2]->X()[isd];
        dist_square += dx*dx;
      }
      double dist = sqrt(dist_square);
      if (dist > max_diameter)
        max_diameter = dist;
    }
  }

  // Todo: tets are not yet supported by this framework!
  minradius_ =  searchradius_fac_*max_diameter;

}//FindSearchRadius

// -------------------------------------------------------------------
// Gmsh Output
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutput()
{
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("std_enriched_map", step_, 30, 0, bgdis_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << "std/enriched/void n\" {\n";
    for (int i=0; i<bgdis_->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis_->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      std::map<int, std::vector<int> >::const_iterator iter = nodeToDof_std_n_.find(actnode->Id());
      std::map<int, std::vector<int> >::const_iterator iteren = nodeToDof_enriched_n_.find(actnode->Id());
      if (iter != nodeToDof_std_n_.end()) kind = 1;//std
      if (iteren != nodeToDof_enriched_n_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "std/enriched/void n+1\" {\n";
    for (int i=0; i<bgdis_->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis_->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      std::map<int, std::vector<int> >::const_iterator iter = nodeToDof_std_np_.find(actnode->Id());
      std::map<int, std::vector<int> >::const_iterator iteren = nodeToDof_enriched_np_.find(actnode->Id());
      if (iter != nodeToDof_std_np_.end()) kind = 1;//std
      if (iteren != nodeToDof_enriched_np_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
}
// -------------------------------------------------------------------
// Gmsh Output for interpolated-Ale FSI-Approach
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutputForInterpolateFSI( const Teuchos::RCP<const Epetra_Vector> & aledispnp,
                                                                    const Teuchos::RCP<const Epetra_Vector> & aledispn)
{
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("emb_element_node_id", 0, 0, 0, bgdis_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    {
      // draw embedded elements with associated gid at the old  position
      gmshfilecontent << "View \" " << "emb Element(old)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        std::vector<double> myolddisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispn, myolddisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + myolddisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + myolddisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + myolddisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded elements with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Element(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actual node
          LINALG::Matrix<3,1> inodepos(pelenodes[inode]->X());
          inodepos(0,0) += mydisp[0+(inode*4)];
          inodepos(1,0) += mydisp[1+(inode*4)];
          inodepos(2,0) += mydisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded nodes with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Node(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(pelenodes[inode]->X());
          inodepos(0,0) += mydisp[0+(inode*4)];
          inodepos(1,0) += mydisp[1+(inode*4)];
          inodepos(2,0) += mydisp[2+(inode*4)];
          IO::GMSH::cellWithScalarToStream(DRT::Element::point1, pelenodes[inode]->Id(), inodepos, gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    }
  }
  gmshfilecontent.close();
}

// -------------------------------------------------------------------
// Subsequent enforcement of incompressibility on background element
// patch after projection by solution of an optimization problem
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::EnforceIncompAfterProjection(
  Teuchos::RCP<GEO::CutWizard>                      wizard,         ///< cut wizard from timestep n+1
  Teuchos::RCP<GEO::CutWizard>                      wizard_n,       ///< cut wizard from timestep n
  Teuchos::RCP<Epetra_Vector> &                     bgstate_velnp,  ///< background fluid velocity at timestep n+1
  Teuchos::RCP<Epetra_Vector> &                     bgstate_veln,   ///< background fluid velocity at timestep n
  Teuchos::RCP<Epetra_Vector> &                     bgstate_velnm,  ///< background fluid velocity at timestep n-1
  const Teuchos::RCP<const LINALG::MapExtractor> &  dbcmaps         ///< background fluid dirichlet map extractor
)
{
  // create an auxiliary patch discretization of the background elements with projected nodal values
  BuildElementPatchForIncompOpt(wizard_n,wizard,dbcmaps);

  // prepare the incompressibility discretization and evaluate
  // incompressibility condition
  PrepareIncompOptDiscret(wizard);
  EvaluateIncompOpt(wizard);

  //solve the optimization problem for the state vectors
  SolveIncompOptProblem(bgstate_velnp);
  SolveIncompOptProblem(bgstate_veln);
  SolveIncompOptProblem(bgstate_velnm);
}
// -----------------------------------------------------------------------
// all cut elements at tn with full dofs are included in incompressibility
// patch. If there  are nodes in projectednodeids which still are not in
// incompressibility patch we add them afterwards. These added Elements should
// have the full dofs again.
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::BuildElementPatchForIncompOpt(
    Teuchos::RCP<GEO::CutWizard>                      wizard_n,
    Teuchos::RCP<GEO::CutWizard>                      wizard_np,
    const Teuchos::RCP<const LINALG::MapExtractor>  & dbcmaps)
{
  //---------------------------------------------
  // find the patch at time t_n

  // delete the elements and nodes of the last time step
  incompnodeids_set_.clear();
  incompelementids_set_.clear();


  // call loop over elements
  const int numele = bgdis_->NumMyRowElements();
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = bgdis_->lRowElement(i);
    int numnodes = actele->NumNode();

    // get the element handle of actele at time tn
    GEO::CUT::ElementHandle * e = wizard_n->GetElement( actele );

    const DRT::Node* const* nodesofele = actele->Nodes();

    // We are looking for elements whose all nodes have dofs at tn+1.
    // So check whether the nodes of this elements have dofs
    int numnodeswithdofs = 0;
    for(int inode = 0; inode < numnodes; inode++)
    {
      if(bgdis_->NumDof(nodesofele[inode]) > 0)
        numnodeswithdofs++;
    }

    // xfem element
    if ( e!=NULL )
    {
      // select all cut elements at t_n, which has complete dofs at t_n+1
      if (e->IsCut() and (numnodeswithdofs == numnodes))
      {
        //insert the element
        incompelementids_set_.insert(actele->Id());

        //get the nodes adjust to this element
        const int* nodeids = actele->NodeIds();

        for(int inode=0; inode<actele->NumNode(); ++inode)
        {
          const int nodegid = nodeids[inode];

          //insert the node
          incompnodeids_set_.insert(nodegid);
        }
      }
    }
  }

  //------------------------------------------------
  //  check if all projected nodes are included
  for(std::set<int>::iterator iter = projectedNodes_.begin(); iter!= projectedNodes_.end();
      iter++)
  {
    std::set<int>::const_iterator iterpatchnodes = incompnodeids_set_.find(*iter);

    if (iterpatchnodes == incompnodeids_set_.end())
    {
      IO::cout << "STOP!! Nodes found which were in projected set but not in incompressibility patch! " << *iter << IO::endl;

      //get all adjacent elements of this node
      int numberOfElements = bgdis_->gNode(*iter)->NumElement();
      const DRT::Element* const* elements = bgdis_->gNode(*iter)->Elements();

      // loop over adjacent elements of this node
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get the adjacent element
        const DRT::Element* ele_adj = elements[ele_current];

        // get vector of pointers of node (for this element)
        const DRT::Node* const* nodesofadjele = ele_adj->Nodes();
        const int numberOfNodes = ele_adj->NumNode();

        int numnodeswithdofs = 0;
        // loop nodes of this element to know whether it has complete dofs
        for(int vec_it = 0; vec_it < numberOfNodes; vec_it++)
        {
          // check if the nodes of this elements have dofs
          if(bgdis_->NumDof(nodesofadjele[vec_it]) > 0)
            numnodeswithdofs++;
        }

        // if all nodes of this elements have dofs add it to incomp patch
        if (numnodeswithdofs == numberOfNodes )
        {
          incompelementids_set_.insert(ele_adj->Id());
          IO::cout << "element found " << ele_adj->Id() << IO::endl;

          //get the nodes of this element
          const int* nodeids = ele_adj->NodeIds();

          for (int inode=0; inode<ele_adj->NumNode(); ++inode)
          {
            const int nodegid = nodeids[inode];

            //insert the node
            incompnodeids_set_.insert(nodegid);
          }
          // break the element loop, we've already found one element for this node
          break;
        }
      }
    }
  }

  //------------------------------------------
  //check it again..
  for(std::set<int>::iterator iter = projectedNodes_.begin(); iter!= projectedNodes_.end();
      iter++)
  {
    std::set<int>::const_iterator iterpatchnodes = incompnodeids_set_.find(*iter);

    if (iterpatchnodes == incompnodeids_set_.end())
      dserror("BUG!! Nodes found which were in projected set but not in incompressibility patch!",*iter);
  }

  //---------------------------------
  // Gmsh debug output
  if (gmsh_debug_out_)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("incom_patch", step_, 5, 0, bgdis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // draw bg elements with associated gid
      gmshfilecontent << "View \" " << "bg Element->Id() \" {\n";
      for (int i=0; i<bgdis_->NumMyColElements(); ++i)
      {
        DRT::Element* actele = bgdis_->lColElement(i);
//         GEO::CUT::ElementHandle * e = wizard_n->GetElement( actele );
        std::set<int>::const_iterator iter = incompelementids_set_.find(actele->Id());
        if ( iter != incompelementids_set_.end())
          IO::GMSH::elementAtInitialPositionToStream(1.0, actele, gmshfilecontent);
        else
          IO::GMSH::elementAtInitialPositionToStream(0.0, actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
  }
}

// -------------------------------------------------------------------
// build an incompressibility discretization
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PrepareIncompOptDiscret(Teuchos::RCP<GEO::CutWizard> wizard_np)
{

  // generate an empty boundary discretisation
  incompdis_ = Teuchos::rcp(new DRT::Discretization((std::string)"incompressibility discretisation",
                                           Teuchos::rcp(bgdis_->Comm().Clone())));

  std::set<int> incompelementids_set_all;
  std::set<int> incompnodeids_set_all;

   // Gather all informations from all processors
  std::vector<int> incompdofsAllproc;
  std::vector<int> incompveldofsAllproc;

  // information how many processors work at all
  std::vector<int> allproc(bgdis_->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<bgdis_->Comm().NumProc(); ++i) allproc[i] = i;

  LINALG::Gather<int>(incompelementids_set_,incompelementids_set_all,(int)bgdis_->Comm().NumProc(),&allproc[0],bgdis_->Comm());
  LINALG::Gather<int>(incompnodeids_set_,incompnodeids_set_all,(int)bgdis_->Comm().NumProc(),&allproc[0],bgdis_->Comm());


  // determine sets of col und row nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;


  // loop all column elements and label all row nodes next to a MHD node
  for (int i=0; i<bgdis_->NumMyColElements(); ++i)
  {
    DRT::Element* actele = bgdis_->lColElement(i);

    // get the node ids of this elements
    const int  numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found=false;

    std::set<int>::const_iterator iter = incompelementids_set_all.find(actele->Id());
    if ( iter != incompelementids_set_all.end()) found=true;

    if(found==true)
    {
      // loop nodeids
      for(int rr=0;rr<numnode;++rr)
      {
        int gid=nodeids[rr];

        if ((bgdis_->NodeRowMap())->LID(gid)>-1)
        {
          adjacent_row.insert(gid);
        }
        adjacent_col.insert(gid);
      }
    }
  }

  // add nodes to incompressibility discretisation
  for(std::set<int>::iterator id = adjacent_row.begin();
      id!=adjacent_row.end(); ++id)
  {
    DRT::Node* actnode=bgdis_->gNode(*id);

    Teuchos::RCP<DRT::Node> incompnode =Teuchos::rcp(actnode->Clone());

    incompdis_->AddNode(incompnode);
  }


  // loop all row elements and add all elements with a MHD node
  for (int i=0; i<bgdis_->NumMyRowElements(); ++i)
  {
    DRT::Element* actele = bgdis_->lRowElement(i);

    bool found=false;

    // check if incompressibility element
    std::set<int>::const_iterator iter = incompelementids_set_all.find(actele->Id());
    if ( iter != incompelementids_set_all.end()) found=true;

    // yes, we have a MHD condition
    if(found==true)
    {
      Teuchos::RCP<DRT::Element> incompele =Teuchos::rcp(actele->Clone());

      incompdis_->AddElement(incompele);
    }
  }

  //incompelementids_set_ needs a full NodeRowMap and a NodeColMap
  Teuchos::RCP<Epetra_Map> newrownodemap;
  Teuchos::RCP<Epetra_Map> newcolnodemap;

  std::vector<int> rownodes;

  // convert std::set to std::vector
  for(std::set<int>::iterator id = adjacent_row.begin();
      id!=adjacent_row.end();
      ++id)
  {
    rownodes.push_back(*id);
  }

  // build noderowmap for new distribution of nodes
  newrownodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     rownodes.size(),
                                     &rownodes[0],
                                     0,
                                     incompdis_->Comm()));


  std::vector<int> colnodes;
  for(std::set<int>::iterator id = adjacent_col.begin();
      id!=adjacent_col.end();
      ++id)
  {
    colnodes.push_back(*id);
  }

  // build nodecolmap for new distribution of nodes
  newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     colnodes.size(),
                                     &colnodes[0],
                                     0,
                                     incompdis_->Comm()));

  //FIXME:
  // incompdis is created from a XFEM discretization, but we don't pass
  // a XFEM dofset here - in case of multiple dofsets per node,
  // evaluation is going to crash
  incompdis_->Redistribute(*newrownodemap,*newcolnodemap,false,false,false);
  Teuchos::RCP<DRT::DofSet> new_dofset = Teuchos::rcp(new XFEM::XFEMTransparentIndependentDofSet(bgdis_,true,wizard_np));
  incompdis_->ReplaceDofSet(new_dofset); // do not call this with true!!
  incompdis_->FillComplete();
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::EvaluateIncompOpt(Teuchos::RCP<GEO::CutWizard>          wizard)
{

  C_ = LINALG::CreateVector(*incompdis_->DofRowMap(),true);

  // Problem definition:
  //
  // Find u_p*, u_p: interpolated vector
  // min || u_p* - u_p ||^2
  //
  // s.t.  _
  //      |
  //      | div u_p* dOmega = 0
  //     _|
  //
  // equal to: cT.u_p* = 0

  //---------------------------------------
  // Find the vector c:
  //
  //   __   _                  __   _                 __   _
  //   \   |  dN_A             \   |  dN_A            \   |  dN_A
  // ( /   |  ---- dOmega_e1,  /   |  ---- dOmega_e1, /   |  ---- dOmega_e1, ...)
  //   -- -   dx               -- -   dy              -- -   dz
  //   A                       A                      A

  DRT::Element::LocationArray la( 1 );

  // loop over column elements of bgdis
  const int numcolele = incompdis_->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = incompdis_->lColElement(i);

    Teuchos::RCP<MAT::Material> mat = actele->Material();
    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

    GEO::CUT::ElementHandle * e = wizard->GetElement( actele );

    Epetra_SerialDenseVector C_elevec;
    // xfem element
    if ( e!=NULL )
    {

      std::vector<GEO::CUT::plain_volumecell_set> cell_sets;
      std::vector<std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, false); //(include_inner=false)

      if (cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
      if (cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      if (cell_sets.size() == 0)
      {
        IO::cout << "Warning: Element " << actele->Id() << " has all it's volume-cells in void. Check if all nodes are "<<
          "included in other elements." << IO::endl;
        continue;
      }

      for( std::vector<GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        const int setpos = s-cell_sets.begin();
        const std::vector<int> & nds = nds_sets[setpos];


        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*incompdis_,nds,la,false);

        // number of dofs for background element
        // ndof contains all dofs of velocity and pressure (we need just the velocity)
        const size_t ndof  = la[0].lm_.size();
        C_elevec.Shape(ndof,1);


        for( unsigned cellcount=0; cellcount != s->size(); cellcount++ )
        {

          if(!has_xfem_integration_rule)
          {
            DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->CalculateContinuityXFEM(
                ele,
                *incompdis_,
                la[0].lm_,
                C_elevec);
          }
          else
          {
            // call element method
            DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->CalculateContinuityXFEM(
                ele,
                *incompdis_,
                la[0].lm_,
                C_elevec,
                intpoints_sets[setpos][cellcount]);
          }

        }

      }

    }
    else
    {
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*incompdis_,la,false);

      const size_t ndof = la[0].lm_.size();
      C_elevec.Reshape(ndof,1);

      DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                                     *incompdis_,
                                                                                                     la[0].lm_,
                                                                                                     C_elevec);

    }
    LINALG::Assemble(*C_, C_elevec, la[0].lm_, la[0].lmowner_);
  }
}

// -------------------------------------------------------------------
// Solve the optimization problem to enforce incompressibility on patch
// -------------------------------------------------------------------
void  XFEM::XFluidFluidTimeIntegration::SolveIncompOptProblem(Teuchos::RCP<Epetra_Vector>   initialvel)
{
  // ----------------------
  // Prepare the C vector:
  // The vector C still includes the pressure degrees of freedom which are
  // zero. So we need to eliminate them from C so that the vector C has just
  // the velocity degrees of freedom -> C_vel
  LINALG::MapExtractor velpressplitter;
  const unsigned numdim = 3;
  FLD::UTILS::SetupFluidSplit(*incompdis_, numdim, 1,velpressplitter);
  Teuchos::RCP<Epetra_Vector> C_vel = velpressplitter.ExtractOtherVector(C_);

  // length of C_vel on this proc
  const int mylength_cvel = C_vel->MyLength();

  // ----------------------
  // Create needed C_vel-dofmaps
  std::vector<int> C_vel_dofids(mylength_cvel);

  for (int i=0; i<mylength_cvel; ++i)
  {
    C_vel_dofids.push_back(C_vel->Map().GID(i));
  }

  // build dofrowmap for velocity dofs
  Teuchos::RCP<Epetra_Map> veldofrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                                                      C_vel_dofids.size(),
                                                                      &C_vel_dofids[0],
                                                                      0,
                                                                      incompdis_->Comm()));

  Teuchos::RCP<Epetra_Vector> C_vel_copy =  LINALG::CreateVector(*veldofrowmap,true);
  C_vel_copy->Update(1.0, *C_vel, 0.0);

  // -----------------------
  // Initialize Q as identity matrix

//  Teuchos::RCP<Epetra_CrsMatrix> Q = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,maxnumberofentries,false));
//  double ones = 1.0;
//  for(int i=0; i<C_vel->MyLength(); ++i)
//  {
//    int myGID = C_vel->Map().GID(i);
//    int err = Q->InsertGlobalValues(myGID,1,&ones,&myGID);
//    if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
//  }
//  Q->FillComplete();
//
//  Teuchos::RCP<LINALG::SparseMatrix> Q_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q));
//  Q_spr->Complete();

  Teuchos::RCP<LINALG::SparseMatrix> Q = LINALG::Eye(*veldofrowmap);

  // ------------------------
  // Transform the constraint cT.u* = 0 to (Qc)T.Qu* = 0
  //
  // With appropiate Givens Rotations we annul all of the entires of c
  // but the last entry. First we start with the first entry (next).
  // The entry "next" will be annulated through the second
  // entry (next_us). We do the same for all entries of c until all entries
  // of c are zero beside the last entry.
  //
  // next := entry to eliminate, next_us := entry to make next to zero
  //
  // The vector we find at last is Qc


  int pair = veldofrowmap->NumMyElements()-1;
  int maxpair;
  incompdis_->Comm().MaxAll(&pair, &maxpair, 1);

  int next_us = 0;
  const double TOL = 1.0e-14;
  const int maxgid = C_vel->Map().MaxMyGID();


  // the loop over all pairs of C_vel
  for (int next=0; next<maxpair; ++next)
  {
    // build the rotation matrix for current next and next_us
    Teuchos::RCP<Epetra_CrsMatrix> Q_i = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,2));
    // first build Q_i as an identity matrix
    double myval = 1.0;
    for (int j=0; j<mylength_cvel; ++j)
    {
      int myGID = C_vel->Map().GID(j);
      Q_i->InsertGlobalValues(myGID,1,&myval,&myGID);
    }

    if (next < ((C_vel->MyLength())-1))
    {
      // if the next values is not zero find next+1 and do the givens,
      // otherwiese Q_i remains an idendity matrix
      if ( abs((*C_vel)[next]) > TOL )
      {
        next_us = next+1;
        while ( abs((*C_vel)[next_us]) < TOL )
        {
          next_us++;
        }

        if ( next_us <= maxgid )
        {
          // std::cout << "next yy "<< next << " " << (*C_vel)[next] << std::endl;
          // std::cout << "next_us yy "<< next_us << " " <<  (*C_vel)[next_us] << std::endl;

          double Nenner = sqrt((*C_vel)[next]*(*C_vel)[next]+
                               (*C_vel)[next_us]*(*C_vel)[next_us]);
          Nenner = 1/Nenner;

          double cphi = Nenner*(*C_vel)[next_us];
          double sphi = Nenner*(*C_vel)[next];

          int nextGID = veldofrowmap->GID(next);
          int next_usGID = veldofrowmap->GID(next_us);

          double sphi_min = -sphi;
          Q_i->ReplaceGlobalValues(nextGID,1,&cphi,&nextGID);
          Q_i->InsertGlobalValues(nextGID,1,&sphi_min,&next_usGID);
          Q_i->InsertGlobalValues(next_usGID,1,&sphi,&nextGID);
          Q_i->ReplaceGlobalValues(next_usGID,1,&cphi,&next_usGID);
        }
      }
    }
    // do nothing or wait until all processors are finished
    else { }

    incompdis_->Comm().Barrier();
    Q_i->FillComplete(*veldofrowmap,*veldofrowmap);

    Teuchos::RCP<LINALG::SparseMatrix> Q_i_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q_i,View));
    Q_i_spr->Complete(*veldofrowmap,*veldofrowmap);

    // Update of C_vel (which is after every update Qc)
    /* _     _
       | s  -c || x_i |    | s*x_i-c*x_j |   |      0      |
       |       ||     |  = |             | = |             |
       | c   s || x_j |    | c*x_i+s*x_j |   | c*x_i+s*x_j |
       -     -
       x_i: next
       x_j: next_us
    */

    Teuchos::RCP<Epetra_Vector> C_vel_test =  LINALG::CreateVector(*veldofrowmap,true);
    (Q_i_spr->EpetraMatrix())->Multiply(false,*C_vel,*C_vel_test);
    C_vel->Update(1.0,*C_vel_test,0.0);

    // build the final Q  = Qn..Q3*Q2*Q1
    Q = Multiply(*Q_i_spr,false,*Q,false,false,false);
  }

  incompdis_->Comm().Barrier();

  //-----------------------------------------------------
  // Communicate the last entry of C_vel of each processor

  // build an allreduced vector of all last entries of C_vel
  std::vector<int> C_vellast;
  std::vector<int> C_vellast_All;

  if ((maxgid>-1) and ((*C_vel)[veldofrowmap->LID(maxgid)]>TOL))
    C_vellast.push_back(maxgid);

  //information how many processors work at all
  std::vector<int> allproc(incompdis_->Comm().NumProc());

  LINALG::Gather<int>(C_vellast,C_vellast_All,(int)incompdis_->Comm().NumProc(),
                      &allproc[0],incompdis_->Comm());

  // build dofrowmap of last entries
  Teuchos::RCP<Epetra_Map> lastdofrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                                             C_vellast_All.size(),
                                                             &C_vellast_All[0],
                                                             0,
                                                             incompdis_->Comm()));

  const Epetra_Map lastallreduced = *LINALG::AllreduceOverlappingEMap(*lastdofrowmap);
  Teuchos::RCP<Epetra_Vector> C_vel_last =  LINALG::CreateVector(lastallreduced,true);
  int mylastGID = lastallreduced.MaxMyGID();
  LINALG::Export(*C_vel,*C_vel_last);

  // std::cout << "C_vellast " << *C_vel_last << std::endl;
  // std::cout << "length " << C_vel_last->MyLength() << std::endl;

  // loop over all processors
  if (C_vel_last->MyLength() > 1)
  {
    for(int pr=0; pr<incompdis_->Comm().NumProc()-1; ++pr)
    {
      // build the rotation matrix for current next and next_us
      Teuchos::RCP<Epetra_CrsMatrix> Q_i = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,2));
      // first build Q_i as an identity matrix
      double myval = 1.0;
      for(int j=0; j<C_vel->MyLength(); ++j)
      {
        int myGID = C_vel->Map().GID(j);
        Q_i->InsertGlobalValues(myGID,1,&myval,&myGID);
      }

      if(C_vel_last->Map().MyLID(pr+1))
      {
        double next_val = (*C_vel_last)[pr];
        double nextus_val = (*C_vel_last)[pr+1];

        // std::cout << "val " << next_val << " " << nextus_val << std::endl;

        double Nenner = sqrt(next_val*next_val+nextus_val*nextus_val);
        Nenner = 1/Nenner;

        double cphi = Nenner*nextus_val;
        double sphi = Nenner*next_val;

        int nextGID = lastallreduced.GID(pr);
        int next_usGID = lastallreduced.GID(pr+1);
        mylastGID = lastallreduced.GID(pr+1);

        double sphi_min = -sphi;
        Q_i->ReplaceGlobalValues(nextGID,1,&cphi,&nextGID);
        Q_i->InsertGlobalValues(nextGID,1,&sphi_min,&next_usGID);
        Q_i->InsertGlobalValues(next_usGID,1,&sphi,&nextGID);
        Q_i->ReplaceGlobalValues(next_usGID,1,&cphi,&next_usGID);
        Q_i->FillComplete(*veldofrowmap,*veldofrowmap);


        Teuchos::RCP<LINALG::SparseMatrix> Q_i_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q_i,View));
        Q_i_spr->Complete(*veldofrowmap,*veldofrowmap);

        Teuchos::RCP<Epetra_Vector> C_vel_test =  LINALG::CreateVector(*veldofrowmap,true);
        (Q_i_spr->EpetraMatrix())->Multiply(false,*C_vel,*C_vel_test);
        C_vel->Update(1.0,*C_vel_test,0.0);

        // build the final Q  = Qn..Q3*Q2*Q1
        Q = Multiply(*Q_i_spr,false,*Q,false,false,false);
      }
    }
  }

  //------------------------------------------------------------
  // Update the initial velocity vector:
  //
  // find u* from Qu* != Qu
  // Note: From (Qc)T.Qu* = 0  we know that Qu*(next) is zero
  //
  // | Qu1 |   | Qu1    |       | Qu1 |   | Qu1    |   |    0   |
  // | Qu2 |   | Qu2    |       | Qu2 |   | Qu2    |   |    0   |
  // |  .  | = |  .     |   <=> |  .  | = |  .     | - |    0   |
  // |  .  |   |  .     |       |  .  |   |  .     |   |    0   |
  // |  0  |   | Qu_next|       |  0  |   | Qu_next|   | Qu_next|
  //
  // =>
  //               |  0     |       |x ...   next1||    0   |
  //               |  0     |       |  .     next2||    0   |
  // u* = Q'Qu - Q'|  0     | = u - |   .     .   ||    0   |
  //               |  0     |       |         .   ||    0   |
  //               | Qu_next|       |        nextn|| Qu_next|
  //

  // the original velocity vector which we want to improve
  Teuchos::RCP<Epetra_Vector> vel_org = LINALG::CreateVector(*veldofrowmap,true);
  LINALG::Export(*initialvel,*vel_org);

  // Qu
  Teuchos::RCP<Epetra_Vector> Qunext = LINALG::CreateVector(*veldofrowmap,true);
  (Q->EpetraMatrix())->Multiply(false,*vel_org,*Qunext);

  // we need the last entry of Qu (Qu_next)
  for (int i=0; i<Qunext->MyLength(); ++i)
  {
     int mygid = Qunext->Map().GID(i);
     if ( mygid != mylastGID)
       (*Qunext)[Qunext->Map().LID(mygid)] = 0.0;
  }

  Teuchos::RCP<Epetra_Vector> QTQu = LINALG::CreateVector(*veldofrowmap,true);
  (Q->EpetraMatrix())->Multiply(true,*Qunext,*QTQu);

  Teuchos::RCP<Epetra_Vector> u_incomp = LINALG::CreateVector(*veldofrowmap,true);
  u_incomp->Update(1.0, *vel_org, 0.0);
  u_incomp->Update(-1.0, *QTQu, 1.0);

  //incompressibility check before solving the optimization problem
  double sum = 0.0;
  vel_org->Dot(*C_vel_copy, &sum);

  double sum_opt = 0.0;
  u_incomp->Dot(*C_vel_copy, &sum_opt);

  if (myrank_ == 0)
  {
    IO::cout << " Incompressibility Check... " << IO::endl;
    IO::cout << " Original:  "  << sum << ",  After solving optimization problem: "  << sum_opt << IO::endl;
  }

  LINALG::Export(*(u_incomp),*(initialvel));

}

