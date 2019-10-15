/*----------------------------------------------------------------------*/
/*! \file

\brief manages the different types of mesh based coupling conditions and thereby builds the bridge
between the xfluid class and the cut-library

\level 2

\maintainer  Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_mesh.H"

#include "xfem_utils.H"
#include "xfem_interface_utils.H"
#include "xfem_discretization_utils.H"
#include "xfem_xfluid_contact_communicator.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_boundary_parent_calc.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

// Needed for Slave Solid XFSI --> should go to xfem_coulpling_mesh_fsi
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_surface.H"
#include "../drt_mat/elasthyper.H"

// Needed for Slave Fluid XFF --> should go to xfem_coulpling_mesh_ff
#include "../drt_mat/newtonianfluid.H"

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
XFEM::MeshCoupling::MeshCoupling(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,               ///< discretization from which the cutter discretization is derived
    const int coupling_id,      ///< id of composite of coupling conditions
    const double time,          ///< time
    const int step,             ///< time step
    const std::string& suffix,  ///< suffix for cutterdisname
    bool marked_geometry)
    : CouplingBase(bg_dis, cond_name, cond_dis, coupling_id, time, step),
      mark_geometry_(marked_geometry),
      h_scaling_(-1.0),
      firstoutputofrun_(true),
      suffix_(suffix)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::SetCutterDiscretization()
{
  // create a cutter discretization from conditioned nodes of the given coupling discretization
  CreateCutterDisFromCondition(suffix_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::SetConditionsToCopy()
{
  // fill list of conditions that will be copied to the new cutter discretization
  conditions_to_copy_.push_back(cond_name_);

  // additional conditions required for the new boundary conditions
  conditions_to_copy_.push_back("FSICoupling");  // for partitioned XFSI

  // additional conditions required for the displacements of the cutter mesh
  conditions_to_copy_.push_back("XFEMSurfDisplacement");

  // for Navier Slip/ Robin boundary conditions
  // this implementation should be reviewed at some point as it requires these conditions
  //  to have a couplingID. In theory this should not be necessary.
  if (cond_name_ == "XFEMSurfNavierSlip" or cond_name_ == "XFEMSurfNavierSlipTwoPhase")
  {
    conditions_to_copy_.push_back("XFEMRobinDirichletSurf");
    conditions_to_copy_.push_back("XFEMRobinNeumannSurf");
  }
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::CreateCutterDisFromCondition(std::string suffix)
{
  // create name string for new cutter discretization (e.g, "boundary_of_struct_1" or
  // "boundary_of_fluid_2")
  std::string cutterdis_name("boundary_of_");
  cutterdis_name += cond_dis_->Name() + suffix;

  std::ostringstream temp;
  temp << coupling_id_;
  cutterdis_name += "_" + temp.str();

  //--------------------------------
  // create the new cutter discretization form the conditioned coupling discretization
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  cutter_dis_ = discreator->CreateMatchingDiscretizationFromCondition(
      *cond_dis_,      ///< discretization with condition
      cond_name_,      ///< name of the condition, by which the derived discretization is identified
      cutterdis_name,  ///< name of the new discretization
      "",
      conditions_to_copy_,  ///< list of conditions that will be copied to the new discretization
      coupling_id_  ///< coupling id, only elements conditioned with this coupling id are considered
  );
  //--------------------------------


  if (cutter_dis_->NumGlobalNodes() == 0)
  {
    dserror("Empty cutter discretization detected. No coupling can be performed...");
  }

  // for parallel jobs we have to call TransparentDofSet with additional flag true
  bool parallel = cond_dis_->Comm().NumProc() > 1;
  Teuchos::RCP<DRT::DofSet> newdofset =
      Teuchos::rcp(new DRT::TransparentIndependentDofSet(cond_dis_, parallel));

  cutter_dis_->ReplaceDofSet(newdofset);  // do not call this with true!!

  // create node and element distribution with elements and nodes ghosted on all processors
  DRT::UTILS::GhostDiscretizationOnAllProcs(cutter_dis_);
  cutter_dis_->FillComplete();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::GmshOutputDiscretization(std::ostream& gmshfilecontent)
{
  // compute the current boundary position
  std::map<int, LINALG::Matrix<3, 1>> currinterfacepositions;

  // output of cutting discretization
  XFEM::UTILS::ExtractNodeVectors(cutter_dis_, currinterfacepositions, idispnp_);
  XFEM::UTILS::PrintDiscretizationToStream(cutter_dis_, cutter_dis_->Name(), true, true, true, true,
      false, false, gmshfilecontent, &currinterfacepositions);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::PrepareCutterOutput()
{
  // -------------------------------------------------------------------
  // prepare output
  // -------------------------------------------------------------------

  if (!mark_geometry_)  // Do not write for marked geometry!
  {
    cutter_dis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(cutter_dis_)));
    cutter_output_ = cutter_dis_->Writer();
    cutter_output_->WriteMesh(0, 0.0);
  }
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::Output(const int step, const double time, const bool write_restart_data)
{
  if (!mark_geometry_)  // Do not write for marked geometry!
  {
    // output for interface
    cutter_output_->NewStep(step, time);

    cutter_output_->WriteVector("ivelnp", ivelnp_);
    cutter_output_->WriteVector("idispnp", idispnp_);

    cutter_output_->WriteElementData(firstoutputofrun_);
    firstoutputofrun_ = false;
  }
  // do not write restart for general case
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::InitStateVectors()
{
  // move state vectors to extra container class!

  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();

  ivelnp_ = LINALG::CreateVector(*cutterdofrowmap, true);
  iveln_ = LINALG::CreateVector(*cutterdofrowmap, true);
  ivelnm_ = LINALG::CreateVector(*cutterdofrowmap, true);

  idispnp_ = LINALG::CreateVector(*cutterdofrowmap, true);
  idispn_ = LINALG::CreateVector(*cutterdofrowmap, true);
  idispnpi_ = LINALG::CreateVector(*cutterdofrowmap, true);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::ClearState() { cutter_dis_->ClearState(); }

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::SetState()
{
  // set general vector values of cutterdis needed by background element evaluate routine
  ClearState();

  cutter_dis_->SetState("ivelnp", ivelnp_);
  cutter_dis_->SetState("iveln", iveln_);
  cutter_dis_->SetState("idispnp", idispnp_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::SetStateDisplacement()
{
  // set general vector values of cutterdis needed by background element evaluate routine
  ClearState();

  cutter_dis_->SetState("idispnp", idispnp_);
  cutter_dis_->SetState("idispn", idispn_);
  cutter_dis_->SetState("idispnpi", idispnpi_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::UpdateStateVectors()
{
  // update velocity n-1
  ivelnm_->Update(1.0, *iveln_, 0.0);

  // update velocity n
  iveln_->Update(1.0, *ivelnp_, 0.0);

  // update displacement n
  idispn_->Update(1.0, *idispnp_, 0.0);

  // update displacement from last increment (also used for combinations of non-monolithic
  // fluidfluid and monolithic xfsi)
  idispnpi_->Update(1.0, *idispnp_, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::UpdateDisplacementIterationVectors()
{
  // update last iteration interface displacements

  // update displacement from last increment (also used for combinations of non-monolithic
  // fluidfluid and monolithic xfsi)
  idispnpi_->Update(1.0, *idispnp_, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> XFEM::MeshCoupling::GetCutterDispCol()
{
  // export cut-discretization mesh displacements
  Teuchos::RCP<Epetra_Vector> idispcol = LINALG::CreateVector(*cutter_dis_->DofColMap(), true);
  LINALG::Export(*idispnp_, *idispcol);

  return idispcol;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::GetCouplingEleLocationVector(const int sid, std::vector<int>& patchlm)
{
  std::vector<int> patchlmstride, patchlmowner;  // dummy
  return coupl_dis_->gElement(sid)->LocationVector(
      *coupl_dis_, patchlm, patchlmowner, patchlmstride);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
XFEM::MeshVolCoupling::MeshVolCoupling(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,              ///< discretization from which cutter discretization can be derived
    const int coupling_id,     ///< id of composite of coupling conditions
    const double time,         ///< time
    const int step,            ///< time step
    const std::string& suffix  ///< suffix for cutterdisname
    )
    : MeshCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step, suffix),
      init_volcoupling_(false),
      traceEstimate_eigenvalue_update_(INPAR::XFEM::Eigenvalue_update_every_iter),
      reset_step_(-1)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::Init()
{
  XFEM::MeshCoupling::Init();

  // do additional redistributino of embedded discretization and create an auxiliary dis
  if (GetAveragingStrategy() != INPAR::XFEM::Xfluid_Sided)
  {
    // Initialize Volume Coupling
    Init_VolCoupling();

    // Todo: create only for Nitsche+EVP & EOS on outer embedded elements
    CreateAuxiliaryDiscretization();

    ele_to_max_eigenvalue_ = Teuchos::rcp(new std::map<int, double>());

    traceEstimate_eigenvalue_update_ =
        DRT::INPUT::IntegralValue<INPAR::XFEM::TraceEstimate_eigenvalue_update>(
            DRT::Problem::Instance()->XFluidDynamicParams().sublist("STABILIZATION"),
            "UPDATE_EIGENVALUE_TRACE_ESTIMATE");
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::Init_VolCoupling()
{
  if (!init_volcoupling_)
  {
    // ghost coupling elements, that contribute to the cutting discretization
    RedistributeEmbeddedDiscretization();

    init_volcoupling_ = true;
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::GetCouplingEleLocationVector(const int sid, std::vector<int>& patchlm)
{
  std::vector<int> patchlmstride, patchlmowner;  // dummy
  DRT::Element* coupl_ele = GetCouplingElement(sid);
  coupl_ele->LocationVector(*coupl_dis_, patchlm, patchlmowner, patchlmstride);
  return;
}

/*--------------------------------------------------------------------------*
 * Ghost Discretization from which the cutter_dis_ was created
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::RedistributeEmbeddedDiscretization()
{
  //#ifdef DEBUG
  //  // collect conditioned nodes and compare to the overall number of nodes in
  //  // the surface discretization
  //  std::vector<DRT::Condition*> cnd;
  //  cond_dis_->GetCondition(cond_name_,cnd);
  //
  //  // get the set of ids of all xfem nodes
  //  std::set<int> cond_nodeset;
  //  {
  //    for (size_t cond = 0; cond< cnd.size(); ++ cond)
  //    {
  //      // conditioned node ids
  //      const std::vector<int>* nodeids_cnd = cnd[cond]->Nodes();
  //      for (std::vector<int>::const_iterator c = nodeids_cnd->begin();
  //           c != nodeids_cnd->end(); ++c)
  //        cond_nodeset.insert(*c);
  //    }
  //  }
  //
  //  if (cond_nodeset.size() != static_cast<size_t>(cutter_dis_->NumGlobalNodes()))
  //    dserror("Got %d %s nodes but have % dnodes in the boundary discretization created from the
  //    condition",
  //        cond_nodeset.size(), cond_name_.c_str(), cutter_dis_->NumGlobalNodes());
  //#endif

  // get gids of elements (and associated notes), that contribute to the fluid-fluid interface
  std::set<int> adj_eles_row;
  std::set<int> adj_ele_nodes_row;

  const int mypid = cond_dis_->Comm().MyPID();

  // STEP 1: Query
  // loop over nodes of cutter discretization (conditioned nodes)
  for (int icondn = 0; icondn != cutter_dis_->NodeRowMap()->NumMyElements(); ++icondn)
  {
    // get node GID
    const int cond_node_gid = cutter_dis_->NodeRowMap()->GID(icondn);

    // node from coupling discretization (is on this proc, as cutter_dis nodes are
    // a subset!)
    const DRT::Node* cond_node = cond_dis_->gNode(cond_node_gid);

    // get associated elements
    const DRT::Element* const* cond_eles = cond_node->Elements();
    const int num_cond_ele = cond_node->NumElement();

    // loop over associated elements
    for (int ie = 0; ie < num_cond_ele; ++ie)
    {
      if (cond_eles[ie]->Owner() == mypid) adj_eles_row.insert(cond_eles[ie]->Id());

      const int* node_ids = cond_eles[ie]->NodeIds();
      for (int in = 0; in < cond_eles[ie]->NumNode(); ++in)
      {
        if (cond_dis_->gNode(node_ids[in])->Owner() == mypid)
          adj_ele_nodes_row.insert(node_ids[in]);
      }
    }
  }

  // STEP 2 : ghost interface-contributing elements from coupl_dis on all proc

  // collect node & element gids from the auxiliary discetization and
  // store in vector full_{nodes;eles}, which will be appended by the standard
  // column elements/nodes of the discretization we couple with

  std::set<int> full_ele_nodes_col(adj_ele_nodes_row);
  std::set<int> full_eles_col(adj_eles_row);

  for (int in = 0; in < cond_dis_->NumMyColNodes(); in++)
  {
    full_ele_nodes_col.insert(cond_dis_->lColNode(in)->Id());
  }
  for (int ie = 0; ie < cond_dis_->NumMyColElements(); ie++)
  {
    full_eles_col.insert(cond_dis_->lColElement(ie)->Id());
  }

  // create the final column maps
  {
    LINALG::GatherAll(full_ele_nodes_col, cond_dis_->Comm());
    LINALG::GatherAll(full_eles_col, cond_dis_->Comm());

    std::vector<int> full_nodes(full_ele_nodes_col.begin(), full_ele_nodes_col.end());
    std::vector<int> full_eles(full_eles_col.begin(), full_eles_col.end());

    Teuchos::RCP<const Epetra_Map> full_nodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, full_nodes.size(), &full_nodes[0], 0, cond_dis_->Comm()));
    Teuchos::RCP<const Epetra_Map> full_elecolmap =
        Teuchos::rcp(new Epetra_Map(-1, full_eles.size(), &full_eles[0], 0, cond_dis_->Comm()));

    // redistribute nodes and elements to column (ghost) map
    cond_dis_->ExportColumnNodes(*full_nodecolmap);
    cond_dis_->ExportColumnElements(*full_elecolmap);

    cond_dis_->FillComplete(true, true, true);
  }

  // STEP 3: reconnect all parentelement pointers in the cutter_dis_ faceelements
  {
    for (int fele_lid = 0; fele_lid < cutter_dis_->NumMyColElements(); fele_lid++)
    {
      DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(
          cutter_dis_->gElement(cutter_dis_->ElementColMap()->GID(fele_lid)));
      if (!fele) dserror("Cast to FaceElement failed!");

      DRT::Element* ele = cond_dis_->gElement(fele->ParentElementId());
      if (!ele) dserror("Couldn't get Parent Element!");

      fele->SetParentMasterElement(ele, fele->FaceParentNumber());
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
double XFEM::MeshVolCoupling::Get_EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele)
{
  if (ele_to_max_eigenvalue_->find(ele->Id()) == ele_to_max_eigenvalue_->end())
    EstimateNitscheTraceMaxEigenvalue(ele);

  return ele_to_max_eigenvalue_->at(ele->Id());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::ResetEvaluatedTraceEstimates()
{
  switch (traceEstimate_eigenvalue_update_)
  {
    case INPAR::XFEM::Eigenvalue_update_every_iter:
    {
      ele_to_max_eigenvalue_->clear();
      break;
    }
    case INPAR::XFEM::Eigenvalue_update_every_timestep:
    {
      if (reset_step_ < step_)
      {
        ele_to_max_eigenvalue_->clear();
        reset_step_ = step_;
      }
      break;
    }
    case INPAR::XFEM::Eigenvalue_update_once:
    {
      break;
    }
    default:
    {
      dserror("Unknown Eigenvalue update strategy!");
      break;
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::CreateAuxiliaryDiscretization()
{
  std::string aux_coup_disname("auxiliary_coupling_");
  aux_coup_disname += cond_dis_->Name();
  aux_coup_dis_ = Teuchos::rcp(
      new DRT::Discretization(aux_coup_disname, Teuchos::rcp(cond_dis_->Comm().Clone())));

  // make the condition known to the auxiliary discretization
  // we use the same nodal ids and therefore we can just copy the conditions
  // get the set of ids of all xfem nodes
  std::vector<DRT::Condition*> xfemcnd;
  cond_dis_->GetCondition(cond_name_, xfemcnd);

  std::set<int> xfemnodeset;

  for (size_t cond = 0; cond < xfemcnd.size(); ++cond)
  {
    aux_coup_dis_->SetCondition(cond_name_, Teuchos::rcp(new DRT::Condition(*xfemcnd[cond])));
    const std::vector<int>* nodeids_cnd = xfemcnd[cond]->Nodes();
    for (std::vector<int>::const_iterator c = nodeids_cnd->begin(); c != nodeids_cnd->end(); ++c)
      xfemnodeset.insert(*c);
  }

  // determine sets of nodes next to xfem nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;

  // loop all column elements and label all row nodes next to a xfem node
  for (int i = 0; i < cond_dis_->NumMyColElements(); ++i)
  {
    DRT::Element* actele = cond_dis_->lColElement(i);

    // get the node ids of this element
    const int numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found = false;

    // loop the element's nodes, check if a xfem condition is active
    for (int n = 0; n < numnode; ++n)
    {
      const int node_gid(nodeids[n]);
      std::set<int>::iterator curr = xfemnodeset.find(node_gid);
      found = (curr != xfemnodeset.end());
      if (found) break;
    }

    if (!found) continue;

    // if at least one of the element's nodes holds a xfem condition,
    // add all node gids to the adjecent node sets
    for (int n = 0; n < numnode; ++n)
    {
      const int node_gid(nodeids[n]);
      // yes, we have a xfem condition:
      // node stored on this proc? add to the set of row nodes!
      if (coupl_dis_->NodeRowMap()->MyGID(node_gid)) adjacent_row.insert(node_gid);

      // always add to set of col nodes
      adjacent_col.insert(node_gid);
    }

    // add the element to the discretization
    if (cond_dis_->ElementRowMap()->MyGID(actele->Id()))
    {
      Teuchos::RCP<DRT::Element> bndele = Teuchos::rcp(actele->Clone());
      aux_coup_dis_->AddElement(bndele);
    }
  }  // end loop over column elements

  // all row nodes next to a xfem node are now added to the auxiliary discretization
  for (std::set<int>::iterator id = adjacent_row.begin(); id != adjacent_row.end(); ++id)
  {
    DRT::Node* actnode = cond_dis_->gNode(*id);
    Teuchos::RCP<DRT::Node> bndnode = Teuchos::rcp(actnode->Clone());
    aux_coup_dis_->AddNode(bndnode);
  }

  // build nodal row & col maps to redistribute the discretization
  Teuchos::RCP<Epetra_Map> newnoderowmap;
  Teuchos::RCP<Epetra_Map> newnodecolmap;

  {
    // copy row/col node gids to std::vector
    // (expected by Epetra_Map ctor)
    std::vector<int> rownodes(adjacent_row.begin(), adjacent_row.end());
    // build noderowmap for new distribution of nodes
    newnoderowmap =
        Teuchos::rcp(new Epetra_Map(-1, rownodes.size(), &rownodes[0], 0, aux_coup_dis_->Comm()));

    std::vector<int> colnodes(adjacent_col.begin(), adjacent_col.end());

    // build nodecolmap for new distribution of nodes
    newnodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, colnodes.size(), &colnodes[0], 0, aux_coup_dis_->Comm()));

    aux_coup_dis_->Redistribute(*newnoderowmap, *newnodecolmap, false, false, false);

    // make auxiliary discretization have the same dofs as the coupling discretization
    Teuchos::RCP<DRT::DofSet> newdofset =
        Teuchos::rcp(new DRT::TransparentIndependentDofSet(cond_dis_, true));
    aux_coup_dis_->ReplaceDofSet(newdofset,
        false);  // do not call this with true (no replacement in static dofsets intended)
    aux_coup_dis_->FillComplete(true, true, true);
  }
}

//! constructor
XFEM::MeshCouplingBC::MeshCouplingBC(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step,         ///< time step
    bool marked_geometry)
    : MeshCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step, "", marked_geometry)
{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::DoConditionSpecificSetup()
{
  // set the initial interface displacements as they are used for initial cut position at the end of
  // Xfluid::Init()
  SetInterfaceDisplacement();

  // set the interface displacements also to idispn
  idispn_->Update(1.0, *idispnp_, 0.0);

  idispnpi_->Update(1.0, *idispnp_, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
bool XFEM::MeshCouplingBC::HasMovingInterface()
{
  // get the first local col(!) node
  if (cutter_dis_->NumMyColNodes() == 0) dserror("no col node on proc %i", myrank_);

  DRT::Node* lnode = cutter_dis_->lColNode(0);

  std::vector<DRT::Condition*> mycond;
  lnode->GetCondition("XFEMSurfDisplacement", mycond);

  DRT::Condition* cond = mycond[0];

  const std::string* evaltype = cond->Get<std::string>("evaltype");

  if (*evaltype == "zero") return false;

  return true;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::EvaluateCondition(Teuchos::RCP<Epetra_Vector> ivec,
    const std::string& condname, const double time, const double dt)
{
  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < cutter_dis_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor local node
    DRT::Node* lnode = cutter_dis_->lRowNode(lnodeid);
    // the set of degrees of freedom associated with the node
    const std::vector<int> nodedofset = cutter_dis_->Dof(lnode);

    const int numdof = nodedofset.size();

    if (numdof == 0) dserror("node has no dofs");
    std::vector<DRT::Condition*> mycond;
    lnode->GetCondition(condname, mycond);

    // filter out the ones with right coupling id
    std::vector<DRT::Condition*> mycond_by_coupid;

    for (size_t i = 0; i < mycond.size(); ++i)
    {
      DRT::Condition* cond = mycond[i];
      if (cond->GetInt("label") == coupling_id_) mycond_by_coupid.push_back(cond);
    }

    // safety check for unique condition
    // actually for WDBC surf conditions are sufficient due to weak enforcement, however, ivelnp
    // nodalwise then not reasonable for displacements one should set

    const int numconds = mycond_by_coupid.size();

    if (numconds > 1)
      std::cout << " !!! WARNING: more than one condition for node with label " << coupling_id_
                << ", think about implementing line and point conditions in XFEM !!!" << std::endl;

    if (numconds == 0) dserror("no condition available!");

    DRT::Condition* cond = mycond_by_coupid[numconds - 1];  // take the last condition

    // initial value for all nodal dofs to zero
    std::vector<double> final_values(numdof, 0.0);

    if (condname == "XFEMSurfDisplacement")
      EvaluateInterfaceDisplacement(final_values, lnode, cond, time);
    else if (condname == "XFEMSurfWeakDirichlet" or condname == "XFEMRobinDirichletSurf")
      EvaluateInterfaceVelocity(final_values, lnode, cond, time, dt);
    else
      dserror("non supported condname for evaluation %s", condname.c_str());


    // set final values to vector
    for (int dof = 0; dof < numdof; ++dof)
    {
      int gid = nodedofset[dof];
      ivec->ReplaceGlobalValues(1, &final_values[dof], &gid);
    }

  }  // loop row nodes
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::EvaluateInterfaceVelocity(std::vector<double>& final_values,
    DRT::Node* node, DRT::Condition* cond, const double time, const double dt)
{
  const std::string* evaltype = cond->Get<std::string>("evaltype");


  if (*evaltype == "zero")
  {
    // take initialized vector with zero values
  }
  else if (*evaltype == "funct_interpolated")
  {
    // evaluate function at node at current time
    EvaluateFunction(final_values, node->X(), cond, time);
  }
  else if (*evaltype == "funct_gausspoint")
  {
    // do nothing, the evaluate routine is called again directly from the Gaussian point
  }
  else if (*evaltype == "displacement_1storder_wo_initfunct" or
           *evaltype == "displacement_2ndorder_wo_initfunct")
  {
    if (step_ == 0)
      return;  // do not compute velocities from displacements at the beginning and do not set

    ComputeInterfaceVelocityFromDisplacement(final_values, node, dt, evaltype);
  }
  else if (*evaltype == "displacement_1storder_with_initfunct" or
           *evaltype == "displacement_2ndorder_with_initfunct")
  {
    if (step_ == 0)  // evaluate initialization function at node at current time
    {
      EvaluateFunction(final_values, node->X(), cond, time);
    }
    else
      ComputeInterfaceVelocityFromDisplacement(final_values, node, dt, evaltype);
  }
  else
    dserror("evaltype not supported %s", evaltype->c_str());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::EvaluateInterfaceDisplacement(
    std::vector<double>& final_values, DRT::Node* node, DRT::Condition* cond, const double time)
{
  const std::string* evaltype = cond->Get<std::string>("evaltype");

  if (*evaltype == "zero")
  {
    // take initialized vector with zero values
  }
  else if (*evaltype == "funct")
  {
    // evaluate function at node at current time
    EvaluateFunction(final_values, node->X(), cond, time);
  }
  else if (*evaltype == "implementation")
  {
    // evaluate implementation
    // TODO: get the function name from the condition!!!
    std::string function_name = "ROTATING_BEAM";
    EvaluateImplementation(final_values, node->X(), cond, time, function_name);
  }
  else
    dserror("evaltype not supported %s", evaltype->c_str());
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::ComputeInterfaceVelocityFromDisplacement(
    std::vector<double>& final_values, DRT::Node* node, const double dt,
    const std::string* evaltype)
{
  if (dt < 1e-14) dserror("zero or negative time step size not allowed!!!");


  double thetaiface = 0.0;

  if (*evaltype == "displacement_1storder_wo_initfunct" or
      *evaltype == "displacement_1storder_with_initfunct")
    thetaiface = 1.0;  // for backward Euler, OST(1.0)
  else if (*evaltype == "displacement_2ndorder_wo_initfunct" or
           *evaltype == "displacement_2ndorder_with_initfunct")
    thetaiface = 0.5;  // for Crank Nicholson, OST(0.5)
  else
    dserror("not supported");


  const std::vector<int> nodedofset = cutter_dis_->Dof(node);
  const int numdof = nodedofset.size();

  // loop dofs of node
  for (int dof = 0; dof < numdof; ++dof)
  {
    int gid = nodedofset[dof];
    int lid = idispnp_->Map().LID(gid);

    const double dispnp = (*idispnp_)[lid];
    const double dispn = (*idispn_)[lid];
    const double veln = (*iveln_)[lid];

    final_values[dof] =
        1.0 / (thetaiface * dt) * (dispnp - dispn) - (1.0 - thetaiface) / thetaiface * veln;
  }  // loop dofs
}

void XFEM::MeshCouplingBC::EvaluateImplementation(std::vector<double>& final_values,
    const double* x, DRT::Condition* cond, const double time, const std::string& function_name)
{
  const int numdof = final_values.size();

  if (function_name != "ROTATING_BEAM")
    dserror("currently only the rotating beam function is available!");

  double t_1 = 1.0;         // ramp the rotation
  double t_2 = t_1 + 1.0;   // reached the constant angle velocity
  double t_3 = t_2 + 12.0;  // decrease the velocity and turn around
  double t_4 = t_3 + 2.0;   // constant negative angle velocity

  double arg = 0.0;  // prescribe a time-dependent angle

  double T = 16.0;  // time period for 2*Pi
  double angle_vel = 2. * M_PI / T;

  if (time <= t_1)
  {
    arg = 0.0;
  }
  else if (time > t_1 and time <= t_2)
  {
    arg = angle_vel / 2.0 * (time - t_1) -
          angle_vel * (t_2 - t_1) / (2.0 * M_PI) * sin(M_PI * (time - t_1) / (t_2 - t_1));
  }
  else if (time > t_2 and time <= t_3)
  {
    arg = angle_vel * (time - t_2) + M_PI / T * (t_2 - t_1);
  }
  else if (time > t_3 and time <= t_4)
  {
    arg = angle_vel * (t_4 - t_3) / (M_PI)*sin(M_PI * (time - t_3) / (t_4 - t_3)) +
          2.0 * M_PI / T * (t_3 - t_2) + M_PI / T * (t_2 - t_1);
  }
  else if (time > t_4)
  {
    arg = -angle_vel * (time - t_4) + M_PI / T * (t_2 - t_1) + 2.0 * M_PI / T * (t_3 - t_2);
  }
  else
    dserror("for that time we did not define an implemented rotation %f", time);


  // rotation with constant angle velocity around point
  LINALG::Matrix<3, 1> center(true);

  center(0) = 0.0;
  center(1) = 0.0;
  center(2) = 0.0;

  LINALG::Matrix<3, 1> diff(true);
  diff(0) = x[0] - center(0);
  diff(1) = x[1] - center(1);
  diff(2) = x[2] - center(2);

  // rotation matrix
  LINALG::Matrix<3, 3> rot(true);

  rot(0, 0) = cos(arg);
  rot(0, 1) = -sin(arg);
  rot(0, 2) = 0.0;
  rot(1, 0) = sin(arg);
  rot(1, 1) = cos(arg);
  rot(1, 2) = 0.0;
  rot(2, 0) = 0.0;
  rot(2, 1) = 0.0;
  rot(2, 2) = 1.0;

  //          double r= diff.Norm2();
  //
  //          rot.Scale(r);

  LINALG::Matrix<3, 1> x_new(true);
  LINALG::Matrix<3, 1> rotated(true);

  rotated.Multiply(rot, diff);

  x_new.Update(1.0, rotated, -1.0, diff);


  for (int dof = 0; dof < numdof; ++dof)
  {
    final_values[dof] = x_new(dof);
  }
}

/*----------------------------------------------------------------------*
 |  set interface displacement at current time             schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::SetInterfaceDisplacement()
{
  if (myrank_ == 0) IO::cout << "\t set interface displacement, time " << time_ << IO::endl;

  std::string condname = "XFEMSurfDisplacement";

  EvaluateCondition(idispnp_, condname, time_);
}

/*----------------------------------------------------------------------*
 |  set interface velocity at current time                 schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::SetInterfaceVelocity()
{
  if (myrank_ == 0) IO::cout << "\t set interface velocity, time " << time_ << IO::endl;

  EvaluateCondition(ivelnp_, cond_name_, time_, dt_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
//! constructor
XFEM::MeshCouplingWeakDirichlet::MeshCouplingWeakDirichlet(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step,         ///< time step
    bool marked_geometry    ///< is this a marked geometry mesh boundary
    )
    : MeshCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step, marked_geometry)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::DoConditionSpecificSetup()
{
  XFEM::MeshCouplingBC::DoConditionSpecificSetup();

  // set the initial interface velocity and possible initialization function
  SetInterfaceVelocity();

  // set the initial interface velocities also to iveln
  iveln_->Update(1.0, *ivelnp_, 0.0);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::EvaluateCouplingConditions(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<3, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_);

  // no interface traction to be evaluated
  itraction.Clear();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::EvaluateCouplingConditionsOldState(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<3, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_ - dt_);

  // no interface traction to be evaluated
  itraction.Clear();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  SetInterfaceDisplacement();

  // set or compute the current prescribed interface velocities, just for XFEM WDBC
  SetInterfaceVelocity();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::SetupConfigurationMap()
{
  // Configuration of Consistency Terms
  configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::UpdateConfigurationMap_GP(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row].second = full_stab;

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::DoConditionSpecificSetup()
{
  // Call Base Class
  XFEM::MeshCouplingBC::DoConditionSpecificSetup();

  // Check if Inflow Stabilisation is active
  if (!cutterele_conds_.size()) dserror("cutterele_conds_.size = 0!");
  DRT::Condition* cond = (cutterele_conds_[0]).second;
  const int inflow_stab = cond->GetInt("InflowStab") + 1;
  for (uint i = 0; i < cutterele_conds_.size(); ++i)
  {
    DRT::Condition* cond = (cutterele_conds_[i]).second;
    if (inflow_stab != cond->GetInt("InflowStab") + 1)
      dserror(
          "You want to stabilized just some of your Neumann Boundaries? - feel free to implement!");
  }

  if (inflow_stab)
  {
    std::cout << "==| MeshCouplingNeumann: Inflow Stabilization active! |==" << std::endl;
    inflow_stab_ = true;
  }
  else
    inflow_stab_ = false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::SetupConfigurationMap()
{
  if (inflow_stab_)
  {
    // Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
    configuration_map_[INPAR::XFEM::F_Adj_Col] =
        std::pair<bool, double>(true, 1.0);  //<-- IMPORTANT!: used for the constraint scaling
  }
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::UpdateConfigurationMap_GP(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  if (inflow_stab_)
  {
    // Configuration of Penalty Terms
    double veln = normal.Dot(vel_m);  // as the normal is the structural body, inflow is positive
    if (veln < 0)
    {
      configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, -density_m * veln);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF1] =
          std::pair<bool, double>(true, -density_m * normal(0));
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF2] =
          std::pair<bool, double>(true, -density_m * normal(1));
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF3] =
          std::pair<bool, double>(true, -density_m * normal(2));
    }
    else
    {
      configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(false, 0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF1] = std::pair<bool, double>(false, 0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF2] = std::pair<bool, double>(false, 0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF3] = std::pair<bool, double>(false, 0);
    }
  }

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::EvaluateCouplingConditions(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<3, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::EvaluateCouplingConditions(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<6, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::EvaluateCouplingConditionsOldState(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<3, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_ - dt_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  SetInterfaceDisplacement();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
//! constructor
XFEM::MeshCouplingNavierSlip::MeshCouplingNavierSlip(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step,         ///< time step
    bool marked_geometry    ///< is this a marked geometry mesh boundary
    )
    : MeshCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step, marked_geometry)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::DoConditionSpecificSetup()
{
  // TODO: actually the same as in weak Dirichlet!!! move to base class...

  XFEM::MeshCouplingBC::DoConditionSpecificSetup();

  // set the initial interface velocity and possible initialization function
  SetInterfaceVelocity();

  // set the initial interface velocities also to iveln
  iveln_->Update(1.0, *ivelnp_, 0.0);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::SetInterfaceVelocity()
{
  if (myrank_ == 0) IO::cout << "\t set interface velocity, time " << time_ << IO::endl;

  //  EvaluateCondition( ivelnp_, cond_name_, time_, dt_);
  EvaluateCondition(ivelnp_, "XFEMRobinDirichletSurf", time_, dt_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::EvaluateCouplingConditions(
    LINALG::Matrix<3, 1>& ivel, LINALG::Matrix<3, 1>& itraction, LINALG::Matrix<3, 3>& proj_matrix,
    const LINALG::Matrix<3, 1>& x, const LINALG::Matrix<3, 1>& normal, const DRT::Condition* cond,
    const bool& eval_dirich_at_gp,
    double& kappa_m,  ///< fluid sided weighting
    double& visc_m,   ///< fluid sided weighting
    double& visc_s    ///< slave sided dynamic viscosity
)
{
  // Create normal projection matrix.
  SetupProjectionMatrix(proj_matrix, normal);

  // help variable
  int robin_id_dirch;

  if (eval_dirich_at_gp)
  {
    // evaluate interface velocity (given by weak Dirichlet condition)
    robin_id_dirch = cond->GetInt("robin_id_dirch");
    // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
    if (!std::signbit(robin_id_dirch))
      EvaluateDirichletFunction(
          ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_);

// Safety checks
#ifdef DEBUG
    if ((conditionsmap_robin_dirch_.find(robin_id_dirch)) == conditionsmap_robin_dirch_.end())
    {
      dserror("Key was not found in this instance!! Fatal error! (conditionsmap_robin_dirch_)");
    }
#endif
  }

  // evaluate interface traction (given by Neumann condition)
  robin_id_dirch = cond->GetInt("robin_id_neumann");
  if (!std::signbit(robin_id_dirch))
  {
    // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the
    // sliplenght twice (also done in UpdateConfigurationMap_GP ... as soon as this gets relevant we
    // should merge this functions)

    // evaluate interface traction (given by Neumann condition)
    // Add this to the veljump!
    double sliplength = 0.0;
    GetSlipCoefficient(sliplength, x, cond);

    if (sliplength < 0.0) dserror("The slip length can not be negative.");

    if (sliplength != 0.0)
    {
      EvaluateNeumannFunction(
          itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_);

      double sl_visc_fac = sliplength / (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      LINALG::Matrix<3, 1> tmp_itraction(true);
      tmp_itraction.MultiplyTN(proj_matrix, itraction);
      // Project this into tangential direction!!!

      ivel.Update(sl_visc_fac, tmp_itraction, 1.0);

      itraction.Clear();
    }
  }

  if (force_tangvel_map_.find(cond->Id())->second)
  {
    LINALG::Matrix<3, 1> tmp_ivel(true);
    tmp_ivel.MultiplyTN(proj_matrix, ivel);  // apply Projection matrix from the right. (u_0 * P^t)
    ivel.Update(1.0, tmp_ivel, 0.0);
  }

// Safety checks
#ifdef DEBUG
  if (!std::signbit(robin_id_dirch))
  {
    if ((conditionsmap_robin_neumann_.find(robin_id_dirch)) == conditionsmap_robin_neumann_.end())
    {
      dserror("Key was not found in this instance!! Fatal error! (conditionsmap_robin_neumann_)");
    }
  }
  std::map<int, bool>::iterator it_bool;
  if ((it_bool = force_tangvel_map_.find(cond->Id())) == force_tangvel_map_.end())
  {
    dserror("Key was not found in this instance!! Fatal error! (force_tangvel_map_)");
  }
#endif
}

void XFEM::MeshCouplingNavierSlip::EvaluateCouplingConditionsOldState(LINALG::Matrix<3, 1>& ivel,
    LINALG::Matrix<3, 1>& itraction, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // TODO: Add parameters as in call above!!!
  //  //Create normal projection matrix.
  //  LINALG::Matrix<3,3> eye(true);
  //  for(unsigned int i =0; i<3; ++i)
  //    eye(i,i)=1;
  //  for(unsigned int i =0; i<3; ++i)
  //  {
  //    for(unsigned int j =0; j<3; ++j)
  //    {
  //      proj_matrix(i,j)          = eye(i,j) - normal(i,0) * normal(j,0);
  //    }
  //  }

  // evaluate interface velocity (given by weak Dirichlet condition)
  int robin_id_dirch = cond->GetInt("robin_id_dirch");
  // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
  if (!std::signbit(robin_id_dirch))
    EvaluateDirichletFunction(
        ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_ - dt_);

  // evaluate interface traction (given by Neumann condition)
  robin_id_dirch = cond->GetInt("robin_id_neumann");
  if (!std::signbit(robin_id_dirch))
    EvaluateNeumannFunction(
        itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_ - dt_);
}

void XFEM::MeshCouplingNavierSlip::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  SetInterfaceDisplacement();

  //  // set the initial interface velocity and possible initialization function
  //  SetInterfaceVelocity();
  if (myrank_ == 0) IO::cout << "\t set interface velocity, time " << time_ << IO::endl;
  EvaluateCondition(ivelnp_, "XFEMRobinDirichletSurf", time_, dt_);
}

void XFEM::MeshCouplingNavierSlip::GetSlipCoefficient(
    double& slipcoeff, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // Extract correct slip length - bool pair for this condition ID.
  std::pair<double, bool>& tmp_pair = sliplength_map_.find(cond->Id())->second;

  if (tmp_pair.second)  // Is slip length constant?
    slipcoeff = tmp_pair.first;
  else  // Otherwise, evaluate function at gausspoint
    EvaluateScalarFunction(slipcoeff, x.A(), tmp_pair.first, cond, time_);
}

void XFEM::MeshCouplingNavierSlip::CreateRobinIdMap(
    const std::vector<DRT::Condition*>& conditions_NS,
    const std::vector<DRT::Condition*>& conditions_robin, const std::string& robin_id_name,
    std::map<int, DRT::Condition*>& conditionsmap_robin)
{
  // Loop over all Navier Slip conditions
  for (unsigned i = 0; i < conditions_NS.size(); ++i)
  {
    // Extract its robin id (either dirichlet or neumann)
    const int tmp_robin_id = conditions_NS[i]->GetInt(robin_id_name);

    // Is this robin id active? I.e. is it not 0 or negative?
    if (!(tmp_robin_id < 0))
    {
      std::vector<DRT::Condition*> mynewcond;
      GetConditionByRobinId(conditions_robin, tmp_robin_id, mynewcond);

      // The robin id should be unique. I.e. For one Coupling ID only a robin id can only exist
      // once.
      if (mynewcond.size() == 1)
      {
        if (!conditionsmap_robin.insert(std::make_pair(tmp_robin_id, mynewcond[0])).second)
          dserror("ID already existing! For conditionsmap_robin.");
      }
      else
      {
        dserror(
            "Active Robin Dirichlet/Neumann condition provided in Navier Slip condition, but not "
            "in input file!!");
      }
    }
  }

  // Is the created map the same size as the existing provided robin conditions?
  if ((conditionsmap_robin).size() != conditions_robin.size())
  {
    std::cout << "conditionsmap_robin.size(): " << (conditionsmap_robin).size() << std::endl;
    std::cout << "conditions_robin.size(): " << conditions_robin.size() << std::endl;
    dserror("More Robin Dirichlet/Neumann Conditions provided than necessary!");
  }
}

void XFEM::MeshCouplingNavierSlip::SetConditionSpecificParameters()
{
  // Build necessary maps to limit getting integers and strings on Gausspoint level.

  // Get conditions based on cutter discretization.
  std::vector<DRT::Condition*> conditions_dirich;
  cutter_dis_->GetCondition("XFEMRobinDirichletSurf", conditions_dirich);

  std::vector<DRT::Condition*> conditions_neumann;
  cutter_dis_->GetCondition("XFEMRobinNeumannSurf", conditions_neumann);

  if (conditions_neumann.size())
  {
    std::cout << "#################################################################################"
                 "########################\n";
    std::cout << "#################################################################################"
                 "########################\n";
    std::cout << "### WARNING:: XFEM::LevelSetCouplingNavierSlip                              The "
                 "traction jump is      ###\n";
    std::cout << "### divided by the dynviscosity on Gausspoint Level, this might be expensed and "
                 "not really necessary! ###\n";
    std::cout << "#################################################################################"
                 "########################\n";
    std::cout << "#################################################################################"
                 "########################"
              << std::endl;
  }

  std::vector<DRT::Condition*> conditions_NS;
  cutter_dis_->GetCondition(cond_name_, conditions_NS);

  // Establishes unique connection between Navier Slip section and Robin Dirichlet Neumann sections
  CreateRobinIdMap(conditions_NS, conditions_dirich, "robin_id_dirch", conditionsmap_robin_dirch_);

  CreateRobinIdMap(
      conditions_NS, conditions_neumann, "robin_id_neumann", conditionsmap_robin_neumann_);

  // Create maps for easy extraction at gausspoint level
  for (std::vector<DRT::Condition*>::iterator i = conditions_NS.begin(); i != conditions_NS.end();
       ++i)
  {
    int cond_int = (*i)->Id();

    double sliplength = (*i)->GetDouble("slipcoeff");

    // Is the slip length constant? Don't call functions at GP-level unnecessary.
    bool slip_bool = ((*i)->GetInt("funct") < 1);

    bool force_tangential = (((*i)->GetInt("force_tang_vel")) == 1);

    if (!sliplength_map_.insert(std::make_pair(cond_int, std::make_pair(sliplength, slip_bool)))
             .second)
      dserror("ID already existing! For sliplength_map_.");

    if (!force_tangvel_map_.insert(std::make_pair(cond_int, force_tangential)).second)
      dserror("ID already existing! For force_tangvel_map_.");
  }

  // Check if eval-type is same in Navier slip section and
  //       Robin Dirichlet section (Safety check! (not beautiful structure but could be worse..))
  for (std::vector<DRT::Condition*>::iterator i = conditions_NS.begin(); i != conditions_NS.end();
       ++i)
  {
    DRT::Condition* tmp_cond = *i;

    const int tmp_robin_id = tmp_cond->GetInt("robin_id_dirch");
    if (!std::signbit(tmp_robin_id))
    {
      if ((*conditionsmap_robin_dirch_.find(tmp_robin_id)->second->Get<std::string>("evaltype"))
              .compare(*(tmp_cond->Get<std::string>("evaltype"))) != 0)
        dserror("Not same function to evaluate in Dirichlet cond as in Main Cond.");
    }
  }
}

void XFEM::MeshCouplingNavierSlip::GetConditionByRobinId(const std::vector<DRT::Condition*>& mycond,
    const int coupling_id, std::vector<DRT::Condition*>& mynewcond)
{
  mynewcond.clear();

  // select the conditions with specified "couplingID"
  for (size_t i = 0; i < mycond.size(); ++i)
  {
    DRT::Condition* cond = mycond[i];
    const int id = cond->GetInt("robin_id");

    if (id == coupling_id) mynewcond.push_back(cond);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::SetupConfigurationMap()
{
  // Configuration of Consistency Terms
  configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::UpdateConfigurationMap_GP(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
  double sliplength = 0.0;
  GetSlipCoefficient(sliplength, x, cond);

  if (sliplength < 0.0) dserror("The slip length can not be negative.");

  if (sliplength != 0.0)
  {
    double stabnit = 0.0;
    double stabadj = 0.0;
    XFEM::UTILS::GetNavierSlipStabilizationParameters(
        visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
    configuration_map_[INPAR::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[INPAR::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);
  }
  else
  {
    configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = visc_stab_tang;
    configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool, double>(false, 0.0);
    configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool, double>(false, 0.0);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = 1.0;
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
  }

  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_n_Row].second =
      visc_stab_tang;  // full_stab <-- to keep results!

  return;
}

//! constructor
XFEM::MeshCouplingFSI::MeshCouplingFSI(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step          ///< time step
    )
    : MeshVolCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step),
      timefac_(-1.0),
      interfacelaw_(INPAR::XFEM::noslip)
{
}



/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::InitStateVectors()
{
  XFEM::MeshCoupling::InitStateVectors();

  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();
  const Epetra_Map* cutterdofcolmap = cutter_dis_->DofColMap();

  itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap, true);
  iforcecol_ = LINALG::CreateVector(*cutterdofcolmap, true);
}


void XFEM::MeshCouplingFSI::CompleteStateVectors()
{
  //-------------------------------------------------------------------------------
  // finalize itrueresidual vector

  // need to export the interface forces
  Epetra_Vector iforce_tmp(itrueresidual_->Map(), true);
  Epetra_Export exporter_iforce(iforcecol_->Map(), iforce_tmp.Map());
  int err1 = iforce_tmp.Export(*iforcecol_, exporter_iforce, Add);
  if (err1) dserror("Export using exporter returned err=%d", err1);

  // scale the interface trueresidual with -1.0 to get the forces acting on structural side (no
  // residual-scaling!)
  itrueresidual_->Update(-1.0, iforce_tmp, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::ZeroStateVectors_FSI()
{
  itrueresidual_->PutScalar(0.0);
  iforcecol_->PutScalar(0.0);
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFSI::ReadRestart(const int step)
{
  if (myrank_) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");
  //  const int    step = boundaryreader.ReadInt("step");

  if (myrank_ == 0)
  {
    IO::cout << "time: " << time << IO::endl;
    IO::cout << "step: " << step << IO::endl;
  }

  boundaryreader.ReadVector(iveln_, "iveln_res");
  boundaryreader.ReadVector(idispn_, "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_, "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");
  boundaryreader.ReadVector(idispnpi_, "idispnpi_res");

  if (not(cutter_dis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnpi_->Map()))
    dserror("Global dof numbering in maps does not match");
}

void XFEM::MeshCouplingFSI::GetSlipCoefficient(
    double& slipcoeff, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond)
{
  // Extract correct slip length - bool pair for this condition ID.
  std::pair<double, bool>& tmp_pair = sliplength_map_.find(cond->Id())->second;

  if (tmp_pair.second)  // Is slip length constant?
    slipcoeff = tmp_pair.first;
  else  // Otherwise, evaluate function at gausspoint
    EvaluateScalarFunction(slipcoeff, x.A(), tmp_pair.first, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::GmshOutput(const std::string& filename_base, const int step,
    const int gmsh_step_diff, const bool gmsh_debug_out_screen)
{
  std::ostringstream filename_base_fsi;
  filename_base_fsi << filename_base << "_force";

  // compute the current boundary position
  std::map<int, LINALG::Matrix<3, 1>> currinterfacepositions;
  XFEM::UTILS::ExtractNodeVectors(cutter_dis_, currinterfacepositions, idispnp_);


  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      filename_base_fsi.str(), step, gmsh_step_diff, gmsh_debug_out_screen, myrank_);

  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "iforce \" {" << std::endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, itrueresidual_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "idispnp \" {" << std::endl;
    // draw vector field 'idispnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, idispnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "ivelnp \" {" << std::endl;
    // draw vector field 'ivelnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, ivelnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::GmshOutputDiscretization(std::ostream& gmshfilecontent)
{
  // print surface discretization
  XFEM::MeshCoupling::GmshOutputDiscretization(gmshfilecontent);

  // compute the current solid and boundary position
  std::map<int, LINALG::Matrix<3, 1>> currsolidpositions;

  // write dis with zero solid displacements here!
  Teuchos::RCP<Epetra_Vector> solid_dispnp = LINALG::CreateVector(*cond_dis_->DofRowMap(), true);

  XFEM::UTILS::ExtractNodeVectors(cond_dis_, currsolidpositions, solid_dispnp);

  XFEM::UTILS::PrintDiscretizationToStream(cond_dis_, cond_dis_->Name(), true, false, true, false,
      false, false, gmshfilecontent, &currsolidpositions);
}

void XFEM::MeshCouplingFSI::Output(const int step, const double time, const bool write_restart_data)
{
  // output for interface
  cutter_output_->NewStep(step, time);

  cutter_output_->WriteVector("ivelnp", ivelnp_);
  cutter_output_->WriteVector("idispnp", idispnp_);
  cutter_output_->WriteVector("itrueresnp", itrueresidual_);

  cutter_output_->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->WriteVector("iveln_res", iveln_);
    cutter_output_->WriteVector("idispn_res", idispn_);
    cutter_output_->WriteVector("ivelnp_res", ivelnp_);
    cutter_output_->WriteVector("idispnp_res", idispnp_);
    cutter_output_->WriteVector("idispnpi_res", idispnpi_);
  }
}

void XFEM::MeshCouplingFSI::SetConditionSpecificParameters()
{
  std::vector<DRT::Condition*> conditions_XFSI;
  cutter_dis_->GetCondition(cond_name_, conditions_XFSI);

  // Create maps for easy extraction at gausspoint level
  for (std::vector<DRT::Condition*>::iterator i = conditions_XFSI.begin();
       i != conditions_XFSI.end(); ++i)
  {
    int cond_int = (*i)->Id();

    double sliplength = (*i)->GetDouble("slipcoeff");

    // Is the slip length constant? Don't call functions at GP-level unnecessary.
    bool slip_bool = ((*i)->GetInt("funct") < 1);

    if (!sliplength_map_.insert(std::make_pair(cond_int, std::make_pair(sliplength, slip_bool)))
             .second)
      dserror("ID already existing! For sliplength_map_.");

    INPAR::XFEM::InterfaceLaw interfacelaw =
        static_cast<INPAR::XFEM::InterfaceLaw>((*i)->GetInt("INTLAW"));
    if (i != conditions_XFSI.begin())
    {
      if (interfacelaw_ != interfacelaw)
        dserror(
            "XFEM::MeshCouplingFSI::SetConditionSpecificParameters: You defined two different FSI "
            "INTLAWS, not supported yet!");
    }
    interfacelaw_ = interfacelaw;
  }

  if (interfacelaw_ == INPAR::XFEM::navierslip_contact)  // compute h
  {
    double hmax = 0.0;
    for (int ele = 0; ele < bg_dis_->NumMyRowElements(); ++ele)
    {
      DRT::Element* fluid_ele = bg_dis_->lRowElement(ele);
      if (fluid_ele->Shape() == DRT::Element::hex8)
      {
        LINALG::Matrix<3, 8> xyze(true);
        GEO::fillInitialPositionArray(fluid_ele, xyze);
        double vol = XFEM::UTILS::EvalElementVolume<DRT::Element::hex8>(xyze);
        hmax = std::max(hmax, XFEM::UTILS::ComputeVolEqDiameter(vol));
      }
      else
        dserror("Element type != hex8, add it here!");
    }
    bg_dis_->Comm().MaxAll(&hmax, &h_scaling_, 1);
    std::cout << "==| XFEM::MeshCouplingFSI: Computed h_scaling for fluidele is: " << h_scaling_
              << "(Proc: " << bg_dis_->Comm().MyPID() << ")! |==" << std::endl;
  }

  std::cout << "==| XFEM::MeshCouplingFSI: Applied interface law is";
  switch (interfacelaw_)
  {
    case INPAR::XFEM::noslip:
    {
      std::cout << " no-slip! |==" << std::endl;
      break;
    }
    case INPAR::XFEM::noslip_splitpen:
    {
      std::cout << " no-slip with splitted normal and tangential penalty contribution! |=="
                << std::endl;
      break;
    }
    case INPAR::XFEM::slip:
    {
      std::cout << " slip! |==" << std::endl;
      break;
    }
    case INPAR::XFEM::navierslip:
    {
      std::cout << " Navier-slip! |==" << std::endl;
      break;
    }
    case INPAR::XFEM::navierslip_contact:
    {
      std::cout << " Navier-slip with Nitsche contact! |==" << std::endl;
      break;
    }
  }

  // Checks
  // if (interfacelaw_ != INPAR::XFEM::slip && interfacelaw_ != INPAR::XFEM::noslip && interfacelaw_
  // != INPAR::XFEM::noslip_splitpen) dserror("Interface law not implemented!");
}

//----------------------------------------------------------------------
// LiftDrag                                                  chfoe 11/07
//----------------------------------------------------------------------
// calculate lift&drag forces
//
// Lift and drag forces are based upon the right hand side true-residual entities
// of the corresponding nodes. The contribution of the end node of a line is entirely
// added to a present L&D force.
/*----------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::LiftDrag(const int step, const double time) const
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> iforcecol =
      DRT::UTILS::GetColVersionOfRowVector(cutter_dis_, itrueresidual_);

  if (myrank_ == 0)
  {
    // compute force components
    const int nsd = 3;
    const Epetra_Map* dofcolmap = cutter_dis_->DofColMap();
    LINALG::Matrix<3, 1> c(true);
    for (int inode = 0; inode < cutter_dis_->NumMyColNodes(); ++inode)
    {
      const DRT::Node* node = cutter_dis_->lColNode(inode);
      const std::vector<int> dof = cutter_dis_->Dof(node);
      for (int isd = 0; isd < nsd; ++isd)
      {
        // [// minus to get correct sign of lift and drag (force acting on the body) ]
        c(isd) += (*iforcecol)[dofcolmap->LID(dof[isd])];
      }
    }

    // print to file
    std::ostringstream s;
    std::ostringstream header;

    header << std::left << std::setw(10) << "Time" << std::right << std::setw(16) << "F_x"
           << std::right << std::setw(16) << "F_y" << std::right << std::setw(16) << "F_z";
    s << std::left << std::setw(10) << std::scientific << time << std::right << std::setw(16)
      << std::scientific << c(0) << std::right << std::setw(16) << std::scientific << c(1)
      << std::right << std::setw(16) << std::scientific << c(2);

    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName() +
                              ".liftdrag." + cond_name_ + ".txt";
    if (step <= 1)
    {
      f.open(fname.c_str(), std::fstream::trunc);
      f << header.str() << std::endl;
    }
    else
    {
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
    }
    f << s.str() << "\n";
    f.close();

    std::cout << header.str() << std::endl << s.str() << std::endl;
  }
}

/*--------------------------------------------------------------------------*
 * first version without possibility to use Navier Slip
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::SetupConfigurationMap()
{
  if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
  {
    if (GetInterfaceLaw() == INPAR::XFEM::slip)
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::F_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::noslip ||
             GetInterfaceLaw() == INPAR::XFEM::noslip_splitpen)
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

      if (GetInterfaceLaw() == INPAR::XFEM::noslip)
      {
        // Configuration of Penalty Terms
        configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);
      }
      else
      {
        // Configuration of Penalty Terms
        configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);

        // to guarantee correct scaling, terms are not evaluted
        configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool, double>(false, 1.0);
      }
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::navierslip ||
             GetInterfaceLaw() == INPAR::XFEM::navierslip_contact)
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_t_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
    }
    else
      dserror("Intlaw not available!");
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided)
  {
    if (GetInterfaceLaw() == INPAR::XFEM::slip)
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::F_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::XS_Con_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::XS_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::noslip)
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::XS_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::XS_Adj_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);
    }
    else
      dserror("Intlaw not available!");
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::invalid)
    dserror("XFEM::MeshCouplingFSI: Averaging Strategy not set!");
  else
    dserror("XFEM::MeshCouplingFSI: You want to initialize another strategy than Xfluid_Sided?");
  return;
}

/*--------------------------------------------------------------------------*
 * first version without possibility to use Navier Slip
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef DEBUG
  if ((kappa_m != 1 && GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided) ||
      (kappa_m != 0 && GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided))
    dserror("XFEM::MeshCouplingFSI::UpdateConfigurationMap_GP: kappa_m == %f", kappa_m);
#endif

  if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
  {
    if (GetInterfaceLaw() == INPAR::XFEM::slip)
    {
      // Configuration of Penalty Terms
      // configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::noslip)
    {
      // configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool, double>(true, full_stab);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::noslip_splitpen)
    {
      // configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
      configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = visc_stab_tang;
      configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, visc_stab_tang);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::navierslip)
    {
      double sliplength = 0.0;
      GetSlipCoefficient(sliplength, x, cond);

      if (sliplength < 0.0) dserror("The slip should not be negative!");

      double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      double stabnit = 0.0;
      double stabadj = 0.0;
      XFEM::UTILS::GetNavierSlipStabilizationParameters(
          visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);

      configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
      configuration_map_[INPAR::XFEM::X_Pen_t_Row].second = stabnit;
      configuration_map_[INPAR::XFEM::F_Con_t_Col] =
          std::pair<bool, double>(true, sliplength / dynvisc);
      configuration_map_[INPAR::XFEM::F_Con_t_Row] =
          std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
      configuration_map_[INPAR::XFEM::X_Con_t_Row] =
          std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
      configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
      configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row].second =
          full_stab;  // full_stab <-- to keep results!
      configuration_map_[INPAR::XFEM::X_Pen_n_Row].second =
          full_stab;  // full_stab <-- to keep results!
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::navierslip_contact)
    {
      UpdateConfigurationMap_GP_Contact(kappa_m, visc_m, visc_s, density_m, visc_stab_tang,
          full_stab, x, cond, ele, bele, funct, derxy, rst_slave, normal, vel_m, fulltraction);
    }
    else
      dserror("Intlaw not available!");
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided)
  {
    if (GetInterfaceLaw() == INPAR::XFEM::slip)
    {
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[INPAR::XFEM::XS_Adj_n_Row] = std::pair<bool,double>(true,-1.0);
    }
    else if (GetInterfaceLaw() == INPAR::XFEM::noslip)
    {
      configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[INPAR::XFEM::XS_Adj_Row] = std::pair<bool,double>(true,-1.0);
    }
    else
      dserror("Intlaw not available!");
  }
  return;
}

/*--------------------------------------------------------------------------*
 * UpdateConfigurationMap_GP_Contact for XFSCI
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::UpdateConfigurationMap_GP_Contact(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef DEBUG
  dsassert(xf_c_comm_ != Teuchos::null,
      "UpdateConfigurationMap_GP_Contact but no Xfluid Contact Communicator assigned!");
#endif

  static const double MAX_sliplength = 1e40;  // large number for slip case
  static const int MAX_h = 1;  // distance from contact zone at which no-slip is prescribed
  static const int MIN_h = 0;  // distance from contact zone at which full-slip is prescribed
  static const double scaling = 1. / (MAX_h - MIN_h);

  LINALG::Matrix<2, 1> xsi(rst_slave.A(), true);  // 3-->2
  bool pure_fsi = true;                           // do we integrate only fsi
  double gap =
      MAX_h * h_scaling_;  // initialize with large value as this should be the default value ...
  pure_fsi = xf_c_comm_->Get_Contact_State(
      bele->Id(), GetName(), xsi, *fulltraction, gap);  // get gap and if contact is integrated


  double sliplength = 0.0;
  GetSlipCoefficient(sliplength, x, cond);  // get reference slip length

  if ((gap - MIN_h * h_scaling_) * (MAX_sliplength + scaling) <
      h_scaling_)  // larger than maximal allows sliplength
  {
    sliplength *= h_scaling_ * MAX_sliplength;
  }
  else if (gap > MAX_h * h_scaling_)  // no-slip case
  {
    sliplength = 0.;
  }
  else  // scaling scase
  {
    sliplength *= h_scaling_ * (h_scaling_ / (gap - MIN_h * h_scaling_) - scaling);
  }

  if (sliplength < 0.0) dserror("The sliplength should not be negative!");

  double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
  double stabnit = 0.0;
  double stabadj = 0.0;
  XFEM::UTILS::GetNavierSlipStabilizationParameters(
      visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);  // sliplength is input for this

#ifdef WRITE_GMSH
  {
    xf_c_comm_->Gmsh_Write(x, *fulltraction, 1);
    xf_c_comm_->Gmsh_Write(x, (double)pure_fsi, 3);
    xf_c_comm_->Gmsh_Write(x, sliplength, 6);
  }
#endif

  if (pure_fsi)  // standard FSI with gernal Navier-slip --> Case I
  {
    xf_c_comm_->Inc_GP(3);
    configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[INPAR::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[INPAR::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[INPAR::XFEM::X_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, stabadj);
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

    // Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
  }
  else  // we evaluate only the fluid terms here --> Case III (see Remark 5) -> the fluid terms, the
        // solid contributions are added in the contact part
  {
    xf_c_comm_->Inc_GP(4);
    configuration_map_[INPAR::XFEM::X_Con_Row] =
        std::pair<bool, double>(false, 0.0);  // no solid consistency
    configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[INPAR::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[INPAR::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    if (sliplength > 1e-40)                       // this is basically the no-slip case ...
      configuration_map_[INPAR::XFEM::X_Con_t_Row] = std::pair<bool, double>(
          true, -stabnit + dynvisc / sliplength);  //+sign for penalty!+tangcons
    else                                           // avoid to evaluate this term ...
      configuration_map_[INPAR::XFEM::X_Con_t_Row] =
          std::pair<bool, double>(false, 0);  //+sign for penalty!+tancons
    configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, stabadj);
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

    // Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    configuration_map_[INPAR::XFEM::X_Pen_n_Row] =
        std::pair<bool, double>(false, 0.0);  // no normal penalty
  }
}

/*--------------------------------------------------------------------------*
 * Evaluate the Structural Cauchy Stress Matrix and it's linearization with respect to the
 *displacements
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::EvaluateStructuralCauchyStress(DRT::Element* coupl_ele,
    LINALG::Matrix<3, 1>& rst_slave, std::vector<double>& eledisp,
    const LINALG::Matrix<3, 1>& normal, std::vector<Epetra_SerialDenseMatrix>& solid_stress)
{
  if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided) return;

  if (coupl_ele->Shape() == DRT::Element::hex8)
  {
    DRT::ELEMENTS::So_hex8* solid_ele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(coupl_ele);
    if (solid_ele == NULL)
      dserror(
          "XFEM::MeshCouplingFSI::EvaluateStructuralCauchyStress: Cast of coupl_ele to solid_ele "
          "failed!");

    solid_stress.resize(5);  // traction,dtdd,d2dddx,d2dddy,d2dddz
    solid_stress[0].Reshape(NUMDIM_SOH8, 1);
    LINALG::Matrix<NUMDIM_SOH8, 1> traction(solid_stress[0].A(), true);

    solid_stress[1].Reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> dtraction_dd_l(solid_stress[1].A(), true);

    traction.Clear();

    static Epetra_SerialDenseMatrix dtraction_dd_i;
    for (int i = 0; i < NUMDIM_SOH8; ++i)
    {
      LINALG::Matrix<NUMDIM_SOH8, 1> ei(true);
      ei(i, 0) = 1.;
      solid_ele->GetCauchyAtXi(rst_slave, eledisp, normal, ei, traction(i, 0), &dtraction_dd_i,
          &solid_stress[2 + i], NULL, NULL, NULL, NULL, NULL, NULL);
      LINALG::Matrix<NUMDOF_SOH8, 1> dtraction_dd_i_l(dtraction_dd_i.A(), true);
      for (int col = 0; col < NUMDOF_SOH8; ++col) dtraction_dd_l(col, i) = dtraction_dd_i_l(col, 0);
    }
    if (timefac_ > 0)
    {
      // Change from linearization w.r.t. displacements to linearization w.r.t. velocities
      // (All other linearizations on the Nitsche Interface are evaluated like this)
      solid_stress[1].Scale(timefac_);
      for (int idx = 2; idx < 5; ++idx) solid_stress[idx].Scale(-timefac_ * timefac_);
    }
    else
      dserror("XFEM::MeshCouplingFSI::EvaluateStructuralCauchyStress: timefac = %f, not set!",
          timefac_);
  }
  else
    dserror(
        "XFEM::MeshCouplingFSI::EvaluateStructuralCauchyStress:: Element type not implemented "
        "yet!");
  return;
}

/*--------------------------------------------------------------------------*
 * get stress tangent of the slave solid
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::GetStressTangentSlave(DRT::Element* coup_ele,  ///< solid ele
    double& e_s)  ///< stress tangent slavesided
{
  //  if (coup_ele->Material()->MaterialType() == INPAR::MAT::m_elasthyper)
  //    e_s = Teuchos::rcp_dynamic_cast<MAT::ElastHyper>(coup_ele->Material())->GetYoung();
  //  else
  //    dserror("GetCouplingSpecificAverageWeights: Slave Material not a Elasthyper material?");

  // this is a temporal hack as we calculate "E/h" directely with the generalized eigenvalue problem
  // ... need to work on the input section to clarify this ...
  e_s = timefac_;

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele)
{
  DRT::ELEMENTS::StructuralSurface* solidfaceele =
      dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(ele);
  dsassert(solidfaceele != NULL, "Cast to StructuralSurface failed!");

  solidfaceele->SetParentMasterElement(
      coupl_dis_->gElement(solidfaceele->ParentElementId()), solidfaceele->FaceParentNumber());

  DRT::Element::LocationArray la(1);
  solidfaceele->ParentElement()->LocationVector(*coupl_dis_, la, false);

  // extract eledisp here
  // parent and boundary displacement at n+1
  std::vector<double> eledisp((la[0].lm_).size());
  Teuchos::RCP<const Epetra_Vector> dispnp = coupl_dis_->GetState("dispnp");
  if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

  DRT::UTILS::ExtractMyValues(*dispnp, eledisp, la[0].lm_);
  (*ele_to_max_eigenvalue_)[ele->Id()] = solidfaceele->EstimateNitscheTraceMaxEigenvalueCombined(
      eledisp);  // this is (E/h) ...basically :-)
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::PrepareSolve()
{
  // Call Base Class
  XFEM::MeshCoupling::PrepareSolve();

  // Estimate Nitsche Trace Max Eigenvalue
  if (GetAveragingStrategy() != INPAR::XFEM::Xfluid_Sided) ResetEvaluatedTraceEstimates();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::RegisterSideProc(int sid)
{
  if (GetInterfaceLaw() == INPAR::XFEM::navierslip_contact)
    Get_Contact_Comm()->RegisterSideProc(sid);
  return;
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
bool XFEM::MeshCouplingFSI::InitializeFluidState(Teuchos::RCP<GEO::CutWizard> cutwizard,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    Teuchos::RCP<XFEM::ConditionManager> condition_manager,
    Teuchos::RCP<Teuchos::ParameterList> fluidparams)
{
  if (GetInterfaceLaw() == INPAR::XFEM::navierslip_contact)
    Get_Contact_Comm()->InitializeFluidState(cutwizard, fluiddis, condition_manager, fluidparams);
  return (GetInterfaceLaw() == INPAR::XFEM::navierslip_contact);
}

XFEM::MeshCouplingFluidFluid::MeshCouplingFluidFluid(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived,
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step          ///< time step
    )
    : MeshVolCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step),
      moving_interface_(false)
{
}


/*--------------------------------------------------------------------------*
 * this function should go finally!
 * (evaluates materials on the xfluid element for slave side which basically is wrong)
 * doesn't matter if you have the same material on both sides ...
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::GetInterfaceSlaveMaterial(
    DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele, mat, GEO::CUT::Point::outside);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::RedistributeForErrorCalculation()
{
  if (GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided ||
      GetAveragingStrategy() == INPAR::XFEM::Mean)
    return;

  // Initialize Volume Coupling
  Init_VolCoupling();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::SetupConfigurationMap()
{
  // Configuration of Consistency Terms
  configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);

  if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
  {
    // Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);

    // Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided)
  {
    // Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::XF_Con_Col] = std::pair<bool, double>(true, 1.0);

    // Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::XF_Adj_Row] = std::pair<bool, double>(true, 1.0);
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::Mean)
  {
    // Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 0.5);
    configuration_map_[INPAR::XFEM::XF_Con_Col] = std::pair<bool, double>(true, 0.5);

    // Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 0.5);
    configuration_map_[INPAR::XFEM::XF_Adj_Row] = std::pair<bool, double>(true, 0.5);
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::invalid)
    dserror("XFEM::MeshCouplingFluidFluid: Averaging Strategy not set!");
  else
    dserror("XFEM::MeshCouplingFluidFluid: You want to initialize another strategy?");

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::UpdateConfigurationMap_GP(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef DEBUG
  if (!(GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided ||
          GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided ||
          GetAveragingStrategy() == INPAR::XFEM::Mean))
    dserror(
        "XFEM::MeshCouplingFluidFluid::UpdateConfigurationMap_GP: Does your Averaging strategy "
        "change the weighing during the simulation or between gausspoints?");
#endif
  // Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row].second = full_stab;
  configuration_map_[INPAR::XFEM::X_Pen_Row].second = full_stab;
  return;
}

/*--------------------------------------------------------------------------*
 * get viscosity of the slave fluid
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::GetViscositySlave(DRT::Element* coup_ele,  ///< xfluid ele
    double& visc_s)  ///< viscosity slavesided
{
  Teuchos::RCP<MAT::Material> mat_s;
  XFEM::UTILS::GetVolumeCellMaterial(coup_ele, mat_s, GEO::CUT::Point::outside);
  if (mat_s->MaterialType() == INPAR::MAT::m_fluid)
    visc_s = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(mat_s)->Viscosity();
  else
    dserror("GetCouplingSpecificAverageWeights: Slave Material not a fluid material?");

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele)
{
  Teuchos::ParameterList params;
  DRT::Element::LocationArray la(1);
  params.set<Teuchos::RCP<std::map<int, double>>>(
      "trace_estimate_max_eigenvalue_map", ele_to_max_eigenvalue_);
  Epetra_SerialDenseMatrix dummyelemat;
  Epetra_SerialDenseVector dummyelevec;
  DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);
  if (!faceele) dserror("Cast to faceele failed!");  // todo change to dsassert

  faceele->LocationVector(*coupl_dis_, la, false);

  DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(faceele)->EstimateNitscheTraceMaxEigenvalue(
      faceele, params, *coupl_dis_, la[0].lm_, dummyelemat, dummyelevec);

  return;
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFluidFluid::ReadRestart(const int step)
{
  // copy from FSI!

  if (myrank_) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");
  //  const int    step = boundaryreader.ReadInt("step");

  if (myrank_ == 0)
  {
    IO::cout << "time: " << time << IO::endl;
    IO::cout << "step: " << step << IO::endl;
  }

  boundaryreader.ReadVector(iveln_, "iveln_res");
  boundaryreader.ReadVector(idispn_, "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_, "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");
  boundaryreader.ReadVector(idispnpi_, "idispnpi_res");

  if (not(cutter_dis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnpi_->Map()))
    dserror("Global dof numbering in maps does not match");
}

void XFEM::MeshCouplingFluidFluid::Output(
    const int step, const double time, const bool write_restart_data)
{
  // copy from FSI without the itrueresidual output!

  // output for interface
  cutter_output_->NewStep(step, time);

  cutter_output_->WriteVector("ivelnp", ivelnp_);
  cutter_output_->WriteVector("idispnp", idispnp_);

  cutter_output_->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->WriteVector("iveln_res", iveln_);
    cutter_output_->WriteVector("idispn_res", idispn_);
    cutter_output_->WriteVector("ivelnp_res", ivelnp_);
    cutter_output_->WriteVector("idispnp_res", idispnp_);
    cutter_output_->WriteVector("idispnpi_res", idispnpi_);
  }
}
