/*----------------------------------------------------------------------*/
/*! \file

\brief manages the different types of mesh based coupling conditions and thereby builds the bridge
between the xfluid class and the cut-library

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_xfem_coupling_mesh.hpp"

#include "4C_fem_dofset_transparent_independent.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_boundary_parent_calc.hpp"
#include "4C_fluid_ele_parameter_xfem.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_3D_ele_calc_lib_nitsche.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_xfem_discretization_utils.hpp"
#include "4C_xfem_interface_utils.hpp"
#include "4C_xfem_utils.hpp"
#include "4C_xfem_xfluid_contact_communicator.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
XFEM::MeshCoupling::MeshCoupling(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
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
void XFEM::MeshCoupling::set_cutter_discretization()
{
  // create a cutter discretization from conditioned nodes of the given coupling discretization
  create_cutter_dis_from_condition(suffix_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::set_conditions_to_copy()
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
void XFEM::MeshCoupling::create_cutter_dis_from_condition(std::string suffix)
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
  Teuchos::RCP<Core::FE::DiscretizationCreatorBase> discreator =
      Teuchos::rcp(new Core::FE::DiscretizationCreatorBase());
  cutter_dis_ = discreator->create_matching_discretization_from_condition(
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
    FOUR_C_THROW("Empty cutter discretization detected. No coupling can be performed...");
  }

  // for parallel jobs we have to call TransparentDofSet with additional flag true
  bool parallel = cond_dis_->Comm().NumProc() > 1;
  Teuchos::RCP<Core::DOFSets::DofSet> newdofset =
      Teuchos::rcp(new Core::DOFSets::TransparentIndependentDofSet(cond_dis_, parallel));

  cutter_dis_->ReplaceDofSet(newdofset);  // do not call this with true!!

  // create node and element distribution with elements and nodes ghosted on all processors
  Core::Rebalance::GhostDiscretizationOnAllProcs(cutter_dis_);
  cutter_dis_->fill_complete();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::gmsh_output_discretization(std::ostream& gmshfilecontent)
{
  // compute the current boundary position
  std::map<int, Core::LinAlg::Matrix<3, 1>> currinterfacepositions;

  // output of cutting discretization
  XFEM::UTILS::extract_node_vectors(cutter_dis_, currinterfacepositions, idispnp_);
  XFEM::UTILS::PrintDiscretizationToStream(cutter_dis_, cutter_dis_->Name(), true, true, true, true,
      false, false, gmshfilecontent, &currinterfacepositions);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::prepare_cutter_output()
{
  // -------------------------------------------------------------------
  // prepare output
  // -------------------------------------------------------------------

  if (!mark_geometry_)  // Do not write for marked geometry!
  {
    cutter_dis_->SetWriter(Teuchos::rcp(new Core::IO::DiscretizationWriter(cutter_dis_,
        Global::Problem::Instance()->OutputControlFile(),
        Global::Problem::Instance()->spatial_approximation_type())));
    cutter_output_ = cutter_dis_->Writer();
    cutter_output_->write_mesh(0, 0.0);
  }
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::output(const int step, const double time, const bool write_restart_data)
{
  if (!mark_geometry_)  // Do not write for marked geometry!
  {
    // output for interface
    cutter_output_->new_step(step, time);

    cutter_output_->write_vector("ivelnp", ivelnp_);
    cutter_output_->write_vector("idispnp", idispnp_);

    cutter_output_->write_element_data(firstoutputofrun_);
    firstoutputofrun_ = false;
  }
  // do not write restart for general case
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::init_state_vectors()
{
  // move state vectors to extra container class!

  const Epetra_Map* cutterdofrowmap = cutter_dis_->dof_row_map();

  ivelnp_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
  iveln_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
  ivelnm_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);

  idispnp_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
  idispn_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
  idispnpi_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::ClearState() { cutter_dis_->ClearState(); }

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::set_state()
{
  // set general vector values of cutterdis needed by background element evaluate routine
  ClearState();

  cutter_dis_->set_state("ivelnp", ivelnp_);
  cutter_dis_->set_state("iveln", iveln_);
  cutter_dis_->set_state("idispnp", idispnp_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::set_state_displacement()
{
  // set general vector values of cutterdis needed by background element evaluate routine
  ClearState();

  cutter_dis_->set_state("idispnp", idispnp_);
  cutter_dis_->set_state("idispn", idispn_);
  cutter_dis_->set_state("idispnpi", idispnpi_);
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
void XFEM::MeshCoupling::update_displacement_iteration_vectors()
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
  Teuchos::RCP<Epetra_Vector> idispcol =
      Core::LinAlg::CreateVector(*cutter_dis_->DofColMap(), true);
  Core::LinAlg::Export(*idispnp_, *idispcol);

  return idispcol;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::get_coupling_ele_location_vector(const int sid, std::vector<int>& patchlm)
{
  std::vector<int> patchlmstride, patchlmowner;  // dummy
  return coupl_dis_->gElement(sid)->LocationVector(
      *coupl_dis_, patchlm, patchlmowner, patchlmstride);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
XFEM::MeshVolCoupling::MeshVolCoupling(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
        cond_dis,              ///< discretization from which cutter discretization can be derived
    const int coupling_id,     ///< id of composite of coupling conditions
    const double time,         ///< time
    const int step,            ///< time step
    const std::string& suffix  ///< suffix for cutterdisname
    )
    : MeshCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step, suffix),
      init_volcoupling_(false),
      trace_estimate_eigenvalue_update_(Inpar::XFEM::Eigenvalue_update_every_iter),
      reset_step_(-1)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::init()
{
  XFEM::MeshCoupling::init();

  // do additional redistributino of embedded discretization and create an auxiliary dis
  if (get_averaging_strategy() != Inpar::XFEM::Xfluid_Sided)
  {
    // Initialize Volume Coupling
    Init_VolCoupling();

    // Todo: create only for Nitsche+EVP & EOS on outer embedded elements
    create_auxiliary_discretization();

    ele_to_max_eigenvalue_ = Teuchos::rcp(new std::map<int, double>());

    trace_estimate_eigenvalue_update_ =
        Core::UTILS::IntegralValue<Inpar::XFEM::TraceEstimateEigenvalueUpdate>(
            Global::Problem::Instance()->XFluidDynamicParams().sublist("STABILIZATION"),
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
    redistribute_embedded_discretization();

    init_volcoupling_ = true;
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::get_coupling_ele_location_vector(
    const int sid, std::vector<int>& patchlm)
{
  std::vector<int> patchlmstride, patchlmowner;  // dummy
  Core::Elements::Element* coupl_ele = GetCouplingElement(sid);
  coupl_ele->LocationVector(*coupl_dis_, patchlm, patchlmowner, patchlmstride);
  return;
}

/*--------------------------------------------------------------------------*
 * Ghost discretization from which the cutter_dis_ was created
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::redistribute_embedded_discretization()
{
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
    const Core::Nodes::Node* cond_node = cond_dis_->gNode(cond_node_gid);

    // get associated elements
    const Core::Elements::Element* const* cond_eles = cond_node->Elements();
    const int num_cond_ele = cond_node->NumElement();

    // loop over associated elements
    for (int ie = 0; ie < num_cond_ele; ++ie)
    {
      if (cond_eles[ie]->Owner() == mypid) adj_eles_row.insert(cond_eles[ie]->Id());

      const int* node_ids = cond_eles[ie]->NodeIds();
      for (int in = 0; in < cond_eles[ie]->num_node(); ++in)
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
    Core::LinAlg::GatherAll(full_ele_nodes_col, cond_dis_->Comm());
    Core::LinAlg::GatherAll(full_eles_col, cond_dis_->Comm());

    std::vector<int> full_nodes(full_ele_nodes_col.begin(), full_ele_nodes_col.end());
    std::vector<int> full_eles(full_eles_col.begin(), full_eles_col.end());

    Teuchos::RCP<const Epetra_Map> full_nodecolmap = Teuchos::rcp(
        new Epetra_Map(-1, full_nodes.size(), full_nodes.data(), 0, cond_dis_->Comm()));
    Teuchos::RCP<const Epetra_Map> full_elecolmap =
        Teuchos::rcp(new Epetra_Map(-1, full_eles.size(), full_eles.data(), 0, cond_dis_->Comm()));

    // redistribute nodes and elements to column (ghost) map
    cond_dis_->ExportColumnNodes(*full_nodecolmap);
    cond_dis_->export_column_elements(*full_elecolmap);

    cond_dis_->fill_complete(true, true, true);
  }

  // STEP 3: reconnect all parentelement pointers in the cutter_dis_ faceelements
  {
    for (int fele_lid = 0; fele_lid < cutter_dis_->NumMyColElements(); fele_lid++)
    {
      Core::Elements::FaceElement* fele = dynamic_cast<Core::Elements::FaceElement*>(
          cutter_dis_->gElement(cutter_dis_->ElementColMap()->GID(fele_lid)));
      if (!fele) FOUR_C_THROW("Cast to FaceElement failed!");

      Core::Elements::Element* ele = cond_dis_->gElement(fele->ParentElementId());
      if (!ele) FOUR_C_THROW("Couldn't get Parent Element!");

      fele->set_parent_master_element(ele, fele->FaceParentNumber());
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
double XFEM::MeshVolCoupling::get_estimate_nitsche_trace_max_eigenvalue(
    Core::Elements::Element* ele)
{
  if (ele_to_max_eigenvalue_->find(ele->Id()) == ele_to_max_eigenvalue_->end())
    estimate_nitsche_trace_max_eigenvalue(ele);

  return ele_to_max_eigenvalue_->at(ele->Id());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::reset_evaluated_trace_estimates()
{
  switch (trace_estimate_eigenvalue_update_)
  {
    case Inpar::XFEM::Eigenvalue_update_every_iter:
    {
      ele_to_max_eigenvalue_->clear();
      break;
    }
    case Inpar::XFEM::Eigenvalue_update_every_timestep:
    {
      if (reset_step_ < step_)
      {
        ele_to_max_eigenvalue_->clear();
        reset_step_ = step_;
      }
      break;
    }
    case Inpar::XFEM::Eigenvalue_update_once:
    {
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown Eigenvalue update strategy!");
      break;
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshVolCoupling::create_auxiliary_discretization()
{
  std::string aux_coup_disname("auxiliary_coupling_");
  aux_coup_disname += cond_dis_->Name();
  aux_coup_dis_ = Teuchos::rcp(new Core::FE::Discretization(aux_coup_disname,
      Teuchos::rcp(cond_dis_->Comm().Clone()), Global::Problem::Instance()->NDim()));

  // make the condition known to the auxiliary discretization
  // we use the same nodal ids and therefore we can just copy the conditions
  // get the set of ids of all xfem nodes
  std::vector<Core::Conditions::Condition*> xfemcnd;
  cond_dis_->GetCondition(cond_name_, xfemcnd);

  std::set<int> xfemnodeset;

  for (size_t cond = 0; cond < xfemcnd.size(); ++cond)
  {
    aux_coup_dis_->SetCondition(cond_name_, xfemcnd[cond]->copy_without_geometry());
    const std::vector<int>* nodeids_cnd = xfemcnd[cond]->GetNodes();
    for (std::vector<int>::const_iterator c = nodeids_cnd->begin(); c != nodeids_cnd->end(); ++c)
      xfemnodeset.insert(*c);
  }

  // determine sets of nodes next to xfem nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;

  // loop all column elements and label all row nodes next to a xfem node
  for (int i = 0; i < cond_dis_->NumMyColElements(); ++i)
  {
    Core::Elements::Element* actele = cond_dis_->lColElement(i);

    // get the node ids of this element
    const int numnode = actele->num_node();
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
      Teuchos::RCP<Core::Elements::Element> bndele = Teuchos::rcp(actele->Clone());
      aux_coup_dis_->add_element(bndele);
    }
  }  // end loop over column elements

  // all row nodes next to a xfem node are now added to the auxiliary discretization
  for (std::set<int>::iterator id = adjacent_row.begin(); id != adjacent_row.end(); ++id)
  {
    Core::Nodes::Node* actnode = cond_dis_->gNode(*id);
    Teuchos::RCP<Core::Nodes::Node> bndnode = Teuchos::rcp(actnode->Clone());
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
    newnoderowmap = Teuchos::rcp(
        new Epetra_Map(-1, rownodes.size(), rownodes.data(), 0, aux_coup_dis_->Comm()));

    std::vector<int> colnodes(adjacent_col.begin(), adjacent_col.end());

    // build nodecolmap for new distribution of nodes
    newnodecolmap = Teuchos::rcp(
        new Epetra_Map(-1, colnodes.size(), colnodes.data(), 0, aux_coup_dis_->Comm()));

    aux_coup_dis_->Redistribute(*newnoderowmap, *newnodecolmap, false, false, false);

    // make auxiliary discretization have the same dofs as the coupling discretization
    Teuchos::RCP<Core::DOFSets::DofSet> newdofset =
        Teuchos::rcp(new Core::DOFSets::TransparentIndependentDofSet(cond_dis_, true));
    aux_coup_dis_->ReplaceDofSet(newdofset,
        false);  // do not call this with true (no replacement in static dofsets intended)
    aux_coup_dis_->fill_complete(true, true, true);
  }
}

//! constructor
XFEM::MeshCouplingBC::MeshCouplingBC(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
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
void XFEM::MeshCouplingBC::do_condition_specific_setup()
{
  // set the initial interface displacements as they are used for initial cut position at the end of
  // Xfluid::init()
  set_interface_displacement();

  // set the interface displacements also to idispn
  idispn_->Update(1.0, *idispnp_, 0.0);

  idispnpi_->Update(1.0, *idispnp_, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
bool XFEM::MeshCouplingBC::HasMovingInterface()
{
  // get the first local col(!) node
  if (cutter_dis_->NumMyColNodes() == 0) FOUR_C_THROW("no col node on proc %i", myrank_);

  Core::Nodes::Node* lnode = cutter_dis_->lColNode(0);

  std::vector<Core::Conditions::Condition*> mycond;
  lnode->GetCondition("XFEMSurfDisplacement", mycond);

  Core::Conditions::Condition* cond = mycond[0];

  const std::string& evaltype = cond->parameters().get<std::string>("evaltype");

  return evaltype != "zero";
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::evaluate_condition(Teuchos::RCP<Epetra_Vector> ivec,
    const std::string& condname, const double time, const double dt)
{
  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < cutter_dis_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor local node
    Core::Nodes::Node* lnode = cutter_dis_->lRowNode(lnodeid);
    // the set of degrees of freedom associated with the node
    const std::vector<int> nodedofset = cutter_dis_->Dof(lnode);

    const int numdof = nodedofset.size();

    if (numdof == 0) FOUR_C_THROW("node has no dofs");
    std::vector<Core::Conditions::Condition*> mycond;
    lnode->GetCondition(condname, mycond);

    // filter out the ones with right coupling id
    std::vector<Core::Conditions::Condition*> mycond_by_coupid;

    for (auto* cond : mycond)
    {
      if (cond->parameters().get<int>("label") == coupling_id_) mycond_by_coupid.push_back(cond);
    }

    // safety check for unique condition
    // actually for WDBC surf conditions are sufficient due to weak enforcement, however, ivelnp
    // nodalwise then not reasonable for displacements one should set

    const int numconds = mycond_by_coupid.size();

    if (numconds > 1)
      std::cout << " !!! WARNING: more than one condition for node with label " << coupling_id_
                << ", think about implementing line and point conditions in XFEM !!!" << std::endl;

    if (numconds == 0) FOUR_C_THROW("no condition available!");

    Core::Conditions::Condition* cond = mycond_by_coupid[numconds - 1];  // take the last condition

    // initial value for all nodal dofs to zero
    std::vector<double> final_values(numdof, 0.0);

    if (condname == "XFEMSurfDisplacement")
      evaluate_interface_displacement(final_values, lnode, cond, time);
    else if (condname == "XFEMSurfWeakDirichlet" or condname == "XFEMRobinDirichletSurf")
      evaluate_interface_velocity(final_values, lnode, cond, time, dt);
    else
      FOUR_C_THROW("non supported condname for evaluation %s", condname.c_str());


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
void XFEM::MeshCouplingBC::evaluate_interface_velocity(std::vector<double>& final_values,
    Core::Nodes::Node* node, Core::Conditions::Condition* cond, const double time, const double dt)
{
  const std::string* evaltype = &cond->parameters().get<std::string>("evaltype");

  if (*evaltype == "zero")
  {
    // take initialized vector with zero values
  }
  else if (*evaltype == "funct_interpolated")
  {
    // evaluate function at node at current time
    evaluate_function(final_values, node->X().data(), cond, time);
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

    compute_interface_velocity_from_displacement(final_values, node, dt, evaltype);
  }
  else if (*evaltype == "displacement_1storder_with_initfunct" or
           *evaltype == "displacement_2ndorder_with_initfunct")
  {
    if (step_ == 0)  // evaluate initialization function at node at current time
    {
      evaluate_function(final_values, node->X().data(), cond, time);
    }
    else
      compute_interface_velocity_from_displacement(final_values, node, dt, evaltype);
  }
  else
    FOUR_C_THROW("evaltype not supported %s", evaltype->c_str());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::evaluate_interface_displacement(std::vector<double>& final_values,
    Core::Nodes::Node* node, Core::Conditions::Condition* cond, const double time)
{
  const std::string& evaltype = cond->parameters().get<std::string>("evaltype");

  if (evaltype == "zero")
  {
    // take initialized vector with zero values
  }
  else if (evaltype == "funct")
  {
    // evaluate function at node at current time
    evaluate_function(final_values, node->X().data(), cond, time);
  }
  else if (evaltype == "implementation")
  {
    // evaluate implementation
    // TODO: get the function name from the condition!!!
    std::string function_name = "ROTATING_BEAM";
    evaluate_implementation(final_values, node->X().data(), cond, time, function_name);
  }
  else
    FOUR_C_THROW("evaltype not supported %s", evaltype.c_str());
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::compute_interface_velocity_from_displacement(
    std::vector<double>& final_values, Core::Nodes::Node* node, const double dt,
    const std::string* evaltype)
{
  if (dt < 1e-14) FOUR_C_THROW("zero or negative time step size not allowed!!!");


  double thetaiface = 0.0;

  if (*evaltype == "displacement_1storder_wo_initfunct" or
      *evaltype == "displacement_1storder_with_initfunct")
    thetaiface = 1.0;  // for backward Euler, OST(1.0)
  else if (*evaltype == "displacement_2ndorder_wo_initfunct" or
           *evaltype == "displacement_2ndorder_with_initfunct")
    thetaiface = 0.5;  // for Crank Nicholson, OST(0.5)
  else
    FOUR_C_THROW("not supported");


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

void XFEM::MeshCouplingBC::evaluate_implementation(std::vector<double>& final_values,
    const double* x, Core::Conditions::Condition* cond, const double time,
    const std::string& function_name)
{
  const int numdof = final_values.size();

  if (function_name != "ROTATING_BEAM")
    FOUR_C_THROW("currently only the rotating beam function is available!");

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
    FOUR_C_THROW("for that time we did not define an implemented rotation %f", time);


  // rotation with constant angle velocity around point
  Core::LinAlg::Matrix<3, 1> center(true);

  center(0) = 0.0;
  center(1) = 0.0;
  center(2) = 0.0;

  Core::LinAlg::Matrix<3, 1> diff(true);
  diff(0) = x[0] - center(0);
  diff(1) = x[1] - center(1);
  diff(2) = x[2] - center(2);

  // rotation matrix
  Core::LinAlg::Matrix<3, 3> rot(true);

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

  Core::LinAlg::Matrix<3, 1> x_new(true);
  Core::LinAlg::Matrix<3, 1> rotated(true);

  rotated.multiply(rot, diff);

  x_new.update(1.0, rotated, -1.0, diff);


  for (int dof = 0; dof < numdof; ++dof)
  {
    final_values[dof] = x_new(dof);
  }
}

/*----------------------------------------------------------------------*
 |  set interface displacement at current time             schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::set_interface_displacement()
{
  if (myrank_ == 0)
    Core::IO::cout << "\t set interface displacement, time " << time_ << Core::IO::endl;

  std::string condname = "XFEMSurfDisplacement";

  evaluate_condition(idispnp_, condname, time_);
}

/*----------------------------------------------------------------------*
 |  set interface velocity at current time                 schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::MeshCouplingBC::set_interface_velocity()
{
  if (myrank_ == 0) Core::IO::cout << "\t set interface velocity, time " << time_ << Core::IO::endl;

  evaluate_condition(ivelnp_, cond_name_, time_, dt_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
//! constructor
XFEM::MeshCouplingWeakDirichlet::MeshCouplingWeakDirichlet(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
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
void XFEM::MeshCouplingWeakDirichlet::do_condition_specific_setup()
{
  XFEM::MeshCouplingBC::do_condition_specific_setup();

  // set the initial interface velocity and possible initialization function
  set_interface_velocity();

  // set the initial interface velocities also to iveln
  iveln_->Update(1.0, *ivelnp_, 0.0);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::evaluate_coupling_conditions(Core::LinAlg::Matrix<3, 1>& ivel,
    Core::LinAlg::Matrix<3, 1>& itraction, const Core::LinAlg::Matrix<3, 1>& x,
    const Core::Conditions::Condition* cond)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  evaluate_dirichlet_function(ivel, x, cond, time_);

  // no interface traction to be evaluated
  itraction.clear();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::evaluate_coupling_conditions_old_state(
    Core::LinAlg::Matrix<3, 1>& ivel, Core::LinAlg::Matrix<3, 1>& itraction,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  evaluate_dirichlet_function(ivel, x, cond, time_ - dt_);

  // no interface traction to be evaluated
  itraction.clear();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  set_interface_displacement();

  // set or compute the current prescribed interface velocities, just for XFEM WDBC
  set_interface_velocity();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::setup_configuration_map()
{
  // Configuration of Consistency Terms
  configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[Inpar::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingWeakDirichlet::update_configuration_map_gp(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_Row].second = full_stab;

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::do_condition_specific_setup()
{
  // Call Base Class
  XFEM::MeshCouplingBC::do_condition_specific_setup();

  // Check if Inflow Stabilisation is active
  if (!cutterele_conds_.size()) FOUR_C_THROW("cutterele_conds_.size = 0!");
  Core::Conditions::Condition* cond = (cutterele_conds_[0]).second;
  auto inflow_stab = cond->parameters().get<bool>("InflowStab");
  for (auto& cutterele_cond : cutterele_conds_)
  {
    Core::Conditions::Condition* cond = cutterele_cond.second;
    auto this_inflow = cond->parameters().get<bool>("InflowStab");
    if (inflow_stab != this_inflow)
      FOUR_C_THROW(
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
void XFEM::MeshCouplingNeumann::setup_configuration_map()
{
  if (inflow_stab_)
  {
    // Configuration of Penalty Terms
    configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
    configuration_map_[Inpar::XFEM::F_Adj_Col] =
        std::pair<bool, double>(true, 1.0);  //<-- IMPORTANT!: used for the constraint scaling
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::update_configuration_map_gp(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  if (inflow_stab_)
  {
    // Configuration of Penalty Terms
    double veln = normal.dot(vel_m);  // as the normal is the structural body, inflow is positive
    if (veln < 0)
    {
      configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, -density_m * veln);
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF1] =
          std::pair<bool, double>(true, -density_m * normal(0));
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF2] =
          std::pair<bool, double>(true, -density_m * normal(1));
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF3] =
          std::pair<bool, double>(true, -density_m * normal(2));
    }
    else
    {
      configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(false, 0);
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF1] = std::pair<bool, double>(false, 0);
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF2] = std::pair<bool, double>(false, 0);
      configuration_map_[Inpar::XFEM::F_Pen_Row_linF3] = std::pair<bool, double>(false, 0);
    }
  }

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::evaluate_coupling_conditions(Core::LinAlg::Matrix<3, 1>& ivel,
    Core::LinAlg::Matrix<3, 1>& itraction, const Core::LinAlg::Matrix<3, 1>& x,
    const Core::Conditions::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.clear();

  // evaluate interface traction (given by Neumann condition)
  evaluate_neumann_function(itraction, x, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::evaluate_coupling_conditions(Core::LinAlg::Matrix<3, 1>& ivel,
    Core::LinAlg::Matrix<6, 1>& itraction, const Core::LinAlg::Matrix<3, 1>& x,
    const Core::Conditions::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.clear();

  // evaluate interface traction (given by Neumann condition)
  evaluate_neumann_function(itraction, x, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::evaluate_coupling_conditions_old_state(
    Core::LinAlg::Matrix<3, 1>& ivel, Core::LinAlg::Matrix<3, 1>& itraction,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond)
{
  // no interface velocity to be evaluated
  ivel.clear();

  // evaluate interface traction (given by Neumann condition)
  evaluate_neumann_function(itraction, x, cond, time_ - dt_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNeumann::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  set_interface_displacement();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
//! constructor
XFEM::MeshCouplingNavierSlip::MeshCouplingNavierSlip(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
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
void XFEM::MeshCouplingNavierSlip::do_condition_specific_setup()
{
  // TODO: actually the same as in weak Dirichlet!!! move to base class...

  XFEM::MeshCouplingBC::do_condition_specific_setup();

  // set the initial interface velocity and possible initialization function
  set_interface_velocity();

  // set the initial interface velocities also to iveln
  iveln_->Update(1.0, *ivelnp_, 0.0);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::set_interface_velocity()
{
  if (myrank_ == 0) Core::IO::cout << "\t set interface velocity, time " << time_ << Core::IO::endl;

  //  evaluate_condition( ivelnp_, cond_name_, time_, dt_);
  evaluate_condition(ivelnp_, "XFEMRobinDirichletSurf", time_, dt_);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::evaluate_coupling_conditions(Core::LinAlg::Matrix<3, 1>& ivel,
    Core::LinAlg::Matrix<3, 1>& itraction, Core::LinAlg::Matrix<3, 3>& proj_matrix,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::LinAlg::Matrix<3, 1>& normal,
    const Core::Conditions::Condition* cond, const bool& eval_dirich_at_gp,
    double& kappa_m,  ///< fluid sided weighting
    double& visc_m,   ///< fluid sided weighting
    double& visc_s    ///< slave sided dynamic viscosity
)
{
  // Create normal projection matrix.
  setup_projection_matrix(proj_matrix, normal);

  // help variable
  int robin_id_dirch;

  if (eval_dirich_at_gp)
  {
    // evaluate interface velocity (given by weak Dirichlet condition)
    robin_id_dirch = cond->parameters().get<int>("robin_id_dirch");
    // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
    if (!std::signbit(static_cast<double>(robin_id_dirch)))
      evaluate_dirichlet_function(
          ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_);

// Safety checks
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if ((conditionsmap_robin_dirch_.find(robin_id_dirch)) == conditionsmap_robin_dirch_.end())
    {
      FOUR_C_THROW(
          "Key was not found in this instance!! Fatal error! (conditionsmap_robin_dirch_)");
    }
#endif
  }

  // evaluate interface traction (given by Neumann condition)
  robin_id_dirch = cond->parameters().get<int>("robin_id_neumann");
  if (!std::signbit(static_cast<double>(robin_id_dirch)))
  {
    // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the
    // sliplenght twice (also done in update_configuration_map_gp ... as soon as this gets relevant
    // we should merge this functions)

    // evaluate interface traction (given by Neumann condition)
    // Add this to the veljump!
    double sliplength = 0.0;
    GetSlipCoefficient(sliplength, x, cond);

    if (sliplength < 0.0) FOUR_C_THROW("The slip length can not be negative.");

    if (sliplength != 0.0)
    {
      evaluate_neumann_function(
          itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_);

      double sl_visc_fac = sliplength / (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      Core::LinAlg::Matrix<3, 1> tmp_itraction(true);
      tmp_itraction.multiply_tn(proj_matrix, itraction);
      // Project this into tangential direction!!!

      ivel.update(sl_visc_fac, tmp_itraction, 1.0);

      itraction.clear();
    }
  }

  if (force_tangvel_map_.find(cond->Id())->second)
  {
    Core::LinAlg::Matrix<3, 1> tmp_ivel(true);
    tmp_ivel.multiply_tn(proj_matrix, ivel);  // apply Projection matrix from the right. (u_0 * P^t)
    ivel.update(1.0, tmp_ivel, 0.0);
  }

// Safety checks
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!std::signbit(static_cast<double>(robin_id_dirch)))
  {
    if ((conditionsmap_robin_neumann_.find(robin_id_dirch)) == conditionsmap_robin_neumann_.end())
    {
      FOUR_C_THROW(
          "Key was not found in this instance!! Fatal error! (conditionsmap_robin_neumann_)");
    }
  }
  std::map<int, bool>::iterator it_bool;
  if ((it_bool = force_tangvel_map_.find(cond->Id())) == force_tangvel_map_.end())
  {
    FOUR_C_THROW("Key was not found in this instance!! Fatal error! (force_tangvel_map_)");
  }
#endif
}

void XFEM::MeshCouplingNavierSlip::evaluate_coupling_conditions_old_state(
    Core::LinAlg::Matrix<3, 1>& ivel, Core::LinAlg::Matrix<3, 1>& itraction,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond)
{
  // TODO: Add parameters as in call above!!!
  //  //Create normal projection matrix.
  //  Core::LinAlg::Matrix<3,3> eye(true);
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
  int robin_id_dirch = cond->parameters().get<int>("robin_id_dirch");
  // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
  if (!std::signbit(static_cast<double>(robin_id_dirch)))
    evaluate_dirichlet_function(
        ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_ - dt_);

  // evaluate interface traction (given by Neumann condition)
  robin_id_dirch = cond->parameters().get<int>("robin_id_neumann");
  if (!std::signbit(static_cast<double>(robin_id_dirch)))
    evaluate_neumann_function(
        itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_ - dt_);
}

void XFEM::MeshCouplingNavierSlip::PrepareSolve()
{
  // set the new interface displacements where DBCs or Neumann BCs have to be evaluted
  set_interface_displacement();

  //  // set the initial interface velocity and possible initialization function
  //  set_interface_velocity();
  if (myrank_ == 0) Core::IO::cout << "\t set interface velocity, time " << time_ << Core::IO::endl;
  evaluate_condition(ivelnp_, "XFEMRobinDirichletSurf", time_, dt_);
}

void XFEM::MeshCouplingNavierSlip::GetSlipCoefficient(
    double& slipcoeff, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond)
{
  // Extract correct slip length - bool pair for this condition ID.
  std::pair<double, bool>& tmp_pair = sliplength_map_.find(cond->Id())->second;

  if (tmp_pair.second)  // Is slip length constant?
    slipcoeff = tmp_pair.first;
  else  // Otherwise, evaluate function at gausspoint
    evaluate_scalar_function(slipcoeff, x.data(), tmp_pair.first, cond, time_);
}

void XFEM::MeshCouplingNavierSlip::create_robin_id_map(
    const std::vector<Core::Conditions::Condition*>& conditions_NS,
    const std::vector<Core::Conditions::Condition*>& conditions_robin,
    const std::string& robin_id_name,
    std::map<int, Core::Conditions::Condition*>& conditionsmap_robin)
{
  // Loop over all Navier Slip conditions
  for (unsigned i = 0; i < conditions_NS.size(); ++i)
  {
    // Extract its robin id (either dirichlet or neumann)
    const int tmp_robin_id = conditions_NS[i]->parameters().get<int>(robin_id_name);

    // Is this robin id active? I.e. is it not 0 or negative?
    if (!(tmp_robin_id < 0))
    {
      std::vector<Core::Conditions::Condition*> mynewcond;
      get_condition_by_robin_id(conditions_robin, tmp_robin_id, mynewcond);

      // The robin id should be unique. I.e. For one Coupling ID only a robin id can only exist
      // once.
      if (mynewcond.size() == 1)
      {
        if (!conditionsmap_robin.insert(std::make_pair(tmp_robin_id, mynewcond[0])).second)
          FOUR_C_THROW("ID already existing! For conditionsmap_robin.");
      }
      else
      {
        FOUR_C_THROW(
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
    FOUR_C_THROW("More Robin Dirichlet/Neumann Conditions provided than necessary!");
  }
}

void XFEM::MeshCouplingNavierSlip::set_condition_specific_parameters()
{
  // Build necessary maps to limit getting integers and strings on Gausspoint level.

  // Get conditions based on cutter discretization.
  std::vector<Core::Conditions::Condition*> conditions_dirich;
  cutter_dis_->GetCondition("XFEMRobinDirichletSurf", conditions_dirich);

  std::vector<Core::Conditions::Condition*> conditions_neumann;
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

  std::vector<Core::Conditions::Condition*> conditions_NS;
  cutter_dis_->GetCondition(cond_name_, conditions_NS);

  // Establishes unique connection between Navier Slip section and Robin Dirichlet Neumann sections
  create_robin_id_map(
      conditions_NS, conditions_dirich, "robin_id_dirch", conditionsmap_robin_dirch_);

  create_robin_id_map(
      conditions_NS, conditions_neumann, "robin_id_neumann", conditionsmap_robin_neumann_);

  // Create maps for easy extraction at gausspoint level
  for (auto* cond : conditions_NS)
  {
    int cond_int = cond->Id();

    double sliplength = cond->parameters().get<double>("slipcoeff");

    // Is the slip length constant? Don't call functions at GP-level unnecessary.
    bool slip_bool = (cond->parameters().get<int>("funct") < 1);

    bool force_tangential = (cond->parameters().get<int>("force_tang_vel") == 1);

    if (!sliplength_map_.insert(std::make_pair(cond_int, std::make_pair(sliplength, slip_bool)))
             .second)
      FOUR_C_THROW("ID already existing! For sliplength_map_.");

    if (!force_tangvel_map_.insert(std::make_pair(cond_int, force_tangential)).second)
      FOUR_C_THROW("ID already existing! For force_tangvel_map_.");
  }

  // Check if eval-type is same in Navier slip section and
  //       Robin Dirichlet section (Safety check! (not beautiful structure but could be worse..))
  for (auto* tmp_cond : conditions_NS)
  {
    const int tmp_robin_id = tmp_cond->parameters().get<int>("robin_id_dirch");
    if (!std::signbit(static_cast<double>(tmp_robin_id)))
    {
      if ((conditionsmap_robin_dirch_.find(tmp_robin_id)
                  ->second->parameters()
                  .get<std::string>("evaltype")) !=
          (tmp_cond->parameters().get<std::string>("evaltype")))
        FOUR_C_THROW("Not same function to evaluate in Dirichlet cond as in Main Cond.");
    }
  }
}

void XFEM::MeshCouplingNavierSlip::get_condition_by_robin_id(
    const std::vector<Core::Conditions::Condition*>& mycond, const int coupling_id,
    std::vector<Core::Conditions::Condition*>& mynewcond)
{
  mynewcond.clear();

  // select the conditions with specified "couplingID"
  for (auto* cond : mycond)
  {
    const int id = cond->parameters().get<int>("robin_id");

    if (id == coupling_id) mynewcond.push_back(cond);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::setup_configuration_map()
{
  // Configuration of Consistency Terms
  configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Con_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Con_t_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[Inpar::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingNavierSlip::update_configuration_map_gp(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
  double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
  double sliplength = 0.0;
  GetSlipCoefficient(sliplength, x, cond);

  if (sliplength < 0.0) FOUR_C_THROW("The slip length can not be negative.");

  if (sliplength != 0.0)
  {
    double stabnit = 0.0;
    double stabadj = 0.0;
    XFEM::UTILS::GetNavierSlipStabilizationParameters(
        visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);
    configuration_map_[Inpar::XFEM::F_Pen_t_Row].second = stabnit;
    configuration_map_[Inpar::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[Inpar::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[Inpar::XFEM::F_Adj_t_Row].second = stabadj;
    configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);
  }
  else
  {
    configuration_map_[Inpar::XFEM::F_Pen_t_Row].second = visc_stab_tang;
    configuration_map_[Inpar::XFEM::F_Con_t_Row] = std::pair<bool, double>(false, 0.0);
    configuration_map_[Inpar::XFEM::F_Con_t_Col] = std::pair<bool, double>(false, 0.0);
    configuration_map_[Inpar::XFEM::F_Adj_t_Row].second = 1.0;
    configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
  }

  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_n_Row].second =
      visc_stab_tang;  // full_stab <-- to keep results!
}

//! constructor
XFEM::MeshCouplingFSI::MeshCouplingFSI(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step          ///< time step
    )
    : MeshVolCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step),
      timefac_(-1.0),
      interfacelaw_(Inpar::XFEM::noslip)
{
}



/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::init_state_vectors()
{
  XFEM::MeshCoupling::init_state_vectors();

  const Epetra_Map* cutterdofrowmap = cutter_dis_->dof_row_map();
  const Epetra_Map* cutterdofcolmap = cutter_dis_->DofColMap();

  itrueresidual_ = Core::LinAlg::CreateVector(*cutterdofrowmap, true);
  iforcecol_ = Core::LinAlg::CreateVector(*cutterdofcolmap, true);
}


void XFEM::MeshCouplingFSI::complete_state_vectors()
{
  //-------------------------------------------------------------------------------
  // finalize itrueresidual vector

  // need to export the interface forces
  Epetra_Vector iforce_tmp(itrueresidual_->Map(), true);
  Epetra_Export exporter_iforce(iforcecol_->Map(), iforce_tmp.Map());
  int err1 = iforce_tmp.Export(*iforcecol_, exporter_iforce, Add);
  if (err1) FOUR_C_THROW("Export using exporter returned err=%d", err1);

  // scale the interface trueresidual with -1.0 to get the forces acting on structural side (no
  // residual-scaling!)
  itrueresidual_->Update(-1.0, iforce_tmp, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::zero_state_vectors_fsi()
{
  itrueresidual_->PutScalar(0.0);
  iforcecol_->PutScalar(0.0);
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFSI::read_restart(const int step)
{
  if (myrank_) Core::IO::cout << "read_restart for boundary discretization " << Core::IO::endl;

  //-------- boundary discretization
  Core::IO::DiscretizationReader boundaryreader(
      cutter_dis_, Global::Problem::Instance()->InputControlFile(), step);

  const double time = boundaryreader.read_double("time");
  //  const int    step = boundaryreader.ReadInt("step");

  if (myrank_ == 0)
  {
    Core::IO::cout << "time: " << time << Core::IO::endl;
    Core::IO::cout << "step: " << step << Core::IO::endl;
  }

  boundaryreader.read_vector(iveln_, "iveln_res");
  boundaryreader.read_vector(idispn_, "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.read_vector(ivelnp_, "ivelnp_res");
  boundaryreader.read_vector(idispnp_, "idispnp_res");
  boundaryreader.read_vector(idispnpi_, "idispnpi_res");

  if (not(cutter_dis_->dof_row_map())->SameAs(ivelnp_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(iveln_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispnp_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispn_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispnpi_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
}

void XFEM::MeshCouplingFSI::GetSlipCoefficient(
    double& slipcoeff, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond)
{
  // Extract correct slip length - bool pair for this condition ID.
  std::pair<double, bool>& tmp_pair = sliplength_map_.find(cond->Id())->second;

  if (tmp_pair.second)  // Is slip length constant?
    slipcoeff = tmp_pair.first;
  else  // Otherwise, evaluate function at gausspoint
    evaluate_scalar_function(slipcoeff, x.data(), tmp_pair.first, cond, time_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::GmshOutput(const std::string& filename_base, const int step,
    const int gmsh_step_diff, const bool gmsh_debug_out_screen)
{
  std::ostringstream filename_base_fsi;
  filename_base_fsi << filename_base << "_force";

  // compute the current boundary position
  std::map<int, Core::LinAlg::Matrix<3, 1>> currinterfacepositions;
  XFEM::UTILS::extract_node_vectors(cutter_dis_, currinterfacepositions, idispnp_);


  const std::string filename = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_fsi.str(), cutter_dis_->Writer()->output()->file_name(), step, gmsh_step_diff,
      gmsh_debug_out_screen, myrank_);

  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "iforce \" {" << std::endl;
    // draw vector field 'force' for every node
    Core::IO::Gmsh::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, itrueresidual_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "idispnp \" {" << std::endl;
    // draw vector field 'idispnp' for every node
    Core::IO::Gmsh::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, idispnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "ivelnp \" {" << std::endl;
    // draw vector field 'ivelnp' for every node
    Core::IO::Gmsh::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, ivelnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::gmsh_output_discretization(std::ostream& gmshfilecontent)
{
  // print surface discretization
  XFEM::MeshCoupling::gmsh_output_discretization(gmshfilecontent);

  // compute the current solid and boundary position
  std::map<int, Core::LinAlg::Matrix<3, 1>> currsolidpositions;

  // write dis with zero solid displacements here!
  Teuchos::RCP<Epetra_Vector> solid_dispnp =
      Core::LinAlg::CreateVector(*cond_dis_->dof_row_map(), true);

  XFEM::UTILS::extract_node_vectors(cond_dis_, currsolidpositions, solid_dispnp);

  XFEM::UTILS::PrintDiscretizationToStream(cond_dis_, cond_dis_->Name(), true, false, true, false,
      false, false, gmshfilecontent, &currsolidpositions);
}

void XFEM::MeshCouplingFSI::output(const int step, const double time, const bool write_restart_data)
{
  // output for interface
  cutter_output_->new_step(step, time);

  cutter_output_->write_vector("ivelnp", ivelnp_);
  cutter_output_->write_vector("idispnp", idispnp_);
  cutter_output_->write_vector("itrueresnp", itrueresidual_);

  cutter_output_->write_element_data(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->write_vector("iveln_res", iveln_);
    cutter_output_->write_vector("idispn_res", idispn_);
    cutter_output_->write_vector("ivelnp_res", ivelnp_);
    cutter_output_->write_vector("idispnp_res", idispnp_);
    cutter_output_->write_vector("idispnpi_res", idispnpi_);
  }
}

void XFEM::MeshCouplingFSI::set_condition_specific_parameters()
{
  std::vector<Core::Conditions::Condition*> conditions_XFSI;
  cutter_dis_->GetCondition(cond_name_, conditions_XFSI);

  // Create maps for easy extraction at gausspoint level
  auto i = conditions_XFSI.begin();
  for (auto* cond : conditions_XFSI)
  {
    int cond_int = cond->Id();

    double sliplength = cond->parameters().get<double>("slipcoeff");

    // Is the slip length constant? Don't call functions at GP-level unnecessary.
    bool slip_bool = (cond->parameters().get<int>("funct") < 1);

    if (!sliplength_map_.insert(std::make_pair(cond_int, std::make_pair(sliplength, slip_bool)))
             .second)
      FOUR_C_THROW("ID already existing! For sliplength_map_.");

    Inpar::XFEM::InterfaceLaw interfacelaw =
        static_cast<Inpar::XFEM::InterfaceLaw>(cond->parameters().get<int>("INTLAW"));
    if (i != conditions_XFSI.begin())
    {
      if (interfacelaw_ != interfacelaw)
        FOUR_C_THROW(
            "XFEM::MeshCouplingFSI::set_condition_specific_parameters: You defined two different "
            "FSI "
            "INTLAWS, not supported yet!");
    }
    interfacelaw_ = interfacelaw;
    i++;
  }

  if (interfacelaw_ == Inpar::XFEM::navierslip_contact)  // compute h
  {
    double hmax = 0.0;
    for (int ele = 0; ele < bg_dis_->NumMyRowElements(); ++ele)
    {
      Core::Elements::Element* fluid_ele = bg_dis_->lRowElement(ele);
      if (fluid_ele->Shape() == Core::FE::CellType::hex8)
      {
        Core::LinAlg::Matrix<3, 8> xyze(true);
        Core::Geo::fillInitialPositionArray(fluid_ele, xyze);
        double vol = XFEM::UTILS::EvalElementVolume<Core::FE::CellType::hex8>(xyze);
        hmax = std::max(hmax, XFEM::UTILS::ComputeVolEqDiameter(vol));
      }
      else
        FOUR_C_THROW("Element type != hex8, add it here!");
    }
    bg_dis_->Comm().MaxAll(&hmax, &h_scaling_, 1);
    std::cout << "==| XFEM::MeshCouplingFSI: Computed h_scaling for fluidele is: " << h_scaling_
              << "(Proc: " << bg_dis_->Comm().MyPID() << ")! |==" << std::endl;
  }

  std::cout << "==| XFEM::MeshCouplingFSI: Applied interface law is";
  switch (interfacelaw_)
  {
    case Inpar::XFEM::noslip:
    {
      std::cout << " no-slip! |==" << std::endl;
      break;
    }
    case Inpar::XFEM::noslip_splitpen:
    {
      std::cout << " no-slip with splitted normal and tangential penalty contribution! |=="
                << std::endl;
      break;
    }
    case Inpar::XFEM::slip:
    {
      std::cout << " slip! |==" << std::endl;
      break;
    }
    case Inpar::XFEM::navierslip:
    {
      std::cout << " Navier-slip! |==" << std::endl;
      break;
    }
    case Inpar::XFEM::navierslip_contact:
    {
      std::cout << " Navier-slip with Nitsche contact! |==" << std::endl;
      break;
    }
  }
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
      Core::Rebalance::GetColVersionOfRowVector(cutter_dis_, itrueresidual_);

  if (myrank_ == 0)
  {
    // compute force components
    const int nsd = 3;
    const Epetra_Map* dofcolmap = cutter_dis_->DofColMap();
    Core::LinAlg::Matrix<3, 1> c(true);
    for (int inode = 0; inode < cutter_dis_->NumMyColNodes(); ++inode)
    {
      const Core::Nodes::Node* node = cutter_dis_->lColNode(inode);
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
    const std::string fname = Global::Problem::Instance()->OutputControlFile()->file_name() +
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
void XFEM::MeshCouplingFSI::setup_configuration_map()
{
  if (get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided)
  {
    if (GetInterfaceLaw() == Inpar::XFEM::slip)
    {
      // Configuration of Consistency Terms
      configuration_map_[Inpar::XFEM::F_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Con_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[Inpar::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::noslip ||
             GetInterfaceLaw() == Inpar::XFEM::noslip_splitpen)
    {
      // Configuration of Consistency Terms
      configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[Inpar::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

      if (GetInterfaceLaw() == Inpar::XFEM::noslip)
      {
        // Configuration of Penalty Terms
        configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);
      }
      else
      {
        // Configuration of Penalty Terms
        configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);

        // to guarantee correct scaling, terms are not evaluted
        configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[Inpar::XFEM::X_Pen_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[Inpar::XFEM::X_Adj_n_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[Inpar::XFEM::F_Adj_t_Col] = std::pair<bool, double>(false, 1.0);
        configuration_map_[Inpar::XFEM::X_Adj_t_Col] = std::pair<bool, double>(false, 1.0);
      }
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::navierslip ||
             GetInterfaceLaw() == Inpar::XFEM::navierslip_contact)
    {
      // Configuration of Consistency Terms
      configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Con_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Con_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_t_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[Inpar::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
    }
    else
      FOUR_C_THROW("Intlaw not available!");
  }
  else if (get_averaging_strategy() == Inpar::XFEM::Embedded_Sided)
  {
    if (GetInterfaceLaw() == Inpar::XFEM::slip)
    {
      // Configuration of Consistency Terms
      configuration_map_[Inpar::XFEM::F_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::XS_Con_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[Inpar::XFEM::XS_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::noslip)
    {
      // Configuration of Consistency Terms
      configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::XS_Con_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

      // Configuration of Adjount Consistency Terms
      configuration_map_[Inpar::XFEM::XS_Adj_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

      // Configuration of Penalty Terms
      configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[Inpar::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);
    }
    else
      FOUR_C_THROW("Intlaw not available!");
  }
  else if (get_averaging_strategy() == Inpar::XFEM::invalid)
    FOUR_C_THROW("XFEM::MeshCouplingFSI: Averaging Strategy not set!");
  else
    FOUR_C_THROW(
        "XFEM::MeshCouplingFSI: You want to initialize another strategy than Xfluid_Sided?");
  return;
}

/*--------------------------------------------------------------------------*
 * first version without possibility to use Navier Slip
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if ((kappa_m != 1 && get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided) ||
      (kappa_m != 0 && get_averaging_strategy() == Inpar::XFEM::Embedded_Sided))
    FOUR_C_THROW("XFEM::MeshCouplingFSI::update_configuration_map_gp: kappa_m == %f", kappa_m);
#endif

  if (get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided)
  {
    if (GetInterfaceLaw() == Inpar::XFEM::slip)
    {
      // Configuration of Penalty Terms
      // configuration_map_[Inpar::XFEM::X_Con_n_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::noslip)
    {
      // configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[Inpar::XFEM::F_Pen_Row].second = full_stab;
      configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(true, full_stab);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::noslip_splitpen)
    {
      // configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[Inpar::XFEM::F_Pen_n_Row].second = full_stab;
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
      configuration_map_[Inpar::XFEM::F_Pen_t_Row].second = visc_stab_tang;
      configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, visc_stab_tang);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::navierslip)
    {
      double sliplength = 0.0;
      GetSlipCoefficient(sliplength, x, cond);

      if (sliplength < 0.0) FOUR_C_THROW("The slip should not be negative!");

      double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      double stabnit = 0.0;
      double stabadj = 0.0;
      XFEM::UTILS::GetNavierSlipStabilizationParameters(
          visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);

      configuration_map_[Inpar::XFEM::F_Pen_t_Row].second = stabnit;
      configuration_map_[Inpar::XFEM::X_Pen_t_Row].second = stabnit;
      configuration_map_[Inpar::XFEM::F_Con_t_Col] =
          std::pair<bool, double>(true, sliplength / dynvisc);
      configuration_map_[Inpar::XFEM::F_Con_t_Row] =
          std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
      configuration_map_[Inpar::XFEM::X_Con_t_Row] =
          std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
      configuration_map_[Inpar::XFEM::F_Adj_t_Row].second = stabadj;
      configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

      // Configuration of Penalty Terms
      configuration_map_[Inpar::XFEM::F_Pen_n_Row].second =
          full_stab;  // full_stab <-- to keep results!
      configuration_map_[Inpar::XFEM::X_Pen_n_Row].second =
          full_stab;  // full_stab <-- to keep results!
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::navierslip_contact)
    {
      update_configuration_map_gp_contact(kappa_m, visc_m, visc_s, density_m, visc_stab_tang,
          full_stab, x, cond, ele, bele, funct, derxy, rst_slave, normal, vel_m, fulltraction);
    }
    else
      FOUR_C_THROW("Intlaw not available!");
  }
  else if (get_averaging_strategy() == Inpar::XFEM::Embedded_Sided)
  {
    if (GetInterfaceLaw() == Inpar::XFEM::slip)
    {
      configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[Inpar::XFEM::X_Con_n_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[Inpar::XFEM::XS_Adj_n_Row] = std::pair<bool,double>(true,-1.0);
    }
    else if (GetInterfaceLaw() == Inpar::XFEM::noslip)
    {
      configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
      configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(true, visc_stab_tang);
      // configuration_map_[Inpar::XFEM::XS_Adj_Row] = std::pair<bool,double>(true,-1.0);
    }
    else
      FOUR_C_THROW("Intlaw not available!");
  }
  return;
}

/*--------------------------------------------------------------------------*
 * update_configuration_map_gp_contact for XFSCI
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::update_configuration_map_gp_contact(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  FOUR_C_ASSERT(xf_c_comm_ != Teuchos::null,
      "update_configuration_map_gp_contact but no Xfluid Contact Communicator assigned!");
#endif

  static const double MAX_sliplength = 1e40;  // large number for slip case
  static const int MAX_h = 1;  // distance from contact zone at which no-slip is prescribed
  static const int MIN_h = 0;  // distance from contact zone at which full-slip is prescribed
  static const double scaling = 1. / (MAX_h - MIN_h);

  Core::LinAlg::Matrix<2, 1> xsi(rst_slave.data(), true);  // 3-->2
  bool pure_fsi = true;                                    // do we integrate only fsi
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

  if (sliplength < 0.0) FOUR_C_THROW("The sliplength should not be negative!");

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
    configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);
    configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[Inpar::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[Inpar::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[Inpar::XFEM::X_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    configuration_map_[Inpar::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, stabadj);
    configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

    // Configuration of Penalty Terms
    configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
  }
  else  // we evaluate only the fluid terms here --> Case III (see Remark 5) -> the fluid terms, the
        // solid contributions are added in the contact part
  {
    xf_c_comm_->Inc_GP(4);
    configuration_map_[Inpar::XFEM::X_Con_Row] =
        std::pair<bool, double>(false, 0.0);  // no solid consistency
    configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, stabnit);
    configuration_map_[Inpar::XFEM::F_Con_t_Col] =
        std::pair<bool, double>(true, sliplength / dynvisc);
    configuration_map_[Inpar::XFEM::F_Con_t_Row] =
        std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
    if (sliplength > 1e-40)                       // this is basically the no-slip case ...
      configuration_map_[Inpar::XFEM::X_Con_t_Row] = std::pair<bool, double>(
          true, -stabnit + dynvisc / sliplength);  //+sign for penalty!+tangcons
    else                                           // avoid to evaluate this term ...
      configuration_map_[Inpar::XFEM::X_Con_t_Row] =
          std::pair<bool, double>(false, 0);  //+sign for penalty!+tancons
    configuration_map_[Inpar::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, stabadj);
    configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);

    // Configuration of Penalty Terms
    configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, full_stab);
    configuration_map_[Inpar::XFEM::X_Pen_n_Row] =
        std::pair<bool, double>(false, 0.0);  // no normal penalty
  }
}

/*--------------------------------------------------------------------------*
 * Evaluate the Structural Cauchy Stress Matrix and it's linearization with respect to the
 *displacements
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::evaluate_structural_cauchy_stress(Core::Elements::Element* coupl_ele,
    Core::LinAlg::Matrix<3, 1>& rst_slave, std::vector<double>& eledisp,
    const Core::LinAlg::Matrix<3, 1>& normal,
    std::vector<Core::LinAlg::SerialDenseMatrix>& solid_stress)
{
  if (get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided) return;

  FOUR_C_THROW_UNLESS(coupl_ele->Shape() == Core::FE::CellType::hex8,
      "XFEM::MeshCouplingFSI::evaluate_structural_cauchy_stress is currently only implemented for "
      "hex8 elements");


  auto evaluate_cauchy_n_dir_and_derivatives = std::invoke(
      [&]() -> std::function<void(const Core::LinAlg::Matrix<NUMDIM_SOH8, 1>&, double&,
                Core::LinAlg::SerialDenseMatrix&, Core::LinAlg::SerialDenseMatrix&)>
      {
        if (auto* solid_ele = dynamic_cast<Discret::ELEMENTS::SoBase*>(coupl_ele);
            solid_ele != nullptr)
        {
          return [&, solid_ele](const Core::LinAlg::Matrix<NUMDIM_SOH8, 1>& dir,
                     double& cauchy_n_dir, Core::LinAlg::SerialDenseMatrix& d_cauchy_d_d,
                     Core::LinAlg::SerialDenseMatrix& d2_cauchy_d_d2)
          {
            solid_ele->get_cauchy_n_dir_and_derivatives_at_xi(rst_slave, eledisp, normal, dir,
                cauchy_n_dir, &d_cauchy_d_d, &d2_cauchy_d_d2, nullptr, nullptr, nullptr, nullptr,
                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
          };
        }
        else if (auto* solid_ele = dynamic_cast<Discret::ELEMENTS::Solid*>(coupl_ele);
                 solid_ele != nullptr)
        {
          return [&, solid_ele](const Core::LinAlg::Matrix<NUMDIM_SOH8, 1>& dir,
                     double& cauchy_n_dir, Core::LinAlg::SerialDenseMatrix& d_cauchy_d_d,
                     Core::LinAlg::SerialDenseMatrix& d2_cauchy_d_d2)
          {
            Discret::ELEMENTS::CauchyNDirLinearizations<3> linearizations{};
            linearizations.d_cauchyndir_dd = &d_cauchy_d_d;
            linearizations.d2_cauchyndir_dd2 = &d2_cauchy_d_d2;

            cauchy_n_dir = solid_ele->get_normal_cauchy_stress_at_xi<3>(
                eledisp, rst_slave, normal, dir, linearizations);
          };
        }
        else
        {
          FOUR_C_THROW("Unknown solid element type");
        }
      });

  solid_stress.resize(5);  // traction,dtdd,d2dddx,d2dddy,d2dddz
  solid_stress[0].reshape(NUMDIM_SOH8, 1);
  Core::LinAlg::Matrix<NUMDIM_SOH8, 1> traction(solid_stress[0].values(), true);

  solid_stress[1].reshape(NUMDOF_SOH8, NUMDIM_SOH8);
  Core::LinAlg::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> dtraction_dd_l(solid_stress[1].values(), true);

  traction.clear();

  static Core::LinAlg::SerialDenseMatrix dtraction_dd_i;
  for (int i = 0; i < NUMDIM_SOH8; ++i)
  {
    Core::LinAlg::Matrix<NUMDIM_SOH8, 1> ei(true);
    ei(i, 0) = 1.;

    evaluate_cauchy_n_dir_and_derivatives(ei, traction(i, 0), dtraction_dd_i, solid_stress[2 + i]);

    Core::LinAlg::Matrix<NUMDOF_SOH8, 1> dtraction_dd_i_l(dtraction_dd_i.values(), true);
    for (int col = 0; col < NUMDOF_SOH8; ++col) dtraction_dd_l(col, i) = dtraction_dd_i_l(col, 0);
  }

  FOUR_C_THROW_UNLESS(timefac_ > 0,
      "XFEM::MeshCouplingFSI::evaluate_structural_cauchy_stress: timefac = %f, not set!", timefac_);

  // Change from linearization w.r.t. displacements to linearization w.r.t. velocities
  // (All other linearizations on the Nitsche Interface are evaluated like this)
  solid_stress[1].scale(timefac_);
  for (int idx = 2; idx < 5; ++idx) solid_stress[idx].scale(-timefac_ * timefac_);
}

/*--------------------------------------------------------------------------*
 * get stress tangent of the slave solid
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::get_stress_tangent_slave(
    Core::Elements::Element* coup_ele,  ///< solid ele
    double& e_s)                        ///< stress tangent slavesided
{
  //  if (coup_ele->Material()->MaterialType() == Core::Materials::m_elasthyper)
  //    e_s = Teuchos::rcp_dynamic_cast<Mat::ElastHyper>(coup_ele->Material())->GetYoung();
  //  else
  //    FOUR_C_THROW("get_coupling_specific_average_weights: Slave Material not a Elasthyper
  //    material?");

  // this is a temporal hack as we calculate "E/h" directely with the generalized eigenvalue problem
  // ... need to work on the input section to clarify this ...
  e_s = timefac_;

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::estimate_nitsche_trace_max_eigenvalue(Core::Elements::Element* ele)
{
  Discret::ELEMENTS::StructuralSurface* solidfaceele =
      dynamic_cast<Discret::ELEMENTS::StructuralSurface*>(ele);
  FOUR_C_ASSERT(solidfaceele != nullptr, "Cast to StructuralSurface failed!");

  solidfaceele->set_parent_master_element(
      coupl_dis_->gElement(solidfaceele->ParentElementId()), solidfaceele->FaceParentNumber());

  Core::Elements::Element::LocationArray la(1);
  solidfaceele->parent_element()->LocationVector(*coupl_dis_, la, false);

  // extract eledisp here
  // parent and boundary displacement at n+1
  std::vector<double> eledisp((la[0].lm_).size());
  Teuchos::RCP<const Epetra_Vector> dispnp = coupl_dis_->GetState("dispnp");
  if (dispnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'dispnp'");

  Core::FE::ExtractMyValues(*dispnp, eledisp, la[0].lm_);
  (*ele_to_max_eigenvalue_)[ele->Id()] =
      solidfaceele->estimate_nitsche_trace_max_eigenvalue_combined(
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
  if (get_averaging_strategy() != Inpar::XFEM::Xfluid_Sided) reset_evaluated_trace_estimates();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::RegisterSideProc(int sid)
{
  if (GetInterfaceLaw() == Inpar::XFEM::navierslip_contact)
    Get_Contact_Comm()->RegisterSideProc(sid);
  return;
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
bool XFEM::MeshCouplingFSI::initialize_fluid_state(Teuchos::RCP<Core::Geo::CutWizard> cutwizard,
    Teuchos::RCP<Core::FE::Discretization> fluiddis,
    Teuchos::RCP<XFEM::ConditionManager> condition_manager,
    Teuchos::RCP<Teuchos::ParameterList> fluidparams)
{
  if (GetInterfaceLaw() == Inpar::XFEM::navierslip_contact)
    Get_Contact_Comm()->initialize_fluid_state(cutwizard, fluiddis, condition_manager, fluidparams);
  return (GetInterfaceLaw() == Inpar::XFEM::navierslip_contact);
}

XFEM::MeshCouplingFluidFluid::MeshCouplingFluidFluid(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
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
void XFEM::MeshCouplingFluidFluid::get_interface_slave_material(
    Core::Elements::Element* actele, Teuchos::RCP<Core::Mat::Material>& mat)
{
  XFEM::UTILS::get_volume_cell_material(actele, mat, Core::Geo::Cut::Point::outside);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::redistribute_for_error_calculation()
{
  if (get_averaging_strategy() == Inpar::XFEM::Embedded_Sided ||
      get_averaging_strategy() == Inpar::XFEM::Mean)
    return;

  // Initialize Volume Coupling
  Init_VolCoupling();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::setup_configuration_map()
{
  // Configuration of Consistency Terms
  configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);

  // Configuration of Adjount Consistency Terms
  configuration_map_[Inpar::XFEM::F_Adj_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::X_Adj_Col] = std::pair<bool, double>(true, 1.0);

  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(true, 1.0);
  configuration_map_[Inpar::XFEM::X_Pen_Col] = std::pair<bool, double>(true, 1.0);

  if (get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided)
  {
    // Configuration of Consistency Terms
    configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);

    // Configuration of Adjount Consistency Terms
    configuration_map_[Inpar::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 1.0);
  }
  else if (get_averaging_strategy() == Inpar::XFEM::Embedded_Sided)
  {
    // Configuration of Consistency Terms
    configuration_map_[Inpar::XFEM::XF_Con_Col] = std::pair<bool, double>(true, 1.0);

    // Configuration of Adjount Consistency Terms
    configuration_map_[Inpar::XFEM::XF_Adj_Row] = std::pair<bool, double>(true, 1.0);
  }
  else if (get_averaging_strategy() == Inpar::XFEM::Mean)
  {
    // Configuration of Consistency Terms
    configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(true, 0.5);
    configuration_map_[Inpar::XFEM::XF_Con_Col] = std::pair<bool, double>(true, 0.5);

    // Configuration of Adjount Consistency Terms
    configuration_map_[Inpar::XFEM::F_Adj_Row] = std::pair<bool, double>(true, 0.5);
    configuration_map_[Inpar::XFEM::XF_Adj_Row] = std::pair<bool, double>(true, 0.5);
  }
  else if (get_averaging_strategy() == Inpar::XFEM::invalid)
    FOUR_C_THROW("XFEM::MeshCouplingFluidFluid: Averaging Strategy not set!");
  else
    FOUR_C_THROW("XFEM::MeshCouplingFluidFluid: You want to initialize another strategy?");

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::update_configuration_map_gp(
    double& kappa_m,         //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond,
    Core::Elements::Element* ele,   //< Element
    Core::Elements::Element* bele,  //< Boundary Element
    double* funct,                  //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    Core::LinAlg::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    Core::LinAlg::Matrix<3, 1>& normal,     //< normal at gp
    Core::LinAlg::Matrix<3, 1>& vel_m,      //< master velocity at gp
    double* fulltraction                    //< precomputed fsi traction (sigmaF n + gamma relvel)
)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!(get_averaging_strategy() == Inpar::XFEM::Xfluid_Sided ||
          get_averaging_strategy() == Inpar::XFEM::Embedded_Sided ||
          get_averaging_strategy() == Inpar::XFEM::Mean))
    FOUR_C_THROW(
        "XFEM::MeshCouplingFluidFluid::update_configuration_map_gp: Does your Averaging strategy "
        "change the weighing during the simulation or between gausspoints?");
#endif
  // Configuration of Penalty Terms
  configuration_map_[Inpar::XFEM::F_Pen_Row].second = full_stab;
  configuration_map_[Inpar::XFEM::X_Pen_Row].second = full_stab;
  return;
}

/*--------------------------------------------------------------------------*
 * get viscosity of the slave fluid
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::GetViscositySlave(
    Core::Elements::Element* coup_ele,  ///< xfluid ele
    double& visc_s)                     ///< viscosity slavesided
{
  Teuchos::RCP<Core::Mat::Material> mat_s;
  XFEM::UTILS::get_volume_cell_material(coup_ele, mat_s, Core::Geo::Cut::Point::outside);
  if (mat_s->MaterialType() == Core::Materials::m_fluid)
    visc_s = Teuchos::rcp_dynamic_cast<Mat::NewtonianFluid>(mat_s)->Viscosity();
  else
    FOUR_C_THROW("get_coupling_specific_average_weights: Slave Material not a fluid material?");

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFluidFluid::estimate_nitsche_trace_max_eigenvalue(
    Core::Elements::Element* ele)
{
  Teuchos::ParameterList params;
  Core::Elements::Element::LocationArray la(1);
  params.set<Teuchos::RCP<std::map<int, double>>>(
      "trace_estimate_max_eigenvalue_map", ele_to_max_eigenvalue_);
  Core::LinAlg::SerialDenseMatrix dummyelemat;
  Core::LinAlg::SerialDenseVector dummyelevec;
  Core::Elements::FaceElement* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);
  if (!faceele) FOUR_C_THROW("Cast to faceele failed!");  // todo change to FOUR_C_ASSERT

  faceele->LocationVector(*coupl_dis_, la, false);

  Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(faceele)
      ->estimate_nitsche_trace_max_eigenvalue(
          faceele, params, *coupl_dis_, la[0].lm_, dummyelemat, dummyelevec);

  return;
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFluidFluid::read_restart(const int step)
{
  // copy from FSI!

  if (myrank_) Core::IO::cout << "read_restart for boundary discretization " << Core::IO::endl;

  //-------- boundary discretization
  Core::IO::DiscretizationReader boundaryreader(
      cutter_dis_, Global::Problem::Instance()->InputControlFile(), step);

  const double time = boundaryreader.read_double("time");
  //  const int    step = boundaryreader.ReadInt("step");

  if (myrank_ == 0)
  {
    Core::IO::cout << "time: " << time << Core::IO::endl;
    Core::IO::cout << "step: " << step << Core::IO::endl;
  }

  boundaryreader.read_vector(iveln_, "iveln_res");
  boundaryreader.read_vector(idispn_, "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.read_vector(ivelnp_, "ivelnp_res");
  boundaryreader.read_vector(idispnp_, "idispnp_res");
  boundaryreader.read_vector(idispnpi_, "idispnpi_res");

  if (not(cutter_dis_->dof_row_map())->SameAs(ivelnp_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(iveln_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispnp_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispn_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
  if (not(cutter_dis_->dof_row_map())->SameAs(idispnpi_->Map()))
    FOUR_C_THROW("Global dof numbering in maps does not match");
}

void XFEM::MeshCouplingFluidFluid::output(
    const int step, const double time, const bool write_restart_data)
{
  // copy from FSI without the itrueresidual output!

  // output for interface
  cutter_output_->new_step(step, time);

  cutter_output_->write_vector("ivelnp", ivelnp_);
  cutter_output_->write_vector("idispnp", idispnp_);

  cutter_output_->write_element_data(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->write_vector("iveln_res", iveln_);
    cutter_output_->write_vector("idispn_res", idispn_);
    cutter_output_->write_vector("ivelnp_res", ivelnp_);
    cutter_output_->write_vector("idispnp_res", idispnp_);
    cutter_output_->write_vector("idispnpi_res", idispnpi_);
  }
}

FOUR_C_NAMESPACE_CLOSE
