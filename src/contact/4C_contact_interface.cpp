/*---------------------------------------------------------------------*/
/*! \file
\brief One contact interface

\level 2


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_interface.hpp"

#include "4C_binstrategy.hpp"
#include "4C_contact_coupling2d.hpp"
#include "4C_contact_coupling3d.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_interpolator.hpp"
#include "4C_contact_line_coupling.hpp"
#include "4C_contact_nitsche_utils.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_selfcontact_binarytree_unbiased.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_binarytree.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_dofset.hpp"
#include "4C_mortar_projector.hpp"
#include "4C_rebalance_graph_based.hpp"
#include "4C_scatra_ele_parameter_boundary.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::InterfaceDataContainer::InterfaceDataContainer()
    : selfcontact_(false),
      friction_(false),
      non_smooth_contact_(false),
      two_half_pass_(false),
      constr_direction_(Inpar::CONTACT::constr_vague),
      activenodes_(Teuchos::null),
      activedofs_(Teuchos::null),
      inactivenodes_(Teuchos::null),
      inactivedofs_(Teuchos::null),
      activen_(Teuchos::null),
      activet_(Teuchos::null),
      slipnodes_(Teuchos::null),
      slipdofs_(Teuchos::null),
      slipt_(Teuchos::null),
      nonsmoothnodes_(Teuchos::null),
      smoothnodes_(Teuchos::null),
      sdof_vertex_rowmap_(Teuchos::null),
      sdof_vertex_colmap_(Teuchos::null),
      sdof_edge_rowmap_(Teuchos::null),
      sdof_edge_colmap_(Teuchos::null),
      sdof_surf_rowmap_(Teuchos::null),
      sdof_surf_colmap_(Teuchos::null),
      nextendedghosting_(Teuchos::null),
      eextendedghosting_(Teuchos::null),
      binarytreeself_(Teuchos::null),
      cn_values_(Teuchos::null),
      ct_values_(Teuchos::null),
      smpairs_(0),
      smintpairs_(0),
      intcells_(0)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Interface> CONTACT::Interface::create(const int id, const Epetra_Comm& comm,
    const int spatialDim, const Teuchos::ParameterList& icontact, const bool selfcontact)
{
  Teuchos::RCP<Mortar::InterfaceDataContainer> interfaceData_ptr =
      Teuchos::rcp(new CONTACT::InterfaceDataContainer());
  return Teuchos::rcp(
      new Interface(interfaceData_ptr, id, comm, spatialDim, icontact, selfcontact));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Interface::Interface(const Teuchos::RCP<CONTACT::InterfaceDataContainer>& interfaceData)
    : Mortar::Interface(interfaceData),
      interface_data_(interfaceData),
      selfcontact_(interface_data_->is_self_contact()),
      friction_(interface_data_->is_friction()),
      nonSmoothContact_(interface_data_->is_non_smooth_contact()),
      two_half_pass_(interface_data_->is_two_half_pass()),
      constr_direction_(interface_data_->constraint_direction()),
      activenodes_(interface_data_->active_nodes()),
      activedofs_(interface_data_->active_dofs()),
      inactivenodes_(interface_data_->in_active_nodes()),
      inactivedofs_(interface_data_->in_active_dofs()),
      activen_(interface_data_->active_n()),
      activet_(interface_data_->active_t()),
      slipnodes_(interface_data_->slip_nodes()),
      slipdofs_(interface_data_->slip_dofs()),
      slipt_(interface_data_->slip_t()),
      nonsmoothnodes_(interface_data_->non_smooth_nodes()),
      smoothnodes_(interface_data_->smooth_nodes()),
      sdofVertexRowmap_(interface_data_->sdof_vertex_rowmap()),
      sdofVertexColmap_(interface_data_->sdof_vertex_colmap()),
      sdofEdgeRowmap_(interface_data_->sdof_edge_rowmap()),
      sdofEdgeColmap_(interface_data_->sdof_edge_colmap()),
      sdofSurfRowmap_(interface_data_->sdof_surf_rowmap()),
      sdofSurfColmap_(interface_data_->sdof_surf_colmap()),
      nextendedghosting_(interface_data_->n_extended_ghosting()),
      eextendedghosting_(interface_data_->e_extended_ghosting()),
      binarytreeself_(interface_data_->binary_tree_self()),
      cnValues_(interface_data_->cn_values()),
      ctValues_(interface_data_->ct_values()),
      smpairs_(interface_data_->sm_int_pairs()),
      smintpairs_(interface_data_->sm_int_pairs()),
      intcells_(interface_data_->int_cells())
{
  /* do nothing */
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData,
    const int id, const Epetra_Comm& comm, const int spatialDim,
    const Teuchos::ParameterList& icontact, bool selfcontact)
    : Mortar::Interface(interfaceData, id, comm, spatialDim, icontact),
      interface_data_(
          Teuchos::rcp_dynamic_cast<CONTACT::InterfaceDataContainer>(interfaceData, true)),
      selfcontact_(interface_data_->is_self_contact()),
      friction_(interface_data_->is_friction()),
      nonSmoothContact_(interface_data_->is_non_smooth_contact()),
      two_half_pass_(interface_data_->is_two_half_pass()),
      constr_direction_(interface_data_->constraint_direction()),
      activenodes_(interface_data_->active_nodes()),
      activedofs_(interface_data_->active_dofs()),
      inactivenodes_(interface_data_->in_active_nodes()),
      inactivedofs_(interface_data_->in_active_dofs()),
      activen_(interface_data_->active_n()),
      activet_(interface_data_->active_t()),
      slipnodes_(interface_data_->slip_nodes()),
      slipdofs_(interface_data_->slip_dofs()),
      slipt_(interface_data_->slip_t()),
      nonsmoothnodes_(interface_data_->non_smooth_nodes()),
      smoothnodes_(interface_data_->smooth_nodes()),
      sdofVertexRowmap_(interface_data_->sdof_vertex_rowmap()),
      sdofVertexColmap_(interface_data_->sdof_vertex_colmap()),
      sdofEdgeRowmap_(interface_data_->sdof_edge_rowmap()),
      sdofEdgeColmap_(interface_data_->sdof_edge_colmap()),
      sdofSurfRowmap_(interface_data_->sdof_surf_rowmap()),
      sdofSurfColmap_(interface_data_->sdof_surf_colmap()),
      nextendedghosting_(interface_data_->n_extended_ghosting()),
      eextendedghosting_(interface_data_->e_extended_ghosting()),
      binarytreeself_(interface_data_->binary_tree_self()),
      cnValues_(interface_data_->cn_values()),
      ctValues_(interface_data_->ct_values()),
      smpairs_(interface_data_->sm_int_pairs()),
      smintpairs_(interface_data_->sm_int_pairs()),
      intcells_(interface_data_->int_cells())
{
  selfcontact_ = selfcontact;
  nonSmoothContact_ = Core::UTILS::IntegralValue<int>(icontact, "NONSMOOTH_GEOMETRIES");
  two_half_pass_ = icontact.get<bool>("Two_half_pass");
  constr_direction_ = Core::UTILS::IntegralValue<Inpar::CONTACT::ConstraintDirection>(
      icontact, "CONSTRAINT_DIRECTIONS");
  smpairs_ = 0;
  smintpairs_ = 0;
  intcells_ = 0;

  // set frictional contact status
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(icontact, "FRICTION");
  if (ftype != Inpar::CONTACT::friction_none) friction_ = true;

  // set poro contact
  if (icontact.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
      icontact.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra ||
      icontact.get<int>("PROBTYPE") == Inpar::CONTACT::fpi)
  {
    set_poro_flag(true);
    set_poro_type(Inpar::Mortar::poroelast);
  }
  if (icontact.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra)
    set_poro_type(Inpar::Mortar::poroscatra);

  // set ehl contact
  if (icontact.get<int>("PROBTYPE") == Inpar::CONTACT::ehl) set_ehl_flag(true);

  // check for redundant slave storage
  // needed for self contact but not wanted for general contact
  // for self contact this is ensured in BuildInterfaces in contact_strategy_factory.cpp
  // so we only print a warning here, as it is possible to have another contact interface with a
  // different ID that does not need to be a self contact interface
  if (!(selfcontact_ or nonSmoothContact_) &&
      interface_data_->get_extend_ghosting() == Inpar::Mortar::ExtendGhosting::redundant_all)
    if (Interface::get_comm().MyPID() == 0)
      std::cout << "\n\nWARNING: We do not want redundant interface storage for contact where not "
                   "needed, as it is very expensive. But we need it e.g. for self contact."
                << std::endl;

  // initialize extended ghosting for RR loop
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::update_master_slave_sets()
{
  // call mortar function
  Mortar::Interface::update_master_slave_sets();

  //********************************************************************
  // DOFS
  //********************************************************************
  // do the same business for dofs
  // (get row and column maps of slave and master dofs seperately)
  if (nonSmoothContact_)
  {
    std::vector<int> sVc;  // slave column map
    std::vector<int> sVr;  // slave row map
    std::vector<int> sEc;  // master column map
    std::vector<int> sEr;  // master row map
    std::vector<int> sSc;  // master column map
    std::vector<int> sSr;  // master row map

    for (int i = 0; i < discret().node_col_map()->NumMyElements(); ++i)
    {
      int gid = discret().node_col_map()->GID(i);
      Core::Nodes::Node* node = discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
      const bool isslave = mrtrnode->is_slave();
      const int numdof = mrtrnode->num_dof();

      if (isslave)
      {
        // vertex
        if (mrtrnode->is_on_corner())
        {
          for (int j = 0; j < numdof; ++j) sVc.push_back(mrtrnode->dofs()[j]);

          if (discret().node_row_map()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sVr.push_back(mrtrnode->dofs()[j]);
        }
        // edge
        else if (mrtrnode->is_on_edge())
        {
          for (int j = 0; j < numdof; ++j) sEc.push_back(mrtrnode->dofs()[j]);

          if (discret().node_row_map()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sEr.push_back(mrtrnode->dofs()[j]);
        }
        // surface
        else if (!mrtrnode->is_on_corner_edge())
        {
          for (int j = 0; j < numdof; ++j) sSc.push_back(mrtrnode->dofs()[j]);

          if (discret().node_row_map()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sSr.push_back(mrtrnode->dofs()[j]);
        }
        else
        {
          FOUR_C_THROW("unknown case!");
        }
      }
    }

    sdofVertexRowmap_ =
        Teuchos::rcp(new Epetra_Map(-1, (int)sVr.size(), sVr.data(), 0, get_comm()));
    sdofVertexColmap_ =
        Teuchos::rcp(new Epetra_Map(-1, (int)sVc.size(), sVc.data(), 0, get_comm()));
    sdofEdgeRowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sEr.size(), sEr.data(), 0, get_comm()));
    sdofEdgeColmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sEc.size(), sEc.data(), 0, get_comm()));
    sdofSurfRowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sSr.size(), sSr.data(), 0, get_comm()));
    sdofSurfColmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sSc.size(), sSc.data(), 0, get_comm()));
  }
}

/*----------------------------------------------------------------------*
 |  create and fill cn vector                                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::set_cn_ct_values(const int& iter)
{
  // get cn from the input file
  const double cn = interface_params().get<double>("SEMI_SMOOTH_CN");
  const double ct = interface_params().get<double>("SEMI_SMOOTH_CT");

  // set all nodal cn-values to the input value
  get_cn() = Core::LinAlg::CreateVector(*slave_row_nodes(), true);
  int err = get_cn()->PutScalar(cn);
  if (err != 0) FOUR_C_THROW("cn definition failed!");

  // set all nodal ct-values to the input value
  if (friction_)
  {
    get_ct() = Core::LinAlg::CreateVector(*slave_row_nodes(), true);
    err = get_ct()->PutScalar(ct);
    if (err != 0) FOUR_C_THROW("cn definition failed!");
  }

  // modification for edge/corner nodes
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %i", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // calculate characteristic edge length:
    Core::Elements::Element* ele = cnode->elements()[0];
    Element* cele = dynamic_cast<Element*>(ele);
    std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
    std::array<double, 3> pos2 = {0.0, 0.0, 0.0};
    std::array<double, 3> vec = {0.0, 0.0, 0.0};

    pos1[0] = dynamic_cast<Node*>(cele->nodes()[0])->x()[0];
    pos1[1] = dynamic_cast<Node*>(cele->nodes()[0])->x()[1];
    pos1[2] = dynamic_cast<Node*>(cele->nodes()[0])->x()[2];

    pos2[0] = dynamic_cast<Node*>(cele->nodes()[1])->x()[0];
    pos2[1] = dynamic_cast<Node*>(cele->nodes()[1])->x()[1];
    pos2[2] = dynamic_cast<Node*>(cele->nodes()[1])->x()[2];

    vec[0] = pos1[0] - pos2[0];
    vec[1] = pos1[1] - pos2[1];
    vec[2] = pos1[2] - pos2[2];

    const double length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (length < 1e-12)
    {
      std::cout << "*** WARNING *** element edge nearly zero" << std::endl;
      continue;
    }

    if (cnode->is_on_edge())
    {
      get_cn_ref()[get_cn_ref().Map().LID(cnode->id())] = cn * (length * length);
      if (friction_) get_ct_ref()[get_ct_ref().Map().LID(cnode->id())] = ct * (length * length);
    }

    if (cnode->is_on_corner())
    {
      get_cn_ref()[get_cn_ref().Map().LID(cnode->id())] = cn * (length * length * length * length);
      if (friction_)
        get_ct_ref()[get_ct_ref().Map().LID(cnode->id())] =
            ct * (length * length * length * length);
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::Interface& interface)
{
  interface.print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::print(std::ostream& os) const
{
  if (get_comm().MyPID() == 0) os << "Contact ";
  Mortar::Interface::print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  add contact node (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_node(Teuchos::RCP<CONTACT::Node> cnode)
{
  idiscret_->add_node(cnode);
  return;
}

/*----------------------------------------------------------------------*
 |  add contact element (public)                             mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_element(Teuchos::RCP<CONTACT::Element> cele)
{
  // check for quadratic 2d slave elements to be modified
  if (cele->is_slave() && (cele->shape() == Core::FE::CellType::line3)) quadslave_ = true;

  // check for quadratic 3d slave elements to be modified
  if (cele->is_slave() &&
      (cele->shape() == Core::FE::CellType::quad8 || cele->shape() == Core::FE::CellType::tri6))
    quadslave_ = true;

  idiscret_->add_element(cele);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::update_parallel_layout_and_data_structures(const bool perform_rebalancing,
    const bool enforce_ghosting_update, const int maxdof, const double meanVelocity)
{
  if (perform_rebalancing)
  {
    redistribute();
    fill_complete_new(false, maxdof);
  }

  if (perform_rebalancing || enforce_ghosting_update)
  {
    // Assure that at least some maps are available
    if (!filled()) fill_complete_new(false, maxdof);

    // Finalize interface maps
    extend_interface_ghosting_safely(meanVelocity);
    fill_complete_new(true, maxdof);
  }

  // print new parallel distribution
  if (perform_rebalancing)
  {
    if (get_comm().MyPID() == 0)
      std::cout << "Interface parallel distribution after rebalancing:" << std::endl;
    print_parallel_distribution();
  }

  if (perform_rebalancing || enforce_ghosting_update) create_search_tree();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::fill_complete_new(const bool isFinalParallelDistribution, const int maxdof)
{
  std::stringstream ss;
  ss << "CONTACT::Interface::fill_complete_new of '" << discret().name() << "'";
  TEUCHOS_FUNC_TIME_MONITOR(ss.str());

  // store maximum global dof ID handed in
  // this ID is later needed when setting up the Lagrange multiplier
  // dof map, which of course must not overlap with existing dof ranges
  maxdofglobal_ = maxdof;

  /* We'd like to call idiscret_.fill_complete(true,false,false) but this will assign all nodes new
   * degrees of freedom which we don't want. We would like to use the degrees of freedom that were
   * stored in the mortar nodes. To do so, we have to create and set our own version of a
   * Mortar::DofSet class before we call fill_complete on the interface discretization. The
   * specialized DofSet class will not assign new dofs but will assign the dofs stored in the nodes.
   */
  {
    Teuchos::RCP<Mortar::DofSet> mrtrdofset = Teuchos::rcp(new Mortar::DofSet());
    discret().replace_dof_set(mrtrdofset);
  }

  // fill_complete the interface discretization
  discret().fill_complete(isFinalParallelDistribution, false, false);

  // check whether crosspoints / edge nodes shall be considered or not
  initialize_cross_points();

  // check for const/linear interpolation of 2D/3D quadratic Lagrange multipliers
  initialize_lag_mult_const();
  initialize_lag_mult_lin();

  // check/init corner/edge modification
  initialize_corner_edge();

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily address them
  update_master_slave_sets();

  // initialize node and element data container
  initialize_data_container();

  // Communicate quadslave status among ALL processors
  communicate_quad_slave_status_among_all_procs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::extend_interface_ghosting_safely(const double meanVelocity)
{
  using Teuchos::RCP;

  if (discret().node_col_map() == nullptr) FOUR_C_THROW("NodeColMap not set.");
  if (discret().element_col_map() == nullptr) FOUR_C_THROW("ElementColMap not set.");

  // later we might export node and element column map to extended or even FULL overlap,
  // thus store the standard column maps first
  {
    // get standard nodal column map (overlap=1)
    oldnodecolmap_ = Teuchos::rcp(new Epetra_Map(*(discret().node_col_map())));

    // get standard element column map (overlap=1)
    oldelecolmap_ = Teuchos::rcp(new Epetra_Map(*(discret().element_col_map())));
  }

  switch (interface_data_->get_extend_ghosting())
  {
    case Inpar::Mortar::ExtendGhosting::redundant_all:
    {
      // to ease our search algorithms we'll afford the luxury to ghost all nodes
      // on all processors. To do so, we'll take the node row map and export it to
      // full overlap. Then we export the discretization to full overlap column map.
      // This way, also all mortar elements will be fully ghosted on all processors.

      // we want to do full ghosting on all procs
      std::vector<int> allproc(get_comm().NumProc());
      for (int i = 0; i < get_comm().NumProc(); ++i) allproc[i] = i;

      // fill my own row node ids
      const Epetra_Map* noderowmap = discret().node_row_map();
      std::vector<int> sdata(noderowmap->NumMyElements());
      for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

      // gather all gids of nodes redundantly
      std::vector<int> rdata;
      Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), get_comm());

      // build completely overlapping map of nodes (on ALL processors)
      Teuchos::RCP<Epetra_Map> newnodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, get_comm()));
      sdata.clear();
      rdata.clear();

      // fill my own row element ids
      const Epetra_Map* elerowmap = discret().element_row_map();
      sdata.resize(elerowmap->NumMyElements());
      for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

      // gather all gids of elements redundantly
      rdata.resize(0);
      Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), get_comm());

      // build complete overlapping map of elements (on ALL processors)
      Teuchos::RCP<Epetra_Map> newelecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, get_comm()));
      sdata.clear();
      rdata.clear();
      allproc.clear();

      // redistribute the discretization of the interface according to the
      // new column layout
      discret().export_column_nodes(*newnodecolmap);
      discret().export_column_elements(*newelecolmap);

      break;
    }
    case Inpar::Mortar::ExtendGhosting::redundant_master:
    {
      // to ease our search algorithms we'll afford the luxury to ghost all master
      // nodes on all processors. To do so, we'll take the master node row map and
      // export it to full overlap. Then we export the discretization to partially
      // full overlap column map. This way, also all master elements will be fully
      // ghosted on all processors.

      // at least for master, we want to do full ghosting on all procs
      std::vector<int> allproc(get_comm().NumProc());
      for (int i = 0; i < get_comm().NumProc(); ++i) allproc[i] = i;

      // fill my own master row node ids
      const Epetra_Map* noderowmap = discret().node_row_map();
      std::vector<int> sdata;
      for (int i = 0; i < noderowmap->NumMyElements(); ++i)
      {
        int gid = noderowmap->GID(i);
        Core::Nodes::Node* node = discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Mortar::Node* mrtrnode = dynamic_cast<Mortar::Node*>(node);
        if (!mrtrnode->is_slave()) sdata.push_back(gid);
      }

      // gather all master row node gids redundantly
      std::vector<int> rdata;
      Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), get_comm());

      // add my own slave column node ids (non-redundant, standard overlap)
      const Epetra_Map* nodecolmap = discret().node_col_map();
      for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
      {
        int gid = nodecolmap->GID(i);
        Core::Nodes::Node* node = discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Mortar::Node* mrtrnode = dynamic_cast<Mortar::Node*>(node);
        if (mrtrnode->is_slave()) rdata.push_back(gid);
      }

      // build new node column map (on ALL processors)
      Teuchos::RCP<Epetra_Map> newnodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, get_comm()));
      sdata.clear();
      rdata.clear();

      // fill my own master row element ids
      const Epetra_Map* elerowmap = discret().element_row_map();
      sdata.resize(0);
      for (int i = 0; i < elerowmap->NumMyElements(); ++i)
      {
        int gid = elerowmap->GID(i);
        Core::Elements::Element* ele = discret().g_element(gid);
        if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
        Mortar::Element* mrtrele = dynamic_cast<Mortar::Element*>(ele);
        if (!mrtrele->is_slave()) sdata.push_back(gid);
      }

      // gather all gids of elements redundantly
      rdata.resize(0);
      Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), get_comm());

      // add my own slave column node ids (non-redundant, standard overlap)
      const Epetra_Map* elecolmap = discret().element_col_map();
      for (int i = 0; i < elecolmap->NumMyElements(); ++i)
      {
        int gid = elecolmap->GID(i);
        Core::Elements::Element* ele = discret().g_element(gid);
        if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
        Mortar::Element* mrtrele = dynamic_cast<Mortar::Element*>(ele);
        if (mrtrele->is_slave()) rdata.push_back(gid);
      }

      // build new element column map (on ALL processors)
      Teuchos::RCP<Epetra_Map> newelecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, get_comm()));
      sdata.clear();
      rdata.clear();
      allproc.clear();

      // redistribute the discretization of the interface according to the
      // new node / element column layout (i.e. master = full overlap)
      discret().export_column_nodes(*newnodecolmap);
      discret().export_column_elements(*newelecolmap);

      break;
    }
    case Inpar::Mortar::ExtendGhosting::roundrobin:
    {
      // Nothing to do in case of Round-Robin
      break;
    }
    case Inpar::Mortar::ExtendGhosting::binning:
    {
      // Extend master column map via binning

      // Create the binning strategy
      RCP<Core::Binstrategy::BinningStrategy> binningstrategy =
          setup_binning_strategy(meanVelocity);

      // fill master and slave elements into bins
      std::map<int, std::set<int>> slavebinelemap;
      binningstrategy->distribute_eles_to_bins(discret(), slavebinelemap, true);
      std::map<int, std::set<int>> masterbinelemap;
      binningstrategy->distribute_eles_to_bins(discret(), masterbinelemap, false);

      // Extend ghosting of the master elements
      std::map<int, std::set<int>> ext_bin_to_ele_map;
      RCP<const Epetra_Map> extendedmastercolmap =
          binningstrategy->extend_element_col_map(slavebinelemap, masterbinelemap,
              ext_bin_to_ele_map, Teuchos::null, Teuchos::null, discret().element_col_map());

      discret().export_column_elements(*extendedmastercolmap);

      // get the node ids of the elements that are to be ghosted and create a proper node column
      // map for their export
      std::set<int> nodes;
      const int numMasterColElements = extendedmastercolmap->NumMyElements();
      for (int lid = 0; lid < numMasterColElements; ++lid)
      {
        Core::Elements::Element* ele = discret().g_element(extendedmastercolmap->GID(lid));
        const int* nodeids = ele->node_ids();
        for (int inode = 0; inode < ele->num_node(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), colnodes.data(), 0, get_comm()));

      discret().export_column_nodes(*nodecolmap);
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "This case of redundant interface storage has not been covered, yet. Implement it!");
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 |  redistribute contact interface (public)                   popp 08/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::redistribute()
{
  const Teuchos::ParameterList& mortarParallelRedistParams =
      interface_params().sublist("PARALLEL REDISTRIBUTION");

  // make sure we are supposed to be here
  if (Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == Inpar::Mortar::ParallelRedist::redist_none)
    FOUR_C_THROW(
        "You are not supposed to be here since you did not enable PARALLEL_REDIST in the "
        "input file. ");

  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Interface::get_comm().Clone());
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();
  Teuchos::Time time("", true);
  std::set<int>::const_iterator iter;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  //**********************************************************************
  // (1) SLAVE splitting in close / non-close parts
  //**********************************************************************
  // perform contact search (still with non-optimal distribution)
  initialize();
  if (search_alg() == Inpar::Mortar::search_bfele)
    evaluate_search_brute_force(search_param());
  else if (search_alg() == Inpar::Mortar::search_binarytree)
    evaluate_search_binarytree();
  else
    FOUR_C_THROW("Invalid search algorithm");

  // split slave element row map and build redundant vector of
  // all close / non-close slave node ids on all procs
  std::vector<int> closeele, noncloseele;
  std::vector<int> localcns, localfns;

  split_into_far_and_close_sets(closeele, noncloseele, localcns, localfns);

  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  const int numMySlaveColElements = slave_col_elements()->NumMyElements();
  for (int i = 0; i < numMySlaveColElements; ++i)
  {
    int gid = slave_col_elements()->GID(i);
    Core::Elements::Element* ele = discret().g_element(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    Mortar::Element* mele = dynamic_cast<Mortar::Element*>(ele);

    mele->mo_data().search_elements().resize(0);
  }

  // we need an arbitrary preliminary element row map
  Teuchos::RCP<Epetra_Map> slaveCloseRowEles = Teuchos::rcp(
      new Epetra_Map(-1, (int)closeele.size(), closeele.data(), 0, Interface::get_comm()));
  Teuchos::RCP<Epetra_Map> slaveNonCloseRowEles = Teuchos::rcp(
      new Epetra_Map(-1, (int)noncloseele.size(), noncloseele.data(), 0, Interface::get_comm()));
  Teuchos::RCP<Epetra_Map> masterRowEles = Teuchos::rcp(new Epetra_Map(*master_row_elements()));

  // check for consistency
  if (slaveCloseRowEles->NumGlobalElements() == 0 && slaveNonCloseRowEles->NumGlobalElements() == 0)
    FOUR_C_THROW("CONTACT Redistribute: Both slave sets (close/non-close) are empty");

  //**********************************************************************
  // (2) SPECIAL CASES and output to screen
  //**********************************************************************
  // print element overview
  if (!myrank)
  {
    int cl = slaveCloseRowEles->NumGlobalElements();
    int ncl = slaveNonCloseRowEles->NumGlobalElements();
    int ma = masterRowEles->NumGlobalElements();
    std::cout << "Element overview: " << cl << " / " << ncl << " / " << ma
              << "  (close-S / non-close-S / M)";
  }

  // print old parallel distribution
  if (myrank == 0)
    std::cout << "\nInterface parallel distribution before rebalancing:" << std::endl;
  print_parallel_distribution();

  // use simple base class method if there are ONLY close or non-close elements
  // (return value TRUE, because redistribution performed)
  if (slaveCloseRowEles->NumGlobalElements() == 0 || slaveNonCloseRowEles->NumGlobalElements() == 0)
  {
    Mortar::Interface::redistribute();
    return;
  }

  //**********************************************************************
  // (3a) PREPARATIONS decide how many procs are used
  //**********************************************************************
  // first we assume that all procs will be used
  int scproc = numproc;
  int sncproc = numproc;
  int mproc = numproc;

  // minimum number of elements per proc
  const int minele = mortarParallelRedistParams.get<int>("MIN_ELEPROC");

  // Max. relative imbalance between subdomain sizes
  const double imbalance_tol = mortarParallelRedistParams.get<double>("IMBALANCE_TOL");

  // calculate real number of procs to be used
  if (minele > 0)
  {
    scproc = static_cast<int>((slaveCloseRowEles->NumGlobalElements()) / minele);
    sncproc = static_cast<int>((slaveNonCloseRowEles->NumGlobalElements()) / minele);
    mproc = static_cast<int>((masterRowEles->NumGlobalElements()) / minele);
    if (slaveCloseRowEles->NumGlobalElements() < 2 * minele) scproc = 1;
    if (slaveNonCloseRowEles->NumGlobalElements() < 2 * minele) sncproc = 1;
    if (masterRowEles->NumGlobalElements() < 2 * minele) mproc = 1;
    if (scproc > numproc) scproc = numproc;
    if (sncproc > numproc) sncproc = numproc;
    if (mproc > numproc) mproc = numproc;
  }

  // print message
  if (!myrank)
  {
    std::cout << "\nRedistributing interface '" << discret().name() << "' .........\n";
    std::cout << "Procs used for redistribution: " << scproc << " / " << sncproc << " / " << mproc
              << " (close-S / non-close-S / M)\n";
  }

  //**********************************************************************
  // (3b) PREPARATIONS build initial node graph
  //**********************************************************************
  // create graph object
  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *slave_row_nodes(), 108, false));

  // loop over all row nodes to fill graph
  const int numMySlaveRowNodes = slave_row_nodes()->NumMyElements();
  for (int k = 0; k < numMySlaveRowNodes; ++k)
  {
    int gid = slave_row_nodes()->GID(k);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

    // find adjacent elements first
    for (int k = 0; k < node->num_element(); ++k)
    {
      // store adjacent nodes
      Core::Elements::Element* ele = node->elements()[k];
      int numnode = ele->num_node();
      std::vector<int> nodeids(numnode);
      for (int n = 0; n < numnode; ++n) nodeids[n] = ele->node_ids()[n];

      int err = graph->InsertGlobalIndices(gid, numnode, nodeids.data());
      if (err < 0) FOUR_C_THROW("graph->InsertGlobalIndices returned %d", err);
      if (err == 1) FOUR_C_THROW("graph->InsertGlobalIndices returned %d", err);
    }
  }

  // fill graph and optimize storage
  graph->FillComplete();
  graph->OptimizeStorage();

  //**********************************************************************
  // (4) CLOSE SLAVE redistribution
  //**********************************************************************

  // build redundant vector of all close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalcns;
  std::set<int> setglobalcns;
  std::vector<int> scnids;
  Core::LinAlg::Gather<int>(localcns, globalcns, numproc, allproc.data(), Interface::get_comm());
  for (int i = 0; i < (int)globalcns.size(); ++i) setglobalcns.insert(globalcns[i]);
  for (iter = setglobalcns.begin(); iter != setglobalcns.end(); ++iter) scnids.push_back(*iter);

  //**********************************************************************
  // call parallel redistribution
  Teuchos::RCP<const Epetra_CrsGraph> slaveCloseNodeGraph =
      Core::Rebalance::BuildGraph(idiscret_, slaveCloseRowEles);

  Teuchos::ParameterList slaveCloseRebalanceParams;
  slaveCloseRebalanceParams.set<std::string>("num parts", std::to_string(scproc));
  slaveCloseRebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance_tol));

  const auto& [slaveCloseRowNodes, slaveCloseColNodes] =
      Core::Rebalance::RebalanceNodeMaps(slaveCloseNodeGraph, slaveCloseRebalanceParams);
  //**********************************************************************

  //**********************************************************************
  // (5) NON-CLOSE SLAVE redistribution
  //**********************************************************************

  // build redundant vector of all non-close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalfns;
  std::set<int> setglobalfns;
  std::vector<int> sncnids;
  Core::LinAlg::Gather<int>(localfns, globalfns, numproc, allproc.data(), Interface::get_comm());
  for (int i = 0; i < (int)globalfns.size(); ++i) setglobalfns.insert(globalfns[i]);
  for (iter = setglobalfns.begin(); iter != setglobalfns.end(); ++iter) sncnids.push_back(*iter);

  //**********************************************************************
  // call parallel redistribution
  Teuchos::RCP<const Epetra_CrsGraph> slaveNonCloseNodeGraph =
      Core::Rebalance::BuildGraph(idiscret_, slaveNonCloseRowEles);

  Teuchos::ParameterList slaveNonCloseRebalanceParams;
  slaveNonCloseRebalanceParams.set<std::string>("num parts", std::to_string(sncproc));
  slaveNonCloseRebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance_tol));

  const auto& [slaveNonCloseRowNodes, snccolnodes] =
      Core::Rebalance::RebalanceNodeMaps(slaveNonCloseNodeGraph, slaveNonCloseRebalanceParams);
  //**********************************************************************

  //**********************************************************************
  // (6) MASTER redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> mrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> mcolnodes = Teuchos::null;

  redistribute_master_side(mrownodes, mcolnodes, masterRowEles, comm, mproc, imbalance_tol);

  //**********************************************************************
  // (7) Merge global interface node row and column map
  //**********************************************************************
  // merge slave node row map from close and non-close parts
  Teuchos::RCP<Epetra_Map> srownodes = Teuchos::null;

  //----------------------------------CASE 1: ONE OR BOTH SLAVE SETS EMPTY
  if (slaveCloseRowNodes == Teuchos::null || slaveNonCloseRowNodes == Teuchos::null)
  {
    FOUR_C_THROW("CONTACT Redistribute: Both slave sets (close/non-close) are empty");
  }
  //-------------------------------------CASE 2: BOTH SLAVE SETS NON-EMPTY
  else
  {
    // find intersection set of close and non-close nodes
    std::set<int> intersec;
    for (iter = setglobalcns.begin(); iter != setglobalcns.end(); ++iter)
    {
      std::set<int>::const_iterator found = setglobalfns.find(*iter);
      if (found != setglobalfns.end()) intersec.insert(*found);
    }

    // build slave node row map
    const int numMySlaveCloseRowNodes = slaveCloseRowNodes->NumMyElements();
    const int numMySlaveNonCloseRowNodes = slaveNonCloseRowNodes->NumMyElements();
    std::vector<int> mygids(numMySlaveCloseRowNodes + numMySlaveNonCloseRowNodes);
    int count = slaveCloseRowNodes->NumMyElements();

    // first get GIDs of input slaveCloseRowNodes
    for (int i = 0; i < count; ++i) mygids[i] = slaveCloseRowNodes->GID(i);

    // then add GIDs of input slaveNonCloseRowNodes (only new ones)
    for (int i = 0; i < numMySlaveNonCloseRowNodes; ++i)
    {
      // check for intersection gid
      // don't do anything for intersection gids (slaveCloseRowNodes dominates!!!)
      std::set<int>::const_iterator found = intersec.find(slaveNonCloseRowNodes->GID(i));
      if (found != intersec.end()) continue;

      // check for overlap
      if (slaveCloseRowNodes->MyGID(slaveNonCloseRowNodes->GID(i)))
        FOUR_C_THROW("Core::LinAlg::MergeMap: Result map is overlapping");

      // add new GIDs to mygids
      mygids[count] = slaveNonCloseRowNodes->GID(i);
      ++count;
    }
    mygids.resize(count);
    sort(mygids.begin(), mygids.end());
    srownodes = Teuchos::rcp(
        new Epetra_Map(-1, (int)mygids.size(), mygids.data(), 0, slaveCloseRowNodes->Comm()));
  }

  // merge interface node row map from slave and master parts
  Teuchos::RCP<Epetra_Map> rownodes = Core::LinAlg::MergeMap(srownodes, mrownodes, false);

  // IMPORTANT NOTE:
  // While merging from the two different slave parts of the discretization
  // (close slave, non-close slave) is feasible for the node row map,
  // this is not possible for the node column map. Some necessary
  // information on ghosting at the transition between close and non-close
  // slave region would always be missed! Thus, we reconstruct a
  // suitable slave node column map "by hand" here. This is quite simply
  // done by exporting the initial node graph to the new distribution
  // and by then asking for its column map.

  // create the output graph (with new slave node row map) and export to it
  Teuchos::RCP<Epetra_CrsGraph> outgraph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *srownodes, 108, false));
  Epetra_Export exporter(graph->RowMap(), *srownodes);
  int err = outgraph->Export(*graph, exporter, Add);
  if (err < 0) FOUR_C_THROW("Graph export returned err=%d", err);

  // trash old graph
  graph = Teuchos::null;

  // call fill complete and optimize storage
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  // get column map from the graph -> build slave node column map
  // (do stupid conversion from Epetra_BlockMap to Epetra_Map)
  const Epetra_BlockMap& bcol = outgraph->ColMap();
  Teuchos::RCP<Epetra_Map> scolnodes = Teuchos::rcp(new Epetra_Map(bcol.NumGlobalElements(),
      bcol.NumMyElements(), bcol.MyGlobalElements(), 0, Interface::get_comm()));

  // trash new graph
  outgraph = Teuchos::null;

  // merge interface node column map from slave and master parts
  Teuchos::RCP<Epetra_Map> colnodes = Core::LinAlg::MergeMap(scolnodes, mcolnodes, false);

  //**********************************************************************
  // (8) Get partitioning information into discretization
  //**********************************************************************
  // build reasonable element maps from the already valid and final node maps
  // (note that nothing is actually redistributed in here)
  const auto& [roweles, coleles] = discret().build_element_row_column(*rownodes, *colnodes);

  // export nodes and elements to the row map
  discret().export_row_nodes(*rownodes);
  discret().export_row_elements(*roweles);

  // export nodes and elements to the column map (create ghosting)
  discret().export_column_nodes(*colnodes);
  discret().export_column_elements(*coleles);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::split_into_far_and_close_sets(std::vector<int>& closeele,
    std::vector<int>& noncloseele, std::vector<int>& localcns, std::vector<int>& localfns) const
{
  const bool performSplitting = Core::UTILS::IntegralValue<bool>(
      interface_params().sublist("PARALLEL REDISTRIBUTION"), "EXPLOIT_PROXIMITY");

  if (performSplitting)
  {
    // loop over all row elements to gather the local information
    for (int i = 0; i < slave_row_elements()->NumMyElements(); ++i)
    {
      // get element
      int gid = slave_row_elements()->GID(i);
      Core::Elements::Element* ele = discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      Mortar::Element* cele = dynamic_cast<Mortar::Element*>(ele);

      // store element id and adjacent node ids
      int close = cele->mo_data().num_search_elements();
      if (close > 0)
      {
        closeele.push_back(gid);
        for (int k = 0; k < cele->num_node(); ++k) localcns.push_back(cele->node_ids()[k]);
      }
      else
      {
        noncloseele.push_back(gid);
        for (int k = 0; k < cele->num_node(); ++k) localfns.push_back(cele->node_ids()[k]);
      }
    }
  }
  else
  {
    // loop over all row elements to gather the local information
    for (int i = 0; i < slave_row_elements()->NumMyElements(); ++i)
    {
      // get element
      int gid = slave_row_elements()->GID(i);
      Core::Elements::Element* ele = discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      Mortar::Element* cele = dynamic_cast<Mortar::Element*>(ele);

      // store element id and adjacent node ids
      noncloseele.push_back(gid);
      for (int k = 0; k < cele->num_node(); ++k) localfns.push_back(cele->node_ids()[k]);
    }
  }
}

/*----------------------------------------------------------------------*
 | collect distribution data (public)                         popp 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::collect_distribution_data(int& numColElements, int& numRowElements)
{
  // loop over proc's column slave elements of the interface
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* slaveElement = dynamic_cast<Element*>(ele1);

    // bool indicating coupling partners
    bool add = (slaveElement->mo_data().num_search_elements() > 0);

    // Check if this element has any coupling partners.
    // Increment element counter if so.
    if (add) ++numColElements;

    // check if - in addition - the active proc owns this element.
    // Increment input variable rowele if so.
    if (add && slaveElement->owner() == get_comm().MyPID()) ++numRowElements;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::create_search_tree()
{
  // warning
#ifdef MORTARGMSHCTN
  if (Dim() == 3 && Comm().MyPID() == 0)
  {
    std::cout << "\n******************************************************************\n";
    std::cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    std::cout << "******************************************************************\n";
  }
#endif

  // binary tree search
  if (search_alg() == Inpar::Mortar::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (self_contact())
    {
      // set state in interface to intialize all kinds of quantities
      Teuchos::RCP<Epetra_Vector> zero = Teuchos::rcp(new Epetra_Vector(*idiscret_->dof_row_map()));
      set_state(Mortar::state_new_displacement, *zero);

      // create fully overlapping map of all contact elements
      Teuchos::RCP<Epetra_Map> elefullmap =
          Core::LinAlg::AllreduceEMap(*idiscret_->element_row_map());

      // create binary tree object for self contact search
      if (!two_half_pass())
      {
        // (NOTE THAT SELF CONTACT SEARCH IS NOT YET FULLY PARALLELIZED!)
        binarytreeself_ = Teuchos::rcp(new CONTACT::SelfBinaryTree(
            discret(), interface_params(), elefullmap, n_dim(), search_param()));
      }
      else
      {
        // if we use the two half pass algorithm, we use the unbiased self binary tree
        // implementation
        binarytreeself_ = Teuchos::rcp(new CONTACT::UnbiasedSelfBinaryTree(
            discret(), interface_params(), elefullmap, n_dim(), search_param()));
      }
      // initialize the self binary tree
      binarytreeself_->init();
    }
    //*****TWO BODY CONTACT*****
    else
    {
      Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
      switch (interface_data_->get_extend_ghosting())
      {
        case Inpar::Mortar::ExtendGhosting::roundrobin:
        case Inpar::Mortar::ExtendGhosting::binning:
        {
          melefullmap = melecolmap_;
          break;
        }
        case Inpar::Mortar::ExtendGhosting::redundant_all:
        case Inpar::Mortar::ExtendGhosting::redundant_master:
        {
          melefullmap = Core::LinAlg::AllreduceEMap(*melerowmap_);
          break;
        }
        default:
        {
          FOUR_C_THROW("Chosen parallel strategy not supported!");
          break;
        }
      }

      {
        // get update type of binary tree
        Inpar::Mortar::BinaryTreeUpdateType updatetype =
            Core::UTILS::IntegralValue<Inpar::Mortar::BinaryTreeUpdateType>(
                interface_params(), "BINARYTREE_UPDATETYPE");

        // create binary tree object for contact search and setup tree
        binarytree_ = Teuchos::rcp(new Mortar::BinaryTree(discret(), selecolmap_, melefullmap,
            n_dim(), search_param(), updatetype, search_use_aux_pos()));
        // initialize the binary tree
        binarytree_->init();
      }
    }
  }

  // no binary tree search
  else
  {
    if (self_contact()) FOUR_C_THROW("Binarytree search needed for self contact");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::initialize_data_container()
{
  // call base class functionality
  Mortar::Interface::initialize_data_container();

  // ==================
  // non-smooth contact:
  // we need this master node data container to create an averaged
  // nodal normal field on the master side for the smoothed cpp
  // normal field!
  if (Core::UTILS::IntegralValue<int>(interface_params(), "CPP_NORMALS") || nonSmoothContact_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = Core::LinAlg::AllreduceEMap(*(master_row_nodes()));

    for (int i = 0; i < masternodes->NumMyElements(); ++i)
    {
      int gid = masternodes->GID(i);
      Core::Nodes::Node* node = discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %i", gid);
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
      mnode->initialize_data_container();
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::initialize()
{
  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i = 0; i < idiscret_->num_my_col_nodes(); ++i)
  {
    CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_->l_col_node(i));

    // reset feasible projection and segmentation status
    node->has_proj() = false;
    node->has_segment() = false;
  }

  // init normal data in master node data container for cpp calculation
  if (Core::UTILS::IntegralValue<int>(interface_params(), "CPP_NORMALS"))
  {
    for (int i = 0; i < master_col_nodes()->NumMyElements(); ++i)
    {
      int gid = master_col_nodes()->GID(i);
      Core::Nodes::Node* node = discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // reset derivative maps of normal vector
      for (int j = 0; j < (int)((cnode->data().get_deriv_n()).size()); ++j)
        (cnode->data().get_deriv_n())[j].clear();
      (cnode->data().get_deriv_n()).resize(0, 0);

      // reset derivative maps of tangent vectors
      for (int j = 0; j < (int)((cnode->data().get_deriv_txi()).size()); ++j)
        (cnode->data().get_deriv_txi())[j].clear();
      (cnode->data().get_deriv_txi()).resize(0, 0);
      for (int j = 0; j < (int)((cnode->data().get_deriv_teta()).size()); ++j)
        (cnode->data().get_deriv_teta())[j].clear();
      (cnode->data().get_deriv_teta()).resize(0, 0);

      for (int j = 0; j < (int)((cnode->data().get_deriv_tangent()).size()); ++j)
        (cnode->data().get_deriv_tangent())[j].clear();
      (cnode->data().get_deriv_tangent()).resize(0, 0);
    }
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i = 0; i < slave_col_nodes_bound()->NumMyElements(); ++i)
  {
    int gid = slave_col_nodes_bound()->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // reset nodal Mortar maps
    // for sts
    cnode->mo_data().get_d().clear();
    cnode->mo_data().get_m().clear();
    cnode->mo_data().get_mmod().clear();
    // for nts
    cnode->mo_data().get_dnts().clear();
    cnode->mo_data().get_mnts().clear();
    // for lts
    cnode->mo_data().get_dlts().clear();
    cnode->mo_data().get_mlts().clear();
    // for ltl
    cnode->mo_data().get_dltl().clear();
    cnode->mo_data().get_mltl().clear();

    // reset derivative maps of normal vector
    for (int j = 0; j < (int)((cnode->data().get_deriv_n()).size()); ++j)
      (cnode->data().get_deriv_n())[j].clear();
    (cnode->data().get_deriv_n()).resize(0, 0);

    // reset derivative maps of tangent vectors
    for (int j = 0; j < (int)((cnode->data().get_deriv_txi()).size()); ++j)
      (cnode->data().get_deriv_txi())[j].clear();
    (cnode->data().get_deriv_txi()).resize(0, 0);
    for (int j = 0; j < (int)((cnode->data().get_deriv_teta()).size()); ++j)
      (cnode->data().get_deriv_teta())[j].clear();
    (cnode->data().get_deriv_teta()).resize(0, 0);

    // reset derivative map of Mortar matrices
    (cnode->data().get_deriv_d()).clear();
    (cnode->data().get_deriv_dlts()).clear();
    (cnode->data().get_deriv_dltl()).clear();
    (cnode->data().get_deriv_m()).clear();
    (cnode->data().get_deriv_mnts()).clear();
    (cnode->data().get_deriv_mlts()).clear();
    (cnode->data().get_deriv_mltl()).clear();

    // reset nodal weighted gap and derivative
    cnode->data().getg() = 1.0e12;
    cnode->data().getgnts() = 1.0e12;
    cnode->data().getglts() = 1.0e12;
    cnode->data().getgltl()[0] = 1.0e12;
    cnode->data().getgltl()[1] = 1.0e12;
    cnode->data().getgltl()[2] = 1.0e12;
    cnode->mo_data().get_dscale() = 0.0;
    (cnode->data().get_deriv_g()).clear();
    (cnode->data().get_deriv_gnts()).clear();
    (cnode->data().get_deriv_glts()).clear();
    for (int j = 0; j < (int)cnode->data().get_deriv_gltl().size(); ++j)
      cnode->data().get_deriv_gltl()[j].clear();
    for (int j = 0; j < (int)cnode->data().get_deriv_jumpltl().size(); ++j)
      cnode->data().get_deriv_jumpltl()[j].clear();
    //    (cnode->Data().GetDerivGltl()).resize(0);

    // reset nodal jump
    cnode->data().getjumpltl()[0] = 1.0e12;
    cnode->data().getjumpltl()[1] = 1.0e12;
    cnode->data().getjumpltl()[2] = 1.0e12;

    // hybrid formulation
    cnode->data().get_alpha_n() = -1.0;
    cnode->data().get_alpha().clear();

    // reset derivative map of lagrange multipliers
    for (int j = 0; j < (int)((cnode->data().get_deriv_z()).size()); ++j)
      (cnode->data().get_deriv_z())[j].clear();
    (cnode->data().get_deriv_z()).resize(0);

    if (friction_)
    {
      FriNode* frinode = dynamic_cast<FriNode*>(cnode);

      // reset SNodes and Mnodes
      frinode->fri_data().get_s_nodes().clear();
      frinode->fri_data().get_m_nodes().clear();

      // for gp slip
      if (Core::UTILS::IntegralValue<int>(interface_params(), "GP_SLIP_INCR") == true)
      {
        // reset jump deriv.
        for (int j = 0; j < (int)((frinode->fri_data().get_deriv_var_jump()).size()); ++j)
          (frinode->fri_data().get_deriv_var_jump())[j].clear();

        (frinode->fri_data().get_deriv_var_jump()).resize(2);

        // reset jumps
        frinode->fri_data().jump_var()[0] = 0.0;
        frinode->fri_data().jump_var()[1] = 0.0;
      }
    }

    // just do poro contact relevant stuff!
    if (interface_data().is_poro())
    {
      cnode->poro_data().getn_coup() = 0.0;
      cnode->poro_data().get_derivn_coup().clear();
      cnode->poro_data().get_vel_derivn_coup().clear();
      cnode->poro_data().get_pres_derivn_coup().clear();
    }

    // just do ehl relevant stuff!
    if (ehl_) cnode->ehl_data().clear();
  }

  //**********************************************************************
  // In general, it is sufficient to reset search candidates only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need to reset the search candidates of
  // all slave elements in the fully overlapping column map there. This
  // is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (self_contact())
  {
    // loop over all elements to reset candidates / search lists
    // (use fully overlapping column map of S+M elements)
    for (int i = 0; i < idiscret_->num_my_col_elements(); ++i)
    {
      Core::Elements::Element* ele = idiscret_->l_col_element(i);
      Mortar::Element* mele = dynamic_cast<Mortar::Element*>(ele);

      mele->mo_data().search_elements().resize(0);

      // dual shape function coefficient matrix
      mele->mo_data().reset_dual_shape();
      mele->mo_data().reset_deriv_dual_shape();
    }
  }
  else
  {
    // loop over all elements to reset candidates / search lists
    // (use standard slave column map)
    for (int i = 0; i < slave_col_elements()->NumMyElements(); ++i)
    {
      int gid = slave_col_elements()->GID(i);
      Core::Elements::Element* ele = discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
      Mortar::Element* mele = dynamic_cast<Mortar::Element*>(ele);

      mele->mo_data().search_elements().resize(0);

      // dual shape function coefficient matrix
      mele->mo_data().reset_dual_shape();
      mele->mo_data().reset_deriv_dual_shape();
    }
  }

  // clear all Nitsche data
  if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM") ==
      Inpar::Mortar::algorithm_gpts)
    for (int e = 0; e < discret().element_col_map()->NumMyElements(); ++e)
      dynamic_cast<Mortar::Element*>(discret().g_element(discret().element_col_map()->GID(e)))
          ->get_nitsche_container()
          .clear();

  // reset s/m pairs and intcell counters
  smpairs_ = 0;
  smintpairs_ = 0;
  intcells_ = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::set_element_areas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (self_contact() or Core::UTILS::IntegralValue<int>(interface_params(), "CPP_NORMALS") or
      nonSmoothContact_)
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i = 0; i < idiscret_->num_my_col_elements(); ++i)
    {
      Mortar::Element* element = dynamic_cast<Mortar::Element*>(idiscret_->l_col_element(i));
      element->initialize_data_container();
      element->mo_data().area() = element->compute_area();
    }
  }
  else
  {
    // refer call back to base class version
    Mortar::Interface::set_element_areas();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  pre evaluate to calc normals                            farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::pre_evaluate(const int& step, const int& iter)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::Interface::pre_evaluate");

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (search_alg() == Inpar::Mortar::search_bfele)
    evaluate_search_brute_force(search_param());
  else if (search_alg() == Inpar::Mortar::search_binarytree)
    evaluate_search_binarytree();
  else
    FOUR_C_THROW("Invalid search algorithm");

    // TODO: maybe we can remove this debug functionality
#ifdef MORTARGMSHCELLS
  // reset integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "w");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \"Integration Cells Proc " << proc << "\" {" << std::endl;
  fprintf(fp, gmshfilecontent.str().c_str());
  fclose(fp);
#endif  // #ifdef MORTARGMSHCELLS

  // set global vector of cn values
  set_cn_ct_values(iter);

  // cpp normals or averaged normal field?
  if (Core::UTILS::IntegralValue<int>(interface_params(), "CPP_NORMALS"))
  {
    // evaluate cpp nodal normals on slave side
    evaluate_cpp_normals();
  }
  else
  {
    // evaluate averaged nodal normals on slave side
    evaluate_nodal_normals();

    // export nodal normals to slave node column map
    // this call is very expensive and the computation
    // time scales directly with the proc number !
    export_nodal_normals();
  }

  // set condition specific parameters needed for evaluation
  set_condition_specific_parameters();

  // compute scaling between coupling types
  //  if(nonSmoothContact_)
  //    compute_scaling();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  store nts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::store_nt_svalues()
{
  // create iterators for data types
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all possibly non smooth nodes
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if integration is done
    if (cnode->mo_data().get_dnts().size() < 1) continue;

    // if nonsmooth contact is activated and the node is no corner node continue
    // if non-smooth contact is not activated go on
    if (!cnode->is_on_corner() and nonSmoothContact_) continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize Pairedvector to nts size
    if ((int)cnode->mo_data().get_d().size() == 0)
      cnode->mo_data().get_d().resize(cnode->mo_data().get_dnts().size());

    for (CI p = cnode->mo_data().get_dnts().begin(); p != cnode->mo_data().get_dnts().end(); ++p)
      cnode->mo_data().get_d()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p = cnode->mo_data().get_mnts().begin(); p != cnode->mo_data().get_mnts().end(); ++p)
      cnode->mo_data().get_m()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->data().getg() = cnode->data().getgnts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p = cnode->data().get_deriv_gnts().begin();
         p != cnode->data().get_deriv_gnts().end(); ++p)
      cnode->data().get_deriv_g()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    // --> No D linearization!

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->data().get_deriv_mnts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->data().get_deriv_mnts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->data().get_deriv_m()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          FOUR_C_THROW("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
  }  // end node loop


  return;
}


/*----------------------------------------------------------------------*
 |  store lts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::store_lt_svalues()
{
  // create iterators for data types
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all possibly non smooth nodes
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    double msum = 0.0;
    double ssum = 0.0;

    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is an edge or a corner and nonsmooth contact is activated
    if (!cnode->is_on_edge() and nonSmoothContact_) continue;
    if (cnode->is_on_corner() and nonSmoothContact_) continue;

    // check if integration is done
    if (cnode->mo_data().get_dlts().size() < 1) continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize Pairedvector to nts size
    if ((int)cnode->mo_data().get_d().size() == 0)
      cnode->mo_data().get_d().resize(
          cnode->mo_data().get_dlts().size() + cnode->mo_data().get_dltl().size());

    for (CI p = cnode->mo_data().get_dlts().begin(); p != cnode->mo_data().get_dlts().end(); ++p)
    {
      cnode->mo_data().get_d()[p->first] += (p->second);
      ssum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p = cnode->mo_data().get_mlts().begin(); p != cnode->mo_data().get_mlts().end(); ++p)
    {
      cnode->mo_data().get_m()[p->first] += (p->second);
      msum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->data().getg() = cnode->data().getglts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p = cnode->data().get_deriv_glts().begin();
         p != cnode->data().get_deriv_glts().end(); ++p)
      cnode->data().get_deriv_g()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->data().get_deriv_dlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->data().get_deriv_dlts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->data().get_deriv_d()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          FOUR_C_THROW("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->data().get_deriv_mlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->data().get_deriv_mlts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->data().get_deriv_m()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          FOUR_C_THROW("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
    //    std::cout << "ssum = " << ssum << "  msum = " << msum << "  balance= " << ssum-msum <<
    //    std::endl;
    if (abs(ssum - msum) > 1e-12) FOUR_C_THROW("no slave master balance!");

  }  // end node loop


  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_lforces_fric(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");
  const double penaltytan = interface_params().get<double>("PENALTYPARAMTAN");
  const double frcoeff = interface_params().get<double>("FRCOEFF");

  std::array<double, 3> oldtraction = {0.0, 0.0, 0.0};

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->fri_data().tractionold_ltl()[0] * cnode->fri_data().tractionold_ltl()[0];
    double y = cnode->fri_data().tractionold_ltl()[1] * cnode->fri_data().tractionold_ltl()[1];
    double z = cnode->fri_data().tractionold_ltl()[2] * cnode->fri_data().tractionold_ltl()[2];
    double tracvalue = sqrt(x + y + z);

    if (tracvalue > 1e-8)
    {
      oldtraction[0] = cnode->fri_data().tractionold_ltl()[0];
      oldtraction[1] = cnode->fri_data().tractionold_ltl()[1];
      oldtraction[2] = cnode->fri_data().tractionold_ltl()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if (cnode->data().getgltl()[0] < 1e8 and cnode->data().getgltl()[1] < 1e8 and
        cnode->data().getgltl()[2] < 1e8)
    {
      // normal force
      std::array<double, 3> fn = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < n_dim(); ++dim) fn[dim] = -penalty * cnode->data().getgltl()[dim];

      // f trial tangential
      std::array<double, 3> ftrial = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < n_dim(); ++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->data().getjumpltl()[dim];

      // trial norm
      double trialnorm =
          sqrt(ftrial[0] * ftrial[0] + ftrial[1] * ftrial[1] + ftrial[2] * ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);

      // real traction
      std::array<double, 3> ftan = {0.0, 0.0, 0.0};

      if (trialnorm - frcoeff * maxtrac <= 0.0)
      {
        for (int dim = 0; dim < n_dim(); ++dim) ftan[dim] = ftrial[dim];
      }
      else
      {
        double coeff = frcoeff * maxtrac / trialnorm;
        for (int dim = 0; dim < n_dim(); ++dim) ftan[dim] = coeff * ftrial[dim];
      }

      // store
      cnode->fri_data().traction()[0] = ftan[0];
      cnode->fri_data().traction()[1] = ftan[1];
      cnode->fri_data().traction()[2] = ftan[2];

      // ASSEMBLE
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dltl()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value = (p->second) * ftan[dim];
            const int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->mo_data().get_mltl()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value = -(p->second) * ftan[dim];
            const int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }

      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_lstiffness_fric(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");
  const double penaltytan = interface_params().get<double>("PENALTYPARAMTAN");
  const double frcoeff = interface_params().get<double>("FRCOEFF");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  std::array<double, 3> oldtraction = {0.0, 0.0, 0.0};

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->fri_data().tractionold_ltl()[0] * cnode->fri_data().tractionold_ltl()[0];
    double y = cnode->fri_data().tractionold_ltl()[1] * cnode->fri_data().tractionold_ltl()[1];
    double z = cnode->fri_data().tractionold_ltl()[2] * cnode->fri_data().tractionold_ltl()[2];
    double tracvalue = sqrt(x + y + z);

    if (tracvalue > 1e-8)
    {
      oldtraction[0] = cnode->fri_data().tractionold_ltl()[0];
      oldtraction[1] = cnode->fri_data().tractionold_ltl()[1];
      oldtraction[2] = cnode->fri_data().tractionold_ltl()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if (cnode->data().getgltl()[0] < 1e8 and cnode->data().getgltl()[1] < 1e8 and
        cnode->data().getgltl()[2] < 1e8)
    {
      // state
      bool stick = true;
      Core::Gen::Pairedvector<int, double> coefflin(100);

      // normal force
      std::array<double, 3> fn = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < n_dim(); ++dim) fn[dim] = -penalty * cnode->data().getgltl()[dim];

      // f trial tangential
      std::array<double, 3> ftrial = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < n_dim(); ++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->data().getjumpltl()[dim];

      double coeff = 0.0;

      // trial norm
      double trialnorm =
          sqrt(ftrial[0] * ftrial[0] + ftrial[1] * ftrial[1] + ftrial[2] * ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);

      // real traction
      std::array<double, 3> ftan = {0.0, 0.0, 0.0};

      if (trialnorm - frcoeff * maxtrac <= 0.0)
      {
        stick = true;
        for (int dim = 0; dim < n_dim(); ++dim) ftan[dim] = ftrial[dim];
      }
      else
      {
        stick = false;
        coeff = frcoeff * maxtrac / trialnorm;
        for (int dim = 0; dim < n_dim(); ++dim) ftan[dim] = coeff * ftrial[dim];

        Core::Gen::Pairedvector<int, double> fn_x(100);
        Core::Gen::Pairedvector<int, double> fn_y(100);
        Core::Gen::Pairedvector<int, double> fn_z(100);

        Core::Gen::Pairedvector<int, double> ft_x(100);
        Core::Gen::Pairedvector<int, double> ft_y(100);
        Core::Gen::Pairedvector<int, double> ft_z(100);

        for (CImap pp = cnode->data().get_deriv_gltl()[0].begin();
             pp != cnode->data().get_deriv_gltl()[0].end(); ++pp)
          fn_x[pp->first] -= penalty * (pp->second);
        for (CImap pp = cnode->data().get_deriv_gltl()[1].begin();
             pp != cnode->data().get_deriv_gltl()[1].end(); ++pp)
          fn_y[pp->first] -= penalty * (pp->second);
        for (CImap pp = cnode->data().get_deriv_gltl()[2].begin();
             pp != cnode->data().get_deriv_gltl()[2].end(); ++pp)
          fn_z[pp->first] -= penalty * (pp->second);

        for (CImap pp = cnode->data().get_deriv_jumpltl()[0].begin();
             pp != cnode->data().get_deriv_jumpltl()[0].end(); ++pp)
          ft_x[pp->first] -= penaltytan * (pp->second);
        for (CImap pp = cnode->data().get_deriv_jumpltl()[1].begin();
             pp != cnode->data().get_deriv_jumpltl()[1].end(); ++pp)
          ft_y[pp->first] -= penaltytan * (pp->second);
        for (CImap pp = cnode->data().get_deriv_jumpltl()[2].begin();
             pp != cnode->data().get_deriv_jumpltl()[2].end(); ++pp)
          ft_z[pp->first] -= penaltytan * (pp->second);

        Core::Gen::Pairedvector<int, double> maxtraclin(100);
        for (CI pp = fn_x.begin(); pp != fn_x.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[0] * pp->second;
        for (CI pp = fn_y.begin(); pp != fn_y.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[1] * pp->second;
        for (CI pp = fn_z.begin(); pp != fn_z.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[2] * pp->second;

        Core::Gen::Pairedvector<int, double> trialnormlin(100);
        for (CI pp = ft_x.begin(); pp != ft_x.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[0] * pp->second;
        for (CI pp = ft_y.begin(); pp != ft_y.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[1] * pp->second;
        for (CI pp = ft_z.begin(); pp != ft_z.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[2] * pp->second;

        for (CI pp = maxtraclin.begin(); pp != maxtraclin.end(); ++pp)
          coefflin[pp->first] += frcoeff * pp->second * (1.0 / trialnorm);

        for (CI pp = trialnormlin.begin(); pp != trialnormlin.end(); ++pp)
          coefflin[pp->first] -= frcoeff * pp->second * (1.0 / (trialnorm * trialnorm)) * maxtrac;
      }


      std::map<int, std::map<int, double>>& dderiv = cnode->data().get_deriv_dltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        Core::Nodes::Node* snode = idiscret_->g_node(sgid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->data().get_deriv_dltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = csnode->dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = ftan[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->fe_assemble(val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->data().get_deriv_mltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        Core::Nodes::Node* mnode = idiscret_->g_node(mgid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->data().get_deriv_mltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = cmnode->dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = ftan[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->fe_assemble(-val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      // ****************************************************************
      // **************************************************** stick state
      // ****************************************************************
      if (stick)
      {
        /**************************************************** D-matrix ******/
        if ((cnode->mo_data().get_dltl()).size() > 0)
        {
          Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CImap pp = cnode->data().get_deriv_jumpltl()[dim].begin();
                   pp != cnode->data().get_deriv_jumpltl()[dim].end(); ++pp)
              {
                double value = penaltytan * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no d matrix entries available for ltlt contact");
        }
        /**************************************************** D-matrix ******/
        if ((cnode->mo_data().get_mltl()).size() > 0)
        {
          std::map<int, double> map = cnode->mo_data().get_mltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CImap pp = cnode->data().get_deriv_jumpltl()[dim].begin();
                   pp != cnode->data().get_deriv_jumpltl()[dim].end(); ++pp)
              {
                double value = -penaltytan * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no m matrix entries available for ltlt contact");
        }
      }
      // ****************************************************************
      // **************************************************** slip state
      // ****************************************************************
      else
      {
        /**************************************************** D-matrix ******/
        if ((cnode->mo_data().get_dltl()).size() > 0)
        {
          Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CImap pp = cnode->data().get_deriv_jumpltl()[dim].begin();
                   pp != cnode->data().get_deriv_jumpltl()[dim].end(); ++pp)
              {
                double value = penaltytan * coeff * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->mo_data().get_dltl()).size() > 0)
        {
          Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CI pp = coefflin.begin(); pp != coefflin.end(); ++pp)
              {
                double value = -penaltytan * ftan[dim] * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->mo_data().get_mltl()).size() > 0)
        {
          std::map<int, double> map = cnode->mo_data().get_mltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CImap pp = cnode->data().get_deriv_jumpltl()[dim].begin();
                   pp != cnode->data().get_deriv_jumpltl()[dim].end(); ++pp)
              {
                double value = -penaltytan * coeff * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no m matrix entries available for ltlt contact");
        }
        if ((cnode->mo_data().get_mltl()).size() > 0)
        {
          std::map<int, double> map = cnode->mo_data().get_mltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            Core::Nodes::Node* snode = idiscret_->g_node(gid3);
            if (!snode) FOUR_C_THROW("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < n_dim(); ++dim)
            {
              for (CI pp = coefflin.begin(); pp != coefflin.end(); ++pp)
              {
                double value = penaltytan * ftan[dim] * (p->second) * (pp->second);
                kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          FOUR_C_THROW("no m matrix entries available for ltlt contact");
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add nts penalty forces master                           farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_nt_sforces_master(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for corners
    if (!cnode->is_on_corner()) continue;

    // is gap is in contact
    if (cnode->data().getgnts() < 1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dnts()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dnts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value =
                penalty * (p->second) * cnode->data().getgnts() * cnode->mo_data().n()[dim];
            int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->mo_data().get_mnts()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mnts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value =
                -penalty * (p->second) * cnode->data().getgnts() * cnode->mo_data().n()[dim];
            int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces master                  farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_sforces_master(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for edges (without corner)
    if (!cnode->is_on_corner_edge()) continue;
    if (cnode->is_on_corner()) continue;

    // scale penalty
    double penaltyLts = penalty * cnode->data().kappa();

    // is gap is in contact
    if (cnode->data().getglts() < 1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dlts()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dlts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value =
                penaltyLts * (p->second) * cnode->data().getglts() * cnode->mo_data().n()[dim];
            int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->mo_data().get_mlts()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mlts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value =
                -penaltyLts * (p->second) * cnode->data().getglts() * cnode->mo_data().n()[dim];
            int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_lforces(Teuchos::RCP<Epetra_FEVector> feff)
{
  // gap = g_n * n
  // D/M = sval/mval
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is valid node
    if (cnode->data().getgltl()[0] < 1e8 and cnode->data().getgltl()[1] < 1e8 and
        cnode->data().getgltl()[2] < 1e8)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dltl()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value = penalty * (p->second) * cnode->data().getgltl()[dim];
            int ltlid = csnode->dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->mo_data().get_mltl()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            double value = -penalty * (p->second) * cnode->data().getgltl()[dim];
            int ltlid = {csnode->dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) FOUR_C_THROW("stop");
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_sstiffness_master(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for edges (without corner)
    if (!cnode->is_on_corner_edge()) continue;
    if (cnode->is_on_corner()) continue;

    // scale penalty
    double penaltyLts = penalty * cnode->data().kappa();

    // is gap is in contact
    if (cnode->data().getglts() < 1e-12)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penaltyLts * cnode->data().getglts() * cnode->mo_data().n()[0];
      lm[1] = penaltyLts * cnode->data().getglts() * cnode->mo_data().n()[1];
      lm[2] = penaltyLts * cnode->data().getglts() * cnode->mo_data().n()[2];

      std::map<int, std::map<int, double>>& dderiv = cnode->data().get_deriv_dlts();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        Core::Nodes::Node* snode = idiscret_->g_node(sgid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->data().get_deriv_dlts()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = csnode->dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->fe_assemble(-val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->data().get_deriv_mlts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        Core::Nodes::Node* mnode = idiscret_->g_node(mgid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->data().get_deriv_mlts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = cmnode->dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->fe_assemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dlts()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dlts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            // gap linearization
            for (CImap pp = cnode->data().get_deriv_glts().begin();
                 pp != cnode->data().get_deriv_glts().end(); ++pp)
            {
              double value = -penaltyLts * (p->second) * (pp->second) * cnode->mo_data().n()[dim];
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->data().get_deriv_n()[dim].begin();
                 pp != cnode->data().get_deriv_n()[dim].end(); ++pp)
            {
              double value = -penaltyLts * (p->second) * (pp->second) * cnode->data().getglts();
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_mlts()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mlts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            for (CImap pp = cnode->data().get_deriv_glts().begin();
                 pp != cnode->data().get_deriv_glts().end(); ++pp)
            {
              double value = penaltyLts * (p->second) * (pp->second) * cnode->mo_data().n()[dim];
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->data().get_deriv_n()[dim].begin();
                 pp != cnode->data().get_deriv_n()[dim].end(); ++pp)
            {
              double value = penaltyLts * (p->second) * (pp->second) * cnode->data().getglts();
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_nt_sstiffness_master(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for corners
    if (!cnode->is_on_corner()) continue;

    // is gap is in contact
    if (cnode->data().getgnts() < 1e-12)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penalty * cnode->data().getgnts() * cnode->mo_data().n()[0];
      lm[1] = penalty * cnode->data().getgnts() * cnode->mo_data().n()[1];
      lm[2] = penalty * cnode->data().getgnts() * cnode->mo_data().n()[2];

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->data().get_deriv_mnts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        Core::Nodes::Node* mnode = idiscret_->g_node(mgid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->data().get_deriv_mnts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = cmnode->dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->fe_assemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dnts()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dnts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            // gap linearization
            for (CImap pp = cnode->data().get_deriv_gnts().begin();
                 pp != cnode->data().get_deriv_gnts().end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second) * cnode->mo_data().n()[dim];
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->data().get_deriv_n()[dim].begin();
                 pp != cnode->data().get_deriv_n()[dim].end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second) * cnode->data().getgnts();
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_mnts()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mnts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            for (CImap pp = cnode->data().get_deriv_gnts().begin();
                 pp != cnode->data().get_deriv_gnts().end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second) * cnode->mo_data().n()[dim];
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->data().get_deriv_n()[dim].begin();
                 pp != cnode->data().get_deriv_n()[dim].end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second) * cnode->data().getgnts();
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::add_lt_lstiffness(Teuchos::RCP<Core::LinAlg::SparseMatrix> kteff)
{
  const double penalty = interface_params().get<double>("PENALTYPARAM");

  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is valid node
    if (cnode->data().getgltl()[0] < 1e8 and cnode->data().getgltl()[1] < 1e8 and
        cnode->data().getgltl()[2] < 1e8)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penalty * cnode->data().getgltl()[0];
      lm[1] = penalty * cnode->data().getgltl()[1];
      lm[2] = penalty * cnode->data().getgltl()[2];

      std::map<int, std::map<int, double>>& dderiv = cnode->data().get_deriv_dltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        Core::Nodes::Node* snode = idiscret_->g_node(sgid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->data().get_deriv_dltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = csnode->dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->fe_assemble(-val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->data().get_deriv_mltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        Core::Nodes::Node* mnode = idiscret_->g_node(mgid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->data().get_deriv_mltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < n_dim(); ++prodj)
        {
          int row = cmnode->dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->fe_assemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_dltl()).size() > 0)
      {
        Core::Gen::Pairedvector<int, double> map = cnode->mo_data().get_dltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            for (CImap pp = cnode->data().get_deriv_gltl()[dim].begin();
                 pp != cnode->data().get_deriv_gltl()[dim].end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second);
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->mo_data().get_mltl()).size() > 0)
      {
        std::map<int, double> map = cnode->mo_data().get_mltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          Core::Nodes::Node* snode = idiscret_->g_node(gid3);
          if (!snode) FOUR_C_THROW("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < n_dim(); ++dim)
          {
            for (CImap pp = cnode->data().get_deriv_gltl()[dim].begin();
                 pp != cnode->data().get_deriv_gltl()[dim].end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second);
              kteff->fe_assemble(value, csnode->dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        FOUR_C_THROW("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  post evaluate to scale calculated terms                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::post_evaluate(const int step, const int iter)
{
  // decide which type of coupling should be evaluated
  Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  switch (algo)
  {
    //*********************************
    // Mortar Coupling (STS)    (2D/3D)
    //*********************************
    case Inpar::Mortar::algorithm_mortar:
    {
      // non-smooth contact
      if (nonSmoothContact_)
      {
        // store lts into mortar data container
        store_lt_svalues();

        // store nts into mortar data container
        store_nt_svalues();
      }
      return;
      break;
    }
    //*********************************
    // Gauss-Point-To-Segment (GPTS)
    //*********************************
    case Inpar::Mortar::algorithm_gpts:
    {
      // already stored
      return;
      break;
    }
    //*********************************
    // Line-to-Segment Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_lts:
    {
      // store lts into mortar data container
      store_lt_svalues();
      break;
    }
    //*********************************
    // Node-to-Segment Coupling (2D/3D)
    //*********************************
    case Inpar::Mortar::algorithm_nts:
    {
      // store nts into mortar data container
      store_nt_svalues();
      break;
    }
    //*********************************
    // line-to-line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_ltl:
    {
      return;
      break;
    }
    //*********************************
    // Node-to-Line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_ntl:
    {
      FOUR_C_THROW("not yet implemented!");
      break;
    }
    //*********************************
    // Segment-to-Line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_stl:
    {
      // store lts into mortar data container
      store_lt_svalues();
      break;
    }
    //*********************************
    // Default case
    //*********************************
    default:
    {
      FOUR_C_THROW("Unknown discr. type for constraints!");
      break;
    }
  }

#ifdef MORTARGMSHCELLS
  // finish integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "a");
  std::stringstream gmshfilecontent2;
  gmshfilecontent2 << "};" << std::endl;
  fprintf(fp, gmshfilecontent2.str().c_str());
  fclose(fp);

  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream newfilename;
  newfilename << "o/gmsh_output/cells_";
  if (step < 10)
    newfilename << 0 << 0 << 0 << 0;
  else if (step < 100)
    newfilename << 0 << 0 << 0;
  else if (step < 1000)
    newfilename << 0 << 0;
  else if (step < 10000)
    newfilename << 0;
  else if (step > 99999)
    FOUR_C_THROW("Gmsh output implemented for a maximum of 99.999 time steps");
  newfilename << step;

  // second index = Newton iteration index
  newfilename << "_";
  if (iter < 10)
    newfilename << 0;
  else if (iter > 99)
    FOUR_C_THROW("Gmsh output implemented for a maximum of 99 iterations");
  newfilename << iter << "_p" << proc << ".pos";

  // rename file
  rename(filename.str().c_str(), newfilename.str().c_str());
#endif  // #ifdef MORTARGMSHCELLS

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-segment coupl          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_sts(
    const Epetra_Map& selecolmap, const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  Mortar::Interface::evaluate_sts(selecolmap, mparams_ptr);
  return;
}

/*----------------------------------------------------------------------*
 |  protected evaluate routine                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_coupling(const Epetra_Map& selecolmap,
    const Epetra_Map* snoderowmap, const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // ask if non-smooth contact is activated!
  if (nonSmoothContact_)
  {
    // 2D: only STS and nts has to be performed
    if (n_dim() == 2)
    {
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_sts(selecolmap, mparams_ptr);

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_nts();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else if (n_dim() == 3)
    {
      //********************************************************************
      // TODO: remove this hack!
      // HACK: LTL is not yet included in nonsmooth contact framework!
      // However, we want to test the LTL code separately. Thus, the "if"-
      // statement is included:
      // decide which type of coupling should be evaluated
      //********************************************************************
      Inpar::Mortar::AlgorithmType algo =
          Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");
      if (algo == Inpar::Mortar::algorithm_ltl)
      {
        evaluate_ltl();
        return;
      }

      //********************************************************************
      // 1) try to project slave nodes onto master elements
      // 2) evaluate shape functions at projected positions
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      evaluate_nts();

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //********************************************************************
      evaluate_ltl();

      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //********************************************************************
      evaluate_lts();

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      evaluate_sts(selecolmap, mparams_ptr);

      //********************************************************************
      // perform LTS steps for master edges
      //********************************************************************
      //      EvaluateLTSMaster();

      //********************************************************************
      // perform NTS steps for master edges
      //********************************************************************
      //      EvaluateNTSMaster();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else
    {
      FOUR_C_THROW("Wrong dimension!");
    }
  }
  else
  {
    //********************************************************************
    // call base routine for standard mortar/nts evaluation
    //********************************************************************
    Mortar::Interface::evaluate_coupling(selecolmap, snoderowmap, mparams_ptr);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Check and initialize corner/edge contact                 farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::initialize_corner_edge()
{
  // return if nonsmooth contact is activated
  if (nonSmoothContact_) return;

  // call base function
  Mortar::Interface::initialize_corner_edge();

  return;
}


/*----------------------------------------------------------------------*
 |  stuff for non-smooth contact geometries                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::detect_non_smooth_geometries()
{
  FOUR_C_THROW("outdated!");

  std::vector<int> nonsmoothnodegids(0);
  std::vector<int> smoothnodegids(0);

  // loop over slave nodes to find nodes and eles
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    if (cnode->owner() != get_comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    if (cnode->num_element() < 2)
    {
      nonsmoothnodegids.push_back(cnode->id());
      continue;
    }

    int loc = 0;
    std::array<double, 3> normalsele0 = {0.0, 0.0, 0.0};
    std::array<double, 3> normalsele1 = {0.0, 0.0, 0.0};

    // build normal at node for 1. ele
    Element* sele0 = dynamic_cast<Element*>(cnode->elements()[0]);
    Core::LinAlg::SerialDenseMatrix elen0(6, 1);
    sele0->build_normal_at_node(cnode->id(), loc, elen0);

    // create 1. unit normal
    normalsele0[0] = elen0(0, 0) / elen0(4, 0);
    normalsele0[1] = elen0(1, 0) / elen0(4, 0);
    normalsele0[2] = elen0(2, 0) / elen0(4, 0);

    // build normal at node for 2. ele
    loc = 0;
    Element* sele1 = dynamic_cast<Element*>(cnode->elements()[1]);
    Core::LinAlg::SerialDenseMatrix elen1(6, 1);
    sele1->build_normal_at_node(cnode->id(), loc, elen1);

    // create 2. unit normal
    normalsele1[0] = elen1(0, 0) / elen1(4, 0);
    normalsele1[1] = elen1(1, 0) / elen1(4, 0);
    normalsele1[2] = elen1(2, 0) / elen1(4, 0);

    double dot = normalsele0[0] * normalsele1[0] + normalsele0[1] * normalsele1[1] +
                 normalsele0[2] * normalsele1[2];
    double R = abs(dot);  // this is the abs ( cos(phi) ) --> 1.0 if parallel; 0.0 if orthogonal

    // critical node
    if (R < 0.8)  // circa 90 - 11.5 = 78.5
      nonsmoothnodegids.push_back(cnode->id());
    else
      smoothnodegids.push_back(cnode->id());
  }

  // create maps
  nonsmoothnodes_ = Core::LinAlg::CreateMap(nonsmoothnodegids, get_comm());
  smoothnodes_ = Core::LinAlg::CreateMap(smoothnodegids, get_comm());

  // debug output
  std::cout << "# nonsmooth nodes: " << nonsmoothnodes_->NumGlobalElements() << std::endl;
  std::cout << "#    smooth nodes: " << smoothnodes_->NumGlobalElements() << std::endl;

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  cpp to edge + Lin                                       farah 11/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::compute_normal_node_to_edge(Mortar::Node& snode, Mortar::Element& mele,
    double* normal, std::vector<Core::Gen::Pairedvector<int, double>>& normaltonodelin)
{
  // define tolerance
  const double tol = 1e-8;
  double dist = 1e12;
  int nrow = mele.num_node();

  Node* node1 = dynamic_cast<Node*>(mele.nodes()[0]);
  Node* node2 = dynamic_cast<Node*>(mele.nodes()[1]);

  double length1 = sqrt(node1->mo_data().edge_tangent()[0] * node1->mo_data().edge_tangent()[0] +
                        node1->mo_data().edge_tangent()[1] * node1->mo_data().edge_tangent()[1] +
                        node1->mo_data().edge_tangent()[2] * node1->mo_data().edge_tangent()[2]);
  double length2 = sqrt(node2->mo_data().edge_tangent()[0] * node2->mo_data().edge_tangent()[0] +
                        node2->mo_data().edge_tangent()[1] * node2->mo_data().edge_tangent()[1] +
                        node2->mo_data().edge_tangent()[2] * node2->mo_data().edge_tangent()[2]);

  if (length1 < 1e-12 or length2 < 1e-12) return dist;

  // calc angle between tangents
  std::array<double, 3> t1 = {0.0, 0.0, 0.0};
  std::array<double, 3> t2 = {0.0, 0.0, 0.0};

  t1[0] = node1->mo_data().edge_tangent()[0];
  t1[1] = node1->mo_data().edge_tangent()[1];
  t1[2] = node1->mo_data().edge_tangent()[2];

  t2[0] = node2->mo_data().edge_tangent()[0];
  t2[1] = node2->mo_data().edge_tangent()[1];
  t2[2] = node2->mo_data().edge_tangent()[2];

  double test = t1[0] * t2[0] + t1[1] * t2[1] + t1[2] * t2[2];
  if (test < tol) FOUR_C_THROW("tangents have wrong direction!");


  double f = 0.0;
  double df = 0.0;
  double xi = 0.0;

  // newton loop
  for (int k = 0; k < MORTARMAXITER; ++k)
  {
    //    std::cout << "k= " << k << std::endl;
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    Core::LinAlg::SerialDenseVector sval(nrow);
    Core::LinAlg::SerialDenseMatrix sderiv(nrow, 1);
    mele.evaluate_shape(&xi, sval, sderiv, nrow);

    // tangent part
    std::array<double, 3> tangent = {0.0, 0.0, 0.0};
    tangent[0] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[0];
    tangent[1] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[1];
    tangent[2] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[2];

    tangent[0] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[0];
    tangent[1] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[1];
    tangent[2] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[2];

    double tangentSlave = 0.0;
    tangentSlave = tangent[0] * snode.xspatial()[0] + tangent[1] * snode.xspatial()[1] +
                   tangent[2] * snode.xspatial()[2];

    // master part
    std::array<double, 3> master = {0.0, 0.0, 0.0};
    master[0] += sval[0] * node1->xspatial()[0];
    master[1] += sval[0] * node1->xspatial()[1];
    master[2] += sval[0] * node1->xspatial()[2];

    master[0] += sval[1] * node2->xspatial()[0];
    master[1] += sval[1] * node2->xspatial()[1];
    master[2] += sval[1] * node2->xspatial()[2];

    double tangentMaster = 0.0;
    tangentMaster = tangent[0] * master[0] + tangent[1] * master[1] + tangent[2] * master[2];

    f = tangentSlave - tangentMaster;
    if (abs(f) <= MORTARCONVTOL) break;
    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // lin tangent part
    std::array<double, 3> lintangent = {0.0, 0.0, 0.0};
    lintangent[0] += sderiv(0, 0) * node1->mo_data().edge_tangent()[0];
    lintangent[1] += sderiv(0, 0) * node1->mo_data().edge_tangent()[1];
    lintangent[2] += sderiv(0, 0) * node1->mo_data().edge_tangent()[2];

    lintangent[0] += sderiv(1, 0) * node2->mo_data().edge_tangent()[0];
    lintangent[1] += sderiv(1, 0) * node2->mo_data().edge_tangent()[1];
    lintangent[2] += sderiv(1, 0) * node2->mo_data().edge_tangent()[2];

    double lintangentSlave = 0.0;
    lintangentSlave = lintangent[0] * snode.xspatial()[0] + lintangent[1] * snode.xspatial()[1] +
                      lintangent[2] * snode.xspatial()[2];

    // lin master part
    std::array<double, 3> linmaster = {0.0, 0.0, 0.0};
    linmaster[0] += sderiv(0, 0) * node1->xspatial()[0];
    linmaster[1] += sderiv(0, 0) * node1->xspatial()[1];
    linmaster[2] += sderiv(0, 0) * node1->xspatial()[2];

    linmaster[0] += sderiv(1, 0) * node2->xspatial()[0];
    linmaster[1] += sderiv(1, 0) * node2->xspatial()[1];
    linmaster[2] += sderiv(1, 0) * node2->xspatial()[2];

    double lintangentMaster = 0.0;
    lintangentMaster =
        lintangent[0] * master[0] + lintangent[1] * master[1] + lintangent[2] * master[2];

    double tangentlinMaster = 0.0;
    tangentlinMaster =
        tangent[0] * linmaster[0] + tangent[1] * linmaster[1] + tangent[2] * linmaster[2];

    df = lintangentSlave - lintangentMaster - tangentlinMaster;
    if (abs(df) < 1e-12) FOUR_C_THROW("df zero");
    xi += -f / df;
  }

  //**********************************************
  //   CHECK XI                                 //
  //**********************************************
  if (-1.0 - tol > xi or xi > 1.0 + tol) return dist;

  //**********************************************
  //   LINEARIZATION   df                       //
  //**********************************************

  Core::LinAlg::SerialDenseVector sval(nrow);
  Core::LinAlg::SerialDenseMatrix sderiv(nrow, 1);
  mele.evaluate_shape(&xi, sval, sderiv, nrow);

  // tangent part
  std::array<double, 3> tangent = {0.0, 0.0, 0.0};
  tangent[0] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[0];
  tangent[1] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[1];
  tangent[2] += sval[0] * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[2];

  tangent[0] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[0];
  tangent[1] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[1];
  tangent[2] += sval[1] * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[2];

  // master part
  std::array<double, 3> master = {0.0, 0.0, 0.0};
  master[0] += sval[0] * node1->xspatial()[0];
  master[1] += sval[0] * node1->xspatial()[1];
  master[2] += sval[0] * node1->xspatial()[2];

  master[0] += sval[1] * node2->xspatial()[0];
  master[1] += sval[1] * node2->xspatial()[1];
  master[2] += sval[1] * node2->xspatial()[2];

  // lin tangent part
  std::array<double, 3> lintangent = {0.0, 0.0, 0.0};
  lintangent[0] += sderiv(0, 0) * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[0];
  lintangent[1] += sderiv(0, 0) * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[1];
  lintangent[2] += sderiv(0, 0) * dynamic_cast<Node*>(mele.nodes()[0])->mo_data().edge_tangent()[2];

  lintangent[0] += sderiv(1, 0) * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[0];
  lintangent[1] += sderiv(1, 0) * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[1];
  lintangent[2] += sderiv(1, 0) * dynamic_cast<Node*>(mele.nodes()[1])->mo_data().edge_tangent()[2];

  // lin master part
  std::array<double, 3> linmaster = {0.0, 0.0, 0.0};
  linmaster[0] += sderiv(0, 0) * node1->xspatial()[0];
  linmaster[1] += sderiv(0, 0) * node1->xspatial()[1];
  linmaster[2] += sderiv(0, 0) * node1->xspatial()[2];

  linmaster[0] += sderiv(1, 0) * node2->xspatial()[0];
  linmaster[1] += sderiv(1, 0) * node2->xspatial()[1];
  linmaster[2] += sderiv(1, 0) * node2->xspatial()[2];


  //**********************************************
  //   LINEARIZATION    f                       //
  //**********************************************
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _CI;

  std::vector<Core::Gen::Pairedvector<int, double>> linT(3, 100);  // added all sizes

  for (_CI p = node1->data().get_deriv_tangent()[0].begin();
       p != node1->data().get_deriv_tangent()[0].end(); ++p)
    linT[0][p->first] += sval[0] * p->second;
  for (_CI p = node1->data().get_deriv_tangent()[1].begin();
       p != node1->data().get_deriv_tangent()[1].end(); ++p)
    linT[1][p->first] += sval[0] * p->second;
  for (_CI p = node1->data().get_deriv_tangent()[2].begin();
       p != node1->data().get_deriv_tangent()[2].end(); ++p)
    linT[2][p->first] += sval[0] * p->second;

  for (_CI p = node2->data().get_deriv_tangent()[0].begin();
       p != node2->data().get_deriv_tangent()[0].end(); ++p)
    linT[0][p->first] += sval[1] * p->second;
  for (_CI p = node2->data().get_deriv_tangent()[1].begin();
       p != node2->data().get_deriv_tangent()[1].end(); ++p)
    linT[1][p->first] += sval[1] * p->second;
  for (_CI p = node2->data().get_deriv_tangent()[2].begin();
       p != node2->data().get_deriv_tangent()[2].end(); ++p)
    linT[2][p->first] += sval[1] * p->second;

  std::vector<Core::Gen::Pairedvector<int, double>> linXsl(3, 100);  // added all sizes
  linXsl[0][snode.dofs()[0]] += 1.0;
  linXsl[1][snode.dofs()[1]] += 1.0;
  linXsl[2][snode.dofs()[2]] += 1.0;

  std::vector<Core::Gen::Pairedvector<int, double>> linXm(3, 100);  // added all sizes
  linXm[0][node1->dofs()[0]] += sval[0];
  linXm[1][node1->dofs()[1]] += sval[0];
  linXm[2][node1->dofs()[2]] += sval[0];

  linXm[0][node2->dofs()[0]] += sval[1];
  linXm[1][node2->dofs()[1]] += sval[1];
  linXm[2][node2->dofs()[2]] += sval[1];

  Core::Gen::Pairedvector<int, double> linf(100);  // added all sizes
  for (_CI p = linT[0].begin(); p != linT[0].end(); ++p)
    linf[p->first] += snode.xspatial()[0] * p->second;
  for (_CI p = linT[1].begin(); p != linT[1].end(); ++p)
    linf[p->first] += snode.xspatial()[1] * p->second;
  for (_CI p = linT[2].begin(); p != linT[2].end(); ++p)
    linf[p->first] += snode.xspatial()[2] * p->second;

  for (_CI p = linXsl[0].begin(); p != linXsl[0].end(); ++p)
    linf[p->first] += tangent[0] * p->second;
  for (_CI p = linXsl[1].begin(); p != linXsl[1].end(); ++p)
    linf[p->first] += tangent[1] * p->second;
  for (_CI p = linXsl[2].begin(); p != linXsl[2].end(); ++p)
    linf[p->first] += tangent[2] * p->second;

  for (_CI p = linT[0].begin(); p != linT[0].end(); ++p) linf[p->first] -= master[0] * p->second;
  for (_CI p = linT[1].begin(); p != linT[1].end(); ++p) linf[p->first] -= master[1] * p->second;
  for (_CI p = linT[2].begin(); p != linT[2].end(); ++p) linf[p->first] -= master[2] * p->second;

  for (_CI p = linXm[0].begin(); p != linXm[0].end(); ++p) linf[p->first] -= tangent[0] * p->second;
  for (_CI p = linXm[1].begin(); p != linXm[1].end(); ++p) linf[p->first] -= tangent[1] * p->second;
  for (_CI p = linXm[2].begin(); p != linXm[2].end(); ++p) linf[p->first] -= tangent[2] * p->second;

  Core::Gen::Pairedvector<int, double> linXi(100);  // added all sizes
  for (_CI p = linf.begin(); p != linf.end(); ++p) linXi[p->first] -= p->second / df;

  //**********************************************
  //   CALC NORMAL                              //
  //**********************************************
  std::array<double, 3> auxnormal = {0.0, 0.0, 0.0};
  auxnormal[0] = snode.xspatial()[0] - master[0];
  auxnormal[1] = snode.xspatial()[1] - master[1];
  auxnormal[2] = snode.xspatial()[2] - master[2];

  // calc distance
  dist =
      sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);

  if (abs(dist) < 1e-12) return 1e12;

  //*******************************************
  // Lin:
  std::vector<Core::Gen::Pairedvector<int, double>> auxlin(3, 100);  // added all sizes

  // xslave
  for (int k = 0; k < 3; ++k) (auxlin[k])[snode.dofs()[k]] += 1.0;

  // xmaster n1
  for (int k = 0; k < 3; ++k) (auxlin[k])[node1->dofs()[k]] -= sval[0];
  // xmaster n2
  for (int k = 0; k < 3; ++k) (auxlin[k])[node2->dofs()[k]] -= sval[1];

  for (_CI p = linXi.begin(); p != linXi.end(); ++p)
  {
    for (int k = 0; k < 3; ++k)
    {
      (auxlin[k])[p->first] -= sderiv(0, 0) * node1->xspatial()[k] * p->second;
      (auxlin[k])[p->first] -= sderiv(1, 0) * node2->xspatial()[k] * p->second;
    }
  }

  normal[0] = auxnormal[0];
  normal[1] = auxnormal[1];
  normal[2] = auxnormal[2];

  //******************************
  // Orientation check:
  std::array<double, 3> slavebasednormal = {0.0, 0.0, 0.0};
  int nseg = snode.num_element();
  Core::Elements::Element** adjeles = snode.elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  Core::LinAlg::SerialDenseMatrix elens(6, nseg);
  Mortar::Element* adjmrtrele = dynamic_cast<Mortar::Element*>(adjeles[0]);

  // build element normal at current node
  // (we have to pass in the index i to be able to store the
  // normal and other information at the right place in elens)
  int i = 0;
  adjmrtrele->build_normal_at_node(snode.id(), i, elens);

  // add (weighted) element normal to nodal normal n
  for (int j = 0; j < 3; ++j) slavebasednormal[j] += elens(j, 0) / elens(4, 0);

  // create unit normal vector
  const double length =
      sqrt(slavebasednormal[0] * slavebasednormal[0] + slavebasednormal[1] * slavebasednormal[1] +
           slavebasednormal[2] * slavebasednormal[2]);
  if (abs(length) < 1e-12)
  {
    FOUR_C_THROW("Nodal normal length 0, node ID %i", snode.id());
  }
  else
  {
    for (int j = 0; j < 3; ++j) slavebasednormal[j] /= length;
  }

  const double dotprod =
      -(normal[0] / dist * slavebasednormal[0] + normal[1] / dist * slavebasednormal[1] +
          normal[2] / dist * slavebasednormal[2]);

  if (dotprod < -1e-12)
  {
    // get the cpp normal
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];

    for (int j = 0; j < 3; ++j)
      for (_CI p = auxlin[j].begin(); p != auxlin[j].end(); ++p)
        (normaltonodelin[j])[p->first] -= (p->second);
  }
  else
  {
    // linearization
    for (int j = 0; j < 3; ++j)
      for (_CI p = auxlin[j].begin(); p != auxlin[j].end(); ++p)
        (normaltonodelin[j])[p->first] += (p->second);
  }

  return dist;
}

/*----------------------------------------------------------------------*
 |  cpp to node + Lin                                       farah 01/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::compute_normal_node_to_node(Mortar::Node& snode, Mortar::Node& mnode,
    double* normal, std::vector<Core::Gen::Pairedvector<int, double>>& normaltonodelin)
{
  const int dim = n_dim();

  // distance between node and surface
  double gdist = 1e12;
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<Core::Gen::Pairedvector<int, double>> glin(3, 1000);
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  double dist = 1e12;
  std::array<double, 3> auxnormal = {0.0, 0.0, 0.0};

  // loop over found master nodes
  std::vector<Core::Gen::Pairedvector<int, double>> auxlin(3, 1000);

  // calc vector
  auxnormal[0] = snode.xspatial()[0] - mnode.xspatial()[0];
  auxnormal[1] = snode.xspatial()[1] - mnode.xspatial()[1];
  auxnormal[2] = snode.xspatial()[2] - mnode.xspatial()[2];

  // remove numerical artifacts
  if (abs(auxnormal[0]) < 1e-12) auxnormal[0] = 0.0;
  if (abs(auxnormal[1]) < 1e-12) auxnormal[1] = 0.0;
  if (abs(auxnormal[2]) < 1e-12) auxnormal[2] = 0.0;

  // calc distance
  dist =
      sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);

  // if nodes lying on each other: continue to next master node
  if (abs(dist) < 1e-12) return dist;

  //*******************************************
  // Lin:
  // xslave
  for (int k = 0; k < dim; ++k) (auxlin[k])[snode.dofs()[k]] += 1.0;

  // xmaster
  for (int k = 0; k < dim; ++k) (auxlin[k])[mnode.dofs()[k]] -= 1.0;

  // get normal
  gdist = dist;

  // normalize vector
  gnormal[0] = auxnormal[0];  /// dist;
  gnormal[1] = auxnormal[1];  /// dist;
  gnormal[2] = auxnormal[2];  /// dist;

  // linearization
  glin = auxlin;

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];

  //******************************
  // Orientation check:
  std::array<double, 3> slavebasednormal = {0.0, 0.0, 0.0};
  int nseg = snode.num_element();
  Core::Elements::Element** adjeles = snode.elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  Core::LinAlg::SerialDenseMatrix elens(6, nseg);
  Mortar::Element* adjmrtrele = dynamic_cast<Mortar::Element*>(adjeles[0]);

  // build element normal at current node
  // (we have to pass in the index i to be able to store the
  // normal and other information at the right place in elens)
  int i = 0;
  adjmrtrele->build_normal_at_node(snode.id(), i, elens);

  // add (weighted) element normal to nodal normal n
  for (int j = 0; j < 3; ++j) slavebasednormal[j] += elens(j, 0) / elens(4, 0);

  // create unit normal vector
  const double length =
      sqrt(slavebasednormal[0] * slavebasednormal[0] + slavebasednormal[1] * slavebasednormal[1] +
           slavebasednormal[2] * slavebasednormal[2]);
  if (abs(length) < 1e-12)
  {
    FOUR_C_THROW("Nodal normal length 0, node ID %i", snode.id());
  }
  else
  {
    for (int j = 0; j < 3; ++j) slavebasednormal[j] /= length;
  }

  const double dotprod = -(normal[0] * slavebasednormal[0] + normal[1] * slavebasednormal[1] +
                           normal[2] * slavebasednormal[2]);

  if (dotprod < -1e-12)
  {
    // get the cpp normal
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];

    for (int j = 0; j < dim; ++j)
      for (CI p = glin[j].begin(); p != glin[j].end(); ++p)
        (normaltonodelin[j])[p->first] -= (p->second);
  }
  else
  {
    // linearization
    for (int j = 0; j < dim; ++j)
      for (CI p = glin[j].begin(); p != glin[j].end(); ++p)
        (normaltonodelin[j])[p->first] += (p->second);
  }


  return gdist;
}

/*----------------------------------------------------------------------*
 |  evaluate closest point normals                          farah 08/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_cpp_normals()
{
  // Build averaged normal field on physically smooth surface
  // loop over proc's master nodes of the interface
  // use row map and export to column map later
  for (int i = 0; i < master_row_nodes()->NumMyElements(); ++i)
  {
    int gid = master_row_nodes()->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each master node
    mrtrnode->build_averaged_normal();

    // build tangent
    if (mrtrnode->is_on_edge()) mrtrnode->build_averaged_edge_tangent();
  }

  // export nodal normals
  export_master_nodal_normals();

  // loop over slave nodes
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Mortar::Node* mrtrnode = dynamic_cast<Mortar::Node*>(node);

    if (mrtrnode->owner() != get_comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    // vector with possible contacting master eles/nodes
    std::vector<Mortar::Element*> meles;
    std::vector<Mortar::Node*> mnodes;

    // fill vector with possibly contacting meles
    find_m_eles(*mrtrnode, meles);

    // fallback solution if no mele is available
    if (meles.size() < 1)  // or !mrtrnode->IsOnCornerEdge())
    {
      Node* cnode = dynamic_cast<Node*>(mrtrnode);
      cnode->build_averaged_normal();
      continue;
    }


    // Here we have all found master elements for one slave node.
    // distance for cpp
    double normaltoline[3] = {0.0, 0.0, 0.0};
    std::vector<Core::Gen::Pairedvector<int, double>> normaltolineLin(3, 1);  // 1 dummy

    // Now, calculate distance between node and master line
    double dist = compute_cpp_normal(*mrtrnode, meles, normaltoline, normaltolineLin);

    // if no projection was posible
    if (dist > 1e11)
    {
      Node* cnode = dynamic_cast<Node*>(mrtrnode);
      cnode->build_averaged_normal();
      continue;
    }

    // set the normal and its lineratization
    set_cpp_normal(*mrtrnode, normaltoline, normaltolineLin);
  }

  // export slave normals
  export_nodal_normals();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  export master nodal normals (protected)                  farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::export_master_nodal_normals()
{
  std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>> triad;

  std::map<int, std::vector<int>> n_x_key;
  std::map<int, std::vector<int>> n_y_key;
  std::map<int, std::vector<int>> n_z_key;
  std::map<int, std::vector<int>> txi_x_key;
  std::map<int, std::vector<int>> txi_y_key;
  std::map<int, std::vector<int>> txi_z_key;
  std::map<int, std::vector<int>> teta_x_key;
  std::map<int, std::vector<int>> teta_y_key;
  std::map<int, std::vector<int>> teta_z_key;

  std::map<int, std::vector<double>> n_x_val;
  std::map<int, std::vector<double>> n_y_val;
  std::map<int, std::vector<double>> n_z_val;
  std::map<int, std::vector<double>> txi_x_val;
  std::map<int, std::vector<double>> txi_y_val;
  std::map<int, std::vector<double>> txi_z_val;
  std::map<int, std::vector<double>> teta_x_val;
  std::map<int, std::vector<double>> teta_y_val;
  std::map<int, std::vector<double>> teta_z_val;

  Core::Gen::Pairedvector<int, double>::iterator iter;

  const Teuchos::RCP<Epetra_Map> masternodes = Core::LinAlg::AllreduceEMap(*(mnoderowmap_));

  // build info on row map
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);

    // fill nodal matrix
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc =
        Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, 3));
    (*loc)(0, 0) = cnode->mo_data().n()[0];
    (*loc)(1, 0) = cnode->mo_data().n()[1];
    (*loc)(2, 0) = cnode->mo_data().n()[2];
    (*loc)(0, 1) = cnode->data().txi()[0];
    (*loc)(1, 1) = cnode->data().txi()[1];
    (*loc)(2, 1) = cnode->data().txi()[2];
    (*loc)(0, 2) = cnode->data().teta()[0];
    (*loc)(1, 2) = cnode->data().teta()[1];
    (*loc)(2, 2) = cnode->data().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<Core::Gen::Pairedvector<int, double>>& derivn = cnode->data().get_deriv_n();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivtxi = cnode->data().get_deriv_txi();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivteta = cnode->data().get_deriv_teta();

    for (iter = derivn[0].begin(); iter != derivn[0].end(); ++iter)
    {
      n_x_key[gid].push_back(iter->first);
      n_x_val[gid].push_back(iter->second);
    }
    for (iter = derivn[1].begin(); iter != derivn[1].end(); ++iter)
    {
      n_y_key[gid].push_back(iter->first);
      n_y_val[gid].push_back(iter->second);
    }
    for (iter = derivn[2].begin(); iter != derivn[2].end(); ++iter)
    {
      n_z_key[gid].push_back(iter->first);
      n_z_val[gid].push_back(iter->second);
    }

    for (iter = derivtxi[0].begin(); iter != derivtxi[0].end(); ++iter)
    {
      txi_x_key[gid].push_back(iter->first);
      txi_x_val[gid].push_back(iter->second);
    }
    for (iter = derivtxi[1].begin(); iter != derivtxi[1].end(); ++iter)
    {
      txi_y_key[gid].push_back(iter->first);
      txi_y_val[gid].push_back(iter->second);
    }
    for (iter = derivtxi[2].begin(); iter != derivtxi[2].end(); ++iter)
    {
      txi_z_key[gid].push_back(iter->first);
      txi_z_val[gid].push_back(iter->second);
    }

    for (iter = derivteta[0].begin(); iter != derivteta[0].end(); ++iter)
    {
      teta_x_key[gid].push_back(iter->first);
      teta_x_val[gid].push_back(iter->second);
    }
    for (iter = derivteta[1].begin(); iter != derivteta[1].end(); ++iter)
    {
      teta_y_key[gid].push_back(iter->first);
      teta_y_val[gid].push_back(iter->second);
    }
    for (iter = derivteta[2].begin(); iter != derivteta[2].end(); ++iter)
    {
      teta_z_key[gid].push_back(iter->first);
      teta_z_val[gid].push_back(iter->second);
    }
  }

  // communicate from master node row to column map
  Core::Communication::Exporter ex(*mnoderowmap_, *masternodes, get_comm());
  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for (int i = 0; i < masternodes->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = masternodes->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    int linsize = cnode->get_linsize() + (int)(n_x_key[gid].size());

    if (cnode->owner() == get_comm().MyPID()) continue;

    // extract info
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc = triad[gid];
    cnode->mo_data().n()[0] = (*loc)(0, 0);
    cnode->mo_data().n()[1] = (*loc)(1, 0);
    cnode->mo_data().n()[2] = (*loc)(2, 0);
    cnode->data().txi()[0] = (*loc)(0, 1);
    cnode->data().txi()[1] = (*loc)(1, 1);
    cnode->data().txi()[2] = (*loc)(2, 1);
    cnode->data().teta()[0] = (*loc)(0, 2);
    cnode->data().teta()[1] = (*loc)(1, 2);
    cnode->data().teta()[2] = (*loc)(2, 2);

    // extract derivative info
    std::vector<Core::Gen::Pairedvector<int, double>>& derivn = cnode->data().get_deriv_n();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivtxi = cnode->data().get_deriv_txi();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivteta = cnode->data().get_deriv_teta();

    for (int k = 0; k < (int)(derivn.size()); ++k) derivn[k].clear();
    derivn.resize(3, linsize);
    for (int k = 0; k < (int)(derivtxi.size()); ++k) derivtxi[k].clear();
    derivtxi.resize(3, linsize);
    for (int k = 0; k < (int)(derivteta.size()); ++k) derivteta[k].clear();
    derivteta.resize(3, linsize);

    cnode->data().get_deriv_n()[0].resize(linsize);
    cnode->data().get_deriv_n()[1].resize(linsize);
    cnode->data().get_deriv_n()[2].resize(linsize);

    cnode->data().get_deriv_txi()[0].resize(linsize);
    cnode->data().get_deriv_txi()[1].resize(linsize);
    cnode->data().get_deriv_txi()[2].resize(linsize);

    cnode->data().get_deriv_teta()[0].resize(linsize);
    cnode->data().get_deriv_teta()[1].resize(linsize);
    cnode->data().get_deriv_teta()[2].resize(linsize);

    for (int k = 0; k < (int)(n_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k = 0; k < (int)(n_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k = 0; k < (int)(n_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k = 0; k < (int)(txi_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k = 0; k < (int)(txi_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k = 0; k < (int)(txi_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k = 0; k < (int)(teta_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k = 0; k < (int)(teta_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k = 0; k < (int)(teta_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nodal normals (public)                          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_averaged_nodal_normals()
{
  FOUR_C_THROW("outdated function!");
  // safety
  if (smoothnodes_ == Teuchos::null) FOUR_C_THROW("map of non smooth nodes is wrong!");

  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < smoothnodes_->NumMyElements(); ++i)
  {
    int gid = smoothnodes_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each slave node
    mrtrnode->build_averaged_normal();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 09/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::compute_scaling()
{
  // ltl scaling only for 3D setting
  if (n_dim() == 2) return;

  // compute ltl scaling for point and line contact
  compute_scaling_ltl();

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 01/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::compute_scaling_ltl()
{
  std::cout << "compute_scaling_ltl" << std::endl;

  // define iterator for linerization
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  // angle in degree
  const double minAngle = interface_params().get<double>("HYBRID_ANGLE_MIN");
  const double maxAngle = interface_params().get<double>("HYBRID_ANGLE_MAX");

  // check
  if (minAngle < 0.0 or maxAngle < 0.0) FOUR_C_THROW("invalid hybrid angle!");
  if (minAngle >= maxAngle) FOUR_C_THROW("invalid hybrid angle!");

  // angle in rad
  const double alphaMin = minAngle * (M_PI / 180.0);  // min angle in rad
  const double alphaMax = maxAngle * (M_PI / 180.0);  // max angle in rad

  //======================================================
  // search slave edges and master edges and store them

  // guarantee uniquness of slave edges
  std::set<std::pair<int, int>> donebeforeS;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElementsS;

    if (selement->shape() == Core::FE::CellType::quad4)
    {
      for (int j = 0; j < 4; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->node_ids()[0];
          nodeIds[1] = selement->node_ids()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->node_ids()[1];
          nodeIds[1] = selement->node_ids()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->node_ids()[2];
          nodeIds[1] = selement->node_ids()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->node_ids()[3];
          nodeIds[1] = selement->node_ids()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[0]])->is_on_edge();
        bool node1Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[1]])->is_on_edge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, selement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<Core::Nodes::Node*, 2> nodes = {
              selement->nodes()[nodeLIds[0]], selement->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes.data());

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else
      FOUR_C_THROW("LTL only for quad4!");

    // guarantee uniquness of master edges
    std::set<std::pair<int, int>> donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->mo_data().num_search_elements(); ++k)
    {
      int gid2 = selement->mo_data().search_elements()[k];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->shape() == Core::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->node_ids()[0];
            nodeIds[1] = melement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->node_ids()[1];
            nodeIds[1] = melement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->node_ids()[2];
            nodeIds[1] = melement->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->node_ids()[3];
            nodeIds[1] = melement->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[0]])->is_on_edge();
          bool node1Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[1]])->is_on_edge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
                j, melement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<Core::Nodes::Node*, 2> nodes = {
                melement->nodes()[nodeLIds[0]], melement->nodes()[nodeLIds[1]]};
            lineEle->build_nodal_pointers(nodes.data());

            // init data container for dual shapes
            lineEle->initialize_data_container();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else
        FOUR_C_THROW("LTL only for quad4!");
    }  // end found mele loop

    // loop over slave edges
    for (int s = 0; s < (int)lineElementsS.size(); ++s)
    {
      double gR = 1e12;
      bool parallel = false;

      // loop over master edges
      for (int m = 0; m < (int)lineElementsM.size(); ++m)
      {
        // 1. search slave master pair with cpp point
        LineToLineCouplingPoint3d coup(
            *idiscret_, 3, interface_params(), lineElementsS[s], lineElementsM[m]);

        parallel = coup.check_parallelity();
        if (parallel)
        {
          continue;
          //          FOUR_C_THROW("edges parallel");
          break;
        }
        // create empty points
        double sxi = 0.0;
        double mxi = 0.0;
        Core::Gen::Pairedvector<int, double> dsxi(
            3 * lineElementsM[m]->num_node() + 3 * lineElementsS[s]->num_node());
        Core::Gen::Pairedvector<int, double> dmxi(
            3 * lineElementsM[m]->num_node() + 3 * lineElementsS[s]->num_node());

        coup.line_intersection(&sxi, &mxi, dsxi, dmxi);
        bool valid = coup.check_intersection(&sxi, &mxi);

        // no valid intersection
        if (!valid)
        {
          continue;
        }
        // found valid intersection
        else
        {
          std::cout << "VALID INTERSECTION" << std::endl;
          // 2. calc current angle for this element pair
          Core::Gen::Pairedvector<int, double> linAngle(
              100 + 3 * lineElementsM[m]->num_node() + 3 * lineElementsS[s]->num_node());
          gR = coup.calc_current_angle(linAngle);

          // get vector of nodes
          std::vector<Node*> cnodes;
          for (int nn = 0; nn < lineElementsS[s]->num_node(); ++nn)
            cnodes.push_back(dynamic_cast<Node*>(lineElementsS[s]->nodes()[nn]));

          //======================================================
          // calc alpha:
          double alpha = 0.0;

          std::cout << "gR= " << gR << std::endl;

          // ltl line contact
          if (gR <= alphaMin)
          {
            std::cout << "LTS contact" << std::endl;
            alpha = 1.0;
          }
          // ltl hybrid
          else if (gR > alphaMin and gR < alphaMax)
          {
            std::cout << "current transition angle in degree = " << gR * (180.0 / M_PI)
                      << std::endl;
            alpha = 0.5 * (1.0 - cos(M_PI * (gR - alphaMax) / (alphaMin - alphaMax)));
            std::cout << "alpha = " << alpha << std::endl;
          }
          // ltl point contact
          else
          {
            std::cout << "LTL contact" << std::endl;
            alpha = 0.0;
          }

          // store to data container
          for (int nn = 0; nn < lineElementsS[s]->num_node(); ++nn)
            cnodes[nn]->data().get_alpha_n() = alpha;

          //======================================================
          // lin alpha:
          if (gR <= alphaMin)
          {
            // nothing: pure ltl line contact
          }
          // hybrid linearization:
          else if (gR > alphaMin and gR < alphaMax)
          {
            // clear old data
            for (int nn = 0; nn < lineElementsS[s]->num_node(); ++nn)
              if (cnodes[nn]->data().get_alpha().size() == 0)
                cnodes[nn]->data().get_alpha().resize(1000);

            const double fac = (M_PI / (2.0 * (alphaMin - alphaMax))) *
                               sin(M_PI * (gR - alphaMax) / (alphaMin - alphaMax));
            for (CI p = linAngle.begin(); p != linAngle.end(); ++p)
            {
              for (int nn = 0; nn < lineElementsS[s]->num_node(); ++nn)
                (cnodes[nn]->data().get_alpha())[p->first] = fac * (p->second);
            }
          }
          else
          {
            // nothing: pure ltl point contact
          }
          break;
        }  // valid projection
      }
    }
  }  // end slave loop
  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::scale_normals()
{
  FOUR_C_THROW("outdated!");

  if (n_dim() == 2)
    scale_normals2_d();
  else if (n_dim() == 3)
    scale_normals3_d();
  else
    FOUR_C_THROW("Wrong dimension!");
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::scale_normals2_d()
{
  // define iterator for paired vector
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  //======================================================
  // loop over non smooth slave nodes
  for (int i = 0; i < nonsmoothnodes_->NumMyElements(); ++i)
  {
    int gid = nonsmoothnodes_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    if (cnode->owner() != get_comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    // get cpp normal
    std::array<double, 3> normalcpp = {0.0, 0.0, 0.0};

    // get cpp normal lin.
    std::vector<Core::Gen::Pairedvector<int, double>> dcppaux(3, 1000);

    dcppaux = cnode->data().get_deriv_n();
    normalcpp[0] = cnode->mo_data().n()[0];
    normalcpp[1] = cnode->mo_data().n()[1];
    normalcpp[2] = cnode->mo_data().n()[2];

    if (cnode->active())
    {
      std::cout << "cnode->Data().GetAlphaN()] = " << cnode->data().get_alpha_n() << std::endl;
      std::cout << "normalcpp[0] = " << normalcpp[0] << std::endl;
      std::cout << "normalcpp[1] = " << normalcpp[1] << std::endl;
      std::cout << "normalcpp[2] = " << normalcpp[2] << std::endl;
    }


    // calc element normal
    const int eles = cnode->num_element();
    if (eles > 2) FOUR_C_THROW("element number to high!");

    double gR = 1e12;
    int localid = -1;

    for (int n = 0; n < eles; ++n)
    {
      // create element normal:
      int loc = 0;
      Element* seleaux = dynamic_cast<Element*>(cnode->elements()[n]);
      Core::LinAlg::SerialDenseMatrix elens(6, 1);
      seleaux->build_normal_at_node(cnode->id(), loc, elens);

      std::array<double, 3> normalsele = {0.0, 0.0, 0.0};
      normalsele[0] = elens(0, 0) / elens(4, 0);
      normalsele[1] = elens(1, 0) / elens(4, 0);
      normalsele[2] = elens(2, 0) / elens(4, 0);

      // calculate angle between cpp and elenormal
      double dotcpp = normalsele[0] * normalcpp[0] + normalsele[1] * normalcpp[1] +
                      normalsele[2] * normalcpp[2];
      double Rl = acos(dotcpp);  // current angle in rad

      if (cnode->active())
      {
        std::cout << "normalsele[0] = " << normalsele[0] << std::endl;
        std::cout << "normalsele[1] = " << normalsele[1] << std::endl;
        std::cout << "normalsele[2] = " << normalsele[2] << std::endl;

        std::cout << "dotcpp= " << dotcpp << std::endl;
        std::cout << "Rl= " << Rl << std::endl;
      }

      if (Rl < 0.0) FOUR_C_THROW("angle less than 0.0");

      // safe smaller angle and linearization
      if (Rl < gR)
      {
        gR = Rl;
        localid = n;
      }
    }  // end ele loop

    // create new normals:
    const double alpha = cnode->data().get_alpha_n();

    // create element normal:
    int loc = 0;
    Element* seleaux = dynamic_cast<Element*>(cnode->elements()[localid]);
    Core::LinAlg::SerialDenseMatrix elens(6, 1);
    seleaux->build_normal_at_node(cnode->id(), loc, elens);

    std::array<double, 3> normalsele = {0.0, 0.0, 0.0};
    normalsele[0] = elens(0, 0) / elens(4, 0);
    normalsele[1] = elens(1, 0) / elens(4, 0);
    normalsele[2] = elens(2, 0) / elens(4, 0);

    cnode->mo_data().n()[0] = alpha * normalsele[0] + (1.0 - alpha) * normalcpp[0];
    cnode->mo_data().n()[1] = alpha * normalsele[1] + (1.0 - alpha) * normalcpp[1];
    cnode->mo_data().n()[2] = alpha * normalsele[2] + (1.0 - alpha) * normalcpp[2];


    // deriv unit normal!
    std::vector<Core::Gen::Pairedvector<int, double>> dnele(3, 1000);
    seleaux->deriv_normal_at_node(cnode->id(), loc, elens, dnele);

    cnode->data().get_deriv_n()[0].clear();
    cnode->data().get_deriv_n()[1].clear();
    cnode->data().get_deriv_n()[2].clear();

    // lin ele normal * alpha
    for (CI p = dnele[0].begin(); p != dnele[0].end(); ++p)
      (cnode->data().get_deriv_n()[0])[p->first] += alpha * (p->second);
    for (CI p = dnele[1].begin(); p != dnele[1].end(); ++p)
      (cnode->data().get_deriv_n()[1])[p->first] += alpha * (p->second);
    for (CI p = dnele[2].begin(); p != dnele[2].end(); ++p)
      (cnode->data().get_deriv_n()[2])[p->first] += alpha * (p->second);

    // lin cpp normal * alpha
    for (CI p = dcppaux[0].begin(); p != dcppaux[0].end(); ++p)
      (cnode->data().get_deriv_n()[0])[p->first] += (1.0 - alpha) * (p->second);
    for (CI p = dcppaux[1].begin(); p != dcppaux[1].end(); ++p)
      (cnode->data().get_deriv_n()[1])[p->first] += (1.0 - alpha) * (p->second);
    for (CI p = dcppaux[2].begin(); p != dcppaux[2].end(); ++p)
      (cnode->data().get_deriv_n()[2])[p->first] += (1.0 - alpha) * (p->second);

    // lin alpha
    for (CI p = cnode->data().get_alpha().begin(); p != cnode->data().get_alpha().end(); ++p)
    {
      (cnode->data().get_deriv_n()[0])[p->first] += (p->second) * (normalsele[0] - normalcpp[0]);
      (cnode->data().get_deriv_n()[1])[p->first] += (p->second) * (normalsele[1] - normalcpp[1]);
      (cnode->data().get_deriv_n()[2])[p->first] += (p->second) * (normalsele[2] - normalcpp[2]);
    }

    // 2D Tangent!
    if (cnode->num_dof() == 2)
    {
      // simple definition for txi
      cnode->data().txi()[0] = -cnode->mo_data().n()[1];
      cnode->data().txi()[1] = cnode->mo_data().n()[0];
      cnode->data().txi()[2] = 0.0;

      // teta is z-axis
      cnode->data().teta()[0] = 0.0;
      cnode->data().teta()[1] = 0.0;
      cnode->data().teta()[2] = 1.0;

      cnode->data().get_deriv_txi()[0].clear();
      cnode->data().get_deriv_txi()[1].clear();
      cnode->data().get_deriv_txi()[2].clear();


      for (CI p = cnode->data().get_deriv_n()[1].begin(); p != cnode->data().get_deriv_n()[1].end();
           ++p)
        (cnode->data().get_deriv_txi()[0])[p->first] -= (p->second);
      for (CI p = cnode->data().get_deriv_n()[0].begin(); p != cnode->data().get_deriv_n()[0].end();
           ++p)
        (cnode->data().get_deriv_txi()[1])[p->first] += (p->second);
    }
    else
      FOUR_C_THROW("only 2D");

    double length = sqrt(cnode->mo_data().n()[0] * cnode->mo_data().n()[0] +
                         cnode->mo_data().n()[1] * cnode->mo_data().n()[1] +
                         cnode->mo_data().n()[2] * cnode->mo_data().n()[2]);

    if (cnode->active())
    {
      std::cout << "cnode->MoData().n()[0]= " << cnode->mo_data().n()[0] << std::endl;
      std::cout << "cnode->MoData().n()[1]= " << cnode->mo_data().n()[1] << std::endl;
      std::cout << "cnode->MoData().n()[2]= " << cnode->mo_data().n()[2] << std::endl;
      std::cout << "length= " << length << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::scale_normals3_d() { FOUR_C_THROW("not yet implemented!"); }


/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::compute_cpp_normal2_d(Mortar::Node& mrtrnode,
    std::vector<Mortar::Element*> meles, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin)
{
  // define tolerance
  const double tol = 1e-8;
  const double validAngle = 20;

  // distance between node and surface
  double gdist = 1e12;  // distance
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<Core::Gen::Pairedvector<int, double>> glin(3, 1);  // 1 dummy
  std::set<int> donebeforeMasterCorner;

  bool nodeOnNode = false;     // flag for node on node (corner on corner) setting
  bool pathdependent = false;  // flag if we have to check path from last converged check
  std::array<double, 3> vect = {0.0, 0.0, 0.0};  // patch

  // calc trajectory and node-to-node distance for corner nodes
  if (mrtrnode.is_on_corner())
  {
    pathdependent = true;
    CONTACT::Node& coNode = dynamic_cast<CONTACT::Node&>(mrtrnode);
    if (coNode.active()) pathdependent = false;

    // calculate path
    // check trajectory from considered node
    std::array<double, 3> Posn = {0.0, 0.0, 0.0};
    std::array<double, 3> Posnp = {0.0, 0.0, 0.0};
    Posn[0] = mrtrnode.x()[0] + mrtrnode.uold()[0];
    Posn[1] = mrtrnode.x()[1] + mrtrnode.uold()[1];
    Posn[2] = mrtrnode.x()[2] + mrtrnode.uold()[2];
    Posnp[0] = mrtrnode.xspatial()[0];
    Posnp[1] = mrtrnode.xspatial()[1];
    Posnp[2] = mrtrnode.xspatial()[2];

    vect[0] = Posnp[0] - Posn[0];
    vect[1] = Posnp[1] - Posn[1];
    vect[2] = Posnp[2] - Posn[2];

    double lvec = sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
    if (lvec < 1e-12)
    {
      pathdependent = false;
    }
    else
    {
      vect[0] /= lvec;
      vect[1] /= lvec;
      vect[2] /= lvec;
    }

    // Compute normal to node:
    // loop over found eles
    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0; i < meles[ele]->num_node(); ++i)
      {
        Core::Nodes::Node* node = meles[ele]->nodes()[i];
        if (!node) FOUR_C_THROW("Cannot find master node");
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
        linsize += mnode->get_linsize();

        // if master node is also corner node
        if (mnode->is_on_corner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end()) continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->id());

          // aux variables
          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<Core::Gen::Pairedvector<int, double>> auxlin(
              3, linsize + 1 + meles[ele]->num_node());

          // compute distance between corners
          dist = compute_normal_node_to_node(mrtrnode, *mnode, auxnormal, auxlin);

          // if nodes lying on each other
          if (abs(dist) < 1e-12)
          {
            nodeOnNode = true;
            continue;
          }

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          // get closest valid distance
          if (abs(dist) < abs(gdist))
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }  // end node loop
    }    // end element loop
  }

  // loop over found eles
  for (size_t ele = 0; ele < meles.size(); ++ele)
  {
    bool cornerele = false;
    // get linsize
    int linsize = 0;
    for (int i = 0; i < meles[ele]->num_node(); ++i)
    {
      Core::Nodes::Node* node = meles[ele]->nodes()[i];
      if (!node) FOUR_C_THROW("Cannot find master node");
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
      linsize += mnode->get_linsize();

      if (mnode->is_on_corner()) cornerele = true;
    }

    double xi[2] = {0.0, 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<Core::Gen::Pairedvector<int, double>> auxlin(
        3, linsize + 1 + meles[ele]->num_node());

    // check for nonsmooth mele
    if (cornerele and !nodeOnNode)
    {
      // perform CPP to find normals based on element normals
      Mortar::Projector::impl(*meles[ele])
          ->project_s_node_by_m_normal_lin(mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);
    }
    // compute normal with averaged nodal normal field from master surface
    else
    {
      // perform CPP to find normals based on averaged nodal normal field
      Mortar::Projector::impl(*meles[ele])
          ->project_s_node_by_m_nodal_normal_lin(
              mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);
    }

    // check if found parameter space coordinate is within element domain
    if (meles[ele]->shape() == Core::FE::CellType::line2 or
        meles[ele]->shape() == Core::FE::CellType::line3)
    {
      if (-1.0 - tol > xi[0] or xi[0] > 1.0 + tol) continue;
    }
    else
    {
      FOUR_C_THROW("Unknown ele type!");
    }

    // angle between trajectory and normal
    if (pathdependent)
    {
      double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1]);
      if (auxl < 1e-12) continue;
      double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl);

      angle = 180 * (angle / 3.14159265359);
      if (abs(angle) > validAngle) continue;
    }

    // get closest valid distance
    if (abs(dist) < abs(gdist))
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }  // end mele loop

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::compute_cpp_normal3_d(Mortar::Node& mrtrnode,
    std::vector<Mortar::Element*> meles, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin)
{
  // define tolerance
  const double tol = 1e-8;

  // distance between node and surface
  bool pathdependent = true;
  const double validAngle = 5.0;
  double gdist = 1e12;  // distance
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<Core::Gen::Pairedvector<int, double>> glin(3, 1);  // 1 dummy

  //******************************************************
  //             CALC TRAJECTORY
  //******************************************************
  std::array<double, 3> Posn = {0.0, 0.0, 0.0};
  std::array<double, 3> Posnp = {0.0, 0.0, 0.0};
  std::array<double, 3> vect = {0.0, 0.0, 0.0};

  Posn[0] = mrtrnode.x()[0] + mrtrnode.uold()[0];
  Posn[1] = mrtrnode.x()[1] + mrtrnode.uold()[1];
  Posn[2] = mrtrnode.x()[2] + mrtrnode.uold()[2];
  Posnp[0] = mrtrnode.xspatial()[0];
  Posnp[1] = mrtrnode.xspatial()[1];
  Posnp[2] = mrtrnode.xspatial()[2];

  vect[0] = Posnp[0] - Posn[0];
  vect[1] = Posnp[1] - Posn[1];
  vect[2] = Posnp[2] - Posn[2];

  double lvec = sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
  if (lvec < 1e-12)
  {
    pathdependent = false;
  }
  else
  {
    vect[0] /= lvec;
    vect[1] /= lvec;
    vect[2] /= lvec;
  }

  if (dynamic_cast<CONTACT::Node&>(mrtrnode).active()) pathdependent = false;

  //******************************************************
  //             COMPUTE NORMAL TO SURFACE
  //******************************************************
  // loop over found eles for all geometrical nodes
  for (size_t ele = 0; ele < meles.size(); ++ele)
  {
    double xi[2] = {0.0, 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<Core::Gen::Pairedvector<int, double>> auxlin(3, 1000);

    // perform CPP to find normals
    bool success = Mortar::Projector::impl(*meles[ele])
                       ->project_s_node_by_m_nodal_normal_lin(
                           mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);

    // newton not converged
    if (!success) continue;

    // check if found parameter space coordinate is within element domain
    if (meles[ele]->shape() == Core::FE::CellType::quad4 or
        meles[ele]->shape() == Core::FE::CellType::quad8 or
        meles[ele]->shape() == Core::FE::CellType::quad9)
    {
      if (-1.0 - tol > xi[0] or xi[0] > 1.0 + tol or -1.0 - tol > xi[1] or xi[1] > 1.0 + tol)
        continue;
    }
    else if (meles[ele]->shape() == Core::FE::CellType::tri3 or
             meles[ele]->shape() == Core::FE::CellType::tri6)
    {
      if (xi[0] < 0.0 - tol or xi[1] < 0.0 - tol or xi[0] > 1.0 + tol or xi[1] > 1.0 + tol or
          xi[0] + xi[1] > 1.0 + 2 * tol)
        continue;
    }
    else
    {
      FOUR_C_THROW("Unknown ele type!");
    }

    // angle between trajectory and normal
    if (pathdependent)
    {
      double auxl = sqrt(
          auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);
      if (auxl < 1e-12) continue;
      double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                          vect[2] * auxnormal[2] / auxl);

      angle = 180 * (angle / 3.14159265359);
      if (abs(angle) > validAngle) continue;
    }

    if (dist < gdist)
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }  // end mele loop
  //******************************************************
  //             COMPUTE NORMAL TO LINE
  //******************************************************
  if (mrtrnode.is_on_corner_edge())  // only for edge or corner nodes possible
  {
    // guarantee uniquness
    std::set<std::pair<int, int>> donebefore;

    // calc
    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // loop over master edges -> match node number for quad4
      for (int j = 0; j < meles[ele]->num_node(); ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (meles[ele]->shape() == Core::FE::CellType::quad4)
        {
          if (j == 0)
          {
            nodeIds[0] = meles[ele]->node_ids()[0];
            nodeIds[1] = meles[ele]->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meles[ele]->node_ids()[1];
            nodeIds[1] = meles[ele]->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meles[ele]->node_ids()[2];
            nodeIds[1] = meles[ele]->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = meles[ele]->node_ids()[3];
            nodeIds[1] = meles[ele]->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            FOUR_C_THROW("loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<Mortar::Node*>(meles[ele]->nodes()[nodeLIds[0]])->is_on_corner_edge();
        bool node1Edge =
            dynamic_cast<Mortar::Node*>(meles[ele]->nodes()[nodeLIds[1]])->is_on_corner_edge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, meles[ele]->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          Core::Nodes::Node* nodes[2] = {
              meles[ele]->nodes()[nodeLIds[0]], meles[ele]->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes);

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // call cpp function for edge to edge

          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<Core::Gen::Pairedvector<int, double>> auxlin(
              3, 100 + 1 + meles[ele]->num_node());

          // compute distance between node and edge
          dist = compute_normal_node_to_edge(mrtrnode, *lineEle, auxnormal, auxlin);

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] +
                               auxnormal[2] * auxnormal[2]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                                vect[2] * auxnormal[2] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          if (dist <= gdist + tol)
          {
            std::cout << "CLOSE TO EDGE!!!" << std::endl;
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }  // end mele node loop
    }
  }

  //******************************************************
  //             COMPUTE NORMAL TO NODE
  //******************************************************
  if (mrtrnode.is_on_corner())  // only for corner nodes possible
  {
    std::set<int> donebeforeMasterCorner;

    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0; i < meles[ele]->num_node(); ++i)
      {
        Core::Nodes::Node* node = meles[ele]->nodes()[i];
        if (!node) FOUR_C_THROW("Cannot find master node");
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
        linsize += mnode->get_linsize();

        // if master node is also corner node
        if (mnode->is_on_corner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end()) continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->id());

          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<Core::Gen::Pairedvector<int, double>> auxlin(
              3, linsize + 1 + meles[ele]->num_node());

          // compute distance between corners
          dist = compute_normal_node_to_node(mrtrnode, *mnode, auxnormal, auxlin);

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] +
                               auxnormal[2] * auxnormal[2]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                                vect[2] * auxnormal[2] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          if (dist < gdist)
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }
    }
  }

  //******************************************************
  //             FINAL STORAGE
  //******************************************************

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  // bye bye
  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 05/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::compute_cpp_normal(Mortar::Node& mrtrnode,
    std::vector<Mortar::Element*> meles, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normaltolineLin)
{
  // define distance
  double gdist = 1e12;

  //===================================================================
  //===================================================================
  //                           2D case
  //===================================================================
  //===================================================================
  if (n_dim() == 2)
  {
    gdist = compute_cpp_normal2_d(mrtrnode, meles, normal, normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           3D case
  //===================================================================
  //===================================================================
  else if (n_dim() == 3)
  {
    gdist = compute_cpp_normal3_d(mrtrnode, meles, normal, normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           Invalid
  //===================================================================
  //===================================================================
  else
  {
    FOUR_C_THROW("invalid dimension!");
  }

  // return distance
  return gdist;
}

/*----------------------------------------------------------------------*
 |  set cpp normal                                           farah 01/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::set_cpp_normal(Mortar::Node& snode, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normallin)
{
  Node& cnode = dynamic_cast<Node&>(snode);

  const double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  if (length < 1e-12) FOUR_C_THROW("normal length is zero!!!");

  // negative sign because it is a master normal!
  cnode.mo_data().n()[0] = -normal[0] / length;
  cnode.mo_data().n()[1] = -normal[1] / length;
  cnode.mo_data().n()[2] = -normal[2] / length;

  //  if (cnode.IsOnEdge())
  //    std::cout << "normal =  " << cnode.MoData().n()[0] << "  "<< cnode.MoData().n()[1] << "  "<<
  //    cnode.MoData().n()[2] << std::endl;

  // prepare nodal storage maps for derivative
  if ((int)cnode.data().get_deriv_n().size() == 0)
    cnode.data().get_deriv_n().resize(3, normallin[0].size() * 3);
  if ((int)cnode.data().get_deriv_txi().size() == 0)
    cnode.data().get_deriv_txi().resize(3, normallin[0].size() * 3);
  if ((int)cnode.data().get_deriv_teta().size() == 0)
    cnode.data().get_deriv_teta().resize(3, normallin[0].size() * 3);

  // init tangent length
  double ltxi = -1.0;

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.num_dof() == 2)
  {
    // simple definition for txi
    cnode.data().txi()[0] = -cnode.mo_data().n()[1];
    cnode.data().txi()[1] = cnode.mo_data().n()[0];
    cnode.data().txi()[2] = 0.0;

    // teta is z-axis
    cnode.data().teta()[0] = 0.0;
    cnode.data().teta()[1] = 0.0;
    cnode.data().teta()[2] = 1.0;
  }
  // 3D Tangent!
  else
  {
    if (abs(cnode.mo_data().n()[0]) > 1.0e-6 || abs(cnode.mo_data().n()[1]) > 1.0e-6)
    {
      cnode.data().txi()[0] = -cnode.mo_data().n()[1];
      cnode.data().txi()[1] = cnode.mo_data().n()[0];
      cnode.data().txi()[2] = 0.0;
    }
    else
    {
      cnode.data().txi()[0] = 0.0;
      cnode.data().txi()[1] = -cnode.mo_data().n()[2];
      cnode.data().txi()[2] = cnode.mo_data().n()[1];
    }

    ltxi = sqrt(cnode.data().txi()[0] * cnode.data().txi()[0] +
                cnode.data().txi()[1] * cnode.data().txi()[1] +
                cnode.data().txi()[2] * cnode.data().txi()[2]);
    if (ltxi < 1e-12) FOUR_C_THROW("tangent txi length is zero!!!");
    for (int j = 0; j < 3; ++j) cnode.data().txi()[j] /= ltxi;

    // teta follows from corkscrew rule (teta = n x txi)
    cnode.data().teta()[0] = cnode.mo_data().n()[1] * cnode.data().txi()[2] -
                             cnode.mo_data().n()[2] * cnode.data().txi()[1];
    cnode.data().teta()[1] = cnode.mo_data().n()[2] * cnode.data().txi()[0] -
                             cnode.mo_data().n()[0] * cnode.data().txi()[2];
    cnode.data().teta()[2] = cnode.mo_data().n()[0] * cnode.data().txi()[1] -
                             cnode.mo_data().n()[1] * cnode.data().txi()[0];
  }


  //------------------------------------------------------------------
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;

  for (CI p = normallin[0].begin(); p != normallin[0].end(); ++p)
    (cnode.data().get_deriv_n()[0])[p->first] -= (p->second);
  for (CI p = normallin[1].begin(); p != normallin[1].end(); ++p)
    (cnode.data().get_deriv_n()[1])[p->first] -= (p->second);
  for (CI p = normallin[2].begin(); p != normallin[2].end(); ++p)
    (cnode.data().get_deriv_n()[2])[p->first] -= (p->second);

  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
  Core::Gen::Pairedvector<int, double>& derivnx = cnode.data().get_deriv_n()[0];
  Core::Gen::Pairedvector<int, double>& derivny = cnode.data().get_deriv_n()[1];
  Core::Gen::Pairedvector<int, double>& derivnz = cnode.data().get_deriv_n()[2];
  Core::Gen::Pairedvector<int, double> cderivnx = cnode.data().get_deriv_n()[0];
  Core::Gen::Pairedvector<int, double> cderivny = cnode.data().get_deriv_n()[1];
  Core::Gen::Pairedvector<int, double> cderivnz = cnode.data().get_deriv_n()[2];
  const double nxnx = cnode.mo_data().n()[0] * cnode.mo_data().n()[0];
  const double nxny = cnode.mo_data().n()[0] * cnode.mo_data().n()[1];
  const double nxnz = cnode.mo_data().n()[0] * cnode.mo_data().n()[2];
  const double nyny = cnode.mo_data().n()[1] * cnode.mo_data().n()[1];
  const double nynz = cnode.mo_data().n()[1] * cnode.mo_data().n()[2];
  const double nznz = cnode.mo_data().n()[2] * cnode.mo_data().n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeysn;
  for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivny.begin(); p != derivny.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }

  // normalize x-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnx[allkeysn[j]];
    derivnx[allkeysn[j]] =
        (val - nxnx * val - nxny * cderivny[allkeysn[j]] - nxnz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize y-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivny[allkeysn[j]];
    derivny[allkeysn[j]] =
        (val - nxny * cderivnx[allkeysn[j]] - nyny * val - nynz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize z-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnz[allkeysn[j]];
    derivnz[allkeysn[j]] =
        (val - nxnz * cderivnx[allkeysn[j]] - nynz * cderivny[allkeysn[j]] - nznz * val) / length;
  }

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.num_dof() == 2)
  {
    for (CI p = cnode.data().get_deriv_n()[1].begin(); p != cnode.data().get_deriv_n()[1].end();
         ++p)
      (cnode.data().get_deriv_txi()[0])[p->first] -= (p->second);
    for (CI p = cnode.data().get_deriv_n()[0].begin(); p != cnode.data().get_deriv_n()[0].end();
         ++p)
      (cnode.data().get_deriv_txi()[1])[p->first] += (p->second);
  }
  // 3D Tangent!
  else
  {
    // unnormalized tangent derivative txi
    // use definitions for txi from BuildAveragedNormal()
    if (abs(cnode.mo_data().n()[0]) > 1.0e-6 || abs(cnode.mo_data().n()[1]) > 1.0e-6)
    {
      Core::Gen::Pairedvector<int, double>& derivtxix = cnode.data().get_deriv_txi()[0];
      Core::Gen::Pairedvector<int, double>& derivtxiy = cnode.data().get_deriv_txi()[1];

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxix[p->first] -= (p->second);

      for (CI p = derivnx.begin(); p != derivnx.end(); ++p) derivtxiy[p->first] += (p->second);
    }
    else
    {
      Core::Gen::Pairedvector<int, double>& derivtxiy = cnode.data().get_deriv_txi()[1];
      Core::Gen::Pairedvector<int, double>& derivtxiz = cnode.data().get_deriv_txi()[2];

      for (CI p = derivnz.begin(); p != derivnz.end(); ++p) derivtxiy[p->first] -= (p->second);

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxiz[p->first] += (p->second);
    }

    if (ltxi < 1e-12) FOUR_C_THROW("tangent txi length is zero!!!");

    // normalize txi directional derivative
    // (identical to normalization of normal derivative)
    typedef Core::Gen::Pairedvector<int, double>::const_iterator CI;
    Core::Gen::Pairedvector<int, double>& derivtxix = cnode.data().get_deriv_txi()[0];
    Core::Gen::Pairedvector<int, double>& derivtxiy = cnode.data().get_deriv_txi()[1];
    Core::Gen::Pairedvector<int, double>& derivtxiz = cnode.data().get_deriv_txi()[2];
    Core::Gen::Pairedvector<int, double> cderivtxix = cnode.data().get_deriv_txi()[0];
    Core::Gen::Pairedvector<int, double> cderivtxiy = cnode.data().get_deriv_txi()[1];
    Core::Gen::Pairedvector<int, double> cderivtxiz = cnode.data().get_deriv_txi()[2];
    const double txtx = cnode.data().txi()[0] * cnode.data().txi()[0];
    const double txty = cnode.data().txi()[0] * cnode.data().txi()[1];
    const double txtz = cnode.data().txi()[0] * cnode.data().txi()[2];
    const double tyty = cnode.data().txi()[1] * cnode.data().txi()[1];
    const double tytz = cnode.data().txi()[1] * cnode.data().txi()[2];
    const double tztz = cnode.data().txi()[2] * cnode.data().txi()[2];

    // build a vector with all keys from x,y,z maps
    // (we need this in order not to miss any entry!)
    std::vector<int> allkeyst;
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }

    // normalize x-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxix[allkeyst[j]];
      derivtxix[allkeyst[j]] =
          (val - txtx * val - txty * cderivtxiy[allkeyst[j]] - txtz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize y-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiy[allkeyst[j]];
      derivtxiy[allkeyst[j]] =
          (val - txty * cderivtxix[allkeyst[j]] - tyty * val - tytz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize z-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiz[allkeyst[j]];
      derivtxiz[allkeyst[j]] =
          (val - txtz * cderivtxix[allkeyst[j]] - tytz * cderivtxiy[allkeyst[j]] - tztz * val) /
          ltxi;
    }

    // get normalized tangent derivative teta
    // use corkscrew rule from BuildAveragedNormal()
    Core::Gen::Pairedvector<int, double>& derivtetax = cnode.data().get_deriv_teta()[0];
    Core::Gen::Pairedvector<int, double>& derivtetay = cnode.data().get_deriv_teta()[1];
    Core::Gen::Pairedvector<int, double>& derivtetaz = cnode.data().get_deriv_teta()[2];

    for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
    {
      derivtetay[p->first] -= cnode.data().txi()[2] * (p->second);
      derivtetaz[p->first] += cnode.data().txi()[1] * (p->second);
    }
    for (CI p = derivny.begin(); p != derivny.end(); ++p)
    {
      derivtetax[p->first] += cnode.data().txi()[2] * (p->second);
      derivtetaz[p->first] -= cnode.data().txi()[0] * (p->second);
    }
    for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
    {
      derivtetax[p->first] -= cnode.data().txi()[1] * (p->second);
      derivtetay[p->first] += cnode.data().txi()[0] * (p->second);
    }
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      derivtetay[p->first] += cnode.mo_data().n()[2] * (p->second);
      derivtetaz[p->first] -= cnode.mo_data().n()[1] * (p->second);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      derivtetax[p->first] -= cnode.mo_data().n()[2] * (p->second);
      derivtetaz[p->first] += cnode.mo_data().n()[0] * (p->second);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      derivtetax[p->first] += cnode.mo_data().n()[1] * (p->second);
      derivtetay[p->first] -= cnode.mo_data().n()[0] * (p->second);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::export_nodal_normals() const
{
  // create empty data objects
  std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>> triad;

  std::map<int, std::vector<int>> n_x_key;
  std::map<int, std::vector<int>> n_y_key;
  std::map<int, std::vector<int>> n_z_key;
  std::map<int, std::vector<int>> txi_x_key;
  std::map<int, std::vector<int>> txi_y_key;
  std::map<int, std::vector<int>> txi_z_key;
  std::map<int, std::vector<int>> teta_x_key;
  std::map<int, std::vector<int>> teta_y_key;
  std::map<int, std::vector<int>> teta_z_key;

  std::map<int, std::vector<double>> n_x_val;
  std::map<int, std::vector<double>> n_y_val;
  std::map<int, std::vector<double>> n_z_val;
  std::map<int, std::vector<double>> txi_x_val;
  std::map<int, std::vector<double>> txi_y_val;
  std::map<int, std::vector<double>> txi_z_val;
  std::map<int, std::vector<double>> teta_x_val;
  std::map<int, std::vector<double>> teta_y_val;
  std::map<int, std::vector<double>> teta_z_val;

  std::map<int, double>::iterator iter;
  Core::Gen::Pairedvector<int, double>::iterator _iter;

  // build info on row map
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // fill nodal matrix
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc =
        Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, 3));
    (*loc)(0, 0) = cnode->mo_data().n()[0];
    (*loc)(1, 0) = cnode->mo_data().n()[1];
    (*loc)(2, 0) = cnode->mo_data().n()[2];
    (*loc)(0, 1) = cnode->data().txi()[0];
    (*loc)(1, 1) = cnode->data().txi()[1];
    (*loc)(2, 1) = cnode->data().txi()[2];
    (*loc)(0, 2) = cnode->data().teta()[0];
    (*loc)(1, 2) = cnode->data().teta()[1];
    (*loc)(2, 2) = cnode->data().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<Core::Gen::Pairedvector<int, double>>& derivn = cnode->data().get_deriv_n();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivtxi = cnode->data().get_deriv_txi();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivteta = cnode->data().get_deriv_teta();

    for (_iter = derivn[0].begin(); _iter != derivn[0].end(); ++_iter)
    {
      n_x_key[gid].push_back(_iter->first);
      n_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivn[1].begin(); _iter != derivn[1].end(); ++_iter)
    {
      n_y_key[gid].push_back(_iter->first);
      n_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivn[2].begin(); _iter != derivn[2].end(); ++_iter)
    {
      n_z_key[gid].push_back(_iter->first);
      n_z_val[gid].push_back(_iter->second);
    }

    for (_iter = derivtxi[0].begin(); _iter != derivtxi[0].end(); ++_iter)
    {
      txi_x_key[gid].push_back(_iter->first);
      txi_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivtxi[1].begin(); _iter != derivtxi[1].end(); ++_iter)
    {
      txi_y_key[gid].push_back(_iter->first);
      txi_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivtxi[2].begin(); _iter != derivtxi[2].end(); ++_iter)
    {
      txi_z_key[gid].push_back(_iter->first);
      txi_z_val[gid].push_back(_iter->second);
    }

    for (_iter = derivteta[0].begin(); _iter != derivteta[0].end(); ++_iter)
    {
      teta_x_key[gid].push_back(_iter->first);
      teta_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivteta[1].begin(); _iter != derivteta[1].end(); ++_iter)
    {
      teta_y_key[gid].push_back(_iter->first);
      teta_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivteta[2].begin(); _iter != derivteta[2].end(); ++_iter)
    {
      teta_z_key[gid].push_back(_iter->first);
      teta_z_val[gid].push_back(_iter->second);
    }
  }


  // communicate from slave node row to column map
  Core::Communication::Exporter& ex = interface_data_->exporter();

  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for (int i = 0; i < snodecolmapbound_->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    int linsize = cnode->get_linsize() + (int)(n_x_key[gid].size());

    // extract info
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc = triad[gid];
    cnode->mo_data().n()[0] = (*loc)(0, 0);
    cnode->mo_data().n()[1] = (*loc)(1, 0);
    cnode->mo_data().n()[2] = (*loc)(2, 0);
    cnode->data().txi()[0] = (*loc)(0, 1);
    cnode->data().txi()[1] = (*loc)(1, 1);
    cnode->data().txi()[2] = (*loc)(2, 1);
    cnode->data().teta()[0] = (*loc)(0, 2);
    cnode->data().teta()[1] = (*loc)(1, 2);
    cnode->data().teta()[2] = (*loc)(2, 2);

    // extract derivative info
    std::vector<Core::Gen::Pairedvector<int, double>>& derivn = cnode->data().get_deriv_n();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivtxi = cnode->data().get_deriv_txi();
    std::vector<Core::Gen::Pairedvector<int, double>>& derivteta = cnode->data().get_deriv_teta();

    for (int k = 0; k < (int)(derivn.size()); ++k) derivn[k].clear();
    derivn.resize(3, linsize);
    for (int k = 0; k < (int)(derivtxi.size()); ++k) derivtxi[k].clear();
    derivtxi.resize(3, linsize);
    for (int k = 0; k < (int)(derivteta.size()); ++k) derivteta[k].clear();
    derivteta.resize(3, linsize);

    for (int k = 0; k < (int)(n_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k = 0; k < (int)(n_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k = 0; k < (int)(n_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_n()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k = 0; k < (int)(txi_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k = 0; k < (int)(txi_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k = 0; k < (int)(txi_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_txi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k = 0; k < (int)(teta_x_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k = 0; k < (int)(teta_y_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k = 0; k < (int)(teta_z_key[gid].size()); ++k)
      (cnode->data().get_deriv_teta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();
  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::evaluate_search_binarytree()
{
  // *********************************************************************
  // self contact:
  // *********************************************************************
  // We call evaluate_search(), which does both, the bottom-up update (on the whole interface) and
  // the search. Then the dynamic master/slave assignment routine update_master_slave_sets() is
  // called and the new slave nodes' data containers are initialized.
  // *********************************************************************
  if (self_contact())
  {
    // evaluate search itself
    binarytreeself_->evaluate_search();

    // update master/slave sets of interface
    update_master_slave_sets();

    // initialize node data container
    // (include slave side boundary nodes / crosspoints)
    for (int i = 0; i < slave_col_nodes_bound()->NumMyElements(); ++i)
    {
      int gid = slave_col_nodes_bound()->GID(i);
      Core::Nodes::Node* node = discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %i", gid);
      Mortar::Node* mnode = dynamic_cast<Mortar::Node*>(node);

      // initialize container if not yet initialized before
      mnode->initialize_data_container();
    }

    // no initialization of element data container as this would
    // possibly destroy the information on search elements again
    // (this was already done in set_element_areas())
  }

  else
  {
    // call mortar routine
    Mortar::Interface::evaluate_search_binarytree();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-line coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_stl()
{
  // check
  if (n_dim() == 2) FOUR_C_THROW("LTS algorithm only for 3D simulations!");

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // guarantee uniquness
    std::set<std::pair<int, int>> donebefore;

    // loop over found meles
    for (int j = 0; j < selement->mo_data().num_search_elements(); ++j)
    {
      int gid2 = selement->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->shape() == Core::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->node_ids()[0];
            nodeIds[1] = melement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->node_ids()[1];
            nodeIds[1] = melement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->node_ids()[2];
            nodeIds[1] = melement->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->node_ids()[3];
            nodeIds[1] = melement->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

          // if not then create ele
          if (iter == donebefore.end() and itertw == donebefore.end())
          {
            // add to set of processed nodes
            donebefore.insert(actIDs);
            donebefore.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
                j, melement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<Core::Nodes::Node*, 2> nodes = {
                melement->nodes()[nodeLIds[0]], melement->nodes()[nodeLIds[1]]};
            lineEle->build_nodal_pointers(nodes.data());

            // init data container for dual shapes
            lineEle->initialize_data_container();

            std::vector<Element*> seleElements;
            seleElements.push_back(selement);

            // create coupling object
            LineToSurfaceCoupling3d coup(*idiscret_, 3, interface_params(), *melement, lineEle,
                seleElements, LineToSurfaceCoupling3d::stl);

            // perform evaluate!
            coup.evaluate_coupling();
          }
        }  // end edge loop
      }
      else
        FOUR_C_THROW("LTS only for quad4!");

    }  // end found mele loop
  }    // end slave ele loop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type ndoe-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_nts_master()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(interface_params(), n_dim()));

  // guarantee uniquness
  std::set<int> donebefore;

  // loop over all slave col elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->zero_sized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->mo_data().num_search_elements(); ++j)
    {
      int gid2 = selement->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->zero_sized()) continue;

      for (int n = 0; n < (int)melement->num_node(); ++n)
      {
        Mortar::Node* cnode = dynamic_cast<Mortar::Node*>(melement->nodes()[n]);

        std::set<int>::iterator iter = donebefore.find(cnode->id());

        // if not evaluated before
        if (iter == donebefore.end())
        {
          std::vector<Mortar::Element*> dummy;
          Mortar::Element* sele = dynamic_cast<Mortar::Element*>(selement);
          dummy.push_back(sele);

          // call interpolation functions
          bool success = interpolator->interpolate(*cnode, dummy);
          if (success) donebefore.insert(cnode->id());
        }

      }  // node loop
    }    // found contact eles
  }      // sele loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_lts_master()
{
  // check
  if (n_dim() == 2) FOUR_C_THROW("LTS algorithm only for 3D simulations!");

  // guarantee uniquness
  std::set<std::pair<int, int>> donebefore;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // ele check
    if (selement->shape() != Core::FE::CellType::quad4 and
        selement->shape() != Core::FE::CellType::tri3)
      FOUR_C_THROW("LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElements;
    std::vector<Element*> meleElements;

    // compute slave normal
    double slaveN[3] = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0};

    Core::FE::CellType dt = selement->shape();
    if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
             dt == Core::FE::CellType::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      FOUR_C_THROW("auxiliary_plane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->compute_unit_normal_at_xi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->mo_data().num_search_elements(); ++j)
    {
      int gid2 = selement->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = {0.0, 0.0};

      Core::FE::CellType dt = melement->shape();
      if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
               dt == Core::FE::CellType::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        FOUR_C_THROW("auxiliary_plane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->compute_unit_normal_at_xi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->min_edge_size();
      const double mminedge = melement->min_edge_size();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if (abs(scaprod) < tol) continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if (meleElements.size() < 1) continue;

    for (int m = 0; m < (int)meleElements.size(); ++m)
    {
      // loop over slave edges -> match node number for tri3/quad4
      for (int j = 0; j < meleElements[m]->num_node(); ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (meleElements[m]->shape() == Core::FE::CellType::quad4)
        {
          if (j == 0)
          {
            nodeIds[0] = meleElements[m]->node_ids()[0];
            nodeIds[1] = meleElements[m]->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meleElements[m]->node_ids()[1];
            nodeIds[1] = meleElements[m]->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meleElements[m]->node_ids()[2];
            nodeIds[1] = meleElements[m]->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = meleElements[m]->node_ids()[3];
            nodeIds[1] = meleElements[m]->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            FOUR_C_THROW("loop counter and edge number do not match!");
        }
        else if (meleElements[m]->shape() == Core::FE::CellType::tri3)
        {
          if (j == 0)
          {
            nodeIds[0] = meleElements[m]->node_ids()[0];
            nodeIds[1] = meleElements[m]->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meleElements[m]->node_ids()[1];
            nodeIds[1] = meleElements[m]->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meleElements[m]->node_ids()[2];
            nodeIds[1] = meleElements[m]->node_ids()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }
          else
            FOUR_C_THROW("loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<Mortar::Node*>(meleElements[m]->nodes()[nodeLIds[0]])->is_on_edge();
        bool node1Edge =
            dynamic_cast<Mortar::Node*>(meleElements[m]->nodes()[nodeLIds[1]])->is_on_edge();

        if (nonSmoothContact_ and (!node0Edge or !node1Edge)) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, meleElements[m]->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<Core::Nodes::Node*, 2> nodes = {
              meleElements[m]->nodes()[nodeLIds[0]], meleElements[m]->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes.data());

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // push back into vector
          lineElements.push_back(lineEle);

          std::vector<Element*> dummy;
          dummy.push_back(selement);

          // create coupling object
          LineToSurfaceCoupling3d coup(*idiscret_, 3, interface_params(), *meleElements[m], lineEle,
              dummy, LineToSurfaceCoupling3d::lts);

          // perform evaluate!
          coup.evaluate_coupling();
        }
      }  // end edge loop
    }
  }  // slave ele loop
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_lts()
{
  // check
  if (n_dim() == 2) FOUR_C_THROW("LTS algorithm only for 3D simulations!");

  // guarantee uniquness
  std::set<std::pair<int, int>> donebefore;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // ele check
    if (selement->shape() != Core::FE::CellType::quad4 and
        selement->shape() != Core::FE::CellType::tri3)
      FOUR_C_THROW("LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElements;
    std::vector<Element*> meleElements;

    // compute slave normal
    double slaveN[3] = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0};

    Core::FE::CellType dt = selement->shape();
    if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
             dt == Core::FE::CellType::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      FOUR_C_THROW("auxiliary_plane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->compute_unit_normal_at_xi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->mo_data().num_search_elements(); ++j)
    {
      int gid2 = selement->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = {0.0, 0.0};

      Core::FE::CellType dt = melement->shape();
      if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
               dt == Core::FE::CellType::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        FOUR_C_THROW("auxiliary_plane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->compute_unit_normal_at_xi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->min_edge_size();
      const double mminedge = melement->min_edge_size();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if (abs(scaprod) < tol) continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if (meleElements.size() < 1) continue;

    // loop over slave edges -> match node number for tri3/quad4
    for (int j = 0; j < selement->num_node(); ++j)
    {
      int nodeIds[2] = {0, 0};
      int nodeLIds[2] = {0, 0};

      if (selement->shape() == Core::FE::CellType::quad4)
      {
        if (j == 0)
        {
          nodeIds[0] = selement->node_ids()[0];
          nodeIds[1] = selement->node_ids()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->node_ids()[1];
          nodeIds[1] = selement->node_ids()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->node_ids()[2];
          nodeIds[1] = selement->node_ids()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->node_ids()[3];
          nodeIds[1] = selement->node_ids()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }
        else
          FOUR_C_THROW("loop counter and edge number do not match!");
      }
      else if (selement->shape() == Core::FE::CellType::tri3)
      {
        if (j == 0)
        {
          nodeIds[0] = selement->node_ids()[0];
          nodeIds[1] = selement->node_ids()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->node_ids()[1];
          nodeIds[1] = selement->node_ids()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->node_ids()[2];
          nodeIds[1] = selement->node_ids()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }
        else
          FOUR_C_THROW("loop counter and edge number do not match!");
      }

      // check if both nodes on edge geometry
      bool node0Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[0]])->is_on_edge();
      bool node1Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[1]])->is_on_edge();

      if (nonSmoothContact_ and (!node0Edge or !node1Edge)) continue;

      // create pair
      std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
      std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

      // check if processed before
      std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
      std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

      // if not then create ele
      if (iter == donebefore.end() and itertw == donebefore.end())
      {
        // add to set of processed nodes
        donebefore.insert(actIDs);
        donebefore.insert(actIDstw);

        // create line ele:
        Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
            j, selement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

        // get nodes
        std::array<Core::Nodes::Node*, 2> nodes = {
            selement->nodes()[nodeLIds[0]], selement->nodes()[nodeLIds[1]]};
        lineEle->build_nodal_pointers(nodes.data());

        // init data container for dual shapes
        lineEle->initialize_data_container();

        // push back into vector
        lineElements.push_back(lineEle);
      }
    }  // end edge loop

    // loop over all created line elements
    for (int l = 0; l < (int)lineElements.size(); ++l)
    {
      // create coupling object
      LineToSurfaceCoupling3d coup(*idiscret_, 3, interface_params(), *selement, lineElements[l],
          meleElements, LineToSurfaceCoupling3d::lts);

      // perform evaluate!
      coup.evaluate_coupling();
    }
  }  // slave ele loop

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-line coupl                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_ltl()
{
  // check
  if (n_dim() == 2) FOUR_C_THROW("LTL algorithm only for 3D simulations!");

  // guarantee uniquness of slave edges
  std::set<std::pair<int, int>> donebeforeS;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElementsS;

    if (selement->shape() == Core::FE::CellType::quad4)
    {
      for (int j = 0; j < 4; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->node_ids()[0];
          nodeIds[1] = selement->node_ids()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->node_ids()[1];
          nodeIds[1] = selement->node_ids()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->node_ids()[2];
          nodeIds[1] = selement->node_ids()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->node_ids()[3];
          nodeIds[1] = selement->node_ids()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[0]])->is_on_edge();
        bool node1Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[1]])->is_on_edge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, selement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<Core::Nodes::Node*, 2> nodes = {
              selement->nodes()[nodeLIds[0]], selement->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes.data());

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else if (selement->shape() == Core::FE::CellType::tri3)
    {
      for (int j = 0; j < 3; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->node_ids()[0];
          nodeIds[1] = selement->node_ids()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->node_ids()[1];
          nodeIds[1] = selement->node_ids()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->node_ids()[2];
          nodeIds[1] = selement->node_ids()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[0]])->is_on_edge();
        bool node1Edge = dynamic_cast<Mortar::Node*>(selement->nodes()[nodeLIds[1]])->is_on_edge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
              j, selement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<Core::Nodes::Node*, 2> nodes = {
              selement->nodes()[nodeLIds[0]], selement->nodes()[nodeLIds[1]]};
          lineEle->build_nodal_pointers(nodes.data());

          // init data container for dual shapes
          lineEle->initialize_data_container();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else
      FOUR_C_THROW("LTL only for quad4 and tri3!");

    // guarantee uniquness of master edges
    std::set<std::pair<int, int>> donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<Mortar::Element>> lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->mo_data().num_search_elements(); ++k)
    {
      int gid2 = selement->mo_data().search_elements()[k];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->shape() == Core::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->node_ids()[0];
            nodeIds[1] = melement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->node_ids()[1];
            nodeIds[1] = melement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->node_ids()[2];
            nodeIds[1] = melement->node_ids()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->node_ids()[3];
            nodeIds[1] = melement->node_ids()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[0]])->is_on_edge();
          bool node1Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[1]])->is_on_edge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
                j, melement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<Core::Nodes::Node*, 2> nodes = {
                melement->nodes()[nodeLIds[0]], melement->nodes()[nodeLIds[1]]};
            lineEle->build_nodal_pointers(nodes.data());

            // init data container for dual shapes
            lineEle->initialize_data_container();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else if (melement->shape() == Core::FE::CellType::tri3)
      {
        for (int j = 0; j < 3; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->node_ids()[0];
            nodeIds[1] = melement->node_ids()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->node_ids()[1];
            nodeIds[1] = melement->node_ids()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->node_ids()[2];
            nodeIds[1] = melement->node_ids()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[0]])->is_on_edge();
          bool node1Edge =
              dynamic_cast<Mortar::Node*>(melement->nodes()[nodeLIds[1]])->is_on_edge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<Mortar::Element> lineEle = Teuchos::rcp(new Mortar::Element(
                j, melement->owner(), Core::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<Core::Nodes::Node*, 2> nodes = {
                melement->nodes()[nodeLIds[0]], melement->nodes()[nodeLIds[1]]};
            lineEle->build_nodal_pointers(nodes.data());

            // init data container for dual shapes
            lineEle->initialize_data_container();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else
        FOUR_C_THROW("LTL only for quad4 and tri3!");
    }  // end found mele loop

    // loop over slave edges
    for (int s = 0; s < (int)lineElementsS.size(); ++s)
    {
      // loop over master edges
      for (int m = 0; m < (int)lineElementsM.size(); ++m)
      {
        // create coupling object
        LineToLineCouplingPoint3d coup(
            *idiscret_, 3, interface_params(), lineElementsS[s], lineElementsM[m]);

        // perform evaluate!
        coup.evaluate_coupling();
      }
    }

  }  // end slave loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type node-to-segment coupl             farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_nts()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(interface_params(), n_dim()));

  // loop over slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Mortar::Node* mrtrnode = dynamic_cast<Mortar::Node*>(node);

    if (!mrtrnode->is_on_corner() and nonSmoothContact_) continue;

    if (mrtrnode->owner() != get_comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    // vector with possible contacting master eles
    std::vector<Mortar::Element*> meles;

    // fill vector with possibly contacting meles
    find_m_eles(*mrtrnode, meles);

    // skip calculation if no meles vector is empty
    if (meles.size() < 1) continue;

    // call interpolation functions
    interpolator->interpolate(*mrtrnode, meles);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlaps     popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::mortar_coupling(Mortar::Element* sele, std::vector<Mortar::Element*> mele,
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // do stuff before the actual coupling is going to be evaluated
  pre_mortar_coupling(sele, mele, mparams_ptr);

  // increase counter of slave/master pairs
  smpairs_ += (int)mele.size();

  // check if quadratic interpolation is involved
  bool quadratic = false;
  if (sele->is_quad()) quadratic = true;
  for (int m = 0; m < (int)mele.size(); ++m)
    if (mele[m]->is_quad()) quadratic = true;

  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (n_dim() == 2)
  {
    // *************************************************** linear 2D ***
    // ************************************************ quadratic 2D ***
    // neither quadratic interpolation nor mixed linear and quadratic
    // interpolation need any special treatment in the 2d case

    // create Coupling2dManager
    CONTACT::Coupling2dManager coup(discret(), n_dim(), quadratic, interface_params(), sele, mele);
    // evaluate
    coup.evaluate_coupling(mparams_ptr);

    // increase counter of slave/master integration pairs and intcells
    smintpairs_ += (int)mele.size();
    intcells_ += (int)mele.size();
  }
  // ************************************************************** 3D ***
  else if (n_dim() == 3)
  {
    // *************************************************** linear 3D ***
    if (!quadratic)
    {
      // create Coupling3dManager
      CONTACT::Coupling3dManager coup(
          discret(), n_dim(), quadratic, interface_params(), sele, mele);
      // evaluate
      coup.evaluate_coupling(mparams_ptr);

      // increase counter of slave/master integration pairs and intcells
      smintpairs_ += (int)mele.size();
      intcells_ += coup.integration_cells();
    }

    // ************************************************** quadratic 3D ***
    else
    {
      // create Coupling3dQuadManager
      CONTACT::Coupling3dQuadManager coup(
          discret(), n_dim(), quadratic, interface_params(), sele, mele);
      // evaluate
      coup.evaluate_coupling(mparams_ptr);
    }  // quadratic
  }    // 3D
  else
    FOUR_C_THROW("Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  // do stuff after the coupling evaluation
  post_mortar_coupling(sele, mele, mparams_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::integrate_kappa_penalty(CONTACT::Element& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.shape() == Core::FE::CellType::tri3 || sele.shape() == Core::FE::CellType::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0;
    sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0;
    sxia[1] = -1.0;
    sxib[0] = 1.0;
    sxib[1] = 1.0;
  }

  // ************************************************** quadratic 3D ***
  if (n_dim() == 3 && sele.is_quad())
  {
    // get LM interpolation and testing type
    Inpar::Mortar::LagMultQuad lmtype =
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(interface_params(), "LM_QUAD");

    // build linear integration elements from quadratic CElements
    std::vector<Teuchos::RCP<Mortar::IntElement>> sauxelements(0);
    split_int_elements(sele, sauxelements);

    // different options for mortar integration
    if (lmtype == Inpar::Mortar::lagmult_quad || lmtype == Inpar::Mortar::lagmult_lin)
    {
      // do the element integration of kappa and store into gap
      int nrow = sele.num_node();
      Teuchos::RCP<Core::LinAlg::SerialDenseVector> gseg =
          Teuchos::rcp(new Core::LinAlg::SerialDenseVector(nrow));

      // create a CONTACT integrator instance with correct num_gp and Dim
      CONTACT::Integrator integrator(imortar_, sele.shape(), get_comm());
      integrator.integrate_kappa_penalty(sele, sxia, sxib, gseg);

      // do the assembly into the slave nodes
      integrator.assemble_g(get_comm(), sele, *gseg);
    }

    else if (lmtype == Inpar::Mortar::lagmult_pwlin)
    {
      // integrate each int element seperately
      for (int i = 0; i < (int)sauxelements.size(); ++i)
      {
        // do the int element integration of kappa and store into gap
        int nrow = sauxelements[i]->num_node();
        Teuchos::RCP<Core::LinAlg::SerialDenseVector> gseg =
            Teuchos::rcp(new Core::LinAlg::SerialDenseVector(nrow));

        // create a CONTACT integrator instance with correct num_gp and Dim
        CONTACT::Integrator integrator(imortar_, sauxelements[i]->shape(), get_comm());
        integrator.integrate_kappa_penalty(sele, *(sauxelements[i]), sxia, sxib, gseg);

        // do the assembly into the slave nodes
        integrator.assemble_g(get_comm(), *(sauxelements[i]), *gseg);
      }
    }

    else
    {
      FOUR_C_THROW("integrate_kappa_penalty: Invalid case for 3D mortar contact LM interpolation");
    }
  }

  // *************************************************** other cases ***
  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.num_node();
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> gseg =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(nrow));

    // create a CONTACT integrator instance with correct num_gp and Dim
    CONTACT::Integrator integrator(imortar_, sele.shape(), get_comm());
    integrator.integrate_kappa_penalty(sele, sxia, sxib, gseg);

    // do the assembly into the slave nodes
    integrator.assemble_g(get_comm(), sele, *gseg);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node     gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_rel_mov(const Teuchos::RCP<Epetra_Vector> xsmod,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix> dmatrixmod,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix> doldmod)
{
  if (friction_ == false)
    FOUR_C_THROW(
        "Error in Interface::evaluate_relative_movement(): Only evaluated for frictional contact");

  // parameters
  double pp = interface_params().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);
    double cn = get_cn_ref()[get_cn_ref().Map().LID(cnode->id())];

    // get some informatiom form the node
    double gap = cnode->data().getg();

    const int dim = cnode->num_dof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k = 0; k < 3; ++k) nz += cnode->mo_data().n()[k] * cnode->mo_data().lm()[k];

    std::vector<double> jump(dim);
    for (int dim = 0; dim < n_dim(); dim++) jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k)
      lmuzawan += cnode->mo_data().lmuzawa()[k] * cnode->mo_data().n()[k];

    double kappa = cnode->data().kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
            interface_params(), "STRATEGY") == Inpar::CONTACT::solution_penalty ||
        Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
            interface_params(), "STRATEGY") == Inpar::CONTACT::solution_multiscale)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
                 interface_params(), "STRATEGY") == Inpar::CONTACT::solution_lagmult and
             Core::UTILS::IntegralValue<int>(interface_params(), "SEMI_SMOOTH_NEWTON") != 1)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
                 interface_params(), "STRATEGY") == Inpar::CONTACT::solution_lagmult and
             Core::UTILS::IntegralValue<int>(interface_params(), "SEMI_SMOOTH_NEWTON") == 1)
    {
      if ((nz - cn * gap > 0) or cnode->active()) activeinfuture = true;
    }
    else if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
                 interface_params(), "STRATEGY") == Inpar::CONTACT::solution_uzawa)
    {
      if (lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
      FOUR_C_THROW(
          "Error in Interface::evaluate_relative_movement(): Solution strategy not known!");

    if (activeinfuture == true)
    {
      Core::Gen::Pairedvector<int, double>& dmap = cnode->mo_data().get_d();
      Core::Gen::Pairedvector<int, double>& dmapold = cnode->fri_data().get_d_old();

      std::set<int> snodes = cnode->fri_data().get_s_nodes();

      // check if there are entries in the old D map
      if (dmapold.size() < 1)
        FOUR_C_THROW("Error in Interface::evaluate_relative_movement(): No old D-Map!");

      std::map<int, double>::iterator colcurr;
      std::set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr = snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        Core::Nodes::Node* snode = idiscret_->g_node(gid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* csnode = dynamic_cast<Node*>(snode);

        double dik = dmap[csnode->id()];
        double dikold = dmapold[csnode->id()];

        std::map<int, double>::iterator mcurr;

        for (int dim = 0; dim < csnode->num_dof(); ++dim)
        {
          int locid = (xsmod->Map()).LID(csnode->dofs()[dim]);
          jump[dim] -= (dik - dikold) * (*xsmod)[locid];
        }
      }  //  loop over adjacent slave nodes

      std::map<int, double>& mmap = cnode->mo_data().get_m();
      std::map<int, double>& mmapold = cnode->fri_data().get_m_old();

      const std::set<int>& mnodescurrent = cnode->fri_data().get_m_nodes();
      const std::set<int>& mnodesold = cnode->fri_data().get_m_nodes_old();

      // check if there are entries in the M map
      if (mmap.size() < 1)
        FOUR_C_THROW("Error in Interface::evaluate_relative_movement(): No M-Map!");

      // check if there are entries in the old M map
      if (mmapold.size() < 1)
        FOUR_C_THROW("Error in Interface::evaluate_relative_movement(): No old M-Map!");

      if (mnodesold.size() < 1)
        FOUR_C_THROW("Error in Interface::evaluate_relative_movement(): No old M-Set!");

      std::set<int> mnodes;
      std::set<int>::iterator mcurr;

      for (mcurr = mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr = mnodesold.begin(); mcurr != mnodesold.end(); mcurr++) mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        Core::Nodes::Node* mnode = idiscret_->g_node(gid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        double mik = mmap[cmnode->id()];
        double mikold = mmapold[cmnode->id()];

        std::map<int, double>::iterator mcurr;

        for (int dim = 0; dim < cnode->num_dof(); ++dim)
        {
          jump[dim] += (mik - mikold) * (cmnode->xspatial()[dim]);
        }
      }  //  loop over master nodes

      // write it to nodes
      for (int dim = 0; dim < n_dim(); dim++) cnode->fri_data().jump()[dim] = jump[dim];

      // linearization of jump vector

      // reset derivative map of jump
      for (int j = 0; j < (int)((cnode->fri_data().get_deriv_jump()).size()); ++j)
        (cnode->fri_data().get_deriv_jump())[j].clear();
      (cnode->fri_data().get_deriv_jump()).resize(0);

      /*** 01  **********************************************************/

      if (dmatrixmod == Teuchos::null)
      {
        // loop over according slave nodes
        for (scurr = snodes.begin(); scurr != snodes.end(); scurr++)
        {
          int gid = *scurr;
          Core::Nodes::Node* snode = idiscret_->g_node(gid);
          if (!snode) FOUR_C_THROW("Cannot find node with gid %", gid);
          Node* csnode = dynamic_cast<Node*>(snode);

          double dik = dmap[csnode->id()];
          double dikold = dmapold[csnode->id()];

          for (int dimrow = 0; dimrow < cnode->num_dof(); ++dimrow)
          {
            int col = csnode->dofs()[dimrow];
            double val = -(dik - dikold);
            if (abs(val) > 1e-14) cnode->add_deriv_jump_value(dimrow, col, val);
          }
        }
      }
      // in the 3D quadratic case, the values are obtained from the
      // global matrices Dmod and Doldmod
      else
      {
        // loop over dimension of the node
        for (int dim = 0; dim < cnode->num_dof(); ++dim)
        {
          int NumEntries = 0;
          int NumEntriesOld = 0;
          std::vector<double> Values((dmatrixmod->epetra_matrix())->MaxNumEntries());
          std::vector<int> Indices((dmatrixmod->epetra_matrix())->MaxNumEntries());
          std::vector<double> ValuesOld((dmatrixmod->epetra_matrix())->MaxNumEntries());
          std::vector<int> IndicesOld((dmatrixmod->epetra_matrix())->MaxNumEntries());

          // row
          int row = cnode->dofs()[dim];

          // extract entries of this row from matrix
          int err = (dmatrixmod->epetra_matrix())
                        ->ExtractGlobalRowCopy(row, (dmatrixmod->epetra_matrix())->MaxNumEntries(),
                            NumEntries, Values.data(), Indices.data());
          if (err) FOUR_C_THROW("ExtractMyRowView failed: err=%d", err);

          int errold = (doldmod->epetra_matrix())
                           ->ExtractGlobalRowCopy(row, (doldmod->epetra_matrix())->MaxNumEntries(),
                               NumEntriesOld, ValuesOld.data(), IndicesOld.data());
          if (errold) FOUR_C_THROW("ExtractMyRowView failed: err=%d", err);

          // loop over entries of this vector
          for (int j = 0; j < NumEntries; ++j)
          {
            double ValueOld = 0;
            bool found = false;

            // find value with the same index in vector of Dold
            for (int k = 0; k < NumEntriesOld; ++k)
            {
              if (Indices[k] == Indices[j])
              {
                ValueOld = ValuesOld[k];
                found = true;
                break;
              }
            }

            if (found == false or abs(ValueOld) < 1e-12)
              FOUR_C_THROW("Error in EvaluareRelMov(): No old D value exists");

            // write to node
            cnode->add_deriv_jump_value(dim, Indices[j], (Values[j] - ValueOld));
          }
        }
      }

      /*** 02  **********************************************************/
      // loop over according master nodes
      for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        Core::Nodes::Node* mnode = idiscret_->g_node(gid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        double mik = mmap[cmnode->id()];
        double mikold = mmapold[cmnode->id()];

        for (int dimrow = 0; dimrow < cnode->num_dof(); ++dimrow)
        {
          int col = cmnode->dofs()[dimrow];
          double val = (mik - mikold);
          if (abs(val) > 1e-14) cnode->add_deriv_jump_value(dimrow, col, val);
        }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
      std::map<int, std::map<int, double>>& ddmap = cnode->data().get_deriv_d();
      std::map<int, std::map<int, double>>::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr = ddmap.begin(); dscurr != ddmap.end(); ++dscurr)
      {
        int gid = dscurr->first;
        Core::Nodes::Node* snode = idiscret_->g_node(gid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // compute entry of the current stick node / slave node pair
        std::map<int, double>& thisdmmap = cnode->data().get_deriv_d(gid);

        // loop over all entries of the current derivative map
        for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for (int dim = 0; dim < cnode->num_dof(); ++dim)
          {
            int locid = (xsmod->Map()).LID(csnode->dofs()[dim]);
            double val = -colcurr->second * (*xsmod)[locid];
            if (abs(val) > 1e-14) cnode->add_deriv_jump_value(dim, col, val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
      std::map<int, std::map<int, double>>& dmmap = cnode->data().get_deriv_m();
      std::map<int, std::map<int, double>>::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr = dmmap.begin(); dmcurr != dmmap.end(); ++dmcurr)
      {
        int gid = dmcurr->first;
        Core::Nodes::Node* mnode = idiscret_->g_node(gid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
        std::map<int, double>& thisdmmap = cnode->data().get_deriv_m(gid);

        // loop over all entries of the current derivative map
        for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for (int dimrow = 0; dimrow < cnode->num_dof(); ++dimrow)
          {
            double val = colcurr->second * mxi[dimrow];
            if (abs(val) > 1e-14) cnode->add_deriv_jump_value(dimrow, col, val);
          }
        }
      }

      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        for (int j = 0; j < n_dim(); j++)
          if (cnode->dbc_dofs()[j] == true)
          {
            cnode->fri_data().jump()[j] = 0.;
            cnode->fri_data().get_deriv_jump()[j].clear();
          }

    }  // active nodes
  }    // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  calculate nodal distances (public)                     pfaller Jan15|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_distances(const Teuchos::RCP<const Epetra_Vector>& vec,
    std::map<int, std::vector<double>>& mynormals,
    std::map<int, std::vector<Core::Gen::Pairedvector<int, double>>>& dmynormals,
    std::map<int, double>& mygap, std::map<int, std::map<int, double>>& dmygap)
{
  set_state(Mortar::state_new_displacement, *vec);
  initialize();

  // interface needs to be complete
  if (!filled() && get_comm().MyPID() == 0)
    FOUR_C_THROW("fill_complete() not called on interface %", id_);

  // create an interpolator instance
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(imortar_, n_dim()));

  // create normals
  pre_evaluate(-1, -1);  // dummy values

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  get_comm().Barrier();

  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->g_element(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    CONTACT::Element* selement = dynamic_cast<CONTACT::Element*>(ele1);

    if (selement->mo_data().num_search_elements() < 1)
    {
      std::cout << "WARNING: No elements found!" << std::endl;
      continue;
    }

    // skip zero-sized nurbs elements (slave)
    if (selement->zero_sized()) continue;

    // empty vector of master element pointers
    std::vector<CONTACT::Element*> melements;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->mo_data().num_search_elements(); ++j)
    {
      int gid2 = selement->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      CONTACT::Element* melement = dynamic_cast<CONTACT::Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->zero_sized()) continue;

      melements.push_back(melement);
    }

    //**************************************************************
    //                loop over all Slave nodes
    //**************************************************************
    for (int snodes = 0; snodes < selement->num_node(); ++snodes)
    {
      CONTACT::Node* mynode = dynamic_cast<CONTACT::Node*>(selement->nodes()[snodes]);

      // skip this node if already considered
      if (mynode->has_proj()) continue;

      //                store node normals
      //**************************************************************
      //      int gid = snoderowmapbound_->GID(snodes);
      int gid = mynode->id();

      int numdofs = mynode->num_dof();
      std::vector<double> temp(numdofs, 0.0);
      for (int kk = 0; kk < numdofs; kk++)
      {
        temp[kk] = mynode->mo_data().n()[kk];
      }
      mynormals.insert(std::pair<int, std::vector<double>>(gid, temp));
      dmynormals.insert(std::pair<int, std::vector<Core::Gen::Pairedvector<int, double>>>(
          gid, mynode->data().get_deriv_n()));

      //**************************************************************
      double sxi[2] = {0.0, 0.0};

      if (selement->shape() == Core::FE::CellType::quad4 or
          selement->shape() == Core::FE::CellType::quad8 or
          selement->shape() == Core::FE::CellType::quad9)
      {
        // TODO (pfaller): switch case
        if (snodes == 0)
        {
          sxi[0] = -1;
          sxi[1] = -1;
        }
        else if (snodes == 1)
        {
          sxi[0] = 1;
          sxi[1] = -1;
        }
        else if (snodes == 2)
        {
          sxi[0] = 1;
          sxi[1] = 1;
        }
        else if (snodes == 3)
        {
          sxi[0] = -1;
          sxi[1] = 1;
        }
        else if (snodes == 4)
        {
          sxi[0] = 0;
          sxi[1] = -1;
        }
        else if (snodes == 5)
        {
          sxi[0] = 1;
          sxi[1] = 0;
        }
        else if (snodes == 6)
        {
          sxi[0] = 0;
          sxi[1] = 1;
        }
        else if (snodes == 7)
        {
          sxi[0] = -1;
          sxi[1] = 0;
        }
        else if (snodes == 8)
        {
          sxi[0] = 0;
          sxi[1] = 0;
        }
        else
          FOUR_C_THROW("ERORR: wrong node LID");
      }
      else if (selement->shape() == Core::FE::CellType::tri3 or
               selement->shape() == Core::FE::CellType::tri6)
      {
        if (snodes == 0)
        {
          sxi[0] = 0;
          sxi[1] = 0;
        }
        else if (snodes == 1)
        {
          sxi[0] = 1;
          sxi[1] = 0;
        }
        else if (snodes == 2)
        {
          sxi[0] = 0;
          sxi[1] = 1;
        }
        else if (snodes == 3)
        {
          sxi[0] = 0.5;
          sxi[1] = 0;
        }
        else if (snodes == 4)
        {
          sxi[0] = 0.5;
          sxi[1] = 0.5;
        }
        else if (snodes == 5)
        {
          sxi[0] = 0;
          sxi[1] = 0.5;
        }
        else
          FOUR_C_THROW("ERORR: wrong node LID");
      }
      else
      {
        FOUR_C_THROW("Chosen element type not supported for NTS!");
      }

      //**************************************************************
      //                loop over all Master Elements
      //**************************************************************
      // create vectors to store projection information for several master elements in case
      // projection is not unique
      std::vector<double> gap_vec;
      std::vector<std::map<int, double>> dgap_vec;

      for (int nummaster = 0; nummaster < (int)melements.size(); ++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        double projalpha = 0.0;
        bool is_projected =
            Mortar::Projector::impl(*selement, *melements[nummaster])
                ->project_gauss_point3_d(*selement, sxi, *melements[nummaster], mxi, projalpha);

        bool is_on_mele = true;

        // check GP projection
        Core::FE::CellType dt = melements[nummaster]->shape();
        const double tol = 1e-8;
        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol ||
              mxi[1] > 1.0 + tol)
          {
            is_on_mele = false;
          }
        }
        else
        {
          if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
              mxi[0] + mxi[1] > 1.0 + 2 * tol)
          {
            is_on_mele = false;
          }
        }

        // node on mele?
        if (is_on_mele && is_projected)
        {
          // store information of projection so that this node is not considered again
          mynode->has_proj() = true;

          int ndof = 3;
          int ncol = melements[nummaster]->num_node();
          Core::LinAlg::SerialDenseVector mval(ncol);
          Core::LinAlg::SerialDenseMatrix mderiv(ncol, 2);
          melements[nummaster]->evaluate_shape(mxi, mval, mderiv, ncol, false);

          //          int linsize    = mynode->GetLinsize();
          int linsize = 100;
          double gpn[3] = {0.0, 0.0, 0.0};
          //**************************************************************

          // evalute the GP slave coordinate derivatives --> no entries
          std::vector<Core::Gen::Pairedvector<int, double>> dsxi(2, 0);
          std::vector<Core::Gen::Pairedvector<int, double>> dmxi(2, 4 * linsize + ncol * ndof);

          (*interpolator)
              .deriv_xi_g_p3_d(*selement, *melements[nummaster], sxi, mxi, dsxi, dmxi, projalpha);
          (*interpolator).nw_gap3_d(*mynode, *melements[nummaster], mval, mderiv, dmxi, gpn);

          // store linearization for node
          std::map<int, double> dgap = mynode->data().get_deriv_gnts();  // (dof,value)

          // store gap information
          gap_vec.push_back(mynode->data().getgnts());
          dgap_vec.push_back(dgap);

          // reset nodal weighted gap and derivative
          mynode->data().getgnts() = 1.0e12;
          (mynode->data().get_deriv_gnts()).clear();
        }  // End hit ele
      }    // End Loop over all Master Elements

      if (gap_vec.size() > 0)
      {
        // find projection with smallest absoluate value of gap
        std::vector<double>::iterator iter_min =
            std::min_element(gap_vec.begin(), gap_vec.end(), abs_compare);
        const int i_min = std::distance(gap_vec.begin(), iter_min);

        // save to map at GID
        mygap.insert(std::pair<int, double>(gid, gap_vec[i_min]));
        dmygap.insert(std::pair<int, std::map<int, double>>(gid, dgap_vec[i_min]));
      }
    }
  }

  get_comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::evaluate_tangent_norm(double& cnormtan)
{
  // friction coefficient
  double frcoeff = interface_params().get<double>("FRCOEFF");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < slave_row_nodes()->NumMyElements(); ++i)
  {
    int gid = slave_row_nodes()->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some information from node
    double* n = cnode->mo_data().n();
    const int dim = cnode->num_dof();

    // tangential plane
    Core::LinAlg::SerialDenseMatrix tanplane(dim, dim);
    if (dim == 3)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);
      tanplane(0, 2) = -(n[0] * n[2]);
      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
      tanplane(1, 2) = -(n[1] * n[2]);

      tanplane(2, 0) = -(n[2] * n[0]);
      tanplane(2, 1) = -(n[2] * n[1]);
      tanplane(2, 2) = 1 - (n[2] * n[2]);
    }
    else if (dim == 2)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);

      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
    }
    else
      FOUR_C_THROW("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    Core::LinAlg::SerialDenseMatrix jumpvec(dim, 1);
    for (int i = 0; i < dim; i++) jumpvec(i, 0) = cnode->fri_data().jump()[i];

    // evaluate jump in tangential direction
    Core::LinAlg::SerialDenseMatrix jumptan(dim, 1);
    Core::LinAlg::multiply(jumptan, tanplane, jumpvec);

    // force vector
    Core::LinAlg::SerialDenseMatrix forcevec(dim, 1);
    for (int i = 0; i < dim; i++) forcevec(i, 0) = cnode->mo_data().lm()[i];

    // evaluate force in normal direction
    double forcen = 0.0;
    for (int k = 0; k < dim; ++k) forcen += forcevec(k, 0) * n[k];

    // norm of constraint violation for stick nodes
    if (cnode->active() == true and cnode->fri_data().slip() == false)
    {
      for (int j = 0; j < dim; j++) cnormtan += jumptan(j, 0) * jumptan(j, 0);
    }
    // norm of constraint violation for slip nodes
    else if (cnode->active() == true and cnode->fri_data().slip() == true)
    {
      double part1 = 0.0;
      double jumpnorm = 0.0;

      for (int j = 0; j < dim; j++)
      {
        jumpnorm += jumptan(j, 0) * jumptan(j, 0);
        part1 += jumptan(j, 0) * forcevec(j, 0);
      }

      jumpnorm = sqrt(jumpnorm);
      cnormtan += (part1 - frcoeff * forcen * jumpnorm) * (part1 - frcoeff * forcen * jumpnorm);
    }
  }  // loop over slave nodes

  // get cnorm from all procs
  double sumcnormtanallprocs = 0.0;
  get_comm().SumAll(&cnormtan, &sumcnormtanallprocs, 1);
  cnormtan = sumcnormtanallprocs;

  return;
}


/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::update_active_set_semi_smooth()
{
  // get input parameter ftype
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(interface_params(), "FRICTION");

  // this is the complementarity parameter we use for the decision.
  // it might be scaled with a mesh-size dependent factor
  double cn = 0.;
  double ct = 0.;

  // assume that active set has converged and check for opposite
  int localcheck = true;

  // loop over all slave nodes on the current interface
  for (int j = 0; j < slave_row_nodes()->NumMyElements(); ++j)
  {
    int gid = slave_row_nodes()->GID(j);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

    Node* cnode = dynamic_cast<Node*>(node);

    cn = get_cn_ref()[get_cn_ref().Map().LID(cnode->id())];
    if (friction_) ct = get_ct_ref()[get_ct_ref().Map().LID(cnode->id())];

    // get weighted gap
    double wgap = cnode->data().getg();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k = 0; k < 3; ++k)
    {
      nz += cnode->mo_data().n()[k] * cnode->mo_data().lm()[k];
    }

    // friction
    std::vector<double> tz(n_dim() - 1, 0);
    std::vector<double> tjump(n_dim() - 1, 0);
    double euclidean = 0.0;

    if (friction_)
    {
      // static cast
      FriNode* frinode = dynamic_cast<FriNode*>(cnode);

      // compute tangential parts and of Lagrange multiplier and incremental jumps
      for (int i = 0; i < n_dim(); ++i)
      {
        tz[0] += frinode->data().txi()[i] * frinode->mo_data().lm()[i];
        if (n_dim() == 3) tz[1] += frinode->data().teta()[i] * frinode->mo_data().lm()[i];

        if (Core::UTILS::IntegralValue<int>(interface_params(), "GP_SLIP_INCR") == false)
        {
          tjump[0] += frinode->data().txi()[i] * frinode->fri_data().jump()[i];
          if (n_dim() == 3) tjump[1] += frinode->data().teta()[i] * frinode->fri_data().jump()[i];
        }
      }

      if (Core::UTILS::IntegralValue<int>(interface_params(), "GP_SLIP_INCR") == true)
      {
        tjump[0] = frinode->fri_data().jump_var()[0];
        if (n_dim() == 3) tjump[1] = frinode->fri_data().jump_var()[1];
      }

      // evaluate euclidean norm |tz+ct.tjump|
      std::vector<double> sum(n_dim() - 1, 0);
      sum[0] = tz[0] + ct * tjump[0];
      if (n_dim() == 3) sum[1] = tz[1] + ct * tjump[1];
      if (n_dim() == 2) euclidean = abs(sum[0]);
      if (n_dim() == 3) euclidean = sqrt(sum[0] * sum[0] + sum[1] * sum[1]);
    }

    // adhesion
    double adhbound = 0.0;
    if (Core::UTILS::IntegralValue<Inpar::CONTACT::AdhesionType>(interface_params(), "ADHESION") ==
        Inpar::CONTACT::adhesion_bound)
      adhbound = interface_params().get<double>("ADHESION_BOUND");

    // check nodes of inactive set *************************************
    if (cnode->active() == false)
    {
      // check for penetration and/or tensile contact forces
      if (nz - cn * wgap >
          0)  // no averaging of Lagrange multipliers
              // if ((0.5*nz+0.5*nzold) - cn*wgap > 0) // averaging of Lagrange multipliers
      {
        cnode->active() = true;
        localcheck = false;

        // friction
        if (friction_)
        {
          // nodes coming into contact
          dynamic_cast<FriNode*>(cnode)->fri_data().slip() = true;
        }
      }
    }

    // check nodes of active set ***************************************
    else
    {
      // adhesion modification
      nz += adhbound;

      // check for tensile contact forces and/or penetration
      if (nz - cn * wgap <=
          0)  // no averaging of Lagrange multipliers
              // if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
      {
        cnode->active() = false;
        localcheck = false;

        // friction
        if (friction_) dynamic_cast<FriNode*>(cnode)->fri_data().slip() = false;
      }

      // only do something for friction
      else
      {
        // friction tresca
        if (ftype == Inpar::CONTACT::friction_tresca)
        {
          FriNode* frinode = dynamic_cast<FriNode*>(cnode);

          // CAREFUL: friction bound is now interface-local (popp 08/2012)
          double frbound = interface_params().get<double>("FRBOUND");

          if (frinode->fri_data().slip() == false)
          {
            // check (euclidean)-frbound <= 0
            if (euclidean - frbound <= 0)
            {
            }
            // do nothing (stick was correct)
            else
            {
              frinode->fri_data().slip() = true;
              localcheck = false;
            }
          }
          else
          {
            // check (euclidean)-frbound > 0
            if (euclidean - frbound > 0)
            {
            }
            // do nothing (slip was correct)
            else
            {
              frinode->fri_data().slip() = false;
              localcheck = false;
            }
          }
        }  // if (fytpe=="tresca")

        // friction coulomb
        if (ftype == Inpar::CONTACT::friction_coulomb)
        {
          FriNode* frinode = dynamic_cast<FriNode*>(cnode);

          // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
          double frcoeff = frinode->fr_coeff(interface_params().get<double>("FRCOEFF"));
          double frbound;
          static const bool regularization =
              Core::UTILS::IntegralValue<int>(interface_params(), "REGULARIZED_NORMAL_CONTACT");
          if (!regularization)
            frbound = frcoeff * (nz - cn * wgap);
          else
          {
            static const double k = 1. / interface_params().get<double>("REGULARIZATION_STIFFNESS");
            static const double gmax = interface_params().get<double>("REGULARIZATION_THICKNESS");
            if (cnode->mo_data().get_d().size() != 1)
              FOUR_C_THROW(
                  "we need to have a D-value for active contact nodes\nAnd exactly one due to "
                  "biorthogonality");
            double dval = cnode->mo_data().get_d().at(cnode->id());
            const double gLM = gmax * (1. - exp(-k / gmax * nz));
            frbound = frcoeff * std::max(0., nz - cn * (wgap + dval * gLM));
          }

          if (frinode->fri_data().slip() == false)
          {
            // check (euclidean)-frbound <= 0
            if (euclidean - frbound <= 1e-10)
            {
            }
            // do nothing (stick was correct)
            else
            {
              frinode->fri_data().slip() = true;
              localcheck = false;
            }
          }
          else
          {
            // check (euclidean)-frbound > 0
            if (euclidean - frbound > -1e-10)
            {
            }
            // do nothing (slip was correct)
            else
            {
              frinode->fri_data().slip() = false;
              localcheck = false;
            }
          }
        }  // if (ftype == Inpar::CONTACT::friction_coulomb)
      }    // if (nz - cn*wgap <= 0)
    }      // if (cnode->Active()==false)
  }        // loop over all slave nodes

  // broadcast convergence status among processors
  int convcheck = 0;
  get_comm().MinAll(&localcheck, &convcheck, 1);

  return convcheck;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::update_active_set_initial_status() const
{
  // List of GIDs of all my slave row nodes
  int* my_slave_row_node_gids = slave_row_nodes()->MyGlobalElements();

  // loop over all slave nodes on the current interface
  for (int j = 0; j < slave_row_nodes()->NumMyElements(); ++j)
  {
    // Grab the current slave node
    const int gid = my_slave_row_node_gids[j];
    Node* cnode = dynamic_cast<Node*>(discret().g_node(gid));
    if (!cnode) FOUR_C_THROW("Cannot find node with gid %", gid);

    set_node_initially_active(*cnode);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::build_active_set(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mynodegidsInactive(0);
  std::vector<int> mydofgids(0);
  std::vector<int> mydofgidsInactive(0);
  std::vector<int> myslipnodegids(0);
  std::vector<int> myslipdofgids(0);
  std::vector<int> mymnodegids(0);
  std::vector<int> mymdofgids(0);

  // loop over all slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->num_dof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the Node member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding Nodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // flag for initialization of init active nodes with nodal gaps
      bool initcontactbygap =
          Core::UTILS::IntegralValue<int>(interface_params(), "INITCONTACTBYGAP");
      // value
      double initcontactval = interface_params().get<double>("INITCONTACTGAPVALUE");

      // Either init contact by definition or by gap
      if (cnode->is_init_active() and initcontactbygap)
        FOUR_C_THROW("Init contact either by definition in condition or by gap!");

      // check if node is initially active or, if initialization with nodal, gap,
      // the gap is smaller than the prescribed value
      if (cnode->is_init_active() or (initcontactbygap and cnode->data().getg() < initcontactval))
      {
        cnode->active() = true;
        mynodegids.push_back(cnode->id());

        for (int j = 0; j < numdof; ++j) mydofgids.push_back(cnode->dofs()[j]);
      }

      // check if frictional node is initially in slip state
      if (friction_)
      {
        // do nothing: we always assume STICK at t=0
      }
    }

    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->active())
      {
        mynodegids.push_back(cnode->id());

        for (int j = 0; j < numdof; ++j) mydofgids.push_back(cnode->dofs()[j]);
      }
      else
      {
        mynodegidsInactive.push_back(cnode->id());

        for (int j = 0; j < numdof; ++j) mydofgidsInactive.push_back(cnode->dofs()[j]);
      }

      // check if frictional node is in slip state
      if (friction_)
      {
        if (dynamic_cast<FriNode*>(cnode)->fri_data().slip())
        {
          myslipnodegids.push_back(cnode->id());

          for (int j = 0; j < numdof; ++j) myslipdofgids.push_back(cnode->dofs()[j]);
        }
      }
    }
  }

  // create active node map and active dof map
  activenodes_ = Core::LinAlg::CreateMap(mynodegids, get_comm());
  activedofs_ = Core::LinAlg::CreateMap(mydofgids, get_comm());
  inactivenodes_ = Core::LinAlg::CreateMap(mynodegidsInactive, get_comm());
  inactivedofs_ = Core::LinAlg::CreateMap(mydofgidsInactive, get_comm());

  if (friction_)
  {
    // create slip node map and slip dof map
    slipnodes_ = Core::LinAlg::CreateMap(myslipnodegids, get_comm());
    slipdofs_ = Core::LinAlg::CreateMap(myslipdofgids, get_comm());
  }

  // split active dofs and slip dofs
  split_active_dofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs         popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::split_active_dofs()
{
  // get out of here if active set is empty
  if (activenodes_ == Teuchos::null)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements() == 0)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    return true;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNgids(activenodes_->NumMyElements());
  std::vector<int> myTgids((n_dim() - 1) * activenodes_->NumMyElements());

  // dimension check
  double dimcheck = (activedofs_->NumGlobalElements()) / (activenodes_->NumGlobalElements());
  if (dimcheck != n_dim()) FOUR_C_THROW("SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->num_dof();

    // add first dof to Nmap
    myNgids[countN] = cnode->dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j = 1; j < numdof; ++j)
    {
      myTgids[countT] = cnode->dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  get_comm().SumAll(&countN, &gcountN, 1);
  get_comm().SumAll(&countT, &gcountT, 1);

  // check global dimensions
  if ((gcountN + gcountT) != activedofs_->NumGlobalElements())
    FOUR_C_THROW("SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = Teuchos::rcp(new Epetra_Map(gcountN, countN, myNgids.data(), 0, get_comm()));
  activet_ = Teuchos::rcp(new Epetra_Map(gcountT, countT, myTgids.data(), 0, get_comm()));

  // *******************************************************************
  // FRICTION - EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if there is no friction
  if (friction_ == false) return true;

  // get out of here if slip set is empty
  if (slipnodes_ == Teuchos::null)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    return true;
  }

  if (slipnodes_->NumGlobalElements() == 0)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
    return true;
  }

  // define local variables
  int countslipT = 0;
  std::vector<int> myslipTgids((n_dim() - 1) * slipnodes_->NumMyElements());

  // dimension check
  dimcheck = (slipdofs_->NumGlobalElements()) / (slipnodes_->NumGlobalElements());
  if (dimcheck != n_dim()) FOUR_C_THROW("SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
  {
    int gid = slipnodes_->GID(i);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->num_dof();

    // add dofs to slipTmap
    for (int j = 1; j < numdof; ++j)
    {
      myslipTgids[countslipT] = cnode->dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  get_comm().SumAll(&countslipT, &gcountslipT, 1);

  // create Tslipmap objects
  slipt_ = Teuchos::rcp(new Epetra_Map(gcountslipT, countslipT, myslipTgids.data(), 0, get_comm()));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::get_force_of_node(Core::LinAlg::Matrix<3, 1>& nodal_force,
    const Epetra_Vector& force, const Core::Nodes::Node& node) const
{
  std::vector<int> dofs;
  idiscret_->dof(&node, dofs);

  // reset nodal force vector
  std::fill(nodal_force.data(), nodal_force.data() + 3, 0.0);
  const double* f_vals = force.Values();

  if (dofs.size() > 3) FOUR_C_THROW("The interface node seems to have more than 3 DOFs!");

  for (unsigned i = 0; i < dofs.size(); ++i)
  {
    const int dof = dofs[i];
    const int f_lid = force.Map().LID(dof);
    if (f_lid == -1)
      FOUR_C_THROW("Couldn't find the nodal DOF %d in the global force vector!", dof);

    nodal_force(i, 0) = f_vals[f_lid];
  }
}

/*----------------------------------------------------------------------*
 |  Calculate angular interface moments                  hiermeier 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::eval_resultant_moment(const Epetra_Vector& fs, const Epetra_Vector& fm,
    Core::LinAlg::SerialDenseMatrix* conservation_data_ptr) const
{
  double lresMoSl[3] = {0.0, 0.0, 0.0};  // local slave moment
  double gresMoSl[3] = {0.0, 0.0, 0.0};  // global slave moment
  double lresMoMa[3] = {0.0, 0.0, 0.0};  // local master momemnt
  double gresMoMa[3] = {0.0, 0.0, 0.0};  // global master moment
  double gbalMo[3] = {0.0, 0.0, 0.0};    // global moment balance

  double lresFSl[3] = {0.0, 0.0, 0.0};  // local slave force
  double gresFSl[3] = {0.0, 0.0, 0.0};  // global slave force
  double lresFMa[3] = {0.0, 0.0, 0.0};  // local master force
  double gresFMa[3] = {0.0, 0.0, 0.0};  // global master force
  double gbalF[3] = {0.0, 0.0, 0.0};    // global force balance

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  Core::LinAlg::Matrix<3, 1> nforce(true);
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Node* snode = dynamic_cast<Node*>(idiscret_->g_node(gid));

    get_force_of_node(nforce, fs, *snode);

    if (n_dim() == 2)
      lresMoSl[2] += snode->xspatial()[0] * nforce(1, 0) - snode->xspatial()[1] * nforce(0, 0);
    else
    {
      lresMoSl[0] += snode->xspatial()[1] * nforce(2, 0) - snode->xspatial()[2] * nforce(1, 0);
      lresMoSl[1] += snode->xspatial()[2] * nforce(0, 0) - snode->xspatial()[0] * nforce(2, 0);
      lresMoSl[2] += snode->xspatial()[0] * nforce(1, 0) - snode->xspatial()[1] * nforce(0, 0);
    }
    for (unsigned k = 0; k < 3; ++k) lresFSl[k] += nforce(k, 0);
  }
  get_comm().SumAll(lresMoSl, gresMoSl, 3);
  get_comm().SumAll(lresFSl, gresFSl, 3);

  // loop over proc's master nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    Node* mnode = dynamic_cast<Node*>(idiscret_->g_node(gid));

    get_force_of_node(nforce, fm, *mnode);

    if (n_dim() == 2)
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1, 0) - mnode->xspatial()[1] * nforce(0, 0);
    else
    {
      lresMoMa[0] += mnode->xspatial()[1] * nforce(2, 0) - mnode->xspatial()[2] * nforce(1, 0);
      lresMoMa[1] += mnode->xspatial()[2] * nforce(0, 0) - mnode->xspatial()[0] * nforce(2, 0);
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1, 0) - mnode->xspatial()[1] * nforce(0, 0);
    }
    for (unsigned k = 0; k < 3; ++k) lresFMa[k] += nforce(k, 0);
  }
  get_comm().SumAll(lresMoMa, gresMoMa, 3);
  get_comm().SumAll(lresFMa, gresFMa, 3);

  for (int d = 0; d < 3; ++d)
  {
    gbalMo[d] = gresMoSl[d] + gresMoMa[d];
    gbalF[d] = gresFSl[d] + gresFMa[d];
  }
  if (get_comm().MyPID() == 0)
  {
    std::cout << "SLAVE:   "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[2]
              << "]" << std::endl;
    std::cout << "Master:  "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[2]
              << "]" << std::endl;
    std::cout << "Balance: "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[2]
              << "]" << std::endl;
  }

  if (conservation_data_ptr)
  {
    Core::LinAlg::SerialDenseMatrix& conservation_data = *conservation_data_ptr;
    conservation_data.putScalar(0.0);
    if (conservation_data.numRows() < 18) FOUR_C_THROW("conservation_data length is too short!");

    std::copy(gresFSl, gresFSl + 3, conservation_data.values());
    std::copy(gresFMa, gresFMa + 3, conservation_data.values() + 3);
    std::copy(gbalF, gbalF + 3, conservation_data.values() + 6);

    std::copy(gresMoSl, gresMoSl + 3, conservation_data.values() + 9);
    std::copy(gresMoMa, gresMoMa + 3, conservation_data.values() + 12);
    std::copy(gbalMo, gbalMo + 3, conservation_data.values() + 15);
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::Interface& CONTACT::Interface::get_ma_sharing_ref_interface() const
{
  return dynamic_cast<const Interface&>(interface_data_->get_ma_sharing_ref_interface());
}


/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::store_to_old(Mortar::StrategyBase::QuantityType type)
{
  // loop over all slave row nodes on the current interface
  for (int j = 0; j < slave_col_nodes()->NumMyElements(); ++j)
  {
    int gid = slave_col_nodes()->GID(j);
    Core::Nodes::Node* node = idiscret_->g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

    switch (type)
    {
      case Mortar::StrategyBase::dm:
      {
        // store D and M entries
        dynamic_cast<FriNode*>(node)->store_dm_old();
        break;
      }
      case Mortar::StrategyBase::pentrac:
      {
        // store penalty tractions to old ones
        dynamic_cast<FriNode*>(node)->store_trac_old();
        break;
      }
      case Mortar::StrategyBase::n_old:
      {
        dynamic_cast<Node*>(node)->store_old_normal();
        break;
      }
      case Mortar::StrategyBase::activeold:
      {
        dynamic_cast<Node*>(node)->data().active_old() = dynamic_cast<Node*>(node)->active();
        break;
      }
      default:
        FOUR_C_THROW("StoreDMToNodes: Unknown state std::string variable!");
        break;
    }  // switch
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::update_self_contact_lag_mult_set(
    const Epetra_Map& gref_lmmap, const Epetra_Map& gref_smmap)
{
  if (gref_lmmap.NumMyElements() != gref_smmap.NumMyElements()) FOUR_C_THROW("Size mismatch!");

  const int num_sgids = sdofrowmap_->NumMyElements();
  const int* sgids = sdofrowmap_->MyGlobalElements();
  const int* ref_lmgids = gref_lmmap.MyGlobalElements();

  std::vector<int> lmdofs;
  lmdofs.reserve(num_sgids);

  for (int i = 0; i < num_sgids; ++i)
  {
    const int sgid = sgids[i];
    const int ref_lid = gref_smmap.LID(sgid);
    if (ref_lid == -1)
      FOUR_C_THROW(
          "Couldn't find the current slave gid #%d in the reference self "
          "contact slave master map.",
          sgid);
    lmdofs.push_back(ref_lmgids[ref_lid]);
  }

  lmdofmap_ = Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(lmdofs.size()), lmdofs.data(), 0, get_comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::set_node_initially_active(CONTACT::Node& cnode) const
{
  static const bool init_contact_by_gap =
      Core::UTILS::IntegralValue<int>(interface_params(), "INITCONTACTBYGAP");

  const bool node_init_active = cnode.is_init_active();

  // Either init contact by definition or by gap
  if (node_init_active and init_contact_by_gap)
    FOUR_C_THROW("Init contact either by definition in condition or by gap!");
  else if (node_init_active)
    cnode.active() = true;
  else if (init_contact_by_gap)
    set_node_initially_active_by_gap(cnode);

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (node_init_active)
    std::cout << "Node #" << std::setw(5) << cnode.id()
              << " is set initially active via the condition line.\n";
  else if (init_contact_by_gap)
    std::cout << "Node #" << std::setw(5) << cnode.id() << " is set initially active by gap.\n";
#endif

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::set_node_initially_active_by_gap(CONTACT::Node& cnode) const
{
  static const double initcontactval = interface_params().get<double>("INITCONTACTGAPVALUE");

  if (cnode.data().getg() < initcontactval) cnode.active() = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::set_condition_specific_parameters()
{
  if (interface_params().isSublist("ContactS2ICoupling"))
  {
    // read interface parameters and set them to the scatra boundary parameter class
    auto& s2icouplinglist = interface_params().sublist("ContactS2ICoupling", true);
    Discret::ELEMENTS::ScaTraEleParameterBoundary::instance("scatra")->set_parameters(
        s2icouplinglist);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::postprocess_quantities(const Teuchos::ParameterList& outputParams)
{
  using Teuchos::RCP;

  // Check if the given parameter list contains all required data to be written to output
  {
    // Vector with names of all required parameter entries
    std::vector<std::string> requiredEntries;
    requiredEntries.push_back("step");
    requiredEntries.push_back("time");
    requiredEntries.push_back("displacement");
    requiredEntries.push_back("interface traction");
    requiredEntries.push_back("slave forces");
    requiredEntries.push_back("master forces");

    check_output_list(outputParams, requiredEntries);
  }

  // Get the discretization writer and get ready for writing
  RCP<Core::IO::DiscretizationWriter> writer = idiscret_->writer();

  // Get output for this time step started
  {
    const int step = outputParams.get<int>("step");
    const double time = outputParams.get<double>("time");

    writer->clear_map_cache();
    writer->write_mesh(step, time);
    writer->new_step(step, time);
  }

  /* Write interface displacement
   *
   * The interface displacement has been handed in via the parameter list outParams.
   * Grab it from there, then use Core::LinAlg::Export() to extract the interface
   * portion from the global displacement vector. Finally, write the interface
   * portion using this interfaces' discretization writer.
   */
  {
    // Get full displacement vector and extract interface displacement
    RCP<const Epetra_Vector> disp = outputParams.get<RCP<const Epetra_Vector>>("displacement");
    RCP<Epetra_Vector> iDisp = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*disp, *iDisp);

    // Write the interface displacement field
    writer->write_vector("displacement", iDisp, Core::IO::VectorType::dofvector);
  }

  // Write Lagrange multiplier field
  {
    // Get full Lagrange multiplier vector and extract values of this interface
    RCP<const Epetra_Vector> lagMult =
        outputParams.get<RCP<const Epetra_Vector>>("interface traction");
    RCP<Epetra_Vector> iLagMult = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*lagMult, *iLagMult);

    // Write this interface's Lagrange multiplier field
    writer->write_vector("interfacetraction", iLagMult, Core::IO::VectorType::dofvector);
  }

  // Write normal contact stress
  {
    // Get values from parameter list and export to interface dof_row_map
    RCP<const Epetra_Vector> normalStresses =
        outputParams.get<RCP<const Epetra_Vector>>("norcontactstress");
    RCP<Epetra_Vector> iNormalStresses = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*normalStresses, *iNormalStresses);

    // Write this interface's normal contact stress field
    writer->write_vector("norcontactstress", iNormalStresses, Core::IO::VectorType::dofvector);
  }

  // Write tangential contact stress
  {
    // Get values from parameter list and export to interface dof_row_map
    RCP<const Epetra_Vector> tangentialStresses =
        outputParams.get<RCP<const Epetra_Vector>>("tancontactstress");
    RCP<Epetra_Vector> iTangentialStresses = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*tangentialStresses, *iTangentialStresses);

    // Write this interface's normal contact stress field
    writer->write_vector("tancontactstress", iTangentialStresses, Core::IO::VectorType::dofvector);
  }

  // Write nodal forces of slave side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> slaveforces =
        outputParams.get<RCP<const Epetra_Vector>>("slave forces");
    RCP<Epetra_Vector> forces = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*slaveforces, *forces);

    // Write to output
    writer->write_vector("slaveforces", forces, Core::IO::VectorType::dofvector);
  }

  // Write nodal forces of master side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> masterforces =
        outputParams.get<RCP<const Epetra_Vector>>("master forces");
    RCP<Epetra_Vector> forces = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*masterforces, *forces);

    // Write to output
    writer->write_vector("masterforces", forces, Core::IO::VectorType::dofvector);
  }


  // Nodes: node-based vector with '0' at slave nodes and '1' at master nodes
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*mnoderowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> nodeRowMap = Core::LinAlg::MergeMap(snoderowmap_, mnoderowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = Core::LinAlg::CreateVector(*nodeRowMap, true);
    Core::LinAlg::Export(*masterVec, *masterSlaveVec);

    writer->write_vector("slavemasternodes", masterSlaveVec, Core::IO::VectorType::nodevector);
  }

  // Write active set
  {
    // evaluate active set and slip set
    RCP<Epetra_Vector> activeset = Teuchos::rcp(new Epetra_Vector(*activenodes_));
    activeset->PutScalar(1.0);

    if (is_friction())
    {
      RCP<Epetra_Vector> slipset = Teuchos::rcp(new Epetra_Vector(*slipnodes_));
      slipset->PutScalar(1.0);
      RCP<Epetra_Vector> slipsetexp = Teuchos::rcp(new Epetra_Vector(*activenodes_));
      Core::LinAlg::Export(*slipset, *slipsetexp);
      activeset->Update(1.0, *slipsetexp, 1.0);
    }

    // export to interface node row map
    RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*(idiscret_->node_row_map())));
    Core::LinAlg::Export(*activeset, *activesetexp);

    writer->write_vector("activeset", activesetexp, Core::IO::VectorType::nodevector);
  }

  // Elements: element-based vector with '0' at slave elements and '1' at master elements
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*melerowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> eleRowMap = Core::LinAlg::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = Core::LinAlg::CreateVector(*eleRowMap, true);
    Core::LinAlg::Export(*masterVec, *masterSlaveVec);

    writer->write_vector(
        "slavemasterelements", masterSlaveVec, Core::IO::VectorType::elementvector);
  }

  // Write element owners
  {
    RCP<const Epetra_Map> eleRowMap = Core::LinAlg::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> owner = Core::LinAlg::CreateVector(*eleRowMap);

    for (int i = 0; i < idiscret_->element_row_map()->NumMyElements(); ++i)
      (*owner)[i] = idiscret_->l_row_element(i)->owner();

    writer->write_vector("Owner", owner, Core::IO::VectorType::elementvector);
  }
}

FOUR_C_NAMESPACE_CLOSE
