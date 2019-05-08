/*----------------------------------------------------------------------------*/
/**
\file xcontact_levelset_algorithm.cpp

\brief level-set algorithm for the extended contact formulation

\maintainer Matthias Mayr

\date Nov 22, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xcontact_levelset_algorithm.H"
#include "xcontact_levelset_intersection.H"
#include "xcontact_levelset_reinit_elliptic.H"

#include "../drt_cut/cut_point.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::Algorithm::Algorithm(const Teuchos::RCP<DRT::Discretization>& dis,
    const Teuchos::RCP<LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& level_set_params,
    const Teuchos::RCP<Teuchos::ParameterList>& scatradynparams,
    const Teuchos::RCP<Teuchos::ParameterList>& extraparams,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
    : SCATRA::ScaTraTimIntImpl(dis, solver, scatradynparams, extraparams, output),
      SCATRA::LevelSetAlgorithm(
          dis, solver, level_set_params, scatradynparams, extraparams, output),
      ls_intersect_(Teuchos::null),
      reinit_strategy_(Teuchos::null),
      active_row_elements_(Teuchos::null),
      active_row_nodes_(Teuchos::null),
      active_col_elements_(Teuchos::null),
      active_dof_map_extractor_(Teuchos::null),
      l2_sys_mat_diagonal_(Teuchos::null)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::Setup()
{
  SCATRA::LevelSetAlgorithm::Setup();

  // --- set the level-set intersection class
  std::vector<GEO::CUT::Point::PointPosition> desired_pos;
  desired_pos.reserve(2);
  desired_pos.push_back(GEO::CUT::Point::inside);
  desired_pos.push_back(GEO::CUT::Point::undecided);

  ls_intersect_ = Teuchos::rcp(new XCONTACT::LEVELSET::Intersection());
  ls_intersect_->SetDesiredPositions(desired_pos);

  // --- set the reinitialization class
  switch (reinitaction_)
  {
    case INPAR::SCATRA::reinitaction_ellipticeq:
    {
      reinit_strategy_ = Teuchos::rcp(new XCONTACT::LEVELSET::REINIT::Elliptic());
      break;
    }
    default:
      dserror("Unsupported reinitialization action!");
      exit(EXIT_FAILURE);
  }

  // initialize and setup the reinitialization object
  ReinitStrategy().Init(this);
  ReinitStrategy().Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::GetInitialVolumeOfMinusDomain(
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    const Teuchos::RCP<const DRT::Discretization>& scatradis, double& volumedomainminus) const
{
  volumedomainminus = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::Reinitialization()
{
  // return if we have nothing to do here
  if (reinitaction_ == INPAR::SCATRA::reinitaction_none or step_ % ReinitInterval() != 0) return;

  CaptureInterface(zero_iso_line_);

  // compute the reinitialization
  ReinitStrategy().Compute(*phinp_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XCONTACT::LEVELSET::Algorithm::ComputeNodalGradientL2Projection(
    const Epetra_Vector& x, const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_system)
{
  // set element parameters
  Teuchos::ParameterList eleparams;
  eleparams.set<bool>("solve reinit eq", true);
  eleparams.set<int>("action", SCATRA::recon_gradients_at_nodes);

  const Teuchos::ParameterList& scatradyn = problem_->ScalarTransportDynamicParams();
  const int l2_proj_solver = scatradyn.get<int>("L2_PROJ_LINEAR_SOLVER");

  const int dim = problem_->NDim();

  // set given state for element evaluation
  discret_->ClearState();
  discret_->SetState("phinp", Teuchos::rcpFromRef(x));

  l2_sys_mat_diagonal_ = Teuchos::rcp(new Epetra_Vector(ActiveRowNodeMap(), true));

  Teuchos::RCP<Epetra_MultiVector> nodalgradient = DRT::UTILS::ComputeNodalL2Projection(*discret_,
      ActiveRowNodeMap(), ActiveColEleMap(), "phinp", dim, eleparams, l2_proj_solver,
      l2_proj_system, NULL, NULL, l2_sys_mat_diagonal_.get());

  // switch to the dof map ( only meaningful, if we have one scalar per node )
  int err = l2_sys_mat_diagonal_->ReplaceMap(x.Map());
  if (err) dserror("ReplaceMap failed! ( err = %d )", err);

  return nodalgradient;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::CaptureInterface(
    std::map<int, GEO::BoundaryIntCellPtrs>& zero_iso_line)
{
  CheckIsSetup();

  zero_iso_line.clear();
  ls_intersect_->CaptureZeroLevelSet(*phinp_, *discret_, zero_iso_line);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::OutputOfLevelSetSpecificValues()
{
  // do nothing
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::CreateActiveMaps(const Epetra_Vector& phinp)
{
  CreateActiveMaps(*Discretization(), phinp, active_row_elements_, active_row_nodes_,
      active_col_elements_, active_dof_map_extractor_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& XCONTACT::LEVELSET::Algorithm::ActiveRowDofMap() const
{
  return *(ActiveDofMapExtractor().Map(active));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::CreateActiveMaps(const DRT::Discretization& dis,
    const Epetra_Vector& phinp, Teuchos::RCP<Epetra_Map>& active_row_elements,
    Teuchos::RCP<Epetra_Map>& active_row_nodes, Teuchos::RCP<Epetra_Map>& active_col_elements,
    Teuchos::RCP<LINALG::MultiMapExtractor>& active_dof_map_extractor) const
{
  // export phi from row to column map
  Epetra_Vector phinp_col = Epetra_Vector(*dis.DofColMap());
  LINALG::Export(phinp, phinp_col);

  std::set<int> rdof_gids;
  std::set<int> rele_gids;
  std::set<int> rnode_gids;
  std::set<int> cele_gids;

  std::vector<int> dof_gids_per_ele;
  std::vector<int> node_gids_per_ele;

  const int num_ele = dis.ElementColMap()->NumMyElements();
  for (int elid = 0; elid < num_ele; ++elid)
  {
    DRT::Element* ele = dis.lColElement(elid);

    const int numnodes = ele->NumNode();

    node_gids_per_ele.resize(numnodes, -1);

    // expect one dof per node
    dof_gids_per_ele.clear();
    dof_gids_per_ele.reserve(numnodes);


    DRT::Node* const* nodes = ele->Nodes();
    for (int nlid = 0; nlid < numnodes; ++nlid)
    {
      node_gids_per_ele[nlid] = nodes[nlid]->Id();

      const int numdofpernode = dis.NumDof(0, nodes[nlid]);

      for (int d = 0; d < numdofpernode; ++d)
        dof_gids_per_ele.push_back(dis.Dof(0, nodes[nlid], d));
    }

    dsassert(dof_gids_per_ele.size() > 0, "There are no DoF's!");

    // check the level-set values of the current element
    const int dof_lid = phinp_col.Map().LID(dof_gids_per_ele[0]);
    if (dof_lid == -1)
      dserror("The GID %d was not found in the given state vector!", dof_gids_per_ele[0]);

    /* TRUE, if the scalar values are constant over the whole element.
     * We are going to skip these elements in that case. */
    bool is_constant_phi = true;

    double phi_0 = phinp_col[dof_lid];
    for (std::vector<int>::const_iterator cit = ++dof_gids_per_ele.begin();
         cit != dof_gids_per_ele.end(); ++cit)
    {
      const int dof_lid = phinp_col.Map().LID(*cit);

      if (dof_lid == -1) dserror("The GID %d was not found in the given state vector!", *cit);

      if (std::abs(phi_0 - phinp_col[dof_lid]) > 1.0e-12)
      {
        is_constant_phi = false;
        break;
      }
    }

    /* skip elements with constant gap values ( this would lead to a zero gradient,
     * and finally to an undefined behavior of the reinitialization procedure ) */
    if (is_constant_phi) continue;

    if (ele->Owner() == dis.Comm().MyPID())
      // store the row element GID's
      rele_gids.insert(ele->Id());

    if (not ele->HasOnlyGhostNodes(dis.Comm().MyPID()))
    {
      // store the row dof GID's
      for (std::vector<int>::const_iterator cit = dof_gids_per_ele.begin();
           cit != dof_gids_per_ele.end(); ++cit)
      {
        if (dis.DofRowMap()->LID(*cit) > -1) rdof_gids.insert(*cit);
      }

      // store the row node GID's
      for (std::vector<int>::const_iterator cit = node_gids_per_ele.begin();
           cit != node_gids_per_ele.end(); ++cit)
      {
        if (dis.NodeRowMap()->LID(*cit) > -1) rnode_gids.insert(*cit);
      }
    }

    // store the column element GID's
    cele_gids.insert(ele->Id());
  }

  std::vector<int> rele_gids_vec(rele_gids.begin(), rele_gids.end());
  CreateNewMap(rele_gids_vec, dis.Comm(), active_row_elements);

  std::vector<int> rnode_gids_vec(rnode_gids.begin(), rnode_gids.end());
  CreateNewMap(rnode_gids_vec, dis.Comm(), active_row_nodes);

  std::vector<int> cele_gids_vec(cele_gids.begin(), cele_gids.end());
  CreateNewMap(cele_gids_vec, dis.Comm(), active_col_elements);

  Teuchos::RCP<Epetra_Map> dof_map_active = Teuchos::null;
  std::vector<int> rdof_gids_vec(rdof_gids.begin(), rdof_gids.end());
  CreateNewMap(rdof_gids_vec, dis.Comm(), dof_map_active);

  std::vector<Teuchos::RCP<const Epetra_Map>> dof_maps(2, Teuchos::null);
  dof_maps[active] = dof_map_active.getConst();
  dof_maps[inactive] = LINALG::SplitMap(*dis.DofRowMap(), *dof_maps[active]);

  active_dof_map_extractor =
      Teuchos::rcp(new LINALG::MultiMapExtractor(*dis.DofRowMap(), dof_maps));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Algorithm::CreateNewMap(const std::vector<int>& my_entries,
    const Epetra_Comm& comm, Teuchos::RCP<Epetra_Map>& new_map) const
{
  const int* my_entries_ptr = NULL;
  if (my_entries.size() > 0) my_entries_ptr = &my_entries[0];

  new_map = Teuchos::null;
  new_map = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, static_cast<int>(my_entries.size()), my_entries_ptr, 0, comm));
}
