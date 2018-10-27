/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_adaptive_cn.cpp

\brief class which contains different strategies to adapt the regularization
       parameter for special scenarios

\maintainer Michael Hiermeier

\date Mar 6, 2018

\level 3
*/
/*----------------------------------------------------------------------------*/


#include "contact_aug_adaptive_cn.H"
#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/epetra_utils.H"
#include "../drt_lib/drt_discret_interface.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/connectivity.H"
#include "../drt_lib/drt_utils_parallel.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::AdaptiveCn::AdaptiveCn(
    Strategy& strat, const plain_interface_set& interfaces, DataContainer& data)
    : strategy_(strat), interfaces_(interfaces), data_(data)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::AdaptiveCn::Execute(
    const AdaptiveCnType adaptive_cn_type, CONTACT::ParamsInterface& cparams, std::ostream* os)
{
  type_ = adaptive_cn_type;

  switch (type_)
  {
    case AdaptiveCnType::init_nbc:
    {
      InitNbc initnbc(*this);
      initnbc(cparams, os);
      break;
    }
    default:
    {
      dserror("Unknown AdaptiveCnType enumerator: %d", adaptive_cn_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::AdaptiveCn::InitNbc::operator()(
    CONTACT::ParamsInterface& cparams, std::ostream* os)
{
  if (not cparams.IsPredictorState()) return;

  dt_ = Teuchos::Time::wallTime();

  Teuchos::RCP<Epetra_Vector> unit_gap_force = UnitGapForce();

  const STR::MODELEVALUATOR::Contact& cmodel = CModel(cparams);
  Teuchos::RCP<const Epetra_Vector> fextincr = cmodel.GetFextIncr();

  nbc_body_set neumann_bodies;
  SplitStructureIntoDistinctNeumannBodies(neumann_bodies, *fextincr, cmodel.StrDiscret());

  const double unit_gap_force_scale = parent_.data_.SaData().GetUnitGapForceScale();
  cn_init_ = 0.0;
  for (Teuchos::RCP<Body>& body : neumann_bodies)
  {
    body->addContactInterfaceMaps(parent_.interfaces_, parent_.strategy_.Comm());
    body->identifyAndVerifyNumDofsPerNode(cmodel.StrDiscret(), parent_.interfaces_);
    body->computeResultants(cmodel.StrDiscret(), parent_.interfaces_, *fextincr, *unit_gap_force);

    const double cn_tmp = body->leastSquaresForceBalance(unit_gap_force_scale);

    if (cn_init_ < cn_tmp) cn_init_ = cn_tmp;
  }

  double current_cn_max = 0.0;
  parent_.data_.Cn().MaxValue(&current_cn_max);
  if (cn_init_ > current_cn_max)
  {
    parent_.data_.Cn().PutScalar(cn_init_);
  }

  dt_ = Teuchos::Time::wallTime() - dt_;
  PrintUpdate(os);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::AdaptiveCn::InitNbc::PrintUpdate(std::ostream* os_ptr) const
{
  if (not os_ptr) return;

  std::ostream& os = *os_ptr;
  os << "InitNbc has been executed in ... " << std::setw(10) << std::scientific
     << std::setprecision(2) << dt_ << "sec.\n";
  os << "The initial regularization parameter estimate is " << cn_init_ << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::AdaptiveCn::InitNbc::UnitGapForce() const
{
  const LINALG::SparseMatrix& bmatrix = parent_.data_.BMatrix();
  Epetra_Vector unit_interface_gap = Epetra_Vector(bmatrix.RangeMap(), true);

  for (auto& interface_ptr : parent_.interfaces_)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<const Interface&>(*interface_ptr);

    Epetra_Vector unit_igap(*interface.SlaveRowNDofs(), true);
    interface.AssembleActiveUnitGap(unit_igap);

    LINALG::AssembleMyVector(1.0, unit_interface_gap, 1.0, unit_igap);
  }

  Teuchos::RCP<Epetra_Vector> unit_gap_force = Teuchos::rcp(new Epetra_Vector(bmatrix.DomainMap()));
  CATCH_EPETRA_ERROR(bmatrix.Multiply(true, unit_interface_gap, *unit_gap_force));

  return unit_gap_force;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Body::leastSquaresForceBalance(const double scale) const
{
  const double cf = rfext_.Dot(rcontact_);
  const double cc = rcontact_.Dot(rcontact_);

  if (cc == 0.0) dserror("The inner product of the contact forces is equal to zero!");

  return (cf / (scale * cc));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Body::identifyAndVerifyNumDofsPerNode(
    const DRT::DiscretizationInterface& str_discret, const plain_interface_set& interfaces)
{
  num_dofs_per_node_ = numDofsPerNode(str_discret, 0);

  // sanity check
  for (int ilid : interface_lids_)
  {
    if (interfaces[ilid].is_null()) dserror("There is no interface with LID %d", ilid);

    const CoInterface& interface = *interfaces[ilid];
    numDofsPerNode(interface.Discret(), 0, &num_dofs_per_node_);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Body::computeResultants(const DRT::DiscretizationInterface& str_discret,
    const plain_interface_set& interfaces, const Epetra_Vector& fext,
    const Epetra_Vector& unit_gap_force)
{
  // compute Neumann force
  {
    const discret_set discrets(1, &str_discret);
    resultantVector(fext, *node_map_, discrets, rfext_);
  }
  rfext_.Print(std::cout);

  // compute contact forces
  {
    discret_set discrets;
    discrets.reserve(interface_lids_.size());
    for (int ilid : interface_lids_) discrets.push_back(&interfaces[ilid]->Discret());

    resultantVector(unit_gap_force, *contact_inodes_, discrets, rcontact_);
  }
  rcontact_.Print(std::cout);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::AdaptiveCn::InitNbc::SplitStructureIntoDistinctNeumannBodies(
    nbc_body_set& neumann_bodies, const Epetra_Vector& fext,
    const DRT::DiscretizationInterface& full_discret) const
{
  std::vector<Teuchos::RCP<Epetra_Map>> body_node_maps;
  DRT::Connectivity::splitDiscretizationIntoDistinctNodeMaps(
      full_discret, DRT::Connectivity::MapFormat::row, body_node_maps);

  std::vector<Teuchos::RCP<Epetra_Map>> body_dof_maps;
  LINALG::ComputeDofMapsFromNodeMaps(0, body_node_maps, full_discret, body_dof_maps);

  const size_t num_bodies = body_node_maps.size();

  neumann_bodies.reserve(num_bodies);

  for (unsigned b = 0; b < num_bodies; ++b)
  {
    // extract neumann force
    Epetra_Vector fext_b(*body_dof_maps[b]);

    LINALG::ExtractMyVector(fext, fext_b);

    // is there any neumann force acting on this body?
    double fext_magnitude = 0.0;
    fext_b.Norm2(&fext_magnitude);

    // if not, go to the next body
    if (fext_magnitude < std::numeric_limits<double>::epsilon()) continue;

    neumann_bodies.push_back(Teuchos::rcp(new Body()));
    Body& body = *neumann_bodies.back();

    body.addNodeMapPtr(body_node_maps[b]);
    body.addDofMapPtr(body_dof_maps[b]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Body::resultantVector(const Epetra_Vector& vec, const Epetra_Map& node_map,
    const discret_set& discrets, LINALG::SerialDenseVector& resultant) const
{
  resultant.Resize(num_dofs_per_node_);
  resultant.Zero();

  LINALG::SerialDenseVector lresultant(num_dofs_per_node_, true);

  const int num_nodes = node_map.NumMyElements();
  const int* ngids = node_map.MyGlobalElements();

  for (int nlid = 0; nlid < num_nodes; ++nlid)
  {
    int ilid = -1;
    const DRT::Node* node = findRowNode(ngids[nlid], discrets, ilid);
    for (int d = 0; d < num_dofs_per_node_; ++d)
    {
      const int dof_gid = discrets[ilid]->Dof(0, node, d);
      const int dof_lid = vec.Map().LID(dof_gid);
      if (dof_lid == -1) dserror("Couldn't find dof GID #%d.", dof_gid);

      lresultant[d] += vec.Values()[dof_lid];
    }
  }

  const Epetra_Comm& comm = discrets[0]->Comm();
  comm.SumAll(lresultant.A(), resultant.A(), num_dofs_per_node_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Node* CONTACT::AUG::Body::findRowNode(
    int ngid, const discret_set& discrets, int& ilid) const
{
  ilid = 0;
  for (const DRT::DiscretizationInterface* discret : discrets)
  {
    if (discret->NodeRowMap()->MyGID(ngid)) return discret->gNode(ngid);
    ++ilid;
  }
  dserror("Couldn't find the GID in any of the passed discretizations!");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::Body::numDofsPerNode(
    const DRT::DiscretizationInterface& discret, int dofset_id, const int* num_dofs_per_node) const
{
  const Epetra_Map* row_nmap = discret.NodeRowMap();
  const int num_rnodes = row_nmap->NumMyElements();

  int first_num_dofs = (num_dofs_per_node ? *num_dofs_per_node : -1);
  DRT::Node* first_node = NULL;
  if (num_rnodes > 0)
  {
    first_node = discret.lRowNode(0);
    if (first_num_dofs < 0)
      first_num_dofs = discret.NumDof(dofset_id, first_node);
    else if (first_num_dofs != discret.NumDof(dofset_id, first_node))
      dserror("Different number of DoFs per node among different discretizations.");
  }

  for (int nlid = 1; nlid < num_rnodes; ++nlid)
  {
    DRT::Node* node = discret.lRowNode(nlid);
    int num_dofs = discret.NumDof(dofset_id, node);

    if (num_dofs != discret.NumDof(dofset_id, node))
      dserror(
          "Different number of dofs in dof-set #%d detected for node "
          "GID #%d [ndofs = %d] and node GID #%d [ndofs = %d].",
          dofset_id, first_node->Id(), first_num_dofs, node->Id(), num_dofs);
  }

  const Epetra_Comm& comm = discret.Comm();
  std::vector<int> all_num_dofs(comm.NumProc(), -1);
  comm.GatherAll(&first_num_dofs, all_num_dofs.data(), 1);

  for (int p = 0; p < comm.NumProc(); ++p)
  {
    if (all_num_dofs[p] >= 0 and first_num_dofs != all_num_dofs[p])
      dserror("Different number of dofs on proc #d and proc #d.", comm.MyPID(), p);

    if (first_num_dofs < 0 and all_num_dofs[p] > 0) first_num_dofs = all_num_dofs[p];
  }

  return first_num_dofs;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Contact& CONTACT::AUG::AdaptiveCn::InitNbc::CModel(
    CONTACT::ParamsInterface& cparams) const
{
  return dynamic_cast<const STR::MODELEVALUATOR::Contact&>(cparams.GetModelEvaluator());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CONTACT::AUG::Body::extractContactInterfaceMap(
    const Epetra_Map& body_map, const Epetra_Map& map, const Epetra_Map* pmap) const
{
  const Epetra_Map* map_ptr = NULL;
  if (pmap)
    map_ptr = pmap;
  else
    map_ptr = &map;

  Teuchos::RCP<Epetra_Map> overlapping_map = LINALG::ExtractMyOverlappingSubMap(*map_ptr, body_map);

  if (overlapping_map->NumGlobalElements() > 0)
  {
    if (pmap)
      overlapping_map = DRT::UTILS::RedistributeInAccordanceWithReference(map, *overlapping_map);

    return overlapping_map;
  }

  return Teuchos::rcp(new Epetra_Map(0, 0, map.Comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Body::addContactInterfaceMaps(
    const plain_interface_set& interfaces, const Epetra_Comm& comm)
{
  contact_idofs_ = Teuchos::rcp(new Epetra_Map(0, 0, comm));
  contact_inodes_ = Teuchos::rcp(new Epetra_Map(0, 0, comm));
  interface_lids_.clear();

  int i_lid = 0;
  for (const Teuchos::RCP<const CoInterface>& iptr : interfaces)
  {
    const CoInterface& interface = *iptr;

    // --- slave maps ---------------------------------------------------------
    {
      Teuchos::RCP<Epetra_Map> overlapping_dofmap = extractContactInterfaceMap(
          *dof_map_, *interface.SlaveRowDofs(), interface.PSlaveRowDofs().get());
      if (overlapping_dofmap->NumGlobalElements() > 0)
        contact_idofs_ = LINALG::MergeMap(*overlapping_dofmap, *contact_idofs_);

      Teuchos::RCP<Epetra_Map> overlapping_nodemap = extractContactInterfaceMap(
          *node_map_, *interface.SlaveRowNodes(), interface.PSlaveRowNodes().get());
      if (overlapping_nodemap->NumGlobalElements() > 0)
      {
        contact_inodes_ = LINALG::MergeMap(*overlapping_nodemap, *contact_inodes_);
        interface_lids_.insert(i_lid);
      }
    }

    // --- master maps ---------------------------------------------------------
    {
      Teuchos::RCP<Epetra_Map> overlapping_dofmap = extractContactInterfaceMap(
          *dof_map_, *interface.MasterRowDofs(), interface.PMasterRowDofs().get());
      if (overlapping_dofmap->NumGlobalElements() > 0)
        contact_idofs_ = LINALG::MergeMap(*overlapping_dofmap, *contact_idofs_);

      Teuchos::RCP<Epetra_Map> overlapping_nodemap = extractContactInterfaceMap(
          *node_map_, *interface.MasterRowNodes(), interface.PMasterRowNodes().get());
      if (overlapping_nodemap->NumGlobalElements() > 0)
      {
        contact_inodes_ = LINALG::MergeMap(*overlapping_nodemap, *contact_inodes_);
        interface_lids_.insert(i_lid);
      }
    }

    ++i_lid;
  }
}
