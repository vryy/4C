/*----------------------------------------------------------------------------*/
/**
\file xcontact_multi_discretization_wrapper.cpp

\brief wrapper of the discretizations involved in a xcontact simulation

\maintainer Matthias Mayr

\date Nov 10, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xcontact_multi_discretization_wrapper.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_mortar/mortar_node.H"

#include <Epetra_Export.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::MultiDiscretizationWrapper::MultiDiscretizationWrapper()
    : icontact_(0),
      iscatra_(0),
      slave_dof_row_maps_(0),
      slave_dof_col_maps_(0),
      slave_ndof_row_maps_(0),
      master_dof_row_maps_(0),
      master_dof_col_maps_(0),
      contact2scatra_exporter_(0),
      scatra2contact_exporter_(0),
      is_exporter_(false)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::AddContactIDiscret(const XDisPairedPtrVector& idiscrets)
{
  is_exporter_ = false;

  const int inum = idiscrets.size();
  icontact_.resize(inum);

  // resize paired vector sizes
  slave_dof_row_maps_.resize(inum);
  slave_ndof_row_maps_.resize(inum);
  slave_dof_col_maps_.resize(inum);
  master_dof_row_maps_.resize(inum);
  master_dof_col_maps_.resize(inum);

  XDisPairedPtrVector::const_iterator cit;
  for (cit = idiscrets.begin(); cit != idiscrets.end(); ++cit) icontact_[cit->first] = cit->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::AddScaTraIDiscret(
    const XFEM::FieldName& field, Teuchos::RCP<DRT::Discretization> idiscret)
{
  is_exporter_ = false;

  const unsigned inum = iscatra_.size() + 1;
  iscatra_.resize(inum);
  iscatra_[field] = idiscret;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const XCONTACT::MultiDiscretizationWrapper::XDisPairedPtrVector&
XCONTACT::MultiDiscretizationWrapper::GetContactIDiscret() const
{
  return icontact_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<DRT::Discretization>& XCONTACT::MultiDiscretizationWrapper::GetScaTraDiscretPtr(
    const enum XFEM::FieldName& field) const
{
  XDisPairedRCPVector::const_iterator cit = iscatra_.find(field);
  if (cit == iscatra_.end())
    dserror(
        "The given field name \"%s\" could not be found in the wrapped "
        "scatra paired vector!",
        XFEM::FieldName2String(field).c_str());

  return cit->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XCONTACT::MultiDiscretizationWrapper::FillComplete(bool assigndegreesoffreedom,
    bool initelements, bool doboundaryconditions, bool buildsystemmaps, bool setupmapextractor)
{
  int ret = XSTR::MultiDiscretizationWrapper::FillComplete(assigndegreesoffreedom, initelements,
      doboundaryconditions, buildsystemmaps, setupmapextractor);

  SetupContactScaTraDofExporter();

  return ret;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::BuildSlaveMasterMaps()
{
  XDisPairedPtrVector::const_iterator cit;

  std::vector<int> my_slave_row_dofs;
  std::vector<int> my_slave_row_ndofs;
  std::vector<int> my_slave_col_dofs;
  std::vector<int> my_master_row_dofs;
  std::vector<int> my_master_col_dofs;

  for (cit = icontact_.begin(); cit != icontact_.end(); ++cit)
  {
    my_slave_row_dofs.clear();
    my_slave_row_ndofs.clear();
    my_slave_col_dofs.clear();
    my_master_row_dofs.clear();
    my_master_col_dofs.clear();

    const DRT::Discretization& contact_discret = *cit->second;

    unsigned num_my_nodes = contact_discret.NumMyColNodes();
    DRT::Node* const* col_nodes = contact_discret.lColNodes();
    for (unsigned elid = 0; elid < num_my_nodes; ++elid)
    {
      MORTAR::MortarNode* col_mortar_node = dynamic_cast<MORTAR::MortarNode*>(col_nodes[elid]);

      if (not col_mortar_node) dserror("dynamic cast failed!");

      if (col_mortar_node->IsSlave())
      {
        // --- all dofs
        const std::vector<int> dofs = contact_discret.Dof(0, col_mortar_node);

        std::copy(dofs.begin(), dofs.end(), std::back_inserter(my_slave_col_dofs));
        if (col_mortar_node->Owner() == Comm().MyPID())
        {
          std::copy(dofs.begin(), dofs.end(), std::back_inserter(my_slave_row_dofs));

          // --- normal dofs
          my_slave_row_ndofs.push_back(dofs[0]);
        }
      }
      else
      {
        // --- all dofs
        const std::vector<int> dofs = contact_discret.Dof(0, col_mortar_node);

        std::copy(dofs.begin(), dofs.end(), std::back_inserter(my_master_col_dofs));
        if (col_mortar_node->Owner() == Comm().MyPID())
          std::copy(dofs.begin(), dofs.end(), std::back_inserter(my_master_row_dofs));
      }
    }

    // build slave dof row map
    CreateNewMap(my_slave_row_dofs, slave_dof_row_maps_[cit->first]);

    // build slave normal dof row map
    CreateNewMap(my_slave_row_ndofs, slave_ndof_row_maps_[cit->first]);

    // build slave dof column map
    CreateNewMap(my_slave_col_dofs, slave_dof_col_maps_[cit->first]);

    // build master dof row map
    CreateNewMap(my_master_row_dofs, master_dof_row_maps_[cit->first]);

    // build master dof column map
    CreateNewMap(my_master_col_dofs, master_dof_col_maps_[cit->first]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::SetupContactScaTraDofExporter()
{
  BuildSlaveMasterMaps();

  XDisPairedPtrVector::const_iterator cit;
  const unsigned inum = icontact_.size();

  // compare the paired vector sizes
  if (inum != iscatra_.size())
    dserror(
        "There was an unequal number of contact interface discretiztations and"
        " scaTra discretizations detected! (num icontact = %d, num iscaTra = %d)",
        icontact_.size(), iscatra_.size());

  // resize the pairedvectors
  contact2scatra_exporter_.resize(inum);
  scatra2contact_exporter_.resize(inum);

  for (cit = icontact_.begin(); cit != icontact_.end(); ++cit)
    SetupContactScaTraDofExporter(cit->first, *cit->second, *iscatra_[cit->first], false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::SetupContactScaTraDofExporter(
    const XFEM::FieldName& field, const DRT::Discretization& icontact,
    const DRT::Discretization& iscatra, const bool& resize)
{
  if (resize)
  {
    // resize if necessary
    if (contact2scatra_exporter_.size() != icontact_.size())
      contact2scatra_exporter_.resize(icontact_.size());
    // resize if necessary
    if (scatra2contact_exporter_.size() != iscatra_.size())
      scatra2contact_exporter_.resize(iscatra_.size());
  }

  /* first we get the node row maps and if they are equal, we skip
   * the exporter construction */
  const Epetra_Map* scatra_node_map = iscatra.NodeRowMap();
  const Epetra_Map* contact_node_map = icontact.NodeRowMap();

  // check if the ScaTra node row map is a subset of the contact
  // node row map ( slave side only )
  const int num_gids = scatra_node_map->NumMyElements();
  const int* gids = scatra_node_map->MyGlobalElements();

  bool is_subset = true;
  for (int i = 0; i < num_gids; ++i)
  {
    is_subset = contact_node_map->MyGID(gids[i]);
    if (not is_subset) break;
  }

  int global_is_subset = 0;
  int local_is_subset = static_cast<int>(is_subset);
  Comm().MinAll(&local_is_subset, &global_is_subset, 1);

  if (global_is_subset)
  {
    contact2scatra_exporter_[field] = Teuchos::null;
    scatra2contact_exporter_[field] = Teuchos::null;
  }
  else
  {
    const Epetra_Map* scatra_dof_map = iscatra.DofRowMap(0);
    const Epetra_Map* contact_slave_normal_dof_map = slave_ndof_row_maps_.at(field).get();
    contact2scatra_exporter_[field] =
        Teuchos::rcp(new Epetra_Export(*contact_slave_normal_dof_map, *scatra_dof_map));
    scatra2contact_exporter_[field] =
        Teuchos::rcp(new Epetra_Export(*scatra_dof_map, *contact_slave_normal_dof_map));
  }

  is_exporter_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XCONTACT::MultiDiscretizationWrapper::Contact2ScaTra(
    const Teuchos::RCP<Epetra_Vector>& contact_vec, const enum XFEM::FieldName& field,
    const bool& do_normalize) const
{
  return Contact2ScaTra(*contact_vec, field, do_normalize);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XCONTACT::MultiDiscretizationWrapper::Contact2ScaTra(
    const Epetra_Vector& contact_vec, const enum XFEM::FieldName& field,
    const bool& do_normalize) const
{
  Teuchos::RCP<Epetra_Vector> scatra_vec =
      Teuchos::rcp(new Epetra_Vector(*iscatra_.at(field)->DofRowMap(0)));

  Contact2ScaTra(contact_vec, *scatra_vec, field, do_normalize);

  return scatra_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XCONTACT::MultiDiscretizationWrapper::Contact2ScaTra(
    const Teuchos::RCP<Epetra_MultiVector>& contact_vec, const enum XFEM::FieldName& field,
    const bool& do_normalize) const
{
  return Contact2ScaTra(*contact_vec, field, do_normalize);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XCONTACT::MultiDiscretizationWrapper::Contact2ScaTra(
    const Epetra_MultiVector& contact_vec, const enum XFEM::FieldName& field,
    const bool& do_normalize) const
{
  const int num_vec = contact_vec.NumVectors();
  Teuchos::RCP<Epetra_MultiVector> scatra_vec =
      Teuchos::rcp(new Epetra_MultiVector(*iscatra_.at(field)->DofRowMap(0), num_vec));

  Contact2ScaTra(contact_vec, *scatra_vec, field, do_normalize);

  return scatra_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::Contact2ScaTra(const Epetra_MultiVector& contact_vec,
    Epetra_MultiVector& scatra_vec, const enum XFEM::FieldName& field,
    const bool& do_normalize) const
{
  if (not IsExporter()) dserror("Call SetupContactScaTraDofExporter() first!");

  Epetra_Export* c2s_exporter = contact2scatra_exporter_.at(field).get();
  if (c2s_exporter != NULL)
  {
    const int err = scatra_vec.Export(contact_vec, *c2s_exporter, Insert);
    if (err) dserror("Export to ScaTra distribution returned err = %d", err);
  }
  else
  {
    int err = scatra_vec.ReplaceMap(contact_vec.Map());
    if (err) dserror("Replacement of the ScaTra vector map returned err = %d", err);
    scatra_vec.Scale(1.0, contact_vec);
    scatra_vec.ReplaceMap(*iscatra_.at(field)->DofRowMap(0));
  }

  // normalize the extracted contact vector for the ScaTra use
  if (do_normalize) NormalizeInitialScaTraValues(scatra_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::NormalizeInitialScaTraValues(
    Epetra_MultiVector& scatra_vec) const
{
  std::vector<int> inactive_lids;
  std::vector<int> active_lids;

  double lmax_value = 0.0;
  // reciprocal global max value
  double gmax_value_rec = 0.0;

  for (int i = 0; i < scatra_vec.Map().NumMyElements(); ++i)
  {
    const double& val = scatra_vec.Values()[i];
    if (scatra_vec.Values()[i] == 1.0e12)
      inactive_lids.push_back(i);
    else
    {
      active_lids.push_back(i);
      if (lmax_value < std::abs(val))
      {
        lmax_value = std::abs(val);
      }
    }
  }

  scatra_vec.Comm().MaxAll(&lmax_value, &gmax_value_rec, 1);
  gmax_value_rec = 1.0 / gmax_value_rec;

  // normalize the active gap values
  for (std::vector<int>::const_iterator clid = active_lids.begin(); clid != active_lids.end();
       ++clid)
  {
    scatra_vec.Values()[*clid] *= gmax_value_rec;
  }

  // set the inactive values to 2.0
  for (std::vector<int>::const_iterator clid = inactive_lids.begin(); clid != inactive_lids.end();
       ++clid)
  {
    scatra_vec.Values()[*clid] = 2.0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XCONTACT::MultiDiscretizationWrapper::ScaTra2Contact(
    const Teuchos::RCP<Epetra_Vector>& scatra_vec, const enum XFEM::FieldName& field) const
{
  return ScaTra2Contact(*scatra_vec, field);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XCONTACT::MultiDiscretizationWrapper::ScaTra2Contact(
    const Epetra_Vector& scatra_vec, const enum XFEM::FieldName& field) const
{
  Teuchos::RCP<Epetra_Vector> contact_vec =
      Teuchos::rcp(new Epetra_Vector(*slave_ndof_row_maps_.at(field)));

  ScaTra2Contact(scatra_vec, *contact_vec, field);

  return contact_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XCONTACT::MultiDiscretizationWrapper::ScaTra2Contact(
    const Teuchos::RCP<Epetra_MultiVector>& scatra_vec, const enum XFEM::FieldName& field) const
{
  return ScaTra2Contact(*scatra_vec, field);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XCONTACT::MultiDiscretizationWrapper::ScaTra2Contact(
    const Epetra_MultiVector& scatra_vec, const enum XFEM::FieldName& field) const
{
  const int num_vec = scatra_vec.NumVectors();
  Teuchos::RCP<Epetra_MultiVector> contact_vec =
      Teuchos::rcp(new Epetra_MultiVector(*slave_ndof_row_maps_.at(field), num_vec));

  ScaTra2Contact(scatra_vec, *contact_vec, field);

  return contact_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::MultiDiscretizationWrapper::ScaTra2Contact(const Epetra_MultiVector& scatra_vec,
    Epetra_MultiVector& contact_vec, const enum XFEM::FieldName& field) const
{
  if (not IsExporter()) dserror("Call SetupContactScaTraDofExporter() first!");

  Epetra_Export* s2c_exporter = scatra2contact_exporter_.at(field).get();
  if (s2c_exporter != NULL)
  {
    const int err = contact_vec.Export(scatra_vec, *s2c_exporter, Insert);
    if (err) dserror("Export to Contact distribution returned err = %d", err);
  }
  else
  {
    int err = contact_vec.ReplaceMap(scatra_vec.Map());
    if (err) dserror("Replacement of the Contact vector map returned err = %d", err);
    contact_vec.Scale(1.0, scatra_vec);
    contact_vec.ReplaceMap(*slave_ndof_row_maps_.at(field));
  }
}
