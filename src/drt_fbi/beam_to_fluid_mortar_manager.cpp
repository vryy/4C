/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between beams and fluid elements.

\level 3

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/


#include "beam_to_fluid_mortar_manager.H"
#include "beam_to_fluid_meshtying_params.H"
#include "fbi_calc_utils.H"

#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"

#include "../drt_inpar/inpar_fbi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_serialdensevector.H"

#include <Epetra_FEVector.h>


/**
 *
 */
BEAMINTERACTION::BeamToFluidMortarManager::BeamToFluidMortarManager(
    Teuchos::RCP<const DRT::Discretization> discretization1,
    Teuchos::RCP<const DRT::Discretization> discretization2,
    Teuchos::RCP<const FBI::BeamToFluidMeshtyingParams> params, int start_value_lambda_gid)
    : is_setup_(false),
      is_local_maps_build_(false),
      is_global_maps_build_(false),
      start_value_lambda_gid_(start_value_lambda_gid),
      discretization_structure_(discretization1),
      discretization_fluid_(discretization2),
      beam_contact_parameters_ptr_(params),
      lambda_dof_rowmap_(Teuchos::null),
      lambda_dof_colmap_(Teuchos::null),
      beam_dof_rowmap_(Teuchos::null),
      fluid_dof_rowmap_(Teuchos::null),
      node_gid_to_lambda_gid_(Teuchos::null),
      element_gid_to_lambda_gid_(Teuchos::null),
      global_D_(Teuchos::null),
      global_M_(Teuchos::null),
      global_kappa_(Teuchos::null),
      global_active_lambda_(Teuchos::null)
{
  // Get the number of Lagrange multiplier DOF on a beam node and on a beam element.
  switch (params->GetMortarShapeFunctionType())
  {
    case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line2:
    {
      n_lambda_node_ = 1 * 3;
      n_lambda_element_ = 0 * 3;
      break;
    }
    case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line3:
    {
      n_lambda_node_ = 1 * 3;
      n_lambda_element_ = 1 * 3;
      break;
    }
    case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line4:
    {
      n_lambda_node_ = 1 * 3;
      n_lambda_element_ = 2 * 3;
      break;
    }
    default:
      dserror("Mortar shape function not implemented!");
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::Setup()
{
  // Get the global ids of all beam centerline nodes on this rank.
  std::vector<int> my_nodes_gid;
  for (int i_node = 0; i_node < discretization_structure_->NodeRowMap()->NumMyElements(); i_node++)
  {
    DRT::Node const& node = *(discretization_structure_->lRowNode(i_node));
    if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node)) my_nodes_gid.push_back(node.Id());
  }

  // Get the global ids of all beam elements on this rank.
  std::vector<int> my_elements_gid;
  for (int i_element = 0; i_element < discretization_structure_->ElementRowMap()->NumMyElements();
       i_element++)
  {
    DRT::Element const& element = *(discretization_structure_->lRowElement(i_element));
    if (BEAMINTERACTION::UTILS::IsBeamElement(element)) my_elements_gid.push_back(element.Id());
  }

  // Calculate the local number of centerline nodes, beam elements and Lagrange multiplier DOF.
  const unsigned int n_nodes = my_nodes_gid.size();
  const unsigned int n_element = my_elements_gid.size();
  const unsigned int n_lambda_dof = n_nodes * n_lambda_node_ + n_element * n_lambda_element_;


  // Tell all other processors how many lambda DOFs this processor has. This information is needed
  // to construct the lambda_dof_rowmap_.
  std::vector<int> lambda_dof_per_rank(discretization_structure_->Comm().NumProc(), 0);
  int temp_my_n_lambda_dof = (int)n_lambda_dof;
  discretization_structure_->Comm().GatherAll(&temp_my_n_lambda_dof, &lambda_dof_per_rank[0], 1);

  // Get the start GID for the lambda DOFs on this processor.
  int my_lambda_gid_start_value = start_value_lambda_gid_;
  for (int pid = 0; pid < discretization_structure_->Comm().MyPID(); pid++)
    my_lambda_gid_start_value += lambda_dof_per_rank[pid];

  // Fill in all GIDs of the lambda DOFs on this processor.
  std::vector<int> my_lambda_gid(n_lambda_dof, 0);
  for (int my_lid = 0; my_lid < (int)n_lambda_dof; my_lid++)
    my_lambda_gid[my_lid] = my_lambda_gid_start_value + my_lid;

  // Rowmap for the additional GIDs used by the mortar contact discretization.
  lambda_dof_rowmap_ = Teuchos::rcp(new Epetra_Map(
      -1, my_lambda_gid.size(), my_lambda_gid.data(), 0, discretization_structure_->Comm()));


  // We need to be able to get the global ids for a Lagrange multiplier DOF from the global id
  // of a node or element. To do so, we 'abuse' the Epetra_MultiVector as map between the
  // global node / element ids and the global Lagrange multiplier DOF ids.
  Teuchos::RCP<Epetra_Map> node_gid_rowmap = Teuchos::rcp(
      new Epetra_Map(-1, n_nodes, &my_nodes_gid[0], 0, discretization_structure_->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_rowmap = Teuchos::rcp(
      new Epetra_Map(-1, n_element, &my_elements_gid[0], 0, discretization_structure_->Comm()));

  // Map from global node / element ids to global lagrange multiplier ids. Only create the
  // multivector if it hase one or more columns.
  node_gid_to_lambda_gid_ = Teuchos::null;
  element_gid_to_lambda_gid_ = Teuchos::null;
  if (n_lambda_node_ > 0)
    node_gid_to_lambda_gid_ =
        Teuchos::rcp(new Epetra_MultiVector(*node_gid_rowmap, n_lambda_node_, true));
  if (n_lambda_element_ > 0)
    element_gid_to_lambda_gid_ =
        Teuchos::rcp(new Epetra_MultiVector(*element_gid_rowmap, n_lambda_element_, true));

  // Fill in the entries in the node / element global id to Lagrange multiplier global id vector.
  int error_code = 0;
  int lagrange_gid = -1;
  if (node_gid_to_lambda_gid_ != Teuchos::null)
  {
    for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
      for (unsigned int i_lambda = 0; i_lambda < n_lambda_node_; i_lambda++)
      {
        // Get the global Lagrange multiplier id for this node.
        lagrange_gid = lambda_dof_rowmap_->GID(i_node * n_lambda_node_ + i_lambda);

        // Set the global Lagrange multiplier id for this node.
        error_code = node_gid_to_lambda_gid_->ReplaceMyValue(i_node, i_lambda, lagrange_gid);
        if (error_code != 0) dserror("Got error code %d!", error_code);
      }
  }
  if (element_gid_to_lambda_gid_ != Teuchos::null)
  {
    for (unsigned int i_element = 0; i_element < n_element; i_element++)
      for (unsigned int i_lambda = 0; i_lambda < n_lambda_element_; i_lambda++)
      {
        // Get the global Lagrange multiplier id for this element.
        lagrange_gid = lambda_dof_rowmap_->GID(
            n_nodes * n_lambda_node_ + i_element * n_lambda_element_ + i_lambda);

        // Set the global Lagrange multiplier id for this element.
        error_code = element_gid_to_lambda_gid_->ReplaceMyValue(i_element, i_lambda, lagrange_gid);
        if (error_code != 0) dserror("Got error code %d!", error_code);
      }
  }

  // Create the global mortar matrices.
  global_D_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *lambda_dof_rowmap_, 30, true, true, LINALG::SparseMatrix::FE_MATRIX));
  global_M_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *lambda_dof_rowmap_, 100, true, true, LINALG::SparseMatrix::FE_MATRIX));
  global_kappa_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  global_active_lambda_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));

  // Create the maps for beam and solid DOFs.
  SetGlobalMaps();

  // Set flag for successful setup.
  is_setup_ = true;
  is_local_maps_build_ = false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::SetGlobalMaps()
{
  // Loop over all nodes on this processor -> we assume all beam and fluid DOFs are based on nodes.
  std::vector<std::vector<int>> field_dofs(2);

  for (int i_node = 0; i_node < discretization_structure_->NodeRowMap()->NumMyElements(); i_node++)
  {
    const DRT::Node* node = discretization_structure_->lRowNode(i_node);
    if (BEAMINTERACTION::UTILS::IsBeamNode(*node))
      discretization_structure_->Dof(node, field_dofs[0]);
    else
      dserror("The given structure element is not a beam element!");
  }
  for (int i_node = 0; i_node < discretization_fluid_->NodeRowMap()->NumMyElements(); i_node++)
  {
    const DRT::Node* node = discretization_fluid_->lRowNode(i_node);
    discretization_fluid_->Dof(node, field_dofs[1]);
  }

  // Create the beam and fluid maps.
  beam_dof_rowmap_ = Teuchos::rcp(new Epetra_Map(
      -1, field_dofs[0].size(), &(field_dofs[0])[0], 0, discretization_structure_->Comm()));
  fluid_dof_rowmap_ = Teuchos::rcp(new Epetra_Map(
      -1, field_dofs[1].size(), &(field_dofs[1])[0], 0, discretization_fluid_->Comm()));

  // Reset the local maps.
  node_gid_to_lambda_gid_map_.clear();
  element_gid_to_lambda_gid_map_.clear();

  // Set flags for global maps.
  is_global_maps_build_ = true;
  is_local_maps_build_ = false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::SetLocalMaps(
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs)
{
  CheckSetup();
  CheckGlobalMaps();

  // At this point the global multi vectors are filled up completely. To get the map for global
  // node element ids to the global lambda ids we need to be able to extract more than the local
  // values on this processor. Therefore we need a new map that contains all rows we want to
  // access in the global multi vector.
  std::vector<int> node_gid_needed;
  std::vector<int> element_gid_needed;

  // Loop over the pairs and get the global node and element indices needed on this rank.
  for (unsigned int i_pair = 0; i_pair < contact_pairs.size(); i_pair++)
  {
    const Teuchos::RCP<BEAMINTERACTION::BeamContactPair>& pair = contact_pairs[i_pair];

    // The first (beam) element should always be on the same processor as the pair.
    if (pair->Element1()->Owner() != discretization_structure_->Comm().MyPID())
      dserror(
          "The current implementation needs the first element of a beam contact pair to be on the "
          "same processor as the pair!");

    // Get the global id of the nodes / elements that the pairs on this rank need.
    if (n_lambda_node_ > 0)
      // There are nodal lambda DOFs, add the gid for the nodes in this element to the vector.
      // The first two nodes are the centerline nodes.
      for (unsigned int i_node = 0; i_node < 2; i_node++)
        node_gid_needed.push_back(pair->Element1()->Nodes()[i_node]->Id());

    if (n_lambda_element_ > 0)
      // There are element lambda DOFs, add the gid for this element to the vector.
      element_gid_needed.push_back(pair->Element1()->Id());
  }

  // Make the entries in the vectors unique.
  std::vector<int>::iterator it;
  std::sort(node_gid_needed.begin(), node_gid_needed.end());
  it = std::unique(node_gid_needed.begin(), node_gid_needed.end());
  node_gid_needed.resize(std::distance(node_gid_needed.begin(), it));
  std::sort(element_gid_needed.begin(), element_gid_needed.end());
  it = std::unique(element_gid_needed.begin(), element_gid_needed.end());
  element_gid_needed.resize(std::distance(element_gid_needed.begin(), it));

  // Create the maps for the extraction of the values.
  Teuchos::RCP<Epetra_Map> node_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, node_gid_needed.size(), &node_gid_needed[0], 0, discretization_structure_->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, element_gid_needed.size(), &element_gid_needed[0], 0, discretization_structure_->Comm()));

  // Create the Multivectors that will be filled with all values needed on this rank.
  Teuchos::RCP<Epetra_MultiVector> node_gid_to_lambda_gid_copy = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> element_gid_to_lambda_gid_copy = Teuchos::null;
  if (node_gid_to_lambda_gid_ != Teuchos::null)
    node_gid_to_lambda_gid_copy = Teuchos::rcp<Epetra_MultiVector>(
        new Epetra_MultiVector(*node_gid_needed_rowmap, n_lambda_node_, true));
  if (element_gid_to_lambda_gid_ != Teuchos::null)
    element_gid_to_lambda_gid_copy = Teuchos::rcp<Epetra_MultiVector>(
        new Epetra_MultiVector(*element_gid_needed_rowmap, n_lambda_element_, true));

  // Export values from the global multi vector to the ones needed on this rank.
  if (node_gid_to_lambda_gid_ != Teuchos::null)
    LINALG::Export(*node_gid_to_lambda_gid_, *node_gid_to_lambda_gid_copy);
  if (element_gid_to_lambda_gid_ != Teuchos::null)
    LINALG::Export(*element_gid_to_lambda_gid_, *element_gid_to_lambda_gid_copy);

  // Fill in the local maps.
  std::vector<int> lambda_gid_for_col_map;
  lambda_gid_for_col_map.clear();
  node_gid_to_lambda_gid_map_.clear();
  element_gid_to_lambda_gid_map_.clear();
  if (node_gid_to_lambda_gid_ != Teuchos::null)
  {
    std::vector<int> temp_node(n_lambda_node_);
    for (int i_node = 0; i_node < node_gid_needed_rowmap->NumMyElements(); i_node++)
    {
      for (unsigned int i_temp = 0; i_temp < n_lambda_node_; i_temp++)
        temp_node[i_temp] = (int)((*(*node_gid_to_lambda_gid_copy)(i_temp))[i_node]);
      node_gid_to_lambda_gid_map_[node_gid_needed_rowmap->GID(i_node)] = temp_node;
      lambda_gid_for_col_map.insert(
          std::end(lambda_gid_for_col_map), std::begin(temp_node), std::end(temp_node));
    }
  }
  if (element_gid_to_lambda_gid_ != Teuchos::null)
  {
    std::vector<int> temp_elements(n_lambda_element_);
    for (int i_element = 0; i_element < element_gid_needed_rowmap->NumMyElements(); i_element++)
    {
      for (unsigned int i_temp = 0; i_temp < n_lambda_element_; i_temp++)
        temp_elements[i_temp] = (int)((*(*element_gid_to_lambda_gid_copy)(i_temp))[i_element]);
      element_gid_to_lambda_gid_map_[element_gid_needed_rowmap->GID(i_element)] = temp_elements;
      lambda_gid_for_col_map.insert(
          std::end(lambda_gid_for_col_map), std::begin(temp_elements), std::end(temp_elements));
    }
  }

  // Create the global lambda col map.
  lambda_dof_colmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, lambda_gid_for_col_map.size(),
      &lambda_gid_for_col_map[0], 0, discretization_structure_->Comm()));

  // Set flags for local maps.
  is_local_maps_build_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::LocationVector(
    const Teuchos::RCP<const BEAMINTERACTION::BeamContactPair>& contact_pair,
    std::vector<int>& lambda_row) const
{
  CheckSetup();
  CheckLocalMaps();

  // Clear the output vectors.
  lambda_row.clear();

  // Get the global DOFs ids of the nodal Lagrange multipliers.
  if (n_lambda_node_ > 0)
  {
    for (int i_node = 0; i_node < contact_pair->Element1()->NumNode(); i_node++)
    {
      const DRT::Node& node = *(contact_pair->Element1()->Nodes()[i_node]);
      if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node))
      {
        // Get the global id of the node.
        int node_id = node.Id();

        // Check if the id is in the map. If it is, add it to the output vector.
        auto search_key_in_map = node_gid_to_lambda_gid_map_.find(node_id);
        if (search_key_in_map == node_gid_to_lambda_gid_map_.end())
          dserror("Global node id %d not in map!", node_id);
        for (auto const& lambda_gid : search_key_in_map->second) lambda_row.push_back(lambda_gid);
      }
    }
  }

  // Get the global DOFs ids of the element Lagrange multipliers.
  if (n_lambda_element_ > 0)
  {
    if (BEAMINTERACTION::UTILS::IsBeamElement(*contact_pair->Element1()))
    {
      // Get the global id of the element.
      int element_id = contact_pair->Element1()->Id();

      // Check if the id is in the map. If it is, add it to the output vector.
      auto search_key_in_map = element_gid_to_lambda_gid_map_.find(element_id);
      if (search_key_in_map == element_gid_to_lambda_gid_map_.end())
        dserror("Global element id %d not in map!", element_id);
      for (auto const& lambda_gid : search_key_in_map->second) lambda_row.push_back(lambda_gid);
    }
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::EvaluateGlobalDM(
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs)
{
  CheckSetup();
  CheckGlobalMaps();

  // Clear the old values of D, M and kappa.
  int linalg_error = 0;
  linalg_error = global_D_->PutScalar(0.);
  if (linalg_error != 0) dserror("Error in PutScalar!");
  linalg_error = global_M_->PutScalar(0.);
  if (linalg_error != 0) dserror("Error in PutScalar!");
  linalg_error = global_kappa_->PutScalar(0.);
  if (linalg_error != 0) dserror("Error in PutScalar!");

  // Local mortar matrices that will be filled up by EvaluateDM.
  LINALG::SerialDenseMatrix local_D_centerlineDOFs;
  LINALG::SerialDenseMatrix local_M;
  LINALG::SerialDenseVector dummy_constraint_offset;
  LINALG::SerialDenseVector local_kappa;

  // For the D matrix we need to assemble the centerline DOF to the element dof. This is done
  // into this matrix.
  LINALG::SerialDenseMatrix local_D_elementDOFs;

  // Flag if pair has a active contribution.
  bool mortar_is_active = false;

  // Loop over elements and assemble the local D and M matrices into the global ones.
  for (auto& elepairptr : contact_pairs)
  {
    // Evaluate the mortar contributions on the pair, if there are some, assemble into the global
    // matrices.
    mortar_is_active = elepairptr->EvaluateDM(
        local_D_centerlineDOFs, local_M, local_kappa, dummy_constraint_offset);

    if (mortar_is_active)
    {
      // We got contributions from the pair. Now we have to assemble into the global matrices. We
      // use the FEAssembly here, since the contact pairs are not ghosted.

      // Assemble the centerline matrix calculated by EvaluateDM into the full element matrix.
      BEAMINTERACTION::UTILS::AssembleCenterlineDofColMatrixIntoElementColMatrix(
          *discretization_structure_, elepairptr->Element1(), local_D_centerlineDOFs,
          local_D_elementDOFs);

      // Get the GIDs of the Lagrange multipliers.
      std::vector<int> lambda_row;
      LocationVector(elepairptr, lambda_row);

      // Get the GIDs of the beam DOF.
      std::vector<int> beam_row;
      std::vector<int> fluid_temp;
      std::vector<int> fluid_row;
      std::vector<int> dummy_1;
      std::vector<int> dummy_2;
      elepairptr->Element1()->LocationVector(
          *discretization_structure_, beam_row, dummy_1, dummy_2);
      elepairptr->Element2()->LocationVector(*discretization_fluid_, fluid_temp, dummy_1, dummy_2);
      for (unsigned int i = 0; i < fluid_temp.size(); i++)
      {
        if ((i + 1) % 4) fluid_row.push_back(fluid_temp[i]);
      }

      // Save check the matrix sizes.
      if (local_D_elementDOFs.RowDim() != (int)lambda_row.size() &&
          local_D_elementDOFs.ColDim() != (int)beam_row.size())
        dserror("Size of local D matrix does not match the GID vectors!");
      if (local_M.RowDim() != (int)lambda_row.size() && local_M.ColDim() != (int)fluid_row.size())
        dserror("Size of local M matrix does not match the GID vectors!");
      if (local_kappa.RowDim() != (int)lambda_row.size() && local_kappa.ColDim() != 1)
        dserror("Size of local kappa vector does not match the GID vector!");

      // Assemble into the global matrices.
      global_D_->FEAssemble(local_D_elementDOFs, lambda_row, beam_row);
      global_M_->FEAssemble(local_M, lambda_row, fluid_row);
      global_kappa_->SumIntoGlobalValues(
          local_kappa.RowDim(), &lambda_row[0], local_kappa.Values());

      // Set all entries in the local kappa vector to 1 and add them to the active vector.
      for (int i_local = 0; i_local < local_kappa.RowDim(); i_local++) local_kappa(i_local) = 1.;
      global_active_lambda_->SumIntoGlobalValues(
          local_kappa.RowDim(), &lambda_row[0], local_kappa.Values());
    }
  }

  // Complete the global mortar matrices.
  global_D_->Complete(*beam_dof_rowmap_, *lambda_dof_rowmap_);
  global_M_->Complete(*fluid_dof_rowmap_, *lambda_dof_rowmap_);

  // Complete the global scaling vector.
  if (global_kappa_->GlobalAssemble(Add, false)) dserror("Error in GlobalAssemble!");
  if (global_active_lambda_->GlobalAssemble(Add, false)) dserror("Error in GlobalAssemble!");
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMortarManager::AddGlobalForceStiffnessContributions(
    Teuchos::RCP<Epetra_FEVector> fluid_force, Teuchos::RCP<Epetra_FEVector> beam_force,
    Teuchos::RCP<LINALG::SparseMatrix> kbb, Teuchos::RCP<LINALG::SparseMatrix> kbf,
    Teuchos::RCP<LINALG::SparseMatrix> kff, Teuchos::RCP<LINALG::SparseMatrix> kfb,
    Teuchos::RCP<const Epetra_Vector> beam_vel, Teuchos::RCP<const Epetra_Vector> fluid_vel) const
{
  CheckSetup();
  CheckGlobalMaps();

  int linalg_error = 0;

  // Scale D and M with kappa^-1.
  Teuchos::RCP<Epetra_Vector> global_kappa_inv = InvertKappa();
  Teuchos::RCP<LINALG::SparseMatrix> kappa_inv_mat =
      Teuchos::rcp(new LINALG::SparseMatrix(*global_kappa_inv));
  kappa_inv_mat->Complete();
  Teuchos::RCP<LINALG::SparseMatrix> global_D_scaled =
      LINALG::MLMultiply(*kappa_inv_mat, false, *global_D_, false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> global_M_scaled =
      LINALG::MLMultiply(*kappa_inv_mat, false, *global_M_, false, false, false, true);

  // Calculate the needed submatrices.
  Teuchos::RCP<LINALG::SparseMatrix> Dt_kappa_D =
      LINALG::MLMultiply(*global_D_, true, *global_D_scaled, false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> Dt_kappa_M =
      LINALG::MLMultiply(*global_D_, true, *global_M_scaled, false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> Mt_kappa_M =
      LINALG::MLMultiply(*global_M_, true, *global_M_scaled, false, false, false, true);

  if (kff != Teuchos::null) kff->Add(*Mt_kappa_M, false, 1.0, 1.0);

  if (kfb != Teuchos::null) kbf->Add(*Dt_kappa_M, true, -1.0, 1.0);

  if (kbf != Teuchos::null) kfb->Add(*Dt_kappa_M, false, -1.0, 1.0);

  if (kbb != Teuchos::null) kbb->Add(*Dt_kappa_D, false, 1.0, 1.0);


  if (fluid_force != Teuchos::null)
  {
    if (fluid_vel == Teuchos::null || beam_vel == Teuchos::null)
      dserror("Force contributions can only be calculated with a given velocity pointer!");

    // Temporary vectors for matrix-vector multiplication and vector-vector additions.
    Teuchos::RCP<Epetra_Vector> fluid_temp = Teuchos::rcp(new Epetra_Vector(*fluid_dof_rowmap_));

    // Set the values in the global force vector to 0.
    linalg_error = fluid_force->PutScalar(0.);
    if (linalg_error != 0) dserror("Error in PutScalar!");

    linalg_error = Dt_kappa_M->Multiply(true, *beam_vel, *fluid_temp);
    if (linalg_error != 0) dserror("Error in Multiply!");
    linalg_error = fluid_force->Update(-1.0, *fluid_temp, 1.0);
    if (linalg_error != 0) dserror("Error in Update!");
  }

  Teuchos::RCP<Epetra_Vector> beam_temp = Teuchos::rcp(new Epetra_Vector(*beam_dof_rowmap_));
  linalg_error = beam_force->PutScalar(0.);
  if (linalg_error != 0) dserror("Error in PutScalar!");
  // Get the force acting on the beam.
  linalg_error = Dt_kappa_D->Multiply(false, *beam_vel, *beam_temp);
  if (linalg_error != 0) dserror("Error in Multiply!");
  linalg_error = beam_force->Update(1.0, *beam_temp, 1.0);
  if (linalg_error != 0) dserror("Error in Update!");
  linalg_error = Dt_kappa_M->Multiply(false, *fluid_vel, *beam_temp);
  if (linalg_error != 0) dserror("Error in Multiply!");
  linalg_error = beam_force->Update(-1.0, *beam_temp, 1.0);
  if (linalg_error != 0) dserror("Error in Update!");
}


/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToFluidMortarManager::GetGlobalLambda(
    Teuchos::RCP<const Epetra_Vector> vel) const
{
  CheckSetup();
  CheckGlobalMaps();

  // Get the velocity of the beams and the fluid.
  Teuchos::RCP<Epetra_Vector> beam_vel = Teuchos::rcp(new Epetra_Vector(*beam_dof_rowmap_));
  Teuchos::RCP<Epetra_Vector> fluid_vel = Teuchos::rcp(new Epetra_Vector(*fluid_dof_rowmap_));
  LINALG::Export(*vel, *beam_vel);
  LINALG::Export(*vel, *fluid_vel);

  // Set up lambda vector;
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));

  // Create a temporary vector and calculate lambda.
  Teuchos::RCP<Epetra_Vector> lambda_temp_1 = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
  Teuchos::RCP<Epetra_Vector> lambda_temp_2 = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
  int linalg_error = global_D_->Multiply(false, *beam_vel, *lambda_temp_2);
  if (linalg_error != 0) dserror("Error in Multiply!");
  linalg_error = lambda_temp_1->Update(1.0, *lambda_temp_2, 0.0);
  if (linalg_error != 0) dserror("Error in Update!");
  linalg_error = global_M_->Multiply(false, *fluid_vel, *lambda_temp_2);
  if (linalg_error != 0) dserror("Error in Multiply!");
  linalg_error = lambda_temp_1->Update(-1.0, *lambda_temp_2, 1.0);
  if (linalg_error != 0) dserror("Error in Multiply!");

  // Scale Lambda with kappa^-1.
  Teuchos::RCP<Epetra_Vector> global_kappa_inv = InvertKappa();
  Teuchos::RCP<LINALG::SparseMatrix> kappa_inv_mat =
      Teuchos::rcp(new LINALG::SparseMatrix(*global_kappa_inv));
  kappa_inv_mat->Complete();
  linalg_error = kappa_inv_mat->Multiply(false, *lambda_temp_1, *lambda);
  if (linalg_error != 0) dserror("Error in Multiply!");

  return lambda;
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToFluidMortarManager::GetGlobalLambdaCol(
    Teuchos::RCP<const Epetra_Vector> disp) const
{
  Teuchos::RCP<Epetra_Vector> lambda_col = Teuchos::rcp(new Epetra_Vector(*lambda_dof_colmap_));
  LINALG::Export(*GetGlobalLambda(disp), *lambda_col);
  return lambda_col;
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToFluidMortarManager::InvertKappa() const
{
  // Create the inverse vector.
  Teuchos::RCP<Epetra_Vector> global_kappa_inv =
      Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));

  // Calculate the local inverse of kappa.
  double local_kappa_inv_value = 0.;
  for (int lid = 0; lid < lambda_dof_rowmap_->NumMyElements(); lid++)
  {
    if (global_active_lambda_->Values()[lid] > 0.1)
      local_kappa_inv_value = 1. / global_kappa_->Values()[lid];
    else
      local_kappa_inv_value = 0.0;

    global_kappa_inv->ReplaceMyValue(lid, 0, local_kappa_inv_value);
  }

  return global_kappa_inv;
}
