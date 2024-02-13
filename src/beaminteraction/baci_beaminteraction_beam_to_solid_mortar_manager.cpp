/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between beams and solids.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_mortar_manager.hpp"

#include "baci_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "baci_beaminteraction_beam_to_solid_utils.hpp"
#include "baci_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "baci_beaminteraction_calc_utils.hpp"
#include "baci_beaminteraction_contact_pair.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_beaminteraction_str_model_evaluator_datastate.hpp"
#include "baci_geometry_pair.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_multiply.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidMortarManager::BeamToSolidMortarManager(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& params,
    int start_value_lambda_gid)
    : is_setup_(false),
      is_local_maps_build_(false),
      is_global_maps_build_(false),
      start_value_lambda_gid_(start_value_lambda_gid),
      n_lambda_node_translational_(0),
      n_lambda_element_translational_(0),
      n_lambda_node_rotational_(0),
      n_lambda_element_rotational_(0),
      discret_(discret),
      beam_to_solid_params_(params),
      lambda_dof_rowmap_translations_(Teuchos::null),
      lambda_dof_rowmap_rotations_(Teuchos::null),
      lambda_dof_rowmap_(Teuchos::null),
      lambda_dof_colmap_(Teuchos::null),
      beam_dof_rowmap_(Teuchos::null),
      solid_dof_rowmap_(Teuchos::null),
      node_gid_to_lambda_gid_(Teuchos::null),
      element_gid_to_lambda_gid_(Teuchos::null),
      global_constraint_(Teuchos::null),
      global_G_B_(Teuchos::null),
      global_G_S_(Teuchos::null),
      global_FB_L_(Teuchos::null),
      global_FS_L_(Teuchos::null),
      global_kappa_(Teuchos::null),
      global_active_lambda_(Teuchos::null),
      contact_pairs_(Teuchos::null)
{
  // Get the number of Lagrange multiplier DOF on a beam node and on a beam element.
  unsigned int n_lambda_node_temp = 0;
  unsigned int n_lambda_element_temp = 0;
  MortarShapeFunctionsToNumberOfLagrangeValues(beam_to_solid_params_->GetMortarShapeFunctionType(),
      n_lambda_node_temp, n_lambda_element_temp);
  n_lambda_node_ = n_lambda_node_temp;
  n_lambda_node_translational_ = n_lambda_node_temp;
  n_lambda_element_ = n_lambda_element_temp;
  n_lambda_element_translational_ = n_lambda_element_temp;

  // Check if the coupling also consists of rotational coupling.
  auto beam_to_volume_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>(
          beam_to_solid_params_);
  auto beam_to_surface_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>(
          beam_to_solid_params_);
  if (beam_to_volume_params != Teuchos::null)
  {
    // Get the number of Lagrange multiplier DOF for rotational coupling on a beam node and on a
    // beam element.
    MortarShapeFunctionsToNumberOfLagrangeValues(
        beam_to_volume_params->GetMortarShapeFunctionRotationType(), n_lambda_node_temp,
        n_lambda_element_temp);
    n_lambda_node_ += n_lambda_node_temp;
    n_lambda_node_rotational_ = n_lambda_node_temp;
    n_lambda_element_ += n_lambda_element_temp;
    n_lambda_element_rotational_ = n_lambda_element_temp;
  }
  else if (beam_to_surface_params != Teuchos::null)
  {
    if (beam_to_surface_params->GetIsRotationalCoupling())
    {
      MortarShapeFunctionsToNumberOfLagrangeValues(
          beam_to_surface_params->GetMortarShapeFunctionType(), n_lambda_node_temp,
          n_lambda_element_temp);
      n_lambda_node_ += n_lambda_node_temp;
      n_lambda_node_rotational_ = n_lambda_node_temp;
      n_lambda_element_ += n_lambda_element_temp;
      n_lambda_element_rotational_ = n_lambda_element_temp;
    }
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::Setup()
{
  // Get the global ids of all beam centerline nodes on this rank.
  std::vector<int> my_nodes_gid;
  for (int i_node = 0; i_node < discret_->NodeRowMap()->NumMyElements(); i_node++)
  {
    DRT::Node const& node = *(discret_->lRowNode(i_node));
    if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node)) my_nodes_gid.push_back(node.Id());
  }

  // Get the global ids of all beam elements on this rank.
  std::vector<int> my_elements_gid;
  for (int i_element = 0; i_element < discret_->ElementRowMap()->NumMyElements(); i_element++)
  {
    DRT::Element const& element = *(discret_->lRowElement(i_element));
    if (BEAMINTERACTION::UTILS::IsBeamElement(element)) my_elements_gid.push_back(element.Id());
  }

  // Calculate the local number of centerline nodes, beam elements and Lagrange multiplier DOF.
  const unsigned int n_nodes = my_nodes_gid.size();
  const unsigned int n_element = my_elements_gid.size();
  const unsigned int n_lambda_dof = n_nodes * n_lambda_node_ + n_element * n_lambda_element_;


  // Tell all other processors how many lambda DOFs this processor has. This information is needed
  // to construct the lambda_dof_rowmap_.
  std::vector<int> lambda_dof_per_rank(discret_->Comm().NumProc(), 0);
  int temp_my_n_lambda_dof = (int)n_lambda_dof;
  discret_->Comm().GatherAll(&temp_my_n_lambda_dof, lambda_dof_per_rank.data(), 1);

  // Get the start GID for the lambda DOFs on this processor.
  int my_lambda_gid_start_value = start_value_lambda_gid_;
  for (int pid = 0; pid < discret_->Comm().MyPID(); pid++)
    my_lambda_gid_start_value += lambda_dof_per_rank[pid];

  // Fill in all GIDs of the lambda DOFs on this processor (for translations and rotations).
  std::vector<int> my_lambda_gid_translational;
  my_lambda_gid_translational.reserve(
      n_nodes * n_lambda_node_translational_ + n_element * n_lambda_element_translational_);
  std::vector<int> my_lambda_gid_rotational;
  my_lambda_gid_rotational.reserve(
      n_nodes * n_lambda_node_rotational_ + n_element * n_lambda_element_rotational_);
  for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
  {
    for (unsigned int i_dof = 0; i_dof < n_lambda_node_translational_; i_dof++)
      my_lambda_gid_translational.push_back(
          my_lambda_gid_start_value + i_dof + i_node * n_lambda_node_);
    for (unsigned int i_dof = 0; i_dof < n_lambda_node_rotational_; i_dof++)
      my_lambda_gid_rotational.push_back(my_lambda_gid_start_value + i_dof +
                                         i_node * n_lambda_node_ + n_lambda_node_translational_);
  }
  my_lambda_gid_start_value += static_cast<int>(n_nodes * n_lambda_node_);
  for (unsigned int i_element = 0; i_element < n_element; i_element++)
  {
    for (unsigned int i_dof = 0; i_dof < n_lambda_element_translational_; i_dof++)
      my_lambda_gid_translational.push_back(
          my_lambda_gid_start_value + i_dof + i_element * n_lambda_element_);
    for (unsigned int i_dof = 0; i_dof < n_lambda_element_rotational_; i_dof++)
      my_lambda_gid_rotational.push_back(my_lambda_gid_start_value + i_dof +
                                         i_element * n_lambda_element_ +
                                         n_lambda_element_translational_);
  }

  // Rowmap for the additional GIDs used by the mortar contact discretization.
  lambda_dof_rowmap_translations_ = Teuchos::rcp(new Epetra_Map(-1,
      my_lambda_gid_translational.size(), my_lambda_gid_translational.data(), 0, discret_->Comm()));
  lambda_dof_rowmap_rotations_ = Teuchos::rcp(new Epetra_Map(
      -1, my_lambda_gid_rotational.size(), my_lambda_gid_rotational.data(), 0, discret_->Comm()));
  lambda_dof_rowmap_ =
      CORE::LINALG::MergeMap(lambda_dof_rowmap_translations_, lambda_dof_rowmap_rotations_, false);

  // We need to be able to get the global ids for a Lagrange multiplier DOF from the global id
  // of a node or element. To do so, we 'abuse' the Epetra_MultiVector as map between the
  // global node / element ids and the global Lagrange multiplier DOF ids.
  Teuchos::RCP<Epetra_Map> node_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_nodes, my_nodes_gid.data(), 0, discret_->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_element, my_elements_gid.data(), 0, discret_->Comm()));

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

  // Create the maps for beam and solid DOFs.
  SetGlobalMaps();

  // Create the global coupling matrices.
  global_constraint_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  global_G_B_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *lambda_dof_rowmap_, 30, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  global_G_S_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *lambda_dof_rowmap_, 100, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  global_FB_L_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *beam_dof_rowmap_, 30, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  global_FS_L_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *solid_dof_rowmap_, 100, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  global_kappa_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  global_active_lambda_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));

  // Set flag for successful setup.
  is_setup_ = true;
  is_local_maps_build_ = false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::SetGlobalMaps()
{
  // Loop over all nodes on this processor -> we assume all beam and solid DOFs are based on nodes.
  std::vector<int> beam_dofs(0);
  std::vector<int> solid_dofs(0);
  for (int i_node = 0; i_node < discret_->NodeRowMap()->NumMyElements(); i_node++)
  {
    const DRT::Node* node = discret_->lRowNode(i_node);
    if (BEAMINTERACTION::UTILS::IsBeamNode(*node))
      discret_->Dof(node, beam_dofs);
    else
      discret_->Dof(node, solid_dofs);
  }

  // Create the beam and solid maps.
  beam_dof_rowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, beam_dofs.size(), beam_dofs.data(), 0, discret_->Comm()));
  solid_dof_rowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, solid_dofs.size(), solid_dofs.data(), 0, discret_->Comm()));

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
void BEAMINTERACTION::BeamToSolidMortarManager::SetLocalMaps(
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs)
{
  CheckSetup();
  CheckGlobalMaps();

  contact_pairs_ = contact_pairs;

  // At this point the global multi vectors are filled up completely. To get the map for global
  // node element ids to the global lambda ids we need to be able to extract more than the local
  // values on this processor. Therefore we need a new map that contains all rows we want to
  // access in the global multi vector.
  std::vector<int> node_gid_needed;
  std::vector<int> element_gid_needed;

  // Loop over the pairs and get the global node and element indices needed on this rank.
  for (unsigned int i_pair = 0; i_pair < contact_pairs_.size(); i_pair++)
  {
    const Teuchos::RCP<BEAMINTERACTION::BeamContactPair>& pair = contact_pairs_[i_pair];

    // The first (beam) element should always be on the same processor as the pair.
    if (pair->Element1()->Owner() != discret_->Comm().MyPID())
      dserror(
          "The current implementation need the first element of a beam contact pair to be on the "
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
  Teuchos::RCP<Epetra_Map> node_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, node_gid_needed.size(), node_gid_needed.data(), 0, discret_->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, element_gid_needed.size(), element_gid_needed.data(), 0, discret_->Comm()));

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
    CORE::LINALG::Export(*node_gid_to_lambda_gid_, *node_gid_to_lambda_gid_copy);
  if (element_gid_to_lambda_gid_ != Teuchos::null)
    CORE::LINALG::Export(*element_gid_to_lambda_gid_, *element_gid_to_lambda_gid_copy);

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
  lambda_dof_colmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, lambda_gid_for_col_map.size(), lambda_gid_for_col_map.data(), 0, discret_->Comm()));

  // Set flags for local maps.
  is_local_maps_build_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::LocationVector(
    const BEAMINTERACTION::BeamContactPair* contact_pair, std::vector<int>& lambda_row) const
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
void BEAMINTERACTION::BeamToSolidMortarManager::EvaluateGlobalCouplingContributions(
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  CheckSetup();
  CheckGlobalMaps();

  // Reset the global data structures.
  global_constraint_->PutScalar(0.);
  global_G_B_->PutScalar(0.);
  global_G_S_->PutScalar(0.);
  global_FB_L_->PutScalar(0.);
  global_FS_L_->PutScalar(0.);
  global_kappa_->PutScalar(0.);
  global_active_lambda_->PutScalar(0.);

  for (auto& elepairptr : contact_pairs_)
  {
    // Evaluate the mortar contributions of the pair and the pair assembles the terms into the
    // global matrices.
    elepairptr->EvaluateAndAssembleMortarContributions(*discret_, this, *global_G_B_, *global_G_S_,
        *global_FB_L_, *global_FS_L_, *global_constraint_, *global_kappa_, *global_active_lambda_,
        displacement_vector);
  }

  // Complete the global mortar matrices.
  global_G_B_->Complete(*beam_dof_rowmap_, *lambda_dof_rowmap_);
  global_G_S_->Complete(*solid_dof_rowmap_, *lambda_dof_rowmap_);
  global_FB_L_->Complete(*lambda_dof_rowmap_, *beam_dof_rowmap_);
  global_FS_L_->Complete(*lambda_dof_rowmap_, *solid_dof_rowmap_);

  // Complete the global scaling vector.
  if (0 != global_kappa_->GlobalAssemble(Add, false)) dserror("Error in GlobalAssemble!");
  if (0 != global_active_lambda_->GlobalAssemble(Add, false)) dserror("Error in GlobalAssemble!");
  if (0 != global_constraint_->GlobalAssemble(Add, false)) dserror("Error in GlobalAssemble!");
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::AddGlobalForceStiffnessPenaltyContributions(
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff, Teuchos::RCP<Epetra_FEVector> force) const
{
  CheckSetup();
  CheckGlobalMaps();

  int linalg_error = 0;

  if (stiff != Teuchos::null)
  {
    // Scale the linearizations of the constraint equations.
    Teuchos::RCP<Epetra_Vector> global_penalty_kappa_inv = PenaltyInvertKappa();
    Teuchos::RCP<CORE::LINALG::SparseMatrix> penalty_kappa_inv_mat =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(*global_penalty_kappa_inv));
    penalty_kappa_inv_mat->Complete();

    Teuchos::RCP<CORE::LINALG::SparseMatrix> global_G_B_scaled = CORE::LINALG::MLMultiply(
        *penalty_kappa_inv_mat, false, *global_G_B_, false, false, false, true);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> global_G_S_scaled = CORE::LINALG::MLMultiply(
        *penalty_kappa_inv_mat, false, *global_G_S_, false, false, false, true);

    // Calculate the needed submatrices.
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FB_L_times_G_B = CORE::LINALG::MLMultiply(
        *global_FB_L_, false, *global_G_B_scaled, false, false, false, true);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FB_L_times_G_S = CORE::LINALG::MLMultiply(
        *global_FB_L_, false, *global_G_S_scaled, false, false, false, true);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FS_L_times_G_B = CORE::LINALG::MLMultiply(
        *global_FS_L_, false, *global_G_B_scaled, false, false, false, true);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FS_L_times_G_S = CORE::LINALG::MLMultiply(
        *global_FS_L_, false, *global_G_S_scaled, false, false, false, true);

    // Add contributions to the global stiffness matrix.
    stiff->Add(*FB_L_times_G_B, false, 1.0, 1.0);
    stiff->Add(*FB_L_times_G_S, false, 1.0, 1.0);
    stiff->Add(*FS_L_times_G_B, false, 1.0, 1.0);
    stiff->Add(*FS_L_times_G_S, false, 1.0, 1.0);
  }

  if (force != Teuchos::null)
  {
    // Factor for right hand side (forces). 1 corresponds to the meshtying forces being added to the
    // right hand side, -1 to the left hand side.
    const double rhs_factor = -1.0;

    // Get the penalty Lagrange multiplier vector.
    Teuchos::RCP<Epetra_Vector> lambda = GetGlobalLambda();

    // Multiply the lambda vector with FB_L and FS_L to get the forces on the beam and solid,
    // respectively.
    Teuchos::RCP<Epetra_Vector> beam_force = Teuchos::rcp(new Epetra_Vector(*beam_dof_rowmap_));
    Teuchos::RCP<Epetra_Vector> solid_force = Teuchos::rcp(new Epetra_Vector(*solid_dof_rowmap_));
    beam_force->PutScalar(0.);
    solid_force->PutScalar(0.);
    linalg_error = global_FB_L_->Multiply(false, *lambda, *beam_force);
    if (linalg_error != 0) dserror("Error in Multiply!");
    linalg_error = global_FS_L_->Multiply(false, *lambda, *solid_force);
    if (linalg_error != 0) dserror("Error in Multiply!");
    Teuchos::RCP<Epetra_Vector> global_temp =
        Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    CORE::LINALG::Export(*beam_force, *global_temp);
    CORE::LINALG::Export(*solid_force, *global_temp);

    // Add force contributions to global vector.
    linalg_error = force->Update(-1.0 * rhs_factor, *global_temp, 1.0);
    if (linalg_error != 0) dserror("Error in Update");
  }

  // Add the force and stiffness contributions that are assembled directly by the pairs.
  Teuchos::RCP<Epetra_Vector> lambda_col = GetGlobalLambdaCol();
  for (const auto& elepairptr : contact_pairs_)
    elepairptr->EvaluateAndAssemble(
        *discret_, this, force, stiff, *lambda_col, *data_state->GetDisColNp());
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::GetGlobalLambda() const
{
  CheckSetup();
  CheckGlobalMaps();

  // Get the inverted kappa matrix.
  Teuchos::RCP<Epetra_Vector> penalty_global_kappa_inv = PenaltyInvertKappa();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> penalty_kappa_inv_mat =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*penalty_global_kappa_inv));
  penalty_kappa_inv_mat->Complete();

  // Multiply the inverted kappa matrix with the constraint equations and scale them with the
  // penalty parameter.
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
  int linalg_error = penalty_kappa_inv_mat->Multiply(false, *global_constraint_, *lambda);
  if (linalg_error != 0) dserror("Error in Multiply!");

  return lambda;
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::GetGlobalLambdaCol() const
{
  Teuchos::RCP<Epetra_Vector> lambda_col = Teuchos::rcp(new Epetra_Vector(*lambda_dof_colmap_));
  CORE::LINALG::Export(*GetGlobalLambda(), *lambda_col);
  return lambda_col;
}

/**
 *
 */
double BEAMINTERACTION::BeamToSolidMortarManager::GetEnergy() const
{
  // Since this value is also computed for the reference configuration, where the global mortar
  // matrices are not build yet we return 0 in this case.
  if (not global_G_B_->Filled() or not global_G_S_->Filled() or not global_FB_L_->Filled() or
      not global_FS_L_->Filled())
    return 0.0;

  // Calculate the penalty potential.
  Teuchos::RCP<Epetra_Vector> lambda = GetGlobalLambda();
  double dot_product = 0.0;
  global_constraint_->Dot(*lambda, &dot_product);

  // Only rank 0 should return the global energy value.
  if (global_constraint_->Comm().MyPID() == 0)
    return 0.5 * dot_product;
  else
    return 0.0;
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::PenaltyInvertKappa() const
{
  // Create the inverse vector.
  Teuchos::RCP<Epetra_Vector> global_kappa_inv =
      Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));

  // Get the penalty parameters.
  const double penalty_translation = beam_to_solid_params_->GetPenaltyParameter();
  double penalty_rotation = 0.0;
  auto beam_to_volume_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>(
          beam_to_solid_params_);
  auto beam_to_surface_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>(
          beam_to_solid_params_);
  if (beam_to_volume_params != Teuchos::null)
    penalty_rotation = beam_to_volume_params->GetRotationalCouplingPenaltyParameter();
  else if (beam_to_surface_params != Teuchos::null)
    penalty_rotation = beam_to_surface_params->GetRotationalCouplingPenaltyParameter();
  else if (lambda_dof_rowmap_rotations_->NumGlobalElements() > 0)
    dserror(
        "Rotational penalty coupling only implemented for beam-to-volume and beam-to-surface "
        "case.");

  // Calculate the local inverse of kappa.
  double penalty = 0.0;
  double local_kappa_inv_value = 0.;
  for (int lid = 0; lid < lambda_dof_rowmap_->NumMyElements(); lid++)
  {
    if (global_active_lambda_->Values()[lid] > 0.1)
    {
      const int gid = lambda_dof_rowmap_->GID(lid);
      if (lambda_dof_rowmap_translations_->LID(gid) != -1)
        penalty = penalty_translation;
      else if (lambda_dof_rowmap_rotations_->LID(gid) != -1)
        penalty = penalty_rotation;
      else
        dserror("Could not find the GID %d in translation or rotation map", gid);

      local_kappa_inv_value = penalty / global_kappa_->Values()[lid];
    }

    else
      // This LID is inactive.
      local_kappa_inv_value = 0.0;

    global_kappa_inv->ReplaceMyValue(lid, 0, local_kappa_inv_value);
  }

  return global_kappa_inv;
}

BACI_NAMESPACE_CLOSE
