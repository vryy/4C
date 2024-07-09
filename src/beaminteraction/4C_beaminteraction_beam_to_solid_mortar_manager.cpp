/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between beams and solids.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_geometry_pair.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidMortarManager::BeamToSolidMortarManager(
    const Teuchos::RCP<const Core::FE::Discretization>& discret,
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
      beam_to_solid_params_(params)
{
  // Get the dimension of the Lagrange multiplier field
  const bool is_contact =
      !(Teuchos::rcp_dynamic_cast<const BeamToSolidSurfaceContactParams>(beam_to_solid_params_)
              .is_null());
  const unsigned int n_dim = is_contact ? 1 : 3;

  // Get the number of Lagrange multiplier DOF on a beam node and on a beam element.
  const auto& [n_lambda_node_pos, n_lambda_element_pos] =
      MortarShapeFunctionsToNumberOfLagrangeValues(
          beam_to_solid_params_->get_mortar_shape_function_type(), n_dim);
  n_lambda_node_ = n_lambda_node_pos;
  n_lambda_node_translational_ = n_lambda_node_pos;
  n_lambda_element_ = n_lambda_element_pos;
  n_lambda_element_translational_ = n_lambda_element_pos;

  if (beam_to_solid_params_->is_rotational_coupling())
  {
    // Get the mortar shape functions for rotational coupling
    auto mortar_shape_function_rotation = Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::none;

    const auto beam_to_volume_params =
        Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>(
            beam_to_solid_params_);
    const auto beam_to_surface_params =
        Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>(
            beam_to_solid_params_);
    if (beam_to_volume_params == beam_to_surface_params)
    {
      FOUR_C_THROW("The params object should be either of beam-to-solid volume or surface type.");
    }
    else if (beam_to_volume_params != Teuchos::null)
    {
      mortar_shape_function_rotation =
          beam_to_volume_params->get_mortar_shape_function_rotation_type();
    }
    else
    {
      mortar_shape_function_rotation = beam_to_surface_params->get_mortar_shape_function_type();
    }

    // Get the number of Lagrange multiplier DOF for rotational coupling on a beam node and on a
    // beam element.
    const auto& [n_lambda_node_rot, n_lambda_element_rot] =
        MortarShapeFunctionsToNumberOfLagrangeValues(mortar_shape_function_rotation, n_dim);
    n_lambda_node_ += n_lambda_node_rot;
    n_lambda_node_rotational_ = n_lambda_node_rot;
    n_lambda_element_ += n_lambda_element_rot;
    n_lambda_element_rotational_ = n_lambda_element_rot;
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::setup()
{
  // Get the global ids of all beam centerline nodes on this rank.
  std::vector<int> my_nodes_gid;
  for (int i_node = 0; i_node < discret_->node_row_map()->NumMyElements(); i_node++)
  {
    Core::Nodes::Node const& node = *(discret_->l_row_node(i_node));
    if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node)) my_nodes_gid.push_back(node.id());
  }

  // Get the global ids of all beam elements on this rank.
  std::vector<int> my_elements_gid;
  for (int i_element = 0; i_element < discret_->element_row_map()->NumMyElements(); i_element++)
  {
    Core::Elements::Element const& element = *(discret_->l_row_element(i_element));
    if (BEAMINTERACTION::UTILS::IsBeamElement(element)) my_elements_gid.push_back(element.id());
  }

  // Calculate the local number of centerline nodes, beam elements and Lagrange multiplier DOF.
  const unsigned int n_nodes = my_nodes_gid.size();
  const unsigned int n_element = my_elements_gid.size();
  const unsigned int n_lambda_dof = n_nodes * n_lambda_node_ + n_element * n_lambda_element_;


  // Tell all other processors how many lambda DOFs this processor has. This information is needed
  // to construct the lambda_dof_rowmap_.
  std::vector<int> lambda_dof_per_rank(discret_->get_comm().NumProc(), 0);
  int temp_my_n_lambda_dof = (int)n_lambda_dof;
  discret_->get_comm().GatherAll(&temp_my_n_lambda_dof, lambda_dof_per_rank.data(), 1);

  // Get the start GID for the lambda DOFs on this processor.
  int my_lambda_gid_start_value = start_value_lambda_gid_;
  for (int pid = 0; pid < discret_->get_comm().MyPID(); pid++)
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
  lambda_dof_rowmap_translations_ =
      Teuchos::rcp(new Epetra_Map(-1, my_lambda_gid_translational.size(),
          my_lambda_gid_translational.data(), 0, discret_->get_comm()));
  lambda_dof_rowmap_rotations_ = Teuchos::rcp(new Epetra_Map(-1, my_lambda_gid_rotational.size(),
      my_lambda_gid_rotational.data(), 0, discret_->get_comm()));
  lambda_dof_rowmap_ =
      Core::LinAlg::MergeMap(lambda_dof_rowmap_translations_, lambda_dof_rowmap_rotations_, false);

  // We need to be able to get the global ids for a Lagrange multiplier DOF from the global id
  // of a node or element. To do so, we 'abuse' the Epetra_MultiVector as map between the
  // global node / element ids and the global Lagrange multiplier DOF ids.
  Teuchos::RCP<Epetra_Map> node_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_nodes, my_nodes_gid.data(), 0, discret_->get_comm()));
  Teuchos::RCP<Epetra_Map> element_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_element, my_elements_gid.data(), 0, discret_->get_comm()));

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
        if (error_code != 0) FOUR_C_THROW("Got error code %d!", error_code);
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
        if (error_code != 0) FOUR_C_THROW("Got error code %d!", error_code);
      }
  }

  // Create the maps for beam and solid DOFs.
  set_global_maps();

  // Create the global coupling matrices.
  constraint_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  constraint_lin_beam_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 30, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  constraint_lin_solid_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  force_beam_lin_lambda_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *beam_dof_rowmap_, 30, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  force_solid_lin_lambda_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *solid_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  kappa_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  kappa_lin_beam_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 30, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  kappa_lin_solid_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  lambda_active_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));

  // Set flag for successful setup.
  is_setup_ = true;
  is_local_maps_build_ = false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::set_global_maps()
{
  // Loop over all nodes on this processor -> we assume all beam and solid DOFs are based on nodes.
  std::vector<int> beam_dofs(0);
  std::vector<int> solid_dofs(0);
  for (int i_node = 0; i_node < discret_->node_row_map()->NumMyElements(); i_node++)
  {
    const Core::Nodes::Node* node = discret_->l_row_node(i_node);
    if (BEAMINTERACTION::UTILS::IsBeamNode(*node))
      discret_->dof(node, beam_dofs);
    else
      discret_->dof(node, solid_dofs);
  }

  // Create the beam and solid maps.
  beam_dof_rowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, beam_dofs.size(), beam_dofs.data(), 0, discret_->get_comm()));
  solid_dof_rowmap_ = Teuchos::rcp(
      new Epetra_Map(-1, solid_dofs.size(), solid_dofs.data(), 0, discret_->get_comm()));

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
void BEAMINTERACTION::BeamToSolidMortarManager::set_local_maps(
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs)
{
  check_setup();
  check_global_maps();

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
    if (pair->element1()->owner() != discret_->get_comm().MyPID())
      FOUR_C_THROW(
          "The current implementation need the first element of a beam contact pair to be on the "
          "same processor as the pair!");

    // Get the global id of the nodes / elements that the pairs on this rank need.
    if (n_lambda_node_ > 0)
      // There are nodal lambda DOFs, add the gid for the nodes in this element to the vector.
      // The first two nodes are the centerline nodes.
      for (unsigned int i_node = 0; i_node < 2; i_node++)
        node_gid_needed.push_back(pair->element1()->nodes()[i_node]->id());

    if (n_lambda_element_ > 0)
      // There are element lambda DOFs, add the gid for this element to the vector.
      element_gid_needed.push_back(pair->element1()->id());
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
      new Epetra_Map(-1, node_gid_needed.size(), node_gid_needed.data(), 0, discret_->get_comm()));
  Teuchos::RCP<Epetra_Map> element_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, element_gid_needed.size(), element_gid_needed.data(), 0, discret_->get_comm()));

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
    Core::LinAlg::Export(*node_gid_to_lambda_gid_, *node_gid_to_lambda_gid_copy);
  if (element_gid_to_lambda_gid_ != Teuchos::null)
    Core::LinAlg::Export(*element_gid_to_lambda_gid_, *element_gid_to_lambda_gid_copy);

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
      -1, lambda_gid_for_col_map.size(), lambda_gid_for_col_map.data(), 0, discret_->get_comm()));

  // Set flags for local maps.
  is_local_maps_build_ = true;
}

/**
 *
 */
std::pair<std::vector<int>, std::vector<int>>
BEAMINTERACTION::BeamToSolidMortarManager::location_vector(
    const BEAMINTERACTION::BeamContactPair& contact_pair) const
{
  check_setup();
  check_local_maps();

  // Create the output vectors
  std::vector<int> lambda_pos_row;
  std::vector<int> lambda_rot_row;

  // Get the global DOFs ids of the nodal Lagrange multipliers.
  if (n_lambda_node_ > 0)
  {
    for (int i_node = 0; i_node < contact_pair.element1()->num_node(); i_node++)
    {
      const Core::Nodes::Node& node = *(contact_pair.element1()->nodes()[i_node]);
      if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node))
      {
        // Get the global id of the node.
        int node_id = node.id();

        // Check if the id is in the map. If it is, add it to the output vector.
        auto search_key_in_map = node_gid_to_lambda_gid_map_.find(node_id);
        if (search_key_in_map == node_gid_to_lambda_gid_map_.end())
          FOUR_C_THROW("Global node id %d not in map!", node_id);
        const auto node_lambda_gid = search_key_in_map->second;
        for (unsigned int i_pos = 0; i_pos < n_lambda_node_translational_; i_pos++)
        {
          lambda_pos_row.push_back(node_lambda_gid[i_pos]);
        }
        for (unsigned int i_rot = 0; i_rot < n_lambda_node_rotational_; i_rot++)
        {
          lambda_rot_row.push_back(node_lambda_gid[n_lambda_node_translational_ + i_rot]);
        }
      }
    }
  }

  // Get the global DOFs ids of the element Lagrange multipliers.
  if (n_lambda_element_ > 0)
  {
    if (BEAMINTERACTION::UTILS::IsBeamElement(*contact_pair.element1()))
    {
      // Get the global id of the element.
      int element_id = contact_pair.element1()->id();

      // Check if the id is in the map. If it is, add it to the output vector.
      auto search_key_in_map = element_gid_to_lambda_gid_map_.find(element_id);
      if (search_key_in_map == element_gid_to_lambda_gid_map_.end())
        FOUR_C_THROW("Global element id %d not in map!", element_id);
      const auto element_lambda_gid = search_key_in_map->second;
      for (unsigned int i_pos = 0; i_pos < n_lambda_element_translational_; i_pos++)
      {
        lambda_pos_row.push_back(element_lambda_gid[i_pos]);
      }
      for (unsigned int i_rot = 0; i_rot < n_lambda_element_rotational_; i_rot++)
      {
        lambda_rot_row.push_back(element_lambda_gid[n_lambda_element_translational_ + i_rot]);
      }
    }
  }

  return {lambda_pos_row, lambda_rot_row};
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::evaluate_force_stiff_penalty_regularization(
    const Teuchos::RCP<const Solid::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff, Teuchos::RCP<Epetra_FEVector> force)
{
  // Evaluate the global coupling terms
  evaluate_and_assemble_global_coupling_contributions(data_state->get_dis_col_np());

  // Add the penalty terms to the global force and stiffness matrix
  add_global_force_stiffness_penalty_contributions(data_state, stiff, force);
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::get_global_lambda() const
{
  auto penalty_regularization = get_penalty_regularization(false);
  return std::get<0>(penalty_regularization);
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::get_global_lambda_col() const
{
  Teuchos::RCP<Epetra_Vector> lambda_col = Teuchos::rcp(new Epetra_Vector(*lambda_dof_colmap_));
  Core::LinAlg::Export(*get_global_lambda(), *lambda_col);
  return lambda_col;
}

/**
 *
 */
double BEAMINTERACTION::BeamToSolidMortarManager::get_energy() const
{
  // Since this value is also computed for the reference configuration, where the global mortar
  // matrices are not build yet we return 0 in this case.
  if (not constraint_lin_beam_->filled() or not constraint_lin_solid_->filled() or
      not force_beam_lin_lambda_->filled() or not force_solid_lin_lambda_->filled())
    return 0.0;

  // Calculate the penalty potential.
  Teuchos::RCP<Epetra_Vector> lambda = get_global_lambda();
  double dot_product = 0.0;
  constraint_->Dot(*lambda, &dot_product);
  return 0.5 * dot_product;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::evaluate_and_assemble_global_coupling_contributions(
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  check_setup();
  check_global_maps();

  // Reset the global data structures.
  constraint_->PutScalar(0.);
  constraint_lin_beam_->put_scalar(0.);
  constraint_lin_solid_->put_scalar(0.);
  force_beam_lin_lambda_->put_scalar(0.);
  force_solid_lin_lambda_->put_scalar(0.);
  kappa_->PutScalar(0.);
  kappa_lin_beam_->put_scalar(0.);
  kappa_lin_solid_->put_scalar(0.);
  lambda_active_->PutScalar(0.);

  for (auto& elepairptr : contact_pairs_)
  {
    // Evaluate the mortar contributions of the pair and the pair assembles the terms into the
    // global matrices.
    elepairptr->evaluate_and_assemble_mortar_contributions(*discret_, this, *constraint_lin_beam_,
        *constraint_lin_solid_, *force_beam_lin_lambda_, *force_solid_lin_lambda_, *constraint_,
        *kappa_, *kappa_lin_beam_, *kappa_lin_solid_, *lambda_active_, displacement_vector);
  }

  // Complete the global mortar matrices.
  constraint_lin_beam_->complete(*beam_dof_rowmap_, *lambda_dof_rowmap_);
  constraint_lin_solid_->complete(*solid_dof_rowmap_, *lambda_dof_rowmap_);
  force_beam_lin_lambda_->complete(*lambda_dof_rowmap_, *beam_dof_rowmap_);
  force_solid_lin_lambda_->complete(*lambda_dof_rowmap_, *solid_dof_rowmap_);
  kappa_lin_beam_->complete(*beam_dof_rowmap_, *lambda_dof_rowmap_);
  kappa_lin_solid_->complete(*solid_dof_rowmap_, *lambda_dof_rowmap_);

  // Complete the global scaling vector.
  if (0 != kappa_->GlobalAssemble(Add, false))
    FOUR_C_THROW("Failed to perform FE assembly of kappa_.");
  if (0 != lambda_active_->GlobalAssemble(Add, false))
    FOUR_C_THROW("Failed to perform FE assembly of lambda_active_.");
  if (0 != constraint_->GlobalAssemble(Add, false))
    FOUR_C_THROW("Failed to perform FE assembly of constraint_.");
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::add_global_force_stiffness_penalty_contributions(
    const Teuchos::RCP<const Solid::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff, Teuchos::RCP<Epetra_FEVector> force) const
{
  check_setup();
  check_global_maps();

  // Get the penalty regularization
  const bool is_stiff = stiff != Teuchos::null;
  auto penalty_regularization = get_penalty_regularization(is_stiff);
  const auto lambda = std::get<0>(penalty_regularization);

  int linalg_error = 0;

  if (is_stiff)
  {
    // Penalty regularization linearized w.r.t. the constraint equations
    auto penalty_regularization_lin_constaint =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*std::get<1>(penalty_regularization)));
    penalty_regularization_lin_constaint->complete();
    auto regularized_constraint_lin_beam =
        Core::LinAlg::MLMultiply(*penalty_regularization_lin_constaint, false,
            *constraint_lin_beam_, false, false, false, true);
    auto regularized_constraint_lin_solid =
        Core::LinAlg::MLMultiply(*penalty_regularization_lin_constaint, false,
            *constraint_lin_solid_, false, false, false, true);

    // Penalty regularization linearized w.r.t. the scaling vector
    if (kappa_lin_beam_->NormInf() > 1e-12 && kappa_lin_solid_->NormInf() > 1e-12)
    {
      auto penalty_regularization_lin_kappa =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*std::get<2>(penalty_regularization)));
      penalty_regularization_lin_kappa->complete();
      const auto kappa_lin_beam_scaled = Core::LinAlg::MLMultiply(
          *penalty_regularization_lin_kappa, false, *kappa_lin_beam_, false, false, false, true);
      const auto kappa_lin_solid_scaled = Core::LinAlg::MLMultiply(
          *penalty_regularization_lin_kappa, false, *kappa_lin_solid_, false, false, false, true);
      regularized_constraint_lin_beam->add(*kappa_lin_beam_scaled, false, 1.0, 1.0);
      regularized_constraint_lin_solid->add(*kappa_lin_solid_scaled, false, 1.0, 1.0);
    }

    // Calculate the needed submatrices
    const auto force_beam_lin_lambda_times_constaint_lin_beam =
        Core::LinAlg::MLMultiply(*force_beam_lin_lambda_, false, *regularized_constraint_lin_beam,
            false, false, false, true);
    const auto force_beam_lin_lambda_times_constaint_lin_solid =
        Core::LinAlg::MLMultiply(*force_beam_lin_lambda_, false, *regularized_constraint_lin_solid,
            false, false, false, true);
    const auto force_solid_lin_lambda_times_constaint_lin_beam =
        Core::LinAlg::MLMultiply(*force_solid_lin_lambda_, false, *regularized_constraint_lin_beam,
            false, false, false, true);
    const auto force_solid_lin_lambda_times_constaint_lin_solid =
        Core::LinAlg::MLMultiply(*force_solid_lin_lambda_, false, *regularized_constraint_lin_solid,
            false, false, false, true);

    // Add contributions to the global stiffness matrix
    stiff->add(*force_beam_lin_lambda_times_constaint_lin_beam, false, 1.0, 1.0);
    stiff->add(*force_beam_lin_lambda_times_constaint_lin_solid, false, 1.0, 1.0);
    stiff->add(*force_solid_lin_lambda_times_constaint_lin_beam, false, 1.0, 1.0);
    stiff->add(*force_solid_lin_lambda_times_constaint_lin_solid, false, 1.0, 1.0);
  }

  if (force != Teuchos::null)
  {
    // Factor for right hand side (forces). 1 corresponds to the mesh-tying forces being added to
    // the right hand side, -1 to the left hand side.
    const double rhs_factor = -1.0;

    // Multiply the lambda vector with FB_L and FS_L to get the forces on the beam and solid,
    // respectively.
    Teuchos::RCP<Epetra_Vector> beam_force = Teuchos::rcp(new Epetra_Vector(*beam_dof_rowmap_));
    Teuchos::RCP<Epetra_Vector> solid_force = Teuchos::rcp(new Epetra_Vector(*solid_dof_rowmap_));
    beam_force->PutScalar(0.);
    solid_force->PutScalar(0.);
    linalg_error = force_beam_lin_lambda_->multiply(false, *lambda, *beam_force);
    if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");
    linalg_error = force_solid_lin_lambda_->multiply(false, *lambda, *solid_force);
    if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");
    Teuchos::RCP<Epetra_Vector> global_temp =
        Teuchos::rcp(new Epetra_Vector(*discret_->dof_row_map()));
    Core::LinAlg::Export(*beam_force, *global_temp);
    Core::LinAlg::Export(*solid_force, *global_temp);

    // Add force contributions to global vector.
    linalg_error = force->Update(-1.0 * rhs_factor, *global_temp, 1.0);
    if (linalg_error != 0) FOUR_C_THROW("Error in Update");
  }

  // Add the force and stiffness contributions that are assembled directly by the pairs.
  auto lambda_col = Teuchos::rcp(new Epetra_Vector(*lambda_dof_colmap_));
  Core::LinAlg::Export(*lambda, *lambda_col);
  for (const auto& elepairptr : contact_pairs_)
    elepairptr->evaluate_and_assemble(
        *discret_, this, force, stiff, *lambda_col, *data_state->get_dis_col_np());
}

/**
 *
 */
std::tuple<Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_Vector>>
BEAMINTERACTION::BeamToSolidMortarManager::get_penalty_regularization(
    const bool compute_linearization) const
{
  check_setup();
  check_global_maps();

  // Get the inverted kappa matrix.
  Teuchos::RCP<Epetra_Vector> penalty_kappa_inv = penalty_invert_kappa();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> penalty_kappa_inv_mat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*penalty_kappa_inv));
  penalty_kappa_inv_mat->complete();

  // Multiply the inverted kappa matrix with the constraint equations.
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
  int linalg_error = penalty_kappa_inv_mat->multiply(false, *constraint_, *lambda);
  if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");

  if (compute_linearization)
  {
    return {lambda, penalty_kappa_inv, Teuchos::null};
  }
  else
  {
    return {lambda, Teuchos::null, Teuchos::null};
  }
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector> BEAMINTERACTION::BeamToSolidMortarManager::penalty_invert_kappa() const
{
  // Create the inverse vector.
  Teuchos::RCP<Epetra_Vector> kappa_inv = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));

  // Get the penalty parameters.
  const double penalty_translation = beam_to_solid_params_->get_penalty_parameter();
  double penalty_rotation = 0.0;
  auto beam_to_volume_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>(
          beam_to_solid_params_);
  auto beam_to_surface_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>(
          beam_to_solid_params_);
  if (beam_to_volume_params != Teuchos::null)
    penalty_rotation = beam_to_volume_params->get_rotational_coupling_penalty_parameter();
  else if (beam_to_surface_params != Teuchos::null)
    penalty_rotation = beam_to_surface_params->get_rotational_coupling_penalty_parameter();
  else if (lambda_dof_rowmap_rotations_->NumGlobalElements() > 0)
    FOUR_C_THROW(
        "Rotational penalty coupling only implemented for beam-to-volume and beam-to-surface "
        "case.");

  // Calculate the local inverse of kappa.
  double penalty = 0.0;
  double local_kappa_inv_value = 0.;
  for (int lid = 0; lid < lambda_dof_rowmap_->NumMyElements(); lid++)
  {
    if (lambda_active_->Values()[lid] > 0.1)
    {
      const int gid = lambda_dof_rowmap_->GID(lid);
      if (lambda_dof_rowmap_translations_->LID(gid) != -1)
        penalty = penalty_translation;
      else if (lambda_dof_rowmap_rotations_->LID(gid) != -1)
        penalty = penalty_rotation;
      else
        FOUR_C_THROW("Could not find the GID %d in translation or rotation map", gid);

      local_kappa_inv_value = penalty / kappa_->Values()[lid];
    }

    else
      // This LID is inactive.
      local_kappa_inv_value = 0.0;

    kappa_inv->ReplaceMyValue(lid, 0, local_kappa_inv_value);
  }

  return kappa_inv;
}

FOUR_C_NAMESPACE_CLOSE
