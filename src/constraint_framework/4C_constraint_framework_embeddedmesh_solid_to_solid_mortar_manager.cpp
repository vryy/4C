/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between solids.

\level 3
*/

#include "4C_constraint_framework_embeddedmesh_solid_to_solid_mortar_manager.hpp"

#include "4C_constraint_framework_embeddedmesh_interaction_pair.hpp"
#include "4C_constraint_framework_embeddedmesh_params.hpp"
#include "4C_constraint_framework_embeddedmesh_solid_to_solid_utils.hpp"
#include "4C_cut_cutwizard.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

#include <unordered_set>

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::SolidToSolidMortarManager(
    Teuchos::RCP<Core::FE::Discretization>& discret, const Epetra_Vector& displacement_vector,
    CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams& embedded_mesh_coupling_params,
    Teuchos::RCP<Core::IO::VisualizationManager> visualization_manager, int start_value_lambda_gid)
    : discret_(discret),
      start_value_lambda_gid_(start_value_lambda_gid),
      embedded_mesh_coupling_params_(embedded_mesh_coupling_params),
      point_visualization_manager_(visualization_manager),
      lambda_visualization_manager_(visualization_manager)
{
  // Initialize cutwizard instance and perform the cut
  Teuchos::RCP<Cut::CutWizard> cutwizard = Teuchos::rcp(new Cut::CutWizard(discret_));
  CONSTRAINTS::EMBEDDEDMESH::prepare_and_perform_cut(
      cutwizard, discret_, embedded_mesh_coupling_params_);

  // Get the coupling pairs and cut elements
  get_coupling_pairs_and_background_elements(cutwizard, embedded_mesh_coupling_params_, discret_,
      embedded_mesh_solid_pairs_, cut_elements_vector_);

  // Change integration rule of elements if they are cut
  CONSTRAINTS::EMBEDDEDMESH::change_gauss_rule_of_cut_elements(cut_elements_vector_, cutwizard);

  // Get the number of Lagrange multiplier DOF on a solid node and on a solid element
  unsigned int n_lambda_node_temp = 0;
  mortar_shape_functions_to_number_of_lagrange_values(
      embedded_mesh_coupling_params_.embedded_mesh_mortar_shape_function_, n_lambda_node_temp);

  n_lambda_node_ = n_lambda_node_temp;

  point_visualization_manager_->register_visualization_data("background_integration_points");
  point_visualization_manager_->register_visualization_data("interface_integration_points");
  point_visualization_manager_->register_visualization_data("cut_element_integration_points");
  lambda_visualization_manager_->register_visualization_data("lagrange_multipliers");

  auto& background_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("background_integration_points");
  background_integration_points_visualization_data.register_point_data<double>("weights", 1);
  background_integration_points_visualization_data.register_point_data<int>(
      "integration_cell_id", 1);

  auto& interface_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("interface_integration_points");
  interface_integration_points_visualization_data.register_point_data<double>("weights", 1);
  interface_integration_points_visualization_data.register_point_data<int>(
      "integration_cell_id", 1);

  auto& cut_element_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("cut_element_integration_points");
  cut_element_integration_points_visualization_data.register_point_data<double>("weights", 1);
  cut_element_integration_points_visualization_data.register_point_data<int>(
      "integration_cell_id", 1);

  auto& lagrange_multipliers_visualization_data =
      lambda_visualization_manager_->get_visualization_data("lagrange_multipliers");
  lagrange_multipliers_visualization_data.register_point_data<double>("lambda", 3);

  // Setup the solid to solid mortar manager
  setup(displacement_vector);
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::set_state(
    const Epetra_Vector& displacement_vector)
{
  // Check if the coupling pairs are empty, if thats the case, return dserror
  if (embedded_mesh_solid_pairs_.size() == 0)
    FOUR_C_THROW("We cannot set the state of the coupling pairs if they are not defined yet.");

  for (auto couplig_pair_iter : embedded_mesh_solid_pairs_)
    couplig_pair_iter->set_current_element_position(*discret_, displacement_vector);
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::setup(
    const Epetra_Vector& displacement_vector)
{
  // Get the global ids of all mesh nodes on this rank
  std::vector<int> my_nodes_gid;
  for (int i_node = 0; i_node < discret_->node_row_map()->NumMyElements(); i_node++)
  {
    Core::Nodes::Node const& node = *(discret_->l_row_node(i_node));
    if (CONSTRAINTS::EMBEDDEDMESH::is_interface_node(node)) my_nodes_gid.push_back(node.id());
  }

  // Calculate the local number of interface nodes
  const unsigned int n_nodes = my_nodes_gid.size();
  const unsigned int n_lambda_dof = n_nodes * n_lambda_node_;

  // Tell all other processors how many lambda DOFs this processor has. This information is needed
  // to construct the lambda_dof_rowmap_.
  std::vector<int> lambda_dof_per_rank(discret_->get_comm().NumProc(), 0);
  int temp_my_n_lambda_dof = (int)n_lambda_dof;
  discret_->get_comm().GatherAll(&temp_my_n_lambda_dof, &lambda_dof_per_rank[0], 1);

  // Get the start GID for the lambda DOFs on this processor.
  int my_lambda_gid_start_value = start_value_lambda_gid_;
  for (int pid = 0; pid < discret_->get_comm().MyPID(); pid++)
    my_lambda_gid_start_value += lambda_dof_per_rank[pid];

  // Fill in all GIDs of the lambda DOFs on this processor
  std::vector<int> my_lambda_gid;
  my_lambda_gid.reserve(n_nodes * n_lambda_node_);
  for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
  {
    for (unsigned int i_dof = 0; i_dof < n_lambda_node_; i_dof++)
      my_lambda_gid.push_back(my_lambda_gid_start_value + i_dof + i_node * n_lambda_node_);
  }

  // Rowmap for the additional GIDs used by the mortar contact discretization.
  lambda_dof_rowmap_ = Teuchos::rcp(
      new Epetra_Map(-1, my_lambda_gid.size(), my_lambda_gid.data(), 0, discret_->get_comm()));

  // We need to be able to get the global ids for a Lagrange multiplier DOF from the global id
  // of a node. To do so, we 'abuse' the Epetra_MultiVector as map between the
  // global node ids and the global Lagrange multiplier DOF ids.
  Teuchos::RCP<Epetra_Map> node_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_nodes, &my_nodes_gid[0], 0, discret_->get_comm()));

  // Map from global node ids to global lagrange multiplier ids. Only create the
  // multivector if it has one or more columns.
  node_gid_to_lambda_gid_ = Teuchos::null;
  if (n_lambda_node_ > 0)
    node_gid_to_lambda_gid_ =
        Teuchos::rcp(new Epetra_MultiVector(*node_gid_rowmap, n_lambda_node_, true));

  // Fill in the entries in the node global id to Lagrange multiplier global id vector.
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

  // Create the maps for boundary layer and background DOFs..
  set_global_maps();

  // Create the global coupling matrices.
  global_constraint_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  global_g_bl_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 30, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  global_g_bg_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lambda_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  global_fbl_l_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*boundary_layer_interface_dof_rowmap_,
      30, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  global_fbg_l_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *background_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  global_kappa_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));
  global_active_lambda_ = Teuchos::rcp(new Epetra_FEVector(*lambda_dof_rowmap_));

  // Set flag for successful setup.
  is_setup_ = true;
  is_local_maps_build_ = false;

  // Set the local maps
  set_local_maps(displacement_vector);
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::set_global_maps()
{
  // Loop over all nodes on this processor
  std::vector<int> boundary_layer_interface_dofs(0);
  std::vector<int> background_dofs(0);
  for (int i_node = 0; i_node < discret_->node_row_map()->NumMyElements(); i_node++)
  {
    const Core::Nodes::Node* node = discret_->l_row_node(i_node);
    if (is_cut_node(*node))
      discret_->dof(node, background_dofs);
    else if (CONSTRAINTS::EMBEDDEDMESH::is_interface_node(*node))
      discret_->dof(node, boundary_layer_interface_dofs);
  }

  // Create the beam and solid maps.
  boundary_layer_interface_dof_rowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, boundary_layer_interface_dofs.size(),
          &boundary_layer_interface_dofs[0], 0, discret_->get_comm()));
  background_dof_rowmap_ = Teuchos::rcp(
      new Epetra_Map(-1, background_dofs.size(), &background_dofs[0], 0, discret_->get_comm()));

  // Reset the local maps.
  node_gid_to_lambda_gid_map_.clear();

  // Set flags for global maps.
  is_global_maps_build_ = true;
  is_local_maps_build_ = false;
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::set_local_maps(
    const Epetra_Vector& displacement_vector)
{
  check_setup();
  check_global_maps();

  // Set the displacement state of the coupling pairs
  set_state(displacement_vector);

  // At this point the global multi vectors are filled up completely. To get the map for global
  // node ids to the global lambda ids we need to be able to extract more than the local
  // values on this processor. Therefore we need a new map that contains all rows we want to
  // access in the global multi vector.
  std::vector<int> node_gid_needed;

  // Loop over the pairs and get the global node and element indices needed on this rank.
  for (unsigned int i_pair = 0; i_pair < embedded_mesh_solid_pairs_.size(); i_pair++)
  {
    const Teuchos::RCP<EMBEDDEDMESH::SolidInteractionPair>& pair =
        embedded_mesh_solid_pairs_[i_pair];

    // The first (beam) element should always be on the same processor as the pair.
    if (pair->element_1().owner() != discret_->get_comm().MyPID())
      FOUR_C_THROW(
          "The current implementation need the first element of a interface coupling pair to be on "
          "the "
          "same processor as the pair!");

    // Get the global id of the nodes that the pairs on this rank need.
    if (n_lambda_node_ > 0)
      for (int i_node = 0; i_node < pair->element_1().num_node(); i_node++)
        node_gid_needed.push_back(pair->element_1().nodes()[i_node]->id());
  }

  std::vector<int> node_gid_needed_copy = node_gid_needed;

  // Make the entries in the vectors unique.
  std::vector<int>::iterator it;
  std::sort(node_gid_needed.begin(), node_gid_needed.end());
  it = std::unique(node_gid_needed.begin(), node_gid_needed.end());
  node_gid_needed.resize(std::distance(node_gid_needed.begin(), it));

  // Create the maps for the extraction of the values.
  Teuchos::RCP<Epetra_Map> node_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, node_gid_needed.size(), &node_gid_needed[0], 0, discret_->get_comm()));

  // Create the Multivectors that will be filled with all values needed on this rank.
  Teuchos::RCP<Epetra_MultiVector> node_gid_to_lambda_gid_copy = Teuchos::null;
  if (node_gid_to_lambda_gid_ != Teuchos::null)
    node_gid_to_lambda_gid_copy = Teuchos::rcp<Epetra_MultiVector>(
        new Epetra_MultiVector(*node_gid_needed_rowmap, n_lambda_node_, true));

  // Export values from the global multi vector to the ones needed on this rank.
  if (node_gid_to_lambda_gid_ != Teuchos::null)
    Core::LinAlg::export_to(*node_gid_to_lambda_gid_, *node_gid_to_lambda_gid_copy);

  // Fill in the local maps.
  std::vector<int> lambda_gid_for_col_map;
  lambda_gid_for_col_map.clear();
  node_gid_to_lambda_gid_map_.clear();
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

  // Create the global lambda col map.
  lambda_dof_colmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, lambda_gid_for_col_map.size(), &lambda_gid_for_col_map[0], 0, discret_->get_comm()));

  // Set flags for local maps.
  is_local_maps_build_ = true;
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::location_vector(
    const CONSTRAINTS::EMBEDDEDMESH::SolidInteractionPair* interaction_pair,
    std::vector<int>& lambda_row) const
{
  check_setup();
  check_local_maps();

  // Clear the output vectors.
  lambda_row.clear();

  // Get the global DOFs ids of the nodal Lagrange multipliers.
  if (n_lambda_node_ > 0)
  {
    for (int i_node = 0; i_node < interaction_pair->element_1().num_node(); i_node++)
    {
      const Core::Nodes::Node& node = *(interaction_pair->element_1().nodes()[i_node]);
      // Get the global id of the node.
      int node_id = node.id();

      // Check if the id is in the map. If it is, add it to the output vector.
      auto search_key_in_map = node_gid_to_lambda_gid_map_.find(node_id);
      if (search_key_in_map == node_gid_to_lambda_gid_map_.end())
        FOUR_C_THROW("Global node id %d not in map!", node_id);
      for (auto const& lambda_gid : search_key_in_map->second) lambda_row.push_back(lambda_gid);
    }
  }
}

/**
 *
 */
void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::evaluate_global_coupling_contributions(
    const Epetra_Vector& displacement_vector)
{
  check_setup();
  check_global_maps();

  set_state(displacement_vector);

  // Reset the global data structures.
  global_constraint_->PutScalar(0.);
  global_g_bl_->put_scalar(0.);
  global_g_bg_->put_scalar(0.);
  global_fbl_l_->put_scalar(0.);
  global_fbg_l_->put_scalar(0.);
  global_kappa_->PutScalar(0.);
  global_active_lambda_->PutScalar(0.);

  for (auto& elepairptr : embedded_mesh_solid_pairs_)
  {
    elepairptr->evaluate_and_assemble_mortar_contributions(*discret_, this, *global_g_bl_,
        *global_g_bg_, *global_fbl_l_, *global_fbg_l_, *global_constraint_, *global_kappa_,
        *global_active_lambda_);
  }

  // Complete the global mortar matrices.
  global_g_bl_->complete(
      *boundary_layer_interface_dof_rowmap_, *lambda_dof_rowmap_);  // is the complete wrong?
  global_g_bg_->complete(*background_dof_rowmap_, *lambda_dof_rowmap_);
  global_fbl_l_->complete(*lambda_dof_rowmap_, *boundary_layer_interface_dof_rowmap_);
  global_fbg_l_->complete(*lambda_dof_rowmap_, *background_dof_rowmap_);

  // Complete the global scaling vector.
  if (0 != global_kappa_->GlobalAssemble(Add, false)) FOUR_C_THROW("Error in GlobalAssemble!");
  if (0 != global_active_lambda_->GlobalAssemble(Add, false))
    FOUR_C_THROW("Error in GlobalAssemble!");
  if (0 != global_constraint_->GlobalAssemble(Add, false)) FOUR_C_THROW("Error in GlobalAssemble!");
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::
    add_global_force_stiffness_penalty_contributions(
        const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& data_state,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff, Teuchos::RCP<Epetra_Vector> force) const
{
  check_setup();
  check_global_maps();

  int linalg_error = 0;

  if (stiff != Teuchos::null)
  {
    // Scale the linearizations of the constraint equations.
    Teuchos::RCP<Epetra_Vector> global_penalty_kappa_inv = penalty_invert_kappa();
    Teuchos::RCP<Core::LinAlg::SparseMatrix> penalty_kappa_inv_mat =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*global_penalty_kappa_inv));
    penalty_kappa_inv_mat->complete();

    Teuchos::RCP<Core::LinAlg::SparseMatrix> global_G_BL_scaled = Core::LinAlg::MLMultiply(
        *penalty_kappa_inv_mat, false, *global_g_bl_, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> global_G_BG_scaled = Core::LinAlg::MLMultiply(
        *penalty_kappa_inv_mat, false, *global_g_bg_, false, false, false, true);

    // Calculate the needed submatrices.
    Teuchos::RCP<Core::LinAlg::SparseMatrix> FBL_L_times_G_BL = Core::LinAlg::MLMultiply(
        *global_fbl_l_, false, *global_G_BL_scaled, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> FBL_L_times_G_BG = Core::LinAlg::MLMultiply(
        *global_fbl_l_, false, *global_G_BG_scaled, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> FBG_L_times_G_BL = Core::LinAlg::MLMultiply(
        *global_fbg_l_, false, *global_G_BL_scaled, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> FBG_L_times_G_BG = Core::LinAlg::MLMultiply(
        *global_fbg_l_, false, *global_G_BG_scaled, false, false, false, true);

    // Add contributions to the global stiffness matrix.
    stiff->add(*FBL_L_times_G_BL, false, 1.0, 1.0);
    stiff->add(*FBL_L_times_G_BG, false, 1.0, 1.0);
    stiff->add(*FBG_L_times_G_BL, false, 1.0, 1.0);
    stiff->add(*FBG_L_times_G_BG, false, 1.0, 1.0);
  }

  if (force != Teuchos::null)
  {
    // Factor for right hand side (forces). 1 corresponds to the meshtying forces being added to the
    // right hand side, -1 to the left hand side.
    const double rhs_factor = 1.0;

    // Get the penalty Lagrange multiplier vector.
    Teuchos::RCP<Epetra_Vector> lambda = get_global_lambda();

    // Multiply the lambda vector with FBL_L and FBG_L to get the forces on the boundary layer and
    // background, respectively.
    Teuchos::RCP<Epetra_Vector> boundary_layer_interface_force =
        Teuchos::rcp(new Epetra_Vector(*boundary_layer_interface_dof_rowmap_));
    Teuchos::RCP<Epetra_Vector> background_force =
        Teuchos::rcp(new Epetra_Vector(*background_dof_rowmap_));
    boundary_layer_interface_force->PutScalar(0.);
    background_force->PutScalar(0.);
    linalg_error = global_fbl_l_->multiply(false, *lambda, *boundary_layer_interface_force);
    if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");
    linalg_error = global_fbg_l_->multiply(false, *lambda, *background_force);
    if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");
    Teuchos::RCP<Epetra_Vector> global_temp =
        Teuchos::rcp(new Epetra_Vector(*discret_->dof_row_map()));
    Core::LinAlg::export_to(*boundary_layer_interface_force, *global_temp);
    Core::LinAlg::export_to(*background_force, *global_temp);

    // Add force contributions to global vector.
    linalg_error = force->Update(rhs_factor, *global_temp, 1.0);
    if (linalg_error != 0) FOUR_C_THROW("Error in Update");
  }
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector>
CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::get_global_lambda() const
{
  check_setup();
  check_global_maps();

  // Get the inverted kappa matrix.
  Teuchos::RCP<Epetra_Vector> penalty_global_kappa_inv = penalty_invert_kappa();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> penalty_kappa_inv_mat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*penalty_global_kappa_inv));
  penalty_kappa_inv_mat->complete();

  // Multiply the inverted kappa matrix with the constraint equations and scale them with the
  // penalty parameter.
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
  int linalg_error = penalty_kappa_inv_mat->multiply(false, *global_constraint_, *lambda);
  if (linalg_error != 0) FOUR_C_THROW("Error in Multiply!");

  return lambda;
}

/**
 *
 */
Teuchos::RCP<Epetra_Vector>
CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::get_global_lambda_col() const
{
  Teuchos::RCP<Epetra_Vector> lambda_col = Teuchos::rcp(new Epetra_Vector(*lambda_dof_colmap_));
  Core::LinAlg::export_to(*get_global_lambda(), *lambda_col);
  return lambda_col;
}

Teuchos::RCP<Epetra_Vector>
CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::penalty_invert_kappa() const
{
  // Create the inverse vector.
  Teuchos::RCP<Epetra_Vector> global_kappa_inv =
      Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));

  // Get the penalty parameters.
  const double penalty_params =
      embedded_mesh_coupling_params_.embedded_mesh_constraint_penalty_parameter_;

  // Calculate the local inverse of kappa.
  double penalty = 0.0;
  double local_kappa_inv_value = 0.;
  for (int lid = 0; lid < lambda_dof_rowmap_->NumMyElements(); lid++)
  {
    if (global_active_lambda_->Values()[lid] > 0.1)
    {
      const int gid = lambda_dof_rowmap_->GID(lid);
      if (lambda_dof_rowmap_->LID(gid) != -1)
        penalty = penalty_params;
      else
        FOUR_C_THROW("Could not find the GID %d in translation map", gid);

      local_kappa_inv_value = penalty / global_kappa_->Values()[lid];
    }
    else
      // This LID is inactive.
      local_kappa_inv_value = 0.0;

    global_kappa_inv->ReplaceMyValue(lid, 0, local_kappa_inv_value);
  }

  return global_kappa_inv;
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::write_output_lagrange_multipliers(
    double time, int timestep_number)
{
  lambda_visualization_manager_->clear_data();

  auto& lagrange_multipliers_visualization_data =
      lambda_visualization_manager_->get_visualization_data("lagrange_multipliers");

  Teuchos::RCP<Epetra_Vector> lambda = get_global_lambda_col();

  Teuchos::RCP<std::unordered_set<int>> interface_tracker =
      Teuchos::rcp(new std::unordered_set<int>());

  // Loop over pairs
  for (auto& elepairptr : embedded_mesh_solid_pairs_)
  {
    elepairptr->get_pair_visualization(
        lagrange_multipliers_visualization_data, lambda, this, interface_tracker);
  }

  lambda_visualization_manager_->write_to_disk(time, timestep_number);
}

void CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::write_output_integration_points(
    double time, int timestep_number)
{
  point_visualization_manager_->clear_data();

  auto& background_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("background_integration_points");

  auto& interface_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("interface_integration_points");

  auto& cut_element_integration_points_visualization_data =
      point_visualization_manager_->get_visualization_data("cut_element_integration_points");

  // Loop over pairs
  for (auto& elepairptr : embedded_mesh_solid_pairs_)
  {
    unsigned int n_segments = elepairptr->get_num_segments();
    for (size_t iter_segments = 0; iter_segments < n_segments; iter_segments++)
    {
      elepairptr->get_projected_gauss_rule_on_interface(iter_segments,
          background_integration_points_visualization_data,
          interface_integration_points_visualization_data);
    }

    elepairptr->get_projected_gauss_rule_in_cut_element(
        cut_element_integration_points_visualization_data);
  }

  point_visualization_manager_->write_to_disk(time, timestep_number);
}

bool CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager::is_cut_node(
    Core::Nodes::Node const& node)
{
  bool is_cut_node = false;

  // Check if the node belongs to an element that is cut
  for (int num_ele = 0; num_ele < node.num_element(); num_ele++)
  {
    bool is_node_in_cut_ele = std::find(cut_elements_vector_.begin(), cut_elements_vector_.end(),
                                  node.elements()[num_ele]) != cut_elements_vector_.end();
    if (is_node_in_cut_ele) is_cut_node = true;
  }

  return is_cut_node;
}

FOUR_C_NAMESPACE_CLOSE
