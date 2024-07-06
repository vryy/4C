/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the beam-to-solid visualization output writers.

\level 3

*/


#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_utils.hpp"

#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
void BEAMINTERACTION::AddBeamInteractionNodalForces(
    const Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>& visualization,
    const Teuchos::RCP<const Core::FE::Discretization>& discret_ptr,
    const Teuchos::RCP<const Epetra_MultiVector>& displacement,
    const Teuchos::RCP<const Epetra_MultiVector>& force, const bool write_unique_ids)
{
  // Add the reference geometry and displacement to the visualization.
  visualization->add_discretization_nodal_reference_position(discret_ptr);
  visualization->add_discretization_nodal_data("displacement", displacement);

  // Create maps with the GIDs of beam and solid nodes.
  std::vector<int> gid_beam_dof;
  std::vector<int> gid_solid_dof;
  std::vector<int> gid_node;
  for (int i_lid = 0; i_lid < discret_ptr->num_my_row_nodes(); i_lid++)
  {
    gid_node.clear();
    Core::Nodes::Node* current_node = discret_ptr->l_row_node(i_lid);
    discret_ptr->dof(current_node, gid_node);
    if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
      for (unsigned int dim = 0; dim < 3; ++dim) gid_beam_dof.push_back(gid_node[dim]);
    else
      for (unsigned int dim = 0; dim < 3; ++dim) gid_solid_dof.push_back(gid_node[dim]);
  }
  Epetra_Map beam_dof_map(-1, gid_beam_dof.size(), gid_beam_dof.data(), 0, discret_ptr->get_comm());
  Epetra_Map solid_dof_map(
      -1, gid_solid_dof.size(), gid_solid_dof.data(), 0, discret_ptr->get_comm());

  // Extract the forces and add them to the discretization.
  Teuchos::RCP<Epetra_Vector> force_beam =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(beam_dof_map, true));
  Core::LinAlg::Export(*force, *force_beam);
  Teuchos::RCP<Epetra_Vector> force_solid =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(solid_dof_map, true));
  Core::LinAlg::Export(*force, *force_solid);
  visualization->add_discretization_nodal_data("force_beam", force_beam);
  visualization->add_discretization_nodal_data("force_solid", force_solid);

  if (write_unique_ids)
  {
    auto& visualization_data = visualization->get_visualization_data();
    std::vector<double>& unique_id = visualization_data.get_point_data<double>("uid_0_node_id");
    for (int i_lid = 0; i_lid < discret_ptr->num_my_row_nodes(); i_lid++)
      unique_id.push_back(discret_ptr->l_row_node(i_lid)->id());
  }
}

/**
 *
 */
void BEAMINTERACTION::AddAveragedNodalNormals(
    const Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>&
        output_writer_base_ptr,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements,
    const int condition_coupling_id, const bool write_unique_ids)
{
  // Get the visualization vectors.
  auto& visualization_data = output_writer_base_ptr->get_visualization_data();
  std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
  std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");
  std::vector<double>& normal_averaged =
      visualization_data.get_point_data<double>("normal_averaged");
  std::vector<double>& normal_element = visualization_data.get_point_data<double>("normal_element");
  std::vector<double>& coupling_id = visualization_data.get_point_data<double>("coupling_id");

  std::vector<double>* face_id = nullptr;
  if (write_unique_ids)
  {
    face_id = &(visualization_data.get_point_data<double>("uid_0_face_id"));
  }

  // Loop over face elements.
  for (auto const& face_element_iterator : face_elements)
  {
    // Only write the output for the faces that are part of a pair, since otherwise there are faces
    // which have empty position arrays since set_state was never called on them.
    if (face_element_iterator.second->is_part_of_pair())
    {
      // Setup variables.
      Core::LinAlg::Matrix<3, 1, double> X, u, r, n, n_averaged;

      // Set the element parameter coordinates.
      Core::LinAlg::Matrix<2, 1, double> xi(true);
      Core::LinAlg::SerialDenseMatrix nodal_coordinates =
          Core::FE::getEleNodeNumbering_nodes_paramspace(
              face_element_iterator.second->get_drt_face_element()->shape());

      // Loop over element nodes.
      for (int i_node = 0;
           i_node < face_element_iterator.second->get_drt_face_element()->num_node(); i_node++)
      {
        for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
          xi(i_dim) = nodal_coordinates(i_dim, i_node);
        face_element_iterator.second->evaluate_face_position_double(xi, r);
        face_element_iterator.second->evaluate_face_position_double(xi, X, true);
        face_element_iterator.second->evaluate_face_normal_double(xi, n, false, false);
        face_element_iterator.second->evaluate_face_normal_double(xi, n_averaged, false, true);

        u = r;
        u -= X;

        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(X(dim));
          displacement.push_back(u(dim));
          normal_averaged.push_back(n_averaged(dim));
          normal_element.push_back(n(dim));
        }
        coupling_id.push_back(condition_coupling_id);

        if (write_unique_ids)
          face_id->push_back(face_element_iterator.second->get_drt_face_element()->id());
      }
    }
  }
}

/**
 *
 */
void BEAMINTERACTION::GetGlobalCouplingForceResultants(const Core::FE::Discretization& discret,
    const Epetra_MultiVector& force, const Epetra_MultiVector& displacement,
    Core::LinAlg::Matrix<3, 2, double>& beam_resultant,
    Core::LinAlg::Matrix<3, 2, double>& solid_resultant)
{
  beam_resultant.clear();
  solid_resultant.clear();

  // Initialize variables.
  std::vector<double> local_force;
  std::vector<double> local_position;
  std::vector<int> gid_node;

  // Loop over the nodes and sum up the forces and moments.
  for (int i_lid = 0; i_lid < discret.num_my_row_nodes(); i_lid++)
  {
    gid_node.clear();
    Core::Nodes::Node* current_node = discret.l_row_node(i_lid);
    discret.dof(current_node, gid_node);

    // Get the local force and displacement values.
    Core::FE::ExtractMyValues(force, local_force, gid_node);
    Core::FE::ExtractMyValues(displacement, local_position, gid_node);
    for (unsigned int dim = 0; dim < 3; ++dim) local_position[dim] += current_node->x()[dim];

    if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
    {
      if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(*current_node))
        GetNodeCouplingForceResultants(local_force, local_position, beam_resultant);
      else
        // Do nothing for non-centerline nodes.
        continue;
    }
    else
      GetNodeCouplingForceResultants(local_force, local_position, solid_resultant);
  }
}

/**
 *
 */
void BEAMINTERACTION::GetNodeCouplingForceResultants(const std::vector<double>& local_force,
    const std::vector<double>& local_position, Core::LinAlg::Matrix<3, 2, double>& resultant)
{
  Core::LinAlg::Matrix<3, 1, double> node_pos(true);
  Core::LinAlg::Matrix<3, 1, double> node_force(true);
  Core::LinAlg::Matrix<3, 1, double> node_moment(true);

  // Add the force values for this node to the resultants and get the local node position.
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    node_pos(dim) = local_position[dim];
    node_force(dim) = local_force[dim];
    resultant(dim, 0) += node_force(dim);
  }

  // Add the moment values for this node.
  node_moment.cross_product(node_pos, node_force);
  for (unsigned int dim = 0; dim < 3; ++dim) resultant(dim, 1) += node_moment(dim);
}

FOUR_C_NAMESPACE_CLOSE
