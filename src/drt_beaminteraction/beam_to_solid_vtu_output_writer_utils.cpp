/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the beam-to-solid vtu output writers.

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_vtu_output_writer_utils.H"

#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beaminteraction_calc_utils.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


/**
 *
 */
void BEAMINTERACTION::AddBeamInteractionNodalForces(
    const Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>& visualization,
    const Teuchos::RCP<const DRT::Discretization>& discret_ptr,
    const Teuchos::RCP<const Epetra_MultiVector>& displacement,
    const Teuchos::RCP<const Epetra_MultiVector>& force)
{
  // Add the reference geometry and displacement to the visualization.
  visualization->AddDiscretizationNodalReferencePosition(discret_ptr);
  visualization->AddDiscretizationNodalData("displacement", displacement);

  // Create maps with the GIDs of beam and solid nodes.
  std::vector<int> gid_beam_dof;
  std::vector<int> gid_solid_dof;
  std::vector<int> gid_node;
  for (int i_lid = 0; i_lid < discret_ptr->NumMyRowNodes(); i_lid++)
  {
    gid_node.clear();
    DRT::Node* current_node = discret_ptr->lRowNode(i_lid);
    discret_ptr->Dof(current_node, gid_node);
    if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
      for (unsigned int dim = 0; dim < 3; ++dim) gid_beam_dof.push_back(gid_node[dim]);
    else
      for (unsigned int dim = 0; dim < 3; ++dim) gid_solid_dof.push_back(gid_node[dim]);
  }
  Epetra_Map beam_dof_map(-1, gid_beam_dof.size(), &gid_beam_dof[0], 0, discret_ptr->Comm());
  Epetra_Map solid_dof_map(-1, gid_solid_dof.size(), &gid_solid_dof[0], 0, discret_ptr->Comm());

  // Extract the forces and add them to the discretization.
  Teuchos::RCP<Epetra_Vector> force_beam =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(beam_dof_map, true));
  LINALG::Export(*force, *force_beam);
  Teuchos::RCP<Epetra_Vector> force_solid =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(solid_dof_map, true));
  LINALG::Export(*force, *force_solid);
  visualization->AddDiscretizationNodalData("force_beam", force_beam);
  visualization->AddDiscretizationNodalData("force_solid", force_solid);
}

/**
 *
 */
void BEAMINTERACTION::AddAveragedNodalNormals(
    const Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>&
        output_writer_base_ptr,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements,
    const int condition_coupling_id)
{
  // Get the visualization vectors.
  std::vector<double>& point_coordinates =
      output_writer_base_ptr->GetMutablePointCoordinateVector();
  std::vector<double>& displacement =
      output_writer_base_ptr->GetMutablePointDataVector("displacement");
  std::vector<double>& normal_averaged =
      output_writer_base_ptr->GetMutablePointDataVector("normal_averaged");
  std::vector<double>& normal_element =
      output_writer_base_ptr->GetMutablePointDataVector("normal_element");
  std::vector<double>& coupling_id =
      output_writer_base_ptr->GetMutablePointDataVector("coupling_id");

  // Loop over face elements.
  for (auto const& face_element_iterator : face_elements)
  {
    // Setup variables.
    LINALG::Matrix<3, 1, double> X, u, r, n, n_averaged;

    // Set the element parameter coordinates.
    LINALG::Matrix<2, 1, double> xi(true);
    LINALG::SerialDenseMatrix nodal_coordinates = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(
        face_element_iterator.second->GetDrtFaceElement()->Shape());

    // Loop over element nodes.
    for (int i_node = 0; i_node < face_element_iterator.second->GetDrtFaceElement()->NumNode();
         i_node++)
    {
      for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, i_node);
      face_element_iterator.second->EvaluateFacePositionDouble(xi, r);
      face_element_iterator.second->EvaluateFacePositionDouble(xi, X, true);
      face_element_iterator.second->EvaluateFaceNormalDouble(xi, n);
      face_element_iterator.second->EvaluateFaceAveragedNormalDouble(xi, n_averaged);

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
    }
  }
}
