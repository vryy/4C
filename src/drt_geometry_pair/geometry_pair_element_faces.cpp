/*----------------------------------------------------------------------*/
/*! \file

\brief Element classes that represent faces, i.e. surface elements.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#include "geometry_pair_element_faces.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


/**
 *
 */
void GEOMETRYPAIR::FaceElement::GetConnectedFaces(
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements,
    std::vector<Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& connected_faces) const
{
  // First we have to find the global IDs of the volume elements connected to the face.
  std::set<int> connected_elements_id;
  for (int i_node = 0; i_node < drt_face_element_->NumNode(); i_node++)
  {
    const DRT::Node* node = drt_face_element_->Nodes()[i_node];
    for (int i_element = 0; i_element < node->NumElement(); i_element++)
    {
      const DRT::Element* node_element = node->Elements()[i_element];
      connected_elements_id.insert(node_element->Id());
    }
  }

  // Check that the global IDs are in the face_element map, i.e. they are in the surface condition.
  connected_faces.clear();
  for (const auto& volume_id : connected_elements_id)
  {
    // This element does not count as a connected face.
    if (volume_id == drt_face_element_->ParentElementId()) continue;

    auto find_in_faces = face_elements.find(volume_id);
    if (find_in_faces != face_elements.end()) connected_faces.push_back(find_in_faces->second);
  }
}

/**
 *
 */
void GEOMETRYPAIR::FaceElement::GetPatchLocalToGlobalIndices(
    const std::vector<Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& connected_faces,
    std::vector<int>& ltg) const
{
  dserror("Has to be implemented!");
  // The fist DOFs will be the ones from this face element.
  ltg.clear();
}

/**
 *
 */
template <typename surface>
void GEOMETRYPAIR::FaceElementTemplate<surface>::Setup()
{
  // Set the reference position from the nodes.
  face_reference_position_.Clear();
  const DRT::Node* const* nodes = drt_face_element_->Nodes();
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      face_reference_position_(i_node * 3 + i_dim) = nodes[i_node]->X()[i_dim];
}

/**
 *
 */
template <typename surface>
void GEOMETRYPAIR::FaceElementTemplate<surface>::SetState(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const Epetra_Vector>& displacement)
{
  // Set the displacement at the current configuration for this element.
  std::vector<int> lmowner, lmstride;
  std::vector<double> element_displacement;
  drt_face_element_->LocationVector(*discret, face_dof_gid_, lmowner, lmstride);
  DRT::UTILS::ExtractMyValues(*displacement, element_displacement, face_dof_gid_);
  for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
    face_position_(i_dof) = face_reference_position_(i_dof) + element_displacement[i_dof];
};

/**
 *
 */
template <typename surface>
void GEOMETRYPAIR::FaceElementTemplate<surface>::CalculateAveragedNormals(
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get the connected face elements.
  std::vector<Teuchos::RCP<GEOMETRYPAIR::FaceElement>> connected_faces;
  GetConnectedFaces(face_elements, connected_faces);

  // From here on we need this face element in the connected_faces vector.
  connected_faces.push_back(face_elements.at(drt_face_element_->ParentElementId()));

  // Set up nodal averaged normal vectors (in reference and current configuration).
  LINALG::Matrix<surface::n_nodes_, 1, int> normal_ids(true);
  LINALG::Matrix<surface::n_nodes_, 1, int> normal_count(true);
  LINALG::Matrix<3, 1, double> reference_normal(true);
  LINALG::Matrix<3, 1, double> current_normal(true);
  LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, double>> reference_normals(true);
  LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, double>> current_normals(true);

  // Fill in the global node IDs.
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    normal_ids(i_node) = drt_face_element_->Nodes()[i_node]->Id();

  // Calculate the normals on the nodes from all connected (including this one) faces.
  for (const auto& connected_face : connected_faces)
  {
    for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    {
      if (connected_face->EvaluateNodalNormal(normal_ids(i_node), reference_normal, current_normal))
      {
        normal_count(i_node) += 1;
        reference_normals(i_node) += reference_normal;
        current_normals(i_node) += current_normal;
      }
    }
  }

  // Average the normals.
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    if (normal_count(i_node) == 0)
      dserror(
          "No normals calculated for node %d in volume element %d. There has to be at least one.",
          i_node, drt_face_element_->ParentElementId());

    reference_normals(i_node).Scale(1.0 / FADUTILS::VectorNorm(reference_normals(i_node)));
    current_normals(i_node).Scale(1.0 / FADUTILS::VectorNorm(current_normals(i_node)));
  }
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
    {
      reference_normals_(3 * i_node + i_dir) = reference_normals(i_node)(i_dir);
      current_normals_(3 * i_node + i_dir) = current_normals(i_node)(i_dir);
    }
  }
}

/**
 *
 */
template <typename surface>
bool GEOMETRYPAIR::FaceElementTemplate<surface>::EvaluateNodalNormal(const int node_gid,
    LINALG::Matrix<3, 1, double>& reference_normal,
    LINALG::Matrix<3, 1, double>& current_normal) const
{
  // Check if the desired node is part of this face.
  int node_lid = -1;
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    if (drt_face_element_->NodeIds()[i_node] == node_gid)
    {
      node_lid = i_node;
      break;
    }
  }
  if (node_lid < 0)
    // The desired node is not in this face.
    return false;

  // Evaluate this the normal at the desired node.
  {
    reference_normal.Clear();
    current_normal.Clear();

    // Set the parameter coordinate.
    LINALG::Matrix<3, 1, double> xi(true);
    LINALG::SerialDenseMatrix nodal_coordinates =
        DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, node_lid);

    // Calculate the normal on the surface.
    EvaluateSurfaceNormal<surface>(
        xi, face_reference_position_, reference_normal, drt_face_element_.get());
    EvaluateSurfaceNormal<surface>(xi, face_position_, current_normal, drt_face_element_.get());

    return true;
  }
}


/**
 *
 */
Teuchos::RCP<GEOMETRYPAIR::FaceElement> GEOMETRYPAIR::FaceElementFactory(
    const Teuchos::RCP<const DRT::FaceElement>& face_element)
{
  switch (face_element->Shape())
  {
    case DRT::Element::tri3:
      return Teuchos::rcp<FaceElementTemplate<t_tri3>>(
          new FaceElementTemplate<t_tri3>(face_element));
    case DRT::Element::tri6:
      return Teuchos::rcp<FaceElementTemplate<t_tri6>>(
          new FaceElementTemplate<t_tri6>(face_element));
    case DRT::Element::quad4:
      return Teuchos::rcp<FaceElementTemplate<t_quad4>>(
          new FaceElementTemplate<t_quad4>(face_element));
    case DRT::Element::quad8:
      return Teuchos::rcp<FaceElementTemplate<t_quad8>>(
          new FaceElementTemplate<t_quad8>(face_element));
    case DRT::Element::quad9:
      return Teuchos::rcp<FaceElementTemplate<t_quad9>>(
          new FaceElementTemplate<t_quad9>(face_element));
    case DRT::Element::nurbs9:
      return Teuchos::rcp<FaceElementTemplate<t_nurbs9>>(
          new FaceElementTemplate<t_nurbs9>(face_element));
    default:
      dserror("Wrong discretization type given.");
  }

  return Teuchos::null;
}
