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
void GEOMETRYPAIR::FaceElement::Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Initialize class variables.
  patch_dof_gid_.clear();
  connected_faces_.clear();

  // Initialize variables for this method.
  std::vector<int> my_node_gid, other_faces_node_gid;
  my_node_gid.clear();
  my_node_gid.reserve(drt_face_element_->NumNode());
  other_faces_node_gid.clear();

  // Temporary variables.
  std::vector<int> temp_node_dof_gid(3);

  // First add the node GIDs of this face.
  for (int i_node = 0; i_node < drt_face_element_->NumNode(); i_node++)
  {
    const DRT::Node* my_node = drt_face_element_->Nodes()[i_node];
    my_node_gid.push_back(my_node->Id());
    discret->Dof(my_node, 0, temp_node_dof_gid);
    for (const auto& value : temp_node_dof_gid) patch_dof_gid_.push_back(value);
  }

  // Add the node GIDs of the connected faces.
  ConnectedFace temp_connected_face;
  for (int i_node = 0; i_node < drt_face_element_->NumNode(); i_node++)
  {
    // Loop over all elements connected to a node of this face.
    const DRT::Node* node = drt_face_element_->Nodes()[i_node];
    for (int i_element = 0; i_element < node->NumElement(); i_element++)
    {
      const DRT::Element* element = node->Elements()[i_element];

      // Do nothing for this element.
      if (element->Id() == drt_face_element_->ParentElementId()) continue;

      // Check if the element was already searched for.
      if (connected_faces_.find(element->Id()) == connected_faces_.end())
      {
        temp_connected_face.node_lid_map_.clear();
        temp_connected_face.my_node_patch_lid_.clear();

        // Check if the element is part of the surface condition.
        auto find_in_faces = face_elements.find(element->Id());
        if (find_in_faces != face_elements.end())
        {
          // Add the node GIDs of this element.
          for (int i_node_connected_element = 0;
               i_node_connected_element < find_in_faces->second->GetDrtFaceElement()->NumNode();
               i_node_connected_element++)
          {
            const DRT::Node* other_node =
                find_in_faces->second->GetDrtFaceElement()->Nodes()[i_node_connected_element];
            const int node_id = other_node->Id();

            // Check if this node is part of this face element. If not and if it is not already in
            // other_faces_node_gid, it is added to that vector.
            auto it = std::find(my_node_gid.begin(), my_node_gid.end(), node_id);
            if (it == my_node_gid.end())
            {
              if (std::find(other_faces_node_gid.begin(), other_faces_node_gid.end(), node_id) ==
                  other_faces_node_gid.end())
              {
                other_faces_node_gid.push_back(node_id);
                discret->Dof(other_node, 0, temp_node_dof_gid);
                for (const auto& value : temp_node_dof_gid) patch_dof_gid_.push_back(value);
              }
              temp_connected_face.my_node_patch_lid_.push_back(
                  my_node_gid.size() + other_faces_node_gid.size() - 1);
            }
            else
            {
              // The node is part of this face element, add to the map entry.
              const int index_my_node = std::distance(my_node_gid.begin(), it);
              temp_connected_face.node_lid_map_[i_node_connected_element] = index_my_node;
              temp_connected_face.my_node_patch_lid_.push_back(index_my_node);
            }
          }

          // Add this element to the already searched connected elements.
          connected_faces_[element->Id()] = temp_connected_face;
        }
      }
    }
  }

  // If we only want to calculate dependencies on the face DOF and not the patch, we do not need the
  // GID of the connected faces in the GID vector of this face. The connectivity to the other patch
  // faces is still required for the calculation of the averaged reference normals.
  if (not evaluate_current_normals_)
    this->patch_dof_gid_.resize(3 * this->drt_face_element_->NumNode());
}


/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Call setup of base class.
  FaceElement::Setup(discret, face_elements);

  // Set the reference position from the nodes connected to this face.
  face_reference_position_.Clear();
  const DRT::Node* const* nodes = drt_face_element_->Nodes();
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      face_reference_position_(i_node * 3 + i_dim) = nodes[i_node]->X()[i_dim];
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::SetState(
    const Teuchos::RCP<const Epetra_Vector>& displacement,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get all displacements for the current patch.
  std::vector<int> lmowner, lmstride;
  std::vector<double> patch_displacement;
  DRT::UTILS::ExtractMyValues(*displacement, patch_displacement, patch_dof_gid_);

  // Create the full length FAD types.
  const unsigned int n_patch_dof = patch_dof_gid_.size();
  std::vector<scalar_type> patch_displacement_fad(n_patch_dof);
  for (unsigned int i_dof = 0; i_dof < n_patch_dof; i_dof++)
  {
    patch_displacement_fad[i_dof] = FADUTILS::HigherOrderFad<scalar_type>::apply(
        n_patch_dof + n_beam_dof_, n_beam_dof_ + i_dof, patch_displacement[i_dof]);
    if (i_dof < surface::n_dof_)
      face_position_(i_dof) = face_reference_position_(i_dof) + patch_displacement_fad[i_dof];
  }

  if (evaluate_current_normals_)
  {
    // Parameter coordinates corresponding to LIDs of nodes.
    LINALG::Matrix<3, 1, double> xi(true);
    LINALG::SerialDenseMatrix nodal_coordinates =
        DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);

    // Loop over the connected faces and evaluate their nodal normals.
    LINALG::Matrix<surface::n_dof_, 1, scalar_type> q_other_face;
    LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, scalar_type>> normals;
    LINALG::Matrix<3, 1, scalar_type> temp_normal;
    for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    {
      for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, i_node);
      EvaluateSurfaceNormal<surface>(xi, face_position_, normals(i_node), drt_face_element_.get());
    }
    for (const auto& value : connected_faces_)
    {
      // Get the connected face element.
      const Teuchos::RCP<const my_type>& face_element =
          Teuchos::rcp_dynamic_cast<const my_type>(face_elements.at(value.first));

      // Get the vector with the dof values for this element.
      for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          // We need to add the reference position of the other element.
          q_other_face(i_dim + 3 * i_node) =
              face_element->GetFaceReferencePosition()(i_dim + 3 * i_node) +
              patch_displacement_fad[i_dim + 3 * value.second.my_node_patch_lid_[i_node]];
        }

      // Evaluate the normals at the shared nodes.
      for (const auto& node_map_iterator : value.second.node_lid_map_)
      {
        for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
          xi(i_dim) = nodal_coordinates(i_dim, node_map_iterator.first);
        EvaluateSurfaceNormal<surface>(
            xi, q_other_face, temp_normal, face_element->GetDrtFaceElement());
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          normals(node_map_iterator.second)(i_dim) += temp_normal(i_dim);
      }
    }
    AverageNodalNormals(normals, this->current_normals_);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::CalculateAveragedReferenceNormals(
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Parameter coordinates corresponding to LIDs of nodes.
  LINALG::Matrix<3, 1, double> xi(true);
  LINALG::SerialDenseMatrix nodal_coordinates =
      DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);

  // Loop over the connected faces and evaluate their nodal normals.
  LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, double>> normals;
  LINALG::Matrix<3, 1, double> temp_normal;
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, i_node);
    EvaluateSurfaceNormal<surface>(
        xi, face_reference_position_, normals(i_node), drt_face_element_.get());
  }
  for (const auto& value : connected_faces_)
  {
    // Get the connected face element.
    const Teuchos::RCP<const my_type>& face_element =
        Teuchos::rcp_dynamic_cast<const my_type>(face_elements.at(value.first));

    // Evaluate the normals at the shared nodes.
    for (const auto& node_map_iterator : value.second.node_lid_map_)
    {
      for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
        xi(i_dim) = nodal_coordinates(i_dim, node_map_iterator.first);
      EvaluateSurfaceNormal<surface>(xi, face_element->face_reference_position_, temp_normal,
          face_element->GetDrtFaceElement());
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        normals(node_map_iterator.second)(i_dim) += temp_normal(i_dim);
    }
  }
  AverageNodalNormals(normals, this->reference_normals_);
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::EvaluateFacePositionDouble(
    const LINALG::Matrix<2, 1, double>& xi, LINALG::Matrix<3, 1, double>& r, bool reference) const
{
  LINALG::Matrix<surface::n_dof_, 1, double> position_double;
  if (reference)
    position_double = FADUTILS::CastToDouble(face_reference_position_);
  else
    position_double = FADUTILS::CastToDouble(face_position_);

  EvaluatePosition<surface>(xi, position_double, r, drt_face_element_.get());
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::EvaluateFaceNormalDouble(
    const LINALG::Matrix<2, 1, double>& xi, LINALG::Matrix<3, 1, double>& n, const bool reference,
    const bool averaged_normal) const
{
  LINALG::Matrix<surface::n_dof_, 1, double> position_double(true);
  LINALG::Matrix<3 * surface::n_nodes_, 1, double> normals_double;
  bool valid_pointer = false;

  if (averaged_normal)
  {
    if (reference)
      VectorPointerToVectorDouble(this->GetReferenceNormals(), normals_double, valid_pointer);
    else
      VectorPointerToVectorDouble(this->GetCurrentNormals(), normals_double, valid_pointer);
  }

  if (valid_pointer && averaged_normal)
  {
    // Return the normal calculated with the averaged normal field.
    EvaluateSurfaceNormal<surface>(
        xi, position_double, n, drt_face_element_.get(), &normals_double);
  }
  else if ((not valid_pointer) && averaged_normal)
  {
    // Averaged normals are desired, but there is no valid pointer to them -> return a zero vector.
    n.PutScalar(0.);
  }
  else
  {
    // Calculate the normals on the face geometry.
    if (reference)
      position_double = FADUTILS::CastToDouble(face_reference_position_);
    else
      position_double = FADUTILS::CastToDouble(face_position_);

    EvaluateSurfaceNormal<surface>(xi, position_double, n, drt_face_element_.get());
  }
}

/**
 *
 */
template <typename surface, typename scalar_type>
template <typename T>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::AverageNodalNormals(
    LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, T>>& normals,
    LINALG::Matrix<3 * surface::n_nodes_, 1, T>& averaged_normals) const
{
  averaged_normals.PutScalar(0.0);
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    normals(i_node).Scale(1.0 / FADUTILS::VectorNorm(normals(i_node)));
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      averaged_normals(i_dim + 3 * i_node) = normals(i_node)(i_dim);
  }
}


/**
 *
 */
Teuchos::RCP<GEOMETRYPAIR::FaceElement> GEOMETRYPAIR::FaceElementFactory(
    const Teuchos::RCP<const DRT::FaceElement>& face_element, const bool is_fad)
{
  if (not is_fad)
  {
    switch (face_element->Shape())
    {
      case DRT::Element::tri3:
        return Teuchos::rcp(new FaceElementTemplate<t_tri3,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_tri3::n_dof_>>(
            face_element, false));
      case DRT::Element::tri6:
        return Teuchos::rcp(new FaceElementTemplate<t_tri6,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_tri6::n_dof_>>(
            face_element, false));
      case DRT::Element::quad4:
        return Teuchos::rcp(new FaceElementTemplate<t_quad4,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_quad4::n_dof_>>(
            face_element, false));
      case DRT::Element::quad8:
        return Teuchos::rcp(new FaceElementTemplate<t_quad8,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_quad8::n_dof_>>(
            face_element, false));
      case DRT::Element::quad9:
        return Teuchos::rcp(new FaceElementTemplate<t_quad9,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_quad9::n_dof_>>(
            face_element, false));
      case DRT::Element::nurbs9:
        return Teuchos::rcp(new FaceElementTemplate<t_nurbs9,
            Sacado::ELRFad::SLFad<double, t_hermite::n_dof_ + t_nurbs9::n_dof_>>(
            face_element, false));
      default:
        dserror("Wrong discretization type given.");
    }
  }
  else
  {
    switch (face_element->Shape())
    {
      case DRT::Element::tri3:
        return Teuchos::rcp(
            new FaceElementTemplate<t_tri3, Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>>(
                face_element, true));
      case DRT::Element::tri6:
        return Teuchos::rcp(
            new FaceElementTemplate<t_tri6, Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>>(
                face_element, true));
      case DRT::Element::quad4:
        return Teuchos::rcp(
            new FaceElementTemplate<t_quad4, Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>>(
                face_element, true));
      case DRT::Element::quad8:
        return Teuchos::rcp(
            new FaceElementTemplate<t_quad8, Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>>(
                face_element, true));
      case DRT::Element::quad9:
        return Teuchos::rcp(
            new FaceElementTemplate<t_quad9, Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>>(
                face_element, true));
      case DRT::Element::nurbs9:
        return Teuchos::rcp(new FaceElementTemplate<t_nurbs9,
            Sacado::ELRFad::SLFad<Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ +
                                                                    GEOMETRYPAIR::t_nurbs9::n_dof_>,
                GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>>(
            face_element, false));
      default:
        dserror("Wrong discretization type given.");
    }
  }

  return Teuchos::null;
}
