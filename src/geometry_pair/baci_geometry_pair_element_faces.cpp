/*----------------------------------------------------------------------*/
/*! \file

\brief Element classes that represent faces, i.e. surface elements.

\level 1
*/
// End doxygen header.


#include "baci_geometry_pair_element_faces.H"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_geometry_pair_scalar_types.H"
#include "baci_inpar_geometry_pair.H"


/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get the DOF GIDs of this face.
  patch_dof_gid_.clear();
  std::vector<int> lmrowowner;
  std::vector<int> lmstride;
  this->GetDrtFaceElement()->LocationVector(*discret, this->patch_dof_gid_, lmrowowner, lmstride);

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
  // Get all displacements for the current face / patch.
  std::vector<double> patch_displacement;
  DRT::UTILS::ExtractMyValues(*displacement, patch_displacement, patch_dof_gid_);

  // Create the full length FAD types.
  const unsigned int n_patch_dof = patch_dof_gid_.size();
  std::vector<scalar_type> patch_displacement_fad(n_patch_dof);
  for (unsigned int i_dof = 0; i_dof < n_patch_dof; i_dof++)
  {
    patch_displacement_fad[i_dof] = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        n_patch_dof + n_beam_dof_, n_beam_dof_ + i_dof, patch_displacement[i_dof]);
    if (i_dof < surface::n_dof_)
      face_position_(i_dof) = face_reference_position_(i_dof) + patch_displacement_fad[i_dof];
  }
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::EvaluateFacePositionDouble(
    const CORE::LINALG::Matrix<2, 1, double>& xi, CORE::LINALG::Matrix<3, 1, double>& r,
    bool reference) const
{
  CORE::LINALG::Matrix<surface::n_dof_, 1, double> position_double;
  if (reference)
    position_double = CORE::FADUTILS::CastToDouble(face_reference_position_);
  else
    position_double = CORE::FADUTILS::CastToDouble(face_position_);

  EvaluatePosition<surface>(xi, position_double, r, drt_face_element_.get());
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>::EvaluateFaceNormalDouble(
    const CORE::LINALG::Matrix<2, 1, double>& xi, CORE::LINALG::Matrix<3, 1, double>& n,
    const bool reference, const bool averaged_normal) const
{
  if (averaged_normal)
  {
    n.PutScalar(0.0);
    return;
  }
  else
  {
    // Calculate the normals on the face geometry.
    CORE::LINALG::Matrix<surface::n_dof_, 1, double> position_double(true);
    if (reference)
      position_double = CORE::FADUTILS::CastToDouble(face_reference_position_);
    else
      position_double = CORE::FADUTILS::CastToDouble(face_position_);

    EvaluateSurfaceNormal<surface>(xi, position_double, n, drt_face_element_.get());
  }
}


/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Call setup of the base class.
  base_class::Setup(discret, face_elements);

  // Initialize class variables.
  connected_faces_.clear();

  // Initialize variables for this method.
  std::vector<int> my_node_gid, other_faces_node_gid;
  my_node_gid.clear();
  my_node_gid.reserve(this->drt_face_element_->NumNode());
  other_faces_node_gid.clear();

  // Temporary variables.
  std::vector<int> temp_node_dof_gid(3);

  // First add the node GIDs of this face.
  for (int i_node = 0; i_node < this->drt_face_element_->NumNode(); i_node++)
  {
    const DRT::Node* my_node = this->drt_face_element_->Nodes()[i_node];
    my_node_gid.push_back(my_node->Id());
    discret->Dof(my_node, 0, temp_node_dof_gid);
  }

  // Add the node GIDs of the connected faces.
  ConnectedFace temp_connected_face;
  for (int i_node = 0; i_node < this->drt_face_element_->NumNode(); i_node++)
  {
    // Loop over all elements connected to a node of this face.
    const DRT::Node* node = this->drt_face_element_->Nodes()[i_node];
    for (int i_element = 0; i_element < node->NumElement(); i_element++)
    {
      const DRT::Element* element = node->Elements()[i_element];

      // Do nothing for this element.
      if (element->Id() == this->drt_face_element_->ParentElementId()) continue;

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
            const int other_node_id = other_node->Id();

            // Check if the other node is part of this face element.
            auto it = std::find(my_node_gid.begin(), my_node_gid.end(), other_node_id);
            if (it == my_node_gid.end())
            {
              // The other node is not part of this face element. Check if this other node was
              // already processed for this patch.
              auto it_other = std::find(
                  other_faces_node_gid.begin(), other_faces_node_gid.end(), other_node_id);
              if (it_other == other_faces_node_gid.end())
              {
                // This node was not processed yet, so add it to other_faces_node_gid.
                other_faces_node_gid.push_back(other_node_id);
                discret->Dof(other_node, 0, temp_node_dof_gid);
                for (const auto& value : temp_node_dof_gid) this->patch_dof_gid_.push_back(value);

                // Add the patch id of this other node to the face element tracker.
                temp_connected_face.my_node_patch_lid_.push_back(
                    my_node_gid.size() + other_faces_node_gid.size() - 1);
              }
              else
              {
                // Get the patch index of the other node and add it to the face element tracker.
                const int index_other_node = std::distance(other_faces_node_gid.begin(), it_other);
                temp_connected_face.my_node_patch_lid_.push_back(
                    my_node_gid.size() + index_other_node);
              }
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
void GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>::SetState(
    const Teuchos::RCP<const Epetra_Vector>& displacement,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get all displacements for the current face / patch.
  std::vector<double> patch_displacement;
  DRT::UTILS::ExtractMyValues(*displacement, patch_displacement, this->patch_dof_gid_);

  // Create the full length FAD types.
  const unsigned int n_patch_dof = this->patch_dof_gid_.size();
  std::vector<scalar_type> patch_displacement_fad(n_patch_dof);
  for (unsigned int i_dof = 0; i_dof < n_patch_dof; i_dof++)
  {
    patch_displacement_fad[i_dof] = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        n_patch_dof + this->n_beam_dof_, this->n_beam_dof_ + i_dof, patch_displacement[i_dof]);
    if (i_dof < surface::n_dof_)
      this->face_position_(i_dof) =
          this->face_reference_position_(i_dof) + patch_displacement_fad[i_dof];
  }

  if (evaluate_current_normals_)
  {
    // Parameter coordinates corresponding to LIDs of nodes.
    CORE::LINALG::Matrix<3, 1, double> xi(true);
    CORE::LINALG::SerialDenseMatrix nodal_coordinates =
        CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);

    // Loop over the connected faces and evaluate their nodal normals.
    CORE::LINALG::Matrix<surface::n_dof_, 1, scalar_type> q_other_face;
    CORE::LINALG::Matrix<surface::n_nodes_, 1, CORE::LINALG::Matrix<3, 1, scalar_type>> normals;
    CORE::LINALG::Matrix<3, 1, scalar_type> temp_normal;
    for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
    {
      for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, i_node);
      EvaluateSurfaceNormal<surface>(
          xi, this->face_position_, normals(i_node), this->drt_face_element_.get());
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
    AverageNodalNormals(normals, current_normals_);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type>
void GEOMETRYPAIR::FaceElementPatchTemplate<surface,
    scalar_type>::CalculateAveragedReferenceNormals(const std::unordered_map<int,
    Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Parameter coordinates corresponding to LIDs of nodes.
  CORE::LINALG::Matrix<3, 1, double> xi(true);
  CORE::LINALG::SerialDenseMatrix nodal_coordinates =
      CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);

  // Loop over the connected faces and evaluate their nodal normals.
  CORE::LINALG::Matrix<surface::n_nodes_, 1, CORE::LINALG::Matrix<3, 1, double>> normals;
  CORE::LINALG::Matrix<3, 1, double> temp_normal;
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = nodal_coordinates(i_dim, i_node);
    EvaluateSurfaceNormal<surface>(
        xi, this->face_reference_position_, normals(i_node), this->drt_face_element_.get());
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
void GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>::EvaluateFaceNormalDouble(
    const CORE::LINALG::Matrix<2, 1, double>& xi, CORE::LINALG::Matrix<3, 1, double>& n,
    const bool reference, const bool averaged_normal) const
{
  if (averaged_normal)
  {
    CORE::LINALG::Matrix<surface::n_dof_, 1, double> position_double(true);
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, double> normals_double;
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, double>* normals_double_ptr;

    if (reference)
      normals_double_ptr = VectorPointerToVectorDouble(this->GetReferenceNormals(), normals_double);
    else
      normals_double_ptr = VectorPointerToVectorDouble(this->GetCurrentNormals(), normals_double);

    if (normals_double_ptr != nullptr)
    {
      // Return the normal calculated with the averaged normal field.
      EvaluateSurfaceNormal<surface>(
          xi, position_double, n, this->drt_face_element_.get(), normals_double_ptr);
    }
    else
    {
      // Averaged normals are desired, but there is no valid pointer to them -> return a zero
      // vector.
      n.PutScalar(0.);
    }
  }
  else
  {
    // If no averaged normals should be calculated we can call the base method here.
    base_class::EvaluateFaceNormalDouble(xi, n, reference, false);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type>
template <typename T>
void GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>::AverageNodalNormals(
    CORE::LINALG::Matrix<surface::n_nodes_, 1, CORE::LINALG::Matrix<3, 1, T>>& normals,
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, T>& averaged_normals) const
{
  averaged_normals.PutScalar(0.0);
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    normals(i_node).Scale(1.0 / CORE::FADUTILS::VectorNorm(normals(i_node)));
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      averaged_normals(i_dim + 3 * i_node) = normals(i_node)(i_dim);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type, typename volume>
void GEOMETRYPAIR::FaceElementTemplateExtendedVolume<surface, scalar_type, volume>::Setup(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get the DOF GIDs of this face and volume element.
  this->patch_dof_gid_.clear();
  std::vector<int> surface_gid;
  std::vector<int> lmrowowner;
  std::vector<int> lmstride;
  this->GetDrtFaceElement()->LocationVector(*discret, surface_gid, lmrowowner, lmstride);
  this->GetDrtFaceElement()->ParentElement()->LocationVector(
      *discret, this->patch_dof_gid_, lmrowowner, lmstride);

  // Safety checks.
  if (this->patch_dof_gid_.size() != volume::n_dof_ or surface_gid.size() != surface::n_dof_)
    dserror("Mismatch in GID sizes!");

  // Calculate the face to volume DOF map.
  for (unsigned i_dof_surf = 0; i_dof_surf < surface::n_dof_; i_dof_surf++)
  {
    auto dof_iterator =
        find(this->patch_dof_gid_.begin(), this->patch_dof_gid_.end(), surface_gid[i_dof_surf]);
    if (dof_iterator != this->patch_dof_gid_.end())
    {
      // Calculating the index.
      surface_dof_lid_map_(i_dof_surf) = dof_iterator - this->patch_dof_gid_.begin();
    }
    else
      dserror("Could not find the surface DOF %d in the volume DOFs", surface_gid[i_dof_surf]);
  }

  // Set the reference position.
  volume_reference_position_.Clear();
  this->face_reference_position_.Clear();
  const DRT::Node* const* nodes = this->drt_face_element_->ParentElement()->Nodes();
  for (unsigned int i_node = 0; i_node < volume::n_nodes_; i_node++)
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      volume_reference_position_(i_node * 3 + i_dim) = nodes[i_node]->X()[i_dim];
  for (unsigned int i_dof_surf = 0; i_dof_surf < surface::n_dof_; i_dof_surf++)
    this->face_reference_position_(i_dof_surf) =
        volume_reference_position_(surface_dof_lid_map_(i_dof_surf));

  // Surface node to volume node map.
  CORE::LINALG::Matrix<surface::n_nodes_, 1, int> surface_node_lid_map;
  for (unsigned i_node_surf = 0; i_node_surf < surface::n_nodes_; i_node_surf++)
    surface_node_lid_map(i_node_surf) = surface_dof_lid_map_(i_node_surf * 3) / 3;

  if (surface_node_lid_map(0) == 0 and surface_node_lid_map(1) == 4 and
      surface_node_lid_map(2) == 7 and surface_node_lid_map(3) == 3)
  {
    face_to_volume_coordinate_axis_map_(0) = 2;
    face_to_volume_coordinate_axis_map_(1) = 1;
    face_to_volume_coordinate_axis_factor_(0) = 1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 0;
    third_direction_factor_ = -1;
  }
  else if (surface_node_lid_map(0) == 1 and surface_node_lid_map(1) == 2 and
           surface_node_lid_map(2) == 6 and surface_node_lid_map(3) == 5)
  {
    face_to_volume_coordinate_axis_map_(0) = 1;
    face_to_volume_coordinate_axis_map_(1) = 2;
    face_to_volume_coordinate_axis_factor_(0) = 1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 0;
    third_direction_factor_ = 1;
  }
  else if (surface_node_lid_map(0) == 0 and surface_node_lid_map(1) == 1 and
           surface_node_lid_map(2) == 5 and surface_node_lid_map(3) == 4)
  {
    face_to_volume_coordinate_axis_map_(0) = 0;
    face_to_volume_coordinate_axis_map_(1) = 2;
    face_to_volume_coordinate_axis_factor_(0) = 1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 1;
    third_direction_factor_ = -1;
  }
  else if (surface_node_lid_map(0) == 2 and surface_node_lid_map(1) == 3 and
           surface_node_lid_map(2) == 7 and surface_node_lid_map(3) == 6)
  {
    face_to_volume_coordinate_axis_map_(0) = 0;
    face_to_volume_coordinate_axis_map_(1) = 2;
    face_to_volume_coordinate_axis_factor_(0) = -1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 1;
    third_direction_factor_ = 1;
  }
  else if (surface_node_lid_map(0) == 0 and surface_node_lid_map(1) == 3 and
           surface_node_lid_map(2) == 2 and surface_node_lid_map(3) == 1)
  {
    face_to_volume_coordinate_axis_map_(0) = 1;
    face_to_volume_coordinate_axis_map_(1) = 0;
    face_to_volume_coordinate_axis_factor_(0) = 1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 2;
    third_direction_factor_ = -1;
  }
  else if (surface_node_lid_map(0) == 4 and surface_node_lid_map(1) == 5 and
           surface_node_lid_map(2) == 6 and surface_node_lid_map(3) == 7)
  {
    face_to_volume_coordinate_axis_map_(0) = 0;
    face_to_volume_coordinate_axis_map_(1) = 1;
    face_to_volume_coordinate_axis_factor_(0) = 1;
    face_to_volume_coordinate_axis_factor_(1) = 1;
    third_direction_ = 2;
    third_direction_factor_ = 1;
  }
  else
    dserror("Could not map face to volume.");

  // Calculate the reference normals.
  CalculateNormals(volume_reference_position_, this->face_reference_position_, reference_normals_);
}


/**
 *
 */
template <typename surface, typename scalar_type, typename volume>
void GEOMETRYPAIR::FaceElementTemplateExtendedVolume<surface, scalar_type, volume>::SetState(
    const Teuchos::RCP<const Epetra_Vector>& displacement,
    const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
{
  // Get all displacements for the current face / volume.
  std::vector<double> volume_displacement;
  DRT::UTILS::ExtractMyValues(*displacement, volume_displacement, this->patch_dof_gid_);

  // Create the full length FAD types.
  std::vector<scalar_type> patch_displacement_fad(volume::n_dof_);

  volume_position_.PutScalar(0.0);
  for (unsigned int i_dof = 0; i_dof < volume::n_dof_; i_dof++)
  {
    volume_position_(i_dof) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        volume::n_dof_ + this->n_beam_dof_, this->n_beam_dof_ + i_dof,
        volume_displacement[i_dof] + volume_reference_position_(i_dof));
  }
  this->face_position_.PutScalar(0.0);
  for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
    this->face_position_(i_dof) = volume_position_(surface_dof_lid_map_(i_dof));

  CalculateNormals(volume_position_, this->face_position_, current_normals_);
}

/**
 *
 */
template <typename surface, typename scalar_type, typename volume>
template <typename scalar_type_normal>
void GEOMETRYPAIR::FaceElementTemplateExtendedVolume<surface, scalar_type,
    volume>::CalculateNormals(const CORE::LINALG::Matrix<volume::n_dof_, 1, scalar_type_normal>&
                                  volume_position,
    const CORE::LINALG::Matrix<surface::n_dof_, 1, scalar_type_normal>& surface_position,
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type_normal>& normals) const
{
  // Parameter coordinates corresponding to LIDs of nodes.
  CORE::LINALG::Matrix<2, 1, double> xi_surface(true);
  CORE::LINALG::SerialDenseMatrix nodal_coordinates =
      CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);

  // Loop over the faces and evaluate the "normals" at the nodes.
  CORE::LINALG::Matrix<3, 1, double> xi_volume(true);
  CORE::LINALG::Matrix<3, 1, scalar_type_normal> r_surface;
  CORE::LINALG::Matrix<3, 1, scalar_type_normal> r_volume;
  CORE::LINALG::Matrix<3, 3, scalar_type_normal> dr_volume;
  for (unsigned int i_node = 0; i_node < surface::n_nodes_; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
      xi_surface(i_dim) = nodal_coordinates(i_dim, i_node);
    XiFaceToXiVolume(xi_surface, xi_volume);

    EvaluatePosition<surface>(xi_surface, surface_position, r_surface);
    EvaluatePosition<volume>(xi_volume, volume_position, r_volume);
    r_volume -= r_surface;
    if (CORE::FADUTILS::VectorNorm(r_volume) > 1e-10)
      dserror("Nodal positions for face and volume do not match.");

    EvaluatePositionDerivative1<volume>(xi_volume, volume_position, dr_volume);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      normals(3 * i_node + i_dim) = third_direction_factor_ * dr_volume(i_dim, third_direction_);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type, typename volume>
void GEOMETRYPAIR::FaceElementTemplateExtendedVolume<surface, scalar_type,
    volume>::EvaluateFaceNormalDouble(const CORE::LINALG::Matrix<2, 1, double>& xi,
    CORE::LINALG::Matrix<3, 1, double>& n, const bool reference, const bool averaged_normal) const
{
  if (averaged_normal)
  {
    CORE::LINALG::Matrix<surface::n_dof_, 1, double> position_double(true);
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, double> normals_double;
    CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, double>* normals_double_ptr;

    if (reference)
      normals_double_ptr = VectorPointerToVectorDouble(this->GetReferenceNormals(), normals_double);
    else
      normals_double_ptr = VectorPointerToVectorDouble(this->GetCurrentNormals(), normals_double);

    if (normals_double_ptr != nullptr)
    {
      // Return the normal calculated with the averaged normal field.
      EvaluateSurfaceNormal<surface>(
          xi, position_double, n, this->drt_face_element_.get(), normals_double_ptr);
    }
    else
    {
      // Averaged normals are desired, but there is no valid pointer to them -> return a zero
      // vector.
      n.PutScalar(0.);
    }
  }
  else
  {
    // If no averaged normals should be calculated we can call the base method here.
    base_class::EvaluateFaceNormalDouble(xi, n, reference, false);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type, typename volume>
void GEOMETRYPAIR::FaceElementTemplateExtendedVolume<surface, scalar_type,
    volume>::XiFaceToXiVolume(const CORE::LINALG::Matrix<2, 1, double>& xi_face,
    CORE::LINALG::Matrix<3, 1, double>& xi_volume) const
{
  xi_volume(face_to_volume_coordinate_axis_map_(0)) =
      xi_face(0) * face_to_volume_coordinate_axis_factor_(0);
  xi_volume(face_to_volume_coordinate_axis_map_(1)) =
      xi_face(1) * face_to_volume_coordinate_axis_factor_(1);
  xi_volume(third_direction_) = third_direction_factor_;
}


/**
 *
 */
Teuchos::RCP<GEOMETRYPAIR::FaceElement> GEOMETRYPAIR::FaceElementFactory(
    const Teuchos::RCP<const DRT::FaceElement>& face_element, const int fad_order,
    const INPAR::GEOMETRYPAIR::SurfaceNormals surface_normal_strategy)
{
  const bool is_fad = fad_order > 0;
  if (not is_fad)
  {
    switch (face_element->Shape())
    {
      case DRT::Element::quad4:
        return Teuchos::rcp(
            new FaceElementPatchTemplate<t_quad4, line_to_surface_scalar_type<t_hermite, t_quad4>>(
                face_element, false));
      case DRT::Element::quad8:
        return Teuchos::rcp(
            new FaceElementPatchTemplate<t_quad8, line_to_surface_scalar_type<t_hermite, t_quad8>>(
                face_element, false));
      case DRT::Element::quad9:
        return Teuchos::rcp(
            new FaceElementPatchTemplate<t_quad9, line_to_surface_scalar_type<t_hermite, t_quad9>>(
                face_element, false));
      case DRT::Element::tri3:
        return Teuchos::rcp(
            new FaceElementPatchTemplate<t_tri3, line_to_surface_scalar_type<t_hermite, t_tri3>>(
                face_element, false));
      case DRT::Element::tri6:
        return Teuchos::rcp(
            new FaceElementPatchTemplate<t_tri6, line_to_surface_scalar_type<t_hermite, t_tri6>>(
                face_element, false));
      case DRT::Element::nurbs9:
        return Teuchos::rcp(
            new FaceElementTemplate<t_nurbs9, line_to_surface_scalar_type<t_hermite, t_nurbs9>>(
                face_element));
      default:
        dserror("Wrong discretization type given.");
    }
  }
  else
  {
    if (surface_normal_strategy == INPAR::GEOMETRYPAIR::SurfaceNormals::standard)
    {
      switch (fad_order)
      {
        case 1:
        {
          switch (face_element->Shape())
          {
            case DRT::Element::quad4:
              return Teuchos::rcp(new FaceElementPatchTemplate<t_quad4,
                  line_to_surface_patch_scalar_type_1st_order>(face_element, true));
            case DRT::Element::quad8:
              return Teuchos::rcp(new FaceElementPatchTemplate<t_quad8,
                  line_to_surface_patch_scalar_type_1st_order>(face_element, true));
            case DRT::Element::quad9:
              return Teuchos::rcp(new FaceElementPatchTemplate<t_quad9,
                  line_to_surface_patch_scalar_type_1st_order>(face_element, true));
            case DRT::Element::tri3:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_tri3, line_to_surface_patch_scalar_type_1st_order>(
                      face_element, true));
            case DRT::Element::tri6:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_tri6, line_to_surface_patch_scalar_type_1st_order>(
                      face_element, true));
            case DRT::Element::nurbs9:
              return Teuchos::rcp(new FaceElementTemplate<t_nurbs9,
                  line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>>(
                  face_element));
            default:
              dserror("Wrong discretization type given.");
          }
          break;
        }
        case 2:
        {
          switch (face_element->Shape())
          {
            case DRT::Element::quad4:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_quad4, line_to_surface_patch_scalar_type>(
                      face_element, true));
            case DRT::Element::quad8:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_quad8, line_to_surface_patch_scalar_type>(
                      face_element, true));
            case DRT::Element::quad9:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_quad9, line_to_surface_patch_scalar_type>(
                      face_element, true));
            case DRT::Element::tri3:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_tri3, line_to_surface_patch_scalar_type>(
                      face_element, true));
            case DRT::Element::tri6:
              return Teuchos::rcp(
                  new FaceElementPatchTemplate<t_tri6, line_to_surface_patch_scalar_type>(
                      face_element, true));
            case DRT::Element::nurbs9:
              return Teuchos::rcp(new FaceElementTemplate<t_nurbs9,
                  line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>>(face_element));
            default:
              dserror("Wrong discretization type given.");
          }
          break;
        }
        default:
          dserror("Got unexpected fad order.");
      }
    }
    else if (surface_normal_strategy == INPAR::GEOMETRYPAIR::SurfaceNormals::extended_volume)
    {
      switch (face_element->Shape())
      {
        case DRT::Element::quad4:
          return Teuchos::rcp(new FaceElementTemplateExtendedVolume<t_quad4,
              line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hex8>(
              face_element));
        case DRT::Element::quad8:
          return Teuchos::rcp(new FaceElementTemplateExtendedVolume<t_quad8,
              line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hex20>(
              face_element));
        case DRT::Element::quad9:
          return Teuchos::rcp(new FaceElementTemplateExtendedVolume<t_quad9,
              line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hex27>(
              face_element));
        default:
          dserror("Got unexpected face type for extended volume coupling.");
      }
    }
    else
      dserror("Surface normal strategy not recognized.");
  }

  dserror("Could not create a face element.");
  return Teuchos::null;
}
