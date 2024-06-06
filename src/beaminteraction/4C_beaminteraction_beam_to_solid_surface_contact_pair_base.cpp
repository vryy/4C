/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for contact between a 3D beam and a surface element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_contact_pair_base.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"
#include "4C_geometry_pair_scalar_types.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam,
    surface>::BeamToSolidSurfaceContactPairBase()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Clean the segments, as they will be re-evaluated in each iteration.
  this->line_to_3D_segments_.clear();

  // Set the current position of the beam element.
  const int n_patch_dof = face_element_->GetPatchGID().size();
  for (unsigned int i = 0; i < beam::n_dof_; i++)
    this->ele1pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + n_patch_dof, i, beam_centerline_dofvec[i]);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::pre_evaluate()
{
  // Call pre_evaluate on the geometry Pair.
  cast_geometry_pair()->pre_evaluate(
      this->ele1pos_, this->face_element_->GetFaceElementData(), this->line_to_3D_segments_);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam,
    surface>::CreateGeometryPair(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  this->geometry_pair_ =
      GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<scalar_type, beam, surface>(
          element1, element2, geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::SetFaceElement(
    Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>>(
          face_element, true);

  // Set the number of (centerline) degrees of freedom for the beam element in the face element
  face_element_->set_number_of_dof_other_element(
      UTILS::GetNumberOfElementCenterlineDof(this->Element1()));

  // The second element in the pair has to be the face element.
  cast_geometry_pair()->SetElement2(face_element_->GetDrtFaceElement());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, beam, surface>>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::cast_geometry_pair()
    const
{
  return Teuchos::rcp_dynamic_cast<
      GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, beam, surface>>(
      this->geometry_pair_, true);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
std::vector<int>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::get_pair_gid(
    const Discret::Discretization& discret) const
{
  // Get the beam centerline GIDs.
  Core::LinAlg::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  UTILS::GetElementCenterlineGIDIndices(discret, this->Element1(), beam_centerline_gid);

  // Get the patch (in this case just the one face element) GIDs.
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();
  std::vector<int> pair_gid;
  pair_gid.resize(beam::n_dof_ + patch_gid.size());

  // Combine beam and solid GIDs into one vector.
  for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
    pair_gid[i_dof_beam] = beam_centerline_gid(i_dof_beam);
  for (unsigned int i_dof_patch = 0; i_dof_patch < patch_gid.size(); i_dof_patch++)
    pair_gid[beam::n_dof_ + i_dof_patch] = patch_gid[i_dof_patch];

  return pair_gid;
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;


  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
