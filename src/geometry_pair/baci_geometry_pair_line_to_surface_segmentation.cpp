/*----------------------------------------------------------------------*/
/*! \file

\brief Line to surface interaction with segmentation at the boundaries.

\level 1
*/


#include "baci_geometry_pair_line_to_surface_segmentation.hpp"

#include "baci_geometry_pair_element_functions.hpp"
#include "baci_geometry_pair_line_projection.hpp"
#include "baci_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "baci_geometry_pair_scalar_types.hpp"

BACI_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
GEOMETRYPAIR::GeometryPairLineToSurfaceSegmentation<scalar_type, line,
    surface>::GeometryPairLineToSurfaceSegmentation(const DRT::Element* element1,
    const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>& line_to_surface_evaluation_data)
    : GeometryPairLineToSurface<scalar_type, line, surface>(
          element1, element2, line_to_surface_evaluation_data)
{
  // Check if a segment tracker exists for this line element. If not a new one is created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::set<LineSegment<double>>>& segment_tracker_map =
      this->line_to_surface_evaluation_data_->GetSegmentTracker();

  if (segment_tracker_map.find(line_element_id) == segment_tracker_map.end())
  {
    std::set<LineSegment<double>> new_tracking_set;
    new_tracking_set.clear();
    segment_tracker_map[line_element_id] = new_tracking_set;
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceSegmentation<scalar_type, line, surface>::Evaluate(
    const CORE::LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const CORE::LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    std::vector<LineSegment<scalar_type>>& segments,
    const CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DSegmentation<GeometryPairLineToSurfaceSegmentation<scalar_type, line, surface>>::Evaluate(
      this, q_line, q_surface, segments, nodal_normals);
}


/**
 * Explicit template initialization of template class.
 */
namespace GEOMETRYPAIR
{
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_line2,
      t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_line2,
      t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_line2,
      t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_line2,
      t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_line2,
      t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<double, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class GeometryPairLineToSurfaceSegmentation<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class GeometryPairLineToSurfaceSegmentation<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE
