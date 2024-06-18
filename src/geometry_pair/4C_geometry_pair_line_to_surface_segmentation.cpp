/*----------------------------------------------------------------------*/
/*! \file

\brief Line to surface interaction with segmentation at the boundaries.

\level 1
*/


#include "4C_geometry_pair_line_to_surface_segmentation.hpp"

#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_projection.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_geometry_pair_scalar_types.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
GEOMETRYPAIR::GeometryPairLineToSurfaceSegmentation<scalar_type, line,
    surface>::GeometryPairLineToSurfaceSegmentation(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
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
void GEOMETRYPAIR::GeometryPairLineToSurfaceSegmentation<scalar_type, line, surface>::evaluate(
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the pre_evaluate method of the general Gauss point projection class.
  LineTo3DSegmentation<GeometryPairLineToSurfaceSegmentation<scalar_type, line, surface>>::evaluate(
      this, element_data_line, element_data_surface, segments);
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

FOUR_C_NAMESPACE_CLOSE
