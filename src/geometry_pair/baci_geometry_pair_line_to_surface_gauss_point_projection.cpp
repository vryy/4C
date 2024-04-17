/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
*/


#include "baci_geometry_pair_line_to_surface_gauss_point_projection.hpp"

#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_line_projection.hpp"
#include "baci_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "baci_geometry_pair_scalar_types.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
GEOMETRYPAIR::GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line,
    surface>::GeometryPairLineToSurfaceGaussPointProjection(const DRT::Element* element1,
    const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>& line_to_surface_evaluation_data)
    : GeometryPairLineToSurface<scalar_type, line, surface>(
          element1, element2, line_to_surface_evaluation_data)
{
  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_surface_evaluation_data_->GetGaussPointProjectionTracker();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points = this->line_to_surface_evaluation_data_->GetNumberOfGaussPoints();
    std::vector<bool> new_tracking_vector;
    new_tracking_vector.resize(n_gauss_points, false);
    projection_tracker[line_element_id] = new_tracking_vector;
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line,
    surface>::PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line, surface>>::PreEvaluate(this,
      element_data_line, element_data_surface, segments);
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line,
    surface>::Evaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line, surface>>::Evaluate(this,
      element_data_line, element_data_surface, segments);
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
std::vector<bool>& GEOMETRYPAIR::GeometryPairLineToSurfaceGaussPointProjection<scalar_type, line,
    surface>::GetLineProjectionVector() const
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_surface_evaluation_data_->GetGaussPointProjectionTracker();
  return projection_tracker[line_element_id];
}


/**
 * Explicit template initialization of template class.
 */
namespace GEOMETRYPAIR
{

  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<double, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceGaussPointProjection<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceGaussPointProjection<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE
