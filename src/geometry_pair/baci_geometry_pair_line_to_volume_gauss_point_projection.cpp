/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 1
*/


#include "baci_geometry_pair_line_to_volume_gauss_point_projection.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_geometry_pair_element.hpp"
#include "baci_geometry_pair_line_projection.hpp"
#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "baci_geometry_pair_utility_classes.hpp"
#include "baci_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename volume>
GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<scalar_type, line,
    volume>::GeometryPairLineToVolumeGaussPointProjection(const DRT::Element* element1,
    const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
    : GeometryPairLineToVolume<scalar_type, line, volume>(element1, element2, evaluation_data)
{
  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTracker();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points = this->line_to_3d_evaluation_data_->GetNumberOfGaussPoints();
    std::vector<bool> new_tracking_vector;
    new_tracking_vector.resize(n_gauss_points, false);
    projection_tracker[line_element_id] = new_tracking_vector;
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<scalar_type, line,
    volume>::PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<volume, scalar_type>& element_data_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>>::PreEvaluate(this,
      element_data_line, element_data_volume, segments);
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<scalar_type, line,
    volume>::Evaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<volume, scalar_type>& element_data_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>>::Evaluate(this,
      element_data_line, element_data_volume, segments);
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>;

FOUR_C_NAMESPACE_CLOSE
