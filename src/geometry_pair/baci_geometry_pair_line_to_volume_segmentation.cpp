/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with full segmentation of the line.

\level 1
*/


#include "baci_geometry_pair_line_to_volume_segmentation.H"

#include "baci_geometry_pair_element_functions.H"
#include "baci_geometry_pair_line_projection.H"
#include "baci_geometry_pair_line_to_3D_evaluation_data.H"
#include "baci_geometry_pair_utility_classes.H"
#include "baci_utils_fad.H"

BACI_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename volume>
GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, line,
    volume>::GeometryPairLineToVolumeSegmentation(const DRT::Element* element1,
    const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
    : GeometryPairLineToVolume<scalar_type, line, volume>(element1, element2, evaluation_data)
{
  // Check if a segment tracker exists for this line element. If not a new one is created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::set<LineSegment<double>>>& segment_tracker_map =
      this->line_to_3d_evaluation_data_->GetSegmentTracker();

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
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>::Evaluate(
    const CORE::LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const CORE::LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DSegmentation<GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>>::Evaluate(
      this, q_line, q_volume, segments);
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27>;

BACI_NAMESPACE_CLOSE
