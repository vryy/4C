/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with full segmentation of the line.

\level 1
*/


#include "geometry_pair_line_to_volume_segmentation.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"
#include "geometry_pair_utility_classes.H"
#include "geometry_pair_line_projection.H"

#include "../headers/FAD_utils.H"


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, line, volume>::Setup();

  // Check if a segment tracker exists for this line element. If not a new one is created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::set<LineSegment<double>>>& segment_tracker_map =
      this->line_to_3d_evaluation_data_->GetSegmentTrackerMutable();

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
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

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
