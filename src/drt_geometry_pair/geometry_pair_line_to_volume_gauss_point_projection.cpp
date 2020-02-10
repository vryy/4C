/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_gauss_point_projection.H"
#include "geometry_pair_element.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"


#include "geometry_pair_line_projection.H"

/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, line, volume>::Setup();

  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTrackerMutable();

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
    volume>::PreEvaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>>::PreEvaluate(this,
      q_line, q_volume, segments);
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<scalar_type, line,
    volume>::Evaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Call the PreEvaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>>::Evaluate(this,
      q_line, q_volume, segments);
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
