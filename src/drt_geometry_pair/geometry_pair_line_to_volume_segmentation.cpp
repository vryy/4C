/*!
\file geometry_pair_line_to_volume_segmentation.cpp

\brief Line to volume interaction with full segmentation of the line.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_segmentation.H"


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
      n_nodes_element_2, n_nodal_values_element_2>::Setup();
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::Evaluate(const LINALG::TMatrix<scalar_type,
                                            3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
                                            q_line,
    const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
        q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, 2, 2, 8, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, 2, 2, 20, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, 2, 2, 27, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, 2, 2, 4, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, 2, 2, 10, 1>;
