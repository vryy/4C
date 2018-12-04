/*!
\file geometry_pair_line_to_volume_segmentation.H

\brief Line to volume interaction with full segmentation of the line.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H


#include "geometry_pair_line_to_volume.H"

#include "../drt_fem_general/drt_utils_integration.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume by segmenting the
   * line when it intersects a volume.
   * @param scalar_type Type that will be used for scalar values.
   * @param n_nodes_element_1 Number of nodes on element 1.
   * @param n_nodal_values_element_1 Number of values per node (x3) on element 1.
   * @param n_nodes_element_2 Number of nodes on element 2.
   * @param n_nodal_values_element_2 Number of values per node (x3) on element 2.
   */
  template <typename scalar_type, unsigned int n_nodes_element_1,
      unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
      unsigned int n_nodal_values_element_2>
  class GeometryPairLineToVolumeSegmentation
      : public GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
            n_nodes_element_2, n_nodal_values_element_2>
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeSegmentation()
        : GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
              n_nodes_element_2, n_nodal_values_element_2>(){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeSegmentation(){};

    /**
     * \brief Setup the segment tracker for this line element.
     */
    void Setup();

    /**
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void Evaluate(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const;

   private:
  };
}  // namespace GEOMETRYPAIR

#endif
