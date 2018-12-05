/*!
\file geometry_pair_line_to_volume_gauss_point_projection.H

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_H


#include "geometry_pair_line_to_volume.H"

#include "../drt_fem_general/drt_utils_integration.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume by projecting
   * Gauss points from the line to the volume. In case a line pokes out of the volumes it is
   * segmented.
   * @param scalar_type Type that will be used for scalar values.
   * @param n_nodes_element_1 Number of nodes on element 1.
   * @param n_nodal_values_element_1 Number of values per node (x3) on element 1.
   * @param n_nodes_element_2 Number of nodes on element 2.
   * @param n_nodal_values_element_2 Number of values per node (x3) on element 2.
   */
  template <typename scalar_type, unsigned int n_nodes_element_1,
      unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
      unsigned int n_nodal_values_element_2>
  class GeometryPairLineToVolumeGaussPointProjection
      : public GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
            n_nodes_element_2, n_nodal_values_element_2>
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjection()
        : GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
              n_nodes_element_2, n_nodal_values_element_2>(){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeGaussPointProjection(){};

    /**
     * \brief Setup the projection tracker for this line element.
     */
    void Setup();

    /**
     * \brief Try to project all Gauss points to the volume. Only points are checked that do not
     * already have a valid projection in the projection tracker of the evaluation data container.
     * Eventually needed segmentation at lines poking out of the volume is done in the Evaluate
     * method.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void PreEvaluate(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const;

    /**
     * \brief Check if a Gauss point projected valid for this pair in PreEvaluate. If so, all Gauss
     * points have to project valid (in the tracker, since some can be valid on other pairs). If not
     * all project, the beam pokes out of the volumes and in this method segmentation is performed.
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
    /**
     * \brief Get the line projection vector for the line element in this pair.
     * @return Reference to line projection vector.
     */
    std::vector<bool>& GetLineProjectionVectorMutable() const;
  };
}  // namespace GEOMETRYPAIR

#endif
