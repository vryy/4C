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
     * \brief This method performs the segmentation of the line with the volume.
     *
     * First every search point on the beam is projected to the volume. For every search point that
     * has a projection (Newton converged, point does not have to be inside the volume), surface
     * intersections are performed. The resulting surface projection points are sorted and only
     * unique points are kept. Between two successive intersection points there is a line segment
     * (either inside or outside of the beam). To check whether the segment is part of this pair, a
     * point between two successive intersection points is projected, and if it projects valid, the
     * segment is assumed to be in this pair.
     *
     * There can be the special case when two consecutive segments are inside the pair, this happens
     * when the solid is exactly along the surface of a volume and a search point lies on that
     * surface.
     *
     * This approach should work for almost all possible cases. One case that can occur, is when
     * none of the search points project, but the line intersects the volume. In this case more
     * search points are the (expensive!) solution.
     *
     * In the case where a line lies on the surface between two volumes, segments on both volumes
     * are valid. To avoid this case a global segment tracker is kept in the evaluation data, that
     * allows a segment only once on all line to volume pairs. This does not work in the very
     * special case, that the volumes nonconforming discretization.
     *
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
     * \brief Get the segment tracking vector for the line element in this pair. The type is double
     * for the segments there, as this is only needed for comparison purposes. The segment values of
     * the tracker should NOT be used for further calculation.
     * @return Reference to segment tracking vector.
     */
    std::set<LineSegment<double>>& GetSegmentTrackingSetMutable() const;
  };
}  // namespace GEOMETRYPAIR

#endif
