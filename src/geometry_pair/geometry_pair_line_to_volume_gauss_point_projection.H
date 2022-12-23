/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_H


#include "geometry_pair_line_to_volume.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume by projecting
   * Gauss points from the line to the volume. In case a line pokes out of the volumes it is
   * segmented.
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjection
      : public GeometryPairLineToVolume<scalar_type, line, volume>
  {
   public:
    //! Public alias for scalar type so that other classes can use this type.
    using t_scalar_type = scalar_type;

    //! Public alias for line type so that other classes can use this type.
    using t_line = line;

    //! Public alias for volume type so that other classes can use this type.
    using t_other = volume;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjection(
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
        : GeometryPairLineToVolume<scalar_type, line, volume>(evaluation_data){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeGaussPointProjection(){};

    /**
     * \brief Setup the projection tracker for this line element.
     */
    void Setup() override;

    /**
     * \brief Try to project all Gauss points to the volume.
     *
     * Only points are checked that do not
     * already have a valid projection in the projection tracker of the evaluation data container.
     * Eventually needed segmentation at lines poking out of the volume is done in the Evaluate
     * method.
     *
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void PreEvaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;

    /**
     * \brief Check if a Gauss point projected valid for this pair in PreEvaluate.
     *
     * If so, all Gauss points have to project valid (in the tracker, since some can be valid on
     * other pairs). If not all project, the beam pokes out of the volumes and in this method
     * segmentation is performed.
     *
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void Evaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;
  };
}  // namespace GEOMETRYPAIR

#endif
