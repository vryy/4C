/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with full segmentation of the line.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H


#include "geometry_pair_line_to_volume.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume by segmenting the
   * line at all points where it intersects the volume.
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeSegmentation
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
    GeometryPairLineToVolumeSegmentation(
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
        : GeometryPairLineToVolume<scalar_type, line, volume>(evaluation_data){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeSegmentation(){};

    /**
     * \brief Setup the segment tracker for this line element.
     */
    void Setup() override;

    /**
     * \brief This method performs the segmentation of the line with the volume.
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
