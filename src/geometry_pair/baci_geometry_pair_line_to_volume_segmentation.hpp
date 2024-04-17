/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with full segmentation of the line.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_HPP


#include "baci_config.hpp"

#include "baci_geometry_pair_line_to_volume.hpp"

FOUR_C_NAMESPACE_OPEN

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
    GeometryPairLineToVolumeSegmentation(const DRT::Element* element1, const DRT::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data);


    /**
     * \brief This method performs the segmentation of the line with the volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void Evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<volume, scalar_type>& element_data_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
