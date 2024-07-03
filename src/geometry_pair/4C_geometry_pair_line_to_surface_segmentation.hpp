/*----------------------------------------------------------------------*/
/*! \file

\brief Line to surface interaction with segmentation at the boundaries.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_SEGMENTATION_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_SEGMENTATION_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_line_to_surface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a surface by calculating
   * the points where a line intersects the surface (including it's normal dimension).
   * @tparam scalar_type Type that will be used for scalar values.
   * @tparam line Type of line element.
   * @tparam volume Type of volume element.
   */
  template <typename ScalarType, typename Line, typename Surface>
  class GeometryPairLineToSurfaceSegmentation
      : public GeometryPairLineToSurface<ScalarType, Line, Surface>
  {
   public:
    //! Public alias for scalar type so that other classes can use this type.
    using t_scalar_type = ScalarType;

    //! Public alias for line type so that other classes can use this type.
    using t_line = Line;

    //! Public alias for surface type so that other classes can use this type.
    using t_other = Surface;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToSurfaceSegmentation(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>&
            line_to_surface_evaluation_data);


    /**
     * \brief This method performs the segmentation of the line with the surface (and its normal
     * direction).
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Surface, ScalarType>& element_data_surface,
        std::vector<LineSegment<ScalarType>>& segments) const override;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
