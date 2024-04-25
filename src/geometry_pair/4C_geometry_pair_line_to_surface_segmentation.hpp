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
  template <typename scalar_type, typename line, typename surface>
  class GeometryPairLineToSurfaceSegmentation
      : public GeometryPairLineToSurface<scalar_type, line, surface>
  {
   public:
    //! Public alias for scalar type so that other classes can use this type.
    using t_scalar_type = scalar_type;

    //! Public alias for line type so that other classes can use this type.
    using t_line = line;

    //! Public alias for surface type so that other classes can use this type.
    using t_other = surface;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToSurfaceSegmentation(const DRT::Element* element1,
        const DRT::Element* element2,
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
    void Evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const override;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
