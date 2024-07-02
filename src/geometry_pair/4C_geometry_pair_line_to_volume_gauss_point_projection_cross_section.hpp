/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with Gauss point projection on the cylinder surface along the
line.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_line_to_volume.hpp"

FOUR_C_NAMESPACE_OPEN

// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename ScalarType>
  class ProjectionPointVolumeToVolume;
}  // namespace GEOMETRYPAIR
namespace LargeRotations
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}  // namespace LargeRotations


namespace GEOMETRYPAIR
{
  /**
   * \brief This geometry pair projects Gauss points from a line to a volume, but the Gauss points
   * are not exactly on the centerline but on a cross section defined by the centerline.
   *
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename ScalarType, typename Line, typename Volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection
      : public GeometryPairLineToVolume<ScalarType, Line, Volume>
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjectionCrossSection(
        const Core::Elements::Element* element1, const Core::Elements::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data);


    /**
     * \brief Try to project the points on the surface of the line to the volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     * @param line_triad_interpolation (in) Triad interpolation along the line.
     */
    void pre_evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments,
        const LargeRotations::TriadInterpolationLocalRotationVectors<3, double>*
            line_triad_interpolation = nullptr) const;

    /**
     * \brief The only purpose of this method is to check that all points on this line element
     * projected valid (not necessarily in this pair) in pre_evaluate.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments) const override;

   private:
    /**
     * \brief Get the line projection vector for the line element in this pair.
     * @return Reference to line projection vector.
     */
    std::vector<bool>& get_line_projection_vector() const;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
