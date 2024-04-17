/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with Gauss point projection on the cylinder surface along the
line.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_HPP


#include "baci_config.hpp"

#include "baci_geometry_pair_line_to_volume.hpp"

FOUR_C_NAMESPACE_OPEN

// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename scalar_type>
  class ProjectionPointVolumeToVolume;
}  // namespace GEOMETRYPAIR
namespace LARGEROTATIONS
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}  // namespace LARGEROTATIONS


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
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection
      : public GeometryPairLineToVolume<scalar_type, line, volume>
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjectionCrossSection(const DRT::Element* element1,
        const DRT::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data);


    /**
     * \brief Try to project the points on the surface of the line to the volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     * @param line_triad_interpolation (in) Triad interpolation along the line.
     */
    void PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<volume, scalar_type>& element_data_volume,
        std::vector<LineSegment<scalar_type>>& segments,
        const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>*
            line_triad_interpolation = nullptr) const;

    /**
     * \brief The only purpose of this method is to check that all points on this line element
     * projected valid (not necessarily in this pair) in PreEvaluate.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void Evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<volume, scalar_type>& element_data_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;

   private:
    /**
     * \brief Get the line projection vector for the line element in this pair.
     * @return Reference to line projection vector.
     */
    std::vector<bool>& GetLineProjectionVector() const;
  };
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
