/*!

\brief Line to volume interaction with Gauss point projection on the cylinder surface along the
line.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CROSS_SECTION_H


#include "geometry_pair_line_to_volume.H"


// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename scalar_type>
  class ProjectionPointVolumeToVolume;
}  // namespace GEOMETRYPAIR


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
    GeometryPairLineToVolumeGaussPointProjectionCrossSection()
        : GeometryPairLineToVolume<scalar_type, line, volume>(){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeGaussPointProjectionCrossSection(){};

    /**
     * \brief Setup the projection tracker for this line element.
     */
    void Setup();

    /**
     * \brief Try to project the points on the surface of the line to the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void PreEvaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;

    /**
     * \brief The only purpose of this method is to check that all points on this line element
     * projected valid (not necessarily in this pair) in PreEvaluate.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void Evaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const override;

   private:
    /**
     * \brief Get the line projection vector for the line element in this pair.
     * @return Reference to line projection vector.
     */
    std::vector<bool>& GetLineProjectionVectorMutable() const;
  };
}  // namespace GEOMETRYPAIR

#endif
