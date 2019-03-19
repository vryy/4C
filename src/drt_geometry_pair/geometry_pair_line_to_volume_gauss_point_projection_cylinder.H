/*!
\file geometry_pair_line_to_volume_gauss_point_projection_cylinder.H

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CYLINDER_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_GAUSS_POINT_PROJECTION_CYLINDER_H


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
   * \brief Class that handles the geometrical interactions of a line and a volume by projecting
   * Gauss points from the line to the volume. In case a line pokes out of the volumes it is
   * segmented.
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjectionCylinder
      : public GeometryPairLineToVolume<scalar_type, line, volume>
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolumeGaussPointProjectionCylinder()
        : GeometryPairLineToVolume<scalar_type, line, volume>(){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolumeGaussPointProjectionCylinder(){};

    /**
     * \brief Setup the projection tracker for this line element.
     */
    void Setup();

    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call PreEvaluate before
     * Evaluate is called on one of them.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void PreEvaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const;

    /**
     * \brief The only purpose of this method is to check that all points projected valid for this
     * beam in PreEvaluate.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void Evaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
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
