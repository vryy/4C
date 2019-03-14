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
     * \brief Try to project all Gauss points to the volume. Only points are checked that do not
     * already have a valid projection in the projection tracker of the evaluation data container.
     * Eventually needed segmentation at lines poking out of the volume is done in the Evaluate
     * method.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void PreEvaluateCylinder(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<GEOMETRYPAIR::ProjectionPointVolumeToVolume<scalar_type>>&
            cylinder_to_volume_points) const;

    /**
     * \brief Check if a Gauss point projected valid for this pair in PreEvaluate. If so, all Gauss
     * points have to project valid (in the tracker, since some can be valid on other pairs). If not
     * all project, the beam pokes out of the volumes and in this method segmentation is performed.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    void EvaluateCylinder(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<GEOMETRYPAIR::ProjectionPointVolumeToVolume<scalar_type>>&
            cylinder_to_volume_points) const;

   private:
    /**
     * \brief Get the line projection vector for the line element in this pair.
     * @return Reference to line projection vector.
     */
    std::vector<bool>& GetLineProjectionVectorMutable() const;
  };
}  // namespace GEOMETRYPAIR

#endif
