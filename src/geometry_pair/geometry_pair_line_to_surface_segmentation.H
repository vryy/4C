/*----------------------------------------------------------------------*/
/*! \file

\brief Line to surface interaction with segmentation at the boundaries.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_SURFACE_SEGMENTATION_H
#define GEOMETRY_PAIR_LINE_TO_SURFACE_SEGMENTATION_H


#include "geometry_pair_line_to_surface.H"


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
    GeometryPairLineToSurfaceSegmentation(
        const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>&
            line_to_surface_evaluation_data)
        : GeometryPairLineToSurface<scalar_type, line, surface>(line_to_surface_evaluation_data){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToSurfaceSegmentation(){};

    /**
     * \brief Setup the segment tracker for this line element.
     */
    void Setup() override;

    /**
     * \brief This method performs the segmentation of the line with the surface (and its normal
     * direction).
     *
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_surface (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    void Evaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        std::vector<LineSegment<scalar_type>>& segments,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals =
            NULL) const override;
  };
}  // namespace GEOMETRYPAIR

#endif
