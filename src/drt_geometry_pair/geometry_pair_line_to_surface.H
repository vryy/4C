/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_SURFACE_H
#define GEOMETRY_PAIR_LINE_TO_SURFACE_H


#include "geometry_pair.H"

#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace GEOMETRYPAIR
{
  class GeometryPairLineToSurface_TestSuite;
  enum class ProjectionResult;

  template <typename scalar_type>
  class ProjectionPoint1DTo3D;

  class LineTo3DEvaluationData;
}  // namespace GEOMETRYPAIR
namespace LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}  // namespace LINALG
template <typename scalar_type>
class LineSegment;



namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line (element 1) and a surface
   * (element 2).
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param surface Type of surface element.
   */
  template <typename scalar_type, typename line, typename surface>
  class GeometryPairLineToSurface : public GeometryPair
  {
    //! Declare the unit test class as a fried class, so private and protected methods can be used.
    friend GEOMETRYPAIR::GeometryPairLineToSurface_TestSuite;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToSurface(
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& line_to_3d_evaluation_data)
        : GeometryPair(), line_to_3d_evaluation_data_(line_to_3d_evaluation_data){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToSurface(){};

    /**
     * \brief Initialize the pair.
     * @param element1 Pointer to the first element.
     * @param element2 Pointer to the second element.
     */
    void Init(const DRT::Element* element1, const DRT::Element* element2) override;

    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call PreEvaluate before
     * Evaluate is called on one of them.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    virtual void PreEvaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Evaluate the geometry interaction of the line and the surface.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    virtual void Evaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        std::vector<LineSegment<scalar_type>>& segments) const {};

   protected:
    /**
     * \brief Project a point in space to the surface element.
     * @param point (in) Point in space.
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param xi (in/out) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction). The given values are the
     * start values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the projection.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    void ProjectPointToSurface(const LINALG::Matrix<3, 1, scalar_type>& point,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

    /**
     * \brief Evaluate a position on the surface.
     *
     * @template scalar_type_evaluate Scalar type used for this function. This type is different
     * from the scalar type in the class, since we calculate derivatives of the this method inside
     * of other methods.
     *
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction).
     * @param r (out) Position on the surface.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    template <typename scalar_type_evaluate>
    void EvaluateSurfacePosition(const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
        LINALG::Matrix<3, 1, scalar_type_evaluate>& r,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

    /**
     * \brief Evaluate a position on the surface and its derivative w.r.t. xi.
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction).
     * @param r (out) Position on the surface.
     * @param dr (out) Derivative of the position on the surface, w.r.t xi.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    void EvaluateSurfacePositionAndDerivative(
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        const LINALG::Matrix<3, 1, scalar_type>& xi, LINALG::Matrix<3, 1, scalar_type>& r,
        LINALG::Matrix<3, 3, scalar_type>& dr,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

    /**
     * \brief Evaluate the normal at a position on the surface.
     *
     * @template scalar_type_evaluate Scalar type used for this function. This type is different
     * from the scalar type in the class, since we calculate derivatives of the this method inside
     * of other methods.
     *
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction).
     * @param normal (out) Normal on the surface.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    template <typename scalar_type_evaluate>
    void EvaluateSurfaceNormal(const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
        LINALG::Matrix<3, 1, scalar_type_evaluate>& normal,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

    /**
     * \brief Get the intersection between the line and an extended boundary of the surface.
     *
     * The extended boundaries are the edges of the surface, extended by the normal vector there.
     *
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_surface (in) Degrees of freedom for the surface.
     * @param fixed_parameter (in) Index of parameter coordinate to be fixed on the surface. In case
     * of triangle elements, fixed_parameter=2 represents the $r+s=1$ curve.
     * @param fixed_value (in) Value of fixed parameter.
     * @param eta (in/out) Parameter coordinate on the line. The given value is the start value
     * for the Newton iteration.
     * @param xi (in/out) Parameter coordinates for the surface. The given values are the start
     * values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the intersection.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    void IntersectLineWithSurfaceEdge(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
        const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
        LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

    /**
     * \brief Intersect a line with all edges of a surface.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_surface (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param eta_start (in) start value for parameter coordinate on line.
     * @param xi_start (in) start values for parameter coordinates in volume.
     * @param nodal_normals (in) Optional - Normals on the nodes.
     */
    void IntersectLineWithSurface(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_volume,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
        const scalar_type& eta_start, const LINALG::Matrix<3, 1, scalar_type>& xi_start,
        const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const;

   protected:
    //! Link to the geometry evaluation container.
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> line_to_3d_evaluation_data_;
  };
}  // namespace GEOMETRYPAIR

#endif
