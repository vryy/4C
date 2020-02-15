/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and other geometry types.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_PROJECTION_H
#define GEOMETRY_PAIR_LINE_PROJECTION_H


#include <vector>
#include <set>


// Forward declarations.
namespace LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}
namespace DRT
{
  namespace UTILS
  {
    struct IntegrationPoints1D;
  }
}  // namespace DRT
namespace GEOMETRYPAIR
{
  enum class ProjectionResult;

  template <typename scalar_type>
  class ProjectionPoint1DTo3D;

  template <typename scalar_type>
  class LineSegment;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
  /**
   * \brief This class contains static methods for common lint-to-xxx algorithms.
   *
   * This class contains common methods for line-to-3D (volume or surface with normal direction)
   * interactions.
   *
   * The class methods have a variable number of input arguments (using the c++11 feature
   * typename...) by doing so we can handle a variable amount of arguments for the methods, e.g.
   * surface elements can have the optional argument for the averaged nodal normals.
   *
   * @tparam pair_type Class of the line-to-xxx pair that the segmentation / Gauss-point-projection
   * will be performed on.
   */
  template <typename pair_type>
  class LineTo3DBase
  {
   private:
    //! Alias for the scalar type.
    using scalar_type = typename pair_type::t_scalar_type;

    //! Alias for the line type.
    using line = typename pair_type::t_line;

    //! Alias for the other geometry type.
    using other = typename pair_type::t_other;

   public:
    /**
     * \brief Project a point on the line to the other geometry element.
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param eta (in) Parameter coordinate on the line.
     * @param xi (in/out) Parameter coordinates in the other geometry. The given values are the
     * start values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the projection.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void ProjectPointOnLineToOther(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other, const scalar_type& eta,
        LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        optional_type... optional_args);

    /**
     * \brief Project multiple points on the line to the other geometry. The value of eta and xi in
     * the projection points is the start value for the iteration.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param projection_points (in/out) Vector with projection points. The given values for eta and
     * xi are the start values for the iteration.
     * @param n_projections_valid (out) Number of valid projections.
     * @param n_projections (out) Number of points, where the nonlinear system could be solved.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void ProjectPointsOnLineToOther(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
        unsigned int& n_projections_valid, unsigned int& n_projections,
        optional_type... optional_args);

    /**
     * \brief Project multiple points on the line to the other geometry. The value of eta and xi in
     * the projection points is the start value for the iteration.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param projection_points (in/out) Vector with projection points. The given values for eta and
     * xi are the start values for the iteration.
     * @param n_projections_valid (out) Number of valid projections.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void ProjectPointsOnLineToOther(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
        unsigned int& n_projections_valid, optional_type... optional_args);

    /**
     * \brief Project Gauss points on the line segment to the other geometry. If not all points can
     * be projected, an error is thrown. Only projections inside the two elements are considered
     * valid.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param gauss_points (in) Gauss points to project.
     * @param segment (in/out) Vector with found projection points.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void ProjectGaussPointsOnSegmentToOther(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        const DRT::UTILS::IntegrationPoints1D& gauss_points, LineSegment<scalar_type>& segment,
        optional_type... optional_args);

    /**
     * \brief Intersect a line with all surfaces of an other geometry. Use default starting values
     * for eta and xi.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void IntersectLineWithOther(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
        optional_type... optional_args);
  };


  /**
   * \brief This class contains static methods for Gauss point projection algorithms.
   * @tparam pair_type Class of the line-to-xxx pair that the segmentation / Gauss-point-projection
   * will be performed on.
   */
  template <typename pair_type>
  class LineTo3DGaussPointProjection : public LineTo3DBase<pair_type>
  {
   private:
    //! Alias for the scalar type.
    using scalar_type = typename pair_type::t_scalar_type;

    //! Alias for the line type.
    using line = typename pair_type::t_line;

    //! Alias for the other geometry type.
    using other = typename pair_type::t_other;

   public:
    /**
     * \brief Try to project all Gauss points to the geometry.
     *
     * Only points are checked that do not
     * already have a valid projection in the projection tracker of the evaluation data container.
     * Eventually needed segmentation at lines poking out of the other geometry is done in the
     * Evaluate method.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the geometry.
     * @param segments (out) Vector with the segments of this line to xxx pair.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void PreEvaluate(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<LineSegment<scalar_type>>& segments, optional_type... optional_args);

    /**
     * \brief Check if a Gauss point projected valid for this pair in PreEvaluate.
     *
     * If so, all Gauss points have to project valid (in the tracker, since some can be valid on
     * other pairs). If not all project, the beam pokes out of the other geometry and in this method
     * segmentation is performed.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param segments (out) Vector with the segments of this line-to-xxx pair.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void Evaluate(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<LineSegment<scalar_type>>& segments, optional_type... optional_args);

   private:
    /**
     * \brief Get the line projection vector for the line element in the pair.
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @return Mutable reference to line projection vector.
     */
    static std::vector<bool>& GetLineProjectionVectorMutable(const pair_type* pair);
  };


  /**
   * \brief This class contains static methods for segmentation algorithms.
   * @tparam pair_type Class of the line-to-xxx pair that the segmentation / Gauss-point-projection
   * will be performed on.
   */
  template <typename pair_type>
  class LineTo3DSegmentation : public LineTo3DBase<pair_type>
  {
   private:
    //! Alias for the scalar type.
    using scalar_type = typename pair_type::t_scalar_type;

    //! Alias for the line type.
    using line = typename pair_type::t_line;

    //! Alias for the other geometry type.
    using other = typename pair_type::t_other;

   public:
    /**
     * \brief This method performs the segmentation of the line with the other geometry.
     *
     * First every search point on the beam is projected to the other geometry. For every search
     * point that has a projection (Newton converged, point does not have to be inside the other
     * geometry), surface intersections are checked for all surfaces of the other geometry (with the
     * search point as start point for the Newton iteration). The resulting surface intersections
     * points are sorted and only unique points are kept (since some intersections are found from
     * multiple starting points). Between two successive intersection points there is a line segment
     * (either inside or outside of the beam). To check whether the segment is part of this pair, a
     * point between two successive intersection points is projected, and if it projects valid, the
     * segment is assumed to be in this pair.
     *
     * There can be the special case when two consecutive segments are inside the pair, this happens
     * when the solid is exactly along the surface of a other geometry and a search point lies on
     * that surface.
     *
     * This approach should work for almost all possible cases. One case that can occur, is when
     * none of the search points project, but the line intersects the other geometry. In this case
     * more search points are the (expensive!) solution.
     *
     * In the case where a line lies on the interface between two geometries (e.g. the surface
     * between two volumes), segments on both other geometries are valid. To avoid this case a
     * global segment tracker is kept in the evaluation data, that allows a segment only once on all
     * line-to-xxx pairs. This does not work in the very special case, that the other geometry
     * have non-conforming discretization.
     *
     * @param pair (in) Pointer to the pair object that is being evaluated.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_other (in) Degrees of freedom for the other geometry.
     * @param segments (out) Vector with the segments of this line to other geometry pair.
     * @param optional_args (in) Optional - Normals on the nodes for surfaces.
     * @tparam optional_type Parameter packs template type for optional_args.
     */
    template <typename... optional_type>
    static void Evaluate(const pair_type* pair,
        const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
        const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
        std::vector<LineSegment<scalar_type>>& segments, optional_type... optional_args);

   private:
    /**
     * \brief Get the segment tracking vector for the line element the pair. The type is double
     * for the segments there, as this is only needed for comparison purposes. The segment values of
     * the tracker should NOT be used for further calculation as AD will not work correctly in this
     * case.
     * @return Reference to segment tracking vector.
     */
    static std::set<LineSegment<double>>& GetSegmentTrackingSetMutable(const pair_type* pair);
  };

}  // namespace GEOMETRYPAIR

#endif
