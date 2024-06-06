/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

namespace
{
  class GeometryPairLineToSurfaceTest;
}

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace GEOMETRYPAIR
{
  enum class ProjectionResult;

  template <typename scalar_type>
  class ProjectionPoint1DTo3D;

  class LineToSurfaceEvaluationData;
}  // namespace GEOMETRYPAIR
namespace Core::LinAlg
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}  // namespace Core::LinAlg
namespace GEOMETRYPAIR
{
  template <typename scalar_type>
  class LineSegment;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line (element 1) and a surface
   * (element 2).
   * @tparam scalar_type Type that will be used for scalar values.
   * @tparam line Type of line element.
   * @tparam surface Type of surface element.
   */
  template <typename scalar_type, typename line, typename surface>
  class GeometryPairLineToSurface : public GeometryPair
  {
    //! Declare the unit test class as a fried class, so private and protected methods can be used.
    friend GeometryPairLineToSurfaceTest;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToSurface(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>&
            line_to_surface_evaluation_data);


    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call pre_evaluate before
     * Evaluate is called on one of them.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    virtual void pre_evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Evaluate the geometry interaction of the line and the surface.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    virtual void Evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Project a point in space to the surface element.
     * @param point (in) Point in space.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param xi (in/out) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction). The given values are the
     * start values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the projection.
     * @param min_one_iteration (in) Flag if at least one NR iteration should be performed, even if
     * the initial residual satisfies the convergence check.
     */
    void ProjectPointToOther(const Core::LinAlg::Matrix<3, 1, scalar_type>& point,
        const ElementData<surface, scalar_type>& element_data_surface,
        Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const bool min_one_iteration = false) const;

    /**
     * \brief Intersect a line with all edges of a surface.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param eta_start (in) start value for parameter coordinate on line.
     * @param xi_start (in) start values for parameter coordinates in volume.
     */
    void intersect_line_with_other(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
        const scalar_type& eta_start,
        const Core::LinAlg::Matrix<3, 1, scalar_type>& xi_start) const;

    /**
     * \brief Return the pointer to the evaluation data of this pair.
     * @return Pointer to the evaluation data.
     */
    inline const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>& GetEvaluationData() const
    {
      return line_to_surface_evaluation_data_;
    }

   protected:
    /**
     * \brief Get the intersection between the line and an extended boundary of the surface.
     *
     * The extended boundaries are the edges of the surface, extended by the normal vector there.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param fixed_parameter (in) Index of parameter coordinate to be fixed on the surface. In case
     * of triangle elements, fixed_parameter=2 represents the $r+s=1$ curve.
     * @param fixed_value (in) Value of fixed parameter.
     * @param eta (in/out) Parameter coordinate on the line. The given value is the start value
     * for the Newton iteration.
     * @param xi (in/out) Parameter coordinates for the surface. The given values are the start
     * values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the intersection.
     * @param min_one_iteration (in) Flag if at least one NR iteration should be performed, even if
     * the initial residual satisfies the convergence check.
     */
    void intersect_line_with_surface_edge(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
        Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const bool min_one_iteration = false) const;

    /**
     * \brief Get an approximate influence direction along the surface normal direction (both
     * outward and inward).
     */
    double get_surface_normal_influence_direction(
        const ElementData<surface, scalar_type>& element_data_surface) const;

    /**
     * \brief Get an approximate size of the surface.
     * @param (in) element_data_surface
     * @return Maximum distance between the first 3 nodes.
     */
    double get_surface_size(const ElementData<surface, scalar_type>& element_data_surface) const;

    /**
     * \brief Get number of faces for this surface and create a vector with the indices of the
     * faces, so all surfaces can be checked for an intersection with the line.
     * @param n_faces (out) Number of faces.
     * @param face_fixed_parameters (out) Index of fixed parameter.
     * @param face_fixed_values (out) Value of fixed parameter.
     */
    void get_face_fixed_parameters(unsigned int& n_faces,
        std::vector<unsigned int>& face_fixed_parameters,
        std::vector<double>& face_fixed_values) const;

   protected:
    //! Link to the geometry evaluation container.
    Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData> line_to_surface_evaluation_data_;

   private:
    //! Flag if the class is executed by unit tests. If this is the case, the line radius for
    //! ValidParameterSurface will be explicitly given.
    bool is_unit_test_;
  };


  /**
   * \brief Wrapper for line to surface pairs with FAD scalar types.
   *
   * The default GeometryPairLineToSurface can also be used for this. This wrapper exists for
   * performance reasons. In the geometry pairs we have a lot of local Newton iterations. If they
   * are performed with the full FAD types, the performance drops significantly. Therefore, we use
   * geometry pairs with double type to evaluate all intersections and then only perform a single
   * Netwon iteration (with the previously calculated result as starting point) to get the correct
   * FAD dependencies.
   *
   * @tparam scalar_type Type that will be used for scalar values.
   * @tparam line Type of line element.
   * @tparam surface Type of surface element.
   */
  template <typename scalar_type, typename line, typename surface>
  class GeometryPairLineToSurfaceFADWrapper
      : public GeometryPairLineToSurface<scalar_type, line, surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = GeometryPairLineToSurface<scalar_type, line, surface>;

   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToSurfaceFADWrapper(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const Teuchos::RCP<GeometryPairLineToSurface<double, line, surface>>& double_geometry_pair)
        : base_class(element1, element2, double_geometry_pair->GetEvaluationData()),
          geometry_pair_double_(double_geometry_pair){};

    /**
     * \brief Set the pointer to the second element.
     *
     * For surface elements the pairs need the face element representing the surface, not the volume
     * element.
     */
    inline void SetElement2(const Core::Elements::Element* element2) override
    {
      base_class::SetElement2(element2);
      geometry_pair_double_->SetElement2(element2);
    };

    /**
     * \brief Wrap the pre_evaluate call.
     *
     * The returned segment will be directly converted from double to the FAD type, i.e., the
     * derivatives will be WRONG. The correct derivatives are set in the Evaluate wrapper.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    void pre_evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const override;

    /**
     * \brief Wrap the pre_evaluate call.
     *
     * The segments will be evaluated as double segments and only the converged projections /
     * intersections will be revaluated to contain the correct FAD derivatives. This process
     * reduces the performance bottleneck introduced by using the FAD types.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    void Evaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const override;

   private:
    //! Pair to evaluate the intersections with the scalar type double.
    Teuchos::RCP<GeometryPairLineToSurface<double, line, surface>> geometry_pair_double_;
  };

  /**
   * \brief Check if a projection to a surface lies in a variable range.
   *
   * The ValidParameter2D function just checks, if the parameter coordinates are within the
   * surface. We need to additionally check that the normal distance is not to large in order to
   * avoid unwanted results. This is done by comparing the normal distance to a user given
   * normal_influence_direction.
   *
   * @param xi (in) Parameter coordinate
   * @param normal_influence_direction (in) Threshold value for influence distance in normal
   * direction. If this is negative no threshold will be applied and the only check performed, is if
   * the projection lies withing the governing 2D parameter space.
   * @return True if normal distance is in a reasonable range, false otherwise.
   */
  template <typename scalar_type, typename surface>
  bool ValidParameterSurface(
      Core::LinAlg::Matrix<3, 1, scalar_type>& xi, const double normal_influence_direction)
  {
    // We only need to check the normal distance if the coordinates are within the surface.
    if (!ValidParameter2D<surface>(xi)) return false;

    if (normal_influence_direction < 0)
    {
      // No additional check in this case
      return true;
    }
    else
    {
      // Check if the normal direction is in a reasonable user given range
      return (-normal_influence_direction < xi(2) and xi(2) < normal_influence_direction);
    }
  }

  /**
   * \brief Project a point in space to the surface element.
   * @param point (in) Point in space.
   * @param element_data_surface (in) Degrees of freedom for the surface.
   * @param xi (in/out) Parameter coordinates in the volume. The given values are the start values
   * for the Newton iteration.
   * @param projection_result (out) Flag for the result of the projection.
   * @param normal_influence_direction (in) Threshold of normal influence direction for the surface.
   * @param min_one_iteration (in) Flag if at least one NR iteration should be performed, even if
   * the initial residual satisfies the convergence check.
   */
  template <typename scalar_type, typename surface>
  void ProjectPointToSurface(const Core::LinAlg::Matrix<3, 1, scalar_type>& point,
      const ElementData<surface, scalar_type>& element_data_surface,
      Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
      const double normal_influence_direction = -1.0, const bool min_one_iteration = false);
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
