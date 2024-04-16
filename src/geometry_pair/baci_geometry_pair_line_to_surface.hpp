/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_HPP


#include "baci_config.hpp"

#include "baci_geometry_pair.hpp"
#include "baci_geometry_pair_element.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

namespace
{
  class GeometryPairLineToSurfaceTest;
}

BACI_NAMESPACE_OPEN

// Forward declarations.
namespace GEOMETRYPAIR
{
  enum class ProjectionResult;

  template <typename scalar_type>
  class ProjectionPoint1DTo3D;

  class LineToSurfaceEvaluationData;
}  // namespace GEOMETRYPAIR
namespace CORE::LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}  // namespace CORE::LINALG
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
    GeometryPairLineToSurface(const DRT::Element* element1, const DRT::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>&
            line_to_surface_evaluation_data);


    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call PreEvaluate before
     * Evaluate is called on one of them.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    virtual void PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
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
     */
    void ProjectPointToOther(const CORE::LINALG::Matrix<3, 1, scalar_type>& point,
        const ElementData<surface, scalar_type>& element_data_surface,
        CORE::LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const bool min_one_iteration = false) const;

    /**
     * \brief Intersect a line with all edges of a surface.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param eta_start (in) start value for parameter coordinate on line.
     * @param xi_start (in) start values for parameter coordinates in volume.
     */
    void IntersectLineWithOther(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
        const scalar_type& eta_start,
        const CORE::LINALG::Matrix<3, 1, scalar_type>& xi_start) const;

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
     * \brief Evaluate a position on the surface and its derivative w.r.t. xi.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
     * parameter coordiantes, the third one is in the normal direction).
     * @param r (out) Position on the surface.
     * @param dr (out) Derivative of the position on the surface, w.r.t xi.
     */
    void EvaluateSurfacePositionAndDerivative(
        const ElementData<surface, scalar_type>& element_data_surface,
        const CORE::LINALG::Matrix<3, 1, scalar_type>& xi,
        CORE::LINALG::Matrix<3, 1, scalar_type>& r,
        CORE::LINALG::Matrix<3, 3, scalar_type>& dr) const;

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
     */
    void IntersectLineWithSurfaceEdge(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
        CORE::LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
        const bool min_one_iteration = false) const;

    /**
     * \brief Check if a projection to a surface lies in a variable range.
     *
     * The ValidParameter2D function just checks, if the parameter coordinates are within the
     * surface. We need to additionally check that the normal distance is not to large in order to
     * avoid unwanted results. This is done by using a rule of thumb, that the normal distance has
     * be smaller than the approximated surface size or 3 times the beam radius.
     *
     * @param xi
     * @param surface_size
     * @param beam_radius
     * @return True if normal distance is in a reasonable range, false otherwise.
     */
    bool ValidParameterSurface(CORE::LINALG::Matrix<3, 1, scalar_type>& xi,
        const scalar_type& surface_size, const double beam_radius) const;

    /**
     * \brief Get an approximate size of the line radius.
     * @return Radius of the line.
     */
    double GetLineRadius() const;

    /**
     * \brief Get an approximate size of the surface.
     * @param (in) element_data_surface
     * @return Maximum distance between the first 3 nodes.
     */
    scalar_type GetSurfaceSize(const ElementData<surface, scalar_type>& element_data_surface) const;

    /**
     * \brief Get number of faces for this surface and create a vector with the indices of the
     * faces, so all surfaces can be checked for an intersection with the line.
     * @param n_faces (out) Number of faces.
     * @param face_fixed_parameters (out) Index of fixed parameter.
     * @param face_fixed_values (out) Value of fixed parameter.
     */
    void GetFaceFixedParameters(unsigned int& n_faces,
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
    GeometryPairLineToSurfaceFADWrapper(const DRT::Element* element1, const DRT::Element* element2,
        const Teuchos::RCP<GeometryPairLineToSurface<double, line, surface>>& double_geometry_pair)
        : base_class(element1, element2, double_geometry_pair->GetEvaluationData()),
          geometry_pair_double_(double_geometry_pair){};

    /**
     * \brief Set the pointer to the second element.
     *
     * For surface elements the pairs need the face element representing the surface, not the volume
     * element.
     */
    inline void SetElement2(const DRT::Element* element2) override
    {
      base_class::SetElement2(element2);
      geometry_pair_double_->SetElement2(element2);
    };

    /**
     * \brief Wrap the PreEvaluate call.
     *
     * The returned segment will be directly converted from double to the FAD type, i.e., the
     * derivatives will be WRONG. The correct derivatives are set in the Evaluate wrapper.
     *
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_surface (in) Degrees of freedom for the surface.
     * @param segments (out) Vector with the segments of this line to surface pair.
     */
    void PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
        const ElementData<surface, scalar_type>& element_data_surface,
        std::vector<LineSegment<scalar_type>>& segments) const override;

    /**
     * \brief Wrap the PreEvaluate call.
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

}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE

#endif
