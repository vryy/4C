/*!
\file geometry_pair_line_to_volume.H

\brief Class for interaction of lines and volumes.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_VOLUME_H
#define GEOMETRY_PAIR_LINE_TO_VOLUME_H


#include "geometry_pair.H"

#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace LINALG
{
  template <class value_type, unsigned int rows, unsigned int cols>
  class TMatrix;
}
namespace DRT
{
  class Element;

  namespace UTILS
  {
    struct IntegrationPoints1D;
  }
}  // namespace DRT
namespace GEOMETRYPAIR
{
  enum class DiscretizationTypeVolume;

  enum class ProjectionResult;

  template <typename scalar_type>
  class ProjectionPointLineToVolume;

  template <typename scalar_type>
  class LineSegment;

  class GeometryEvaluationDataGlobal;
}  // namespace GEOMETRYPAIR


// Constant values.
const unsigned int LOCALNEWTON_MAX_ITER = 20;
const double LOCALNEWTON_RES_TOL = 1e-10;
const double LOCALNEWTON_DET_TOL = 1e-10;
const double PROJECTION_XI_ETA_TOL = 1e-7;


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line and a volume.
   * @param scalar_type Type that will be used for scalar values.
   * @param n_nodes_element_1 Number of nodes on element 1.
   * @param n_nodal_values_element_1 Number of values per node (x3) on element 1.
   * @param n_nodes_element_2 Number of nodes on element 2.
   * @param n_nodal_values_element_2 Number of values per node (x3) on element 2.
   */
  template <typename scalar_type, unsigned int n_nodes_element_1,
      unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
      unsigned int n_nodal_values_element_2>
  class GeometryPairLineToVolume : public GeometryPair
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolume() : GeometryPair(){};

    /**
     * \brief Destructor.
     */
    ~GeometryPairLineToVolume(){};

    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call PreEvaluate before
     * Evaluate is called on one of them.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void PreEvaluate(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Evaluate the geometry interaction of the line and the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void Evaluate(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Get the position along the first (line) element. This function is templated in case
     * one wants to get the position with a different type than scalar_type.
     * @param xi (in) Parameter coordinate on the line.
     * @param q (in) Degrees of freedom for the line.
     * @param r (out) Position on line.
     */
    template <typename scalar_type_get_pos>
    void GetElement1Position(const scalar_type& xi,
        const LINALG::TMatrix<scalar_type_get_pos, 3 * n_nodes_element_1 * n_nodal_values_element_1,
            1>& q,
        LINALG::TMatrix<scalar_type_get_pos, 3, 1>& r) const;

    /**
     * \brief Get the derivative of the position w.r.t xi along the first (line) element.
     * @param xi (in) Parameter coordinate on the line.
     * @param q (in) Degrees of freedom for the line.
     * @param dr (out) Derivative of the position on the line.
     */
    void GetElement1PositionDerivative(const scalar_type& xi,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>& q,
        LINALG::TMatrix<scalar_type, 3, 1>& dr) const;

    /**
     * \brief Get the position in the second (volume) element. This function is templated in case
     * one wants to get the position with a different type than scalar_type.
     * @param xi (in) Parameter coordinates in the volume.
     * @param q (in) Degrees of freedom for the volume.
     * @param r (out) Position in the volume.
     */
    template <typename scalar_type_get_pos>
    void GetElement2Position(const LINALG::TMatrix<scalar_type, 3, 1>& xi,
        const LINALG::TMatrix<scalar_type_get_pos, 3 * n_nodes_element_2 * n_nodal_values_element_2,
            1>& q,
        LINALG::TMatrix<scalar_type_get_pos, 3, 1>& r) const;

    /**
     * \brief Get the derivative of the position w.r.t xi in the second (volume) element.
     * @param xi (in) Parameter coordinates in the volume.
     * @param q (in) Degrees of freedom for the volume.
     * @param dr (out) Derivative of the position in the volume.
     */
    void GetElement2PositionDerivative(const LINALG::TMatrix<scalar_type, 3, 1>& xi,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>& q,
        LINALG::TMatrix<scalar_type, 3, 3>& dr) const;

   protected:
    /**
     * \brief Project a point on the line to the volume element.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param eta (in) Parameter coordinate on the line.
     * @param xi (in/out) Parameter coordinates in the volume. The given values are the start values
     * for the Newton iteration.
     * @param projection_result (out) Flag for the result of the projection.
     * @return
     */
    void ProjectPointOnLineToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const scalar_type& eta, LINALG::TMatrix<scalar_type, 3, 1>& xi,
        ProjectionResult& projection_result) const;

    /**
     * \brief Project multiple points on the line to the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param eta_points (in) Vector with values for eta.
     * @param projection_points (out) Vector with found projection points.
     * @param allow_not_valid_projections (in) If this is true, point that lie outside of the valid
     * element range will also be returned.
     */
    void ProjectPointsOnLineToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const std::vector<scalar_type>& eta_points,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& projection_points,
        bool allow_not_valid_projections) const;

    /**
     * \brief Project multiple (equidistant) points on the line segment to the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param n_points (in) Number of points.
     * @param segment (in/out) Vector with found projection points.
     * @param allow_not_valid_projections (in) If this is true, point that lie outside of the valid
     * element range will also be returned.
     */
    void ProjectPointsOnSegmentToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        unsigned int n_points, LineSegment<scalar_type>& segment,
        bool allow_not_valid_projections) const;

    /**
     * \brief Project Gauss points on the line segment to the volume. If not all points can be
     * projected, an error is thrown. Only projections inside the two elements are considered valid.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param gauss_points (in) Gauss points to project.
     * @param segment (in/out) Vector with found projection points.
     */
    void ProjectPointsOnSegmentToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const DRT::UTILS::IntegrationPoints1D& gauss_points,
        LineSegment<scalar_type>& segment) const;

    /**
     * \brief Get the intersection between the line and a surface in the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param fixed_parameter (in) Index of parameter coordinate to be fixed on solid. In case
     * of tetraeder elements, fixed_parameter=3 represents the $r+s+t=1$ surface.
     * @param fixed_value (in) Value of fixed parameter.
     * @param eta (in/out) Parameter coordinate on the line. The given value is the start value
     * for the Newton iteration.
     * @param xi (in/out) Parameter coordinates in the volume. The given values are the start
     * values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the intersection.
     */
    void IntersectLineWithSurface(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
        LINALG::TMatrix<scalar_type, 3, 1>& xi, ProjectionResult& projection_result) const;

    /**
     * \brief Intersect a line with all surfaces of a volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param eta_start (in) start value for parameter coordinate on line.
     * @param xi_start (in) start values for parameter coordinates in volume.
     */
    void IntersectLineWithVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points,
        const scalar_type& eta_start, const LINALG::TMatrix<scalar_type, 3, 1>& xi_start) const;

    /**
     * \brief Intersect a line with all surfaces of a volume. Use default starting values for eta
     * and xi.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     */
    void IntersectLineWithVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points) const;

    /**
     * \brief Get the shape of the volume element.
     * @return Shape type of volume element.
     */
    DiscretizationTypeVolume GetVolumeType() const;

    /**
     * \brief Check if the parameter coordinate eta is in the valid range from -1 to 1.
     * @param eta Parameter coordiante on the line.
     * @return True if eta is between -1 and 1
     */
    bool ValidParameterElement1(const scalar_type& eta) const;

    /**
     * \brief Check if the parameter coordinates xi are in the valid range.
     * @param xi Parameter coordiantes in the volume.
     * @return True if xi is valid.
     */
    bool ValidParameterElement2(const LINALG::TMatrix<scalar_type, 3, 1>& xi) const;

    /**
     * \brief Set the parameter coordinate on the line to default starting value for a Newton
     * iteration.
     * @param eta Parameter coordiante on the line.
     */
    void SetStartValuesElement1(scalar_type& eta) const;

    /**
     * \brief Set the parameter coordinates in the volume to default starting values for a Newton
     * iteration.
     * @param xi Parameter coordiantes in the volume.
     */
    void SetStartValuesElement2(LINALG::TMatrix<scalar_type, 3, 1>& xi) const;
  };
}  // namespace GEOMETRYPAIR

#endif
