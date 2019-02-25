/*!
\file geometry_pair_line_to_volume.H

\brief Class for interaction of lines and volumes.

\level 1
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


namespace GEOMETRYPAIR
{
  /**
   * \brief Class that handles the geometrical interactions of a line (element 1) and a volume
   * (element 2).
   * @param scalar_type Type that will be used for scalar values.
   * @param line Type of line element.
   * @param volume Type of volume element.
   */
  template <typename scalar_type, typename line, typename volume>
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
     * \brief Initialize the pair.
     * @param evaluation_data_ptr Pointer to the evaluation data.
     * @param element1 Pointer to the first element.
     * @param element2 Pointer to the second element.
     */
    void Init(Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_ptr,
        const DRT::Element* element1, const DRT::Element* element2) override;

    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call PreEvaluate before
     * Evaluate is called on one of them.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void PreEvaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const {};

    /**
     * \brief Evaluate the geometry interaction of the line and the volume.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void Evaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<LineSegment<scalar_type>>& segments) const {};

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
    void ProjectPointOnLineToVolume(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume, const scalar_type& eta,
        LINALG::TMatrix<scalar_type, 3, 1>& xi, ProjectionResult& projection_result) const;

    /**
     * \brief Project multiple points on the line to the volume. The value of eta and xi in the
     * projection points is the start value for the iteration.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param projection_points (in/out) Vector with projection points. The given values for eta and
     * xi are the start values for the iteration.
     * @param n_projections_valid (out) Number of valid projections.
     * @param n_projections (out) Number of points, where the nonlinear system could be solved.
     */
    void ProjectPointsOnLineToVolume(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& projection_points,
        unsigned int& n_projections_valid, unsigned int& n_projections) const;

    /**
     * \brief Project multiple points on the line to the volume. The value of eta and xi in the
     * projection points is the start value for the iteration.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param projection_points (in/out) Vector with projection points. The given values for eta and
     * xi are the start values for the iteration.
     * @param n_projections_valid (out) Number of valid projections.
     */
    void ProjectPointsOnLineToVolume(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& projection_points,
        unsigned int& n_projections_valid) const;

    /**
     * \brief Project Gauss points on the line segment to the volume. If not all points can be
     * projected, an error is thrown. Only projections inside the two elements are considered valid.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param gauss_points (in) Gauss points to project.
     * @param segment (in/out) Vector with found projection points.
     */
    void ProjectGaussPointsOnSegmentToVolume(
        const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
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
    void IntersectLineWithSurface(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
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
    void IntersectLineWithVolume(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points,
        const scalar_type& eta_start, const LINALG::TMatrix<scalar_type, 3, 1>& xi_start) const;

    /**
     * \brief Intersect a line with all surfaces of a volume. Use default starting values for eta
     * and xi.
     * @param q_line (in) Degrees of freedom for the line.
     * @param q_volume (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     */
    void IntersectLineWithVolume(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
        const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points) const;

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
