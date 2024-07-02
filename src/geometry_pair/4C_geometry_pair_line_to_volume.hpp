/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and volumes.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_VOLUME_HPP

#include "4C_config.hpp"

#include "4C_geometry_pair.hpp"
#include "4C_geometry_pair_element.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace Core::LinAlg
{
  template <unsigned int rows, unsigned int cols, class ValueType>
  class Matrix;
}
namespace Discret
{
  namespace UTILS
  {
    struct IntegrationPoints1D;
  }
}  // namespace Discret

namespace Core::Elements
{
  class Element;
}

namespace GEOMETRYPAIR
{
  enum class DiscretizationTypeVolume;

  enum class ProjectionResult;

  template <typename ScalarType>
  class ProjectionPoint1DTo3D;

  template <typename ScalarType>
  class LineSegment;

  class LineTo3DEvaluationData;
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
  template <typename ScalarType, typename Line, typename Volume>
  class GeometryPairLineToVolume : public GeometryPair
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPairLineToVolume(const Core::Elements::Element* element1,
        const Core::Elements::Element* element2,
        const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& line_to_3d_evaluation_data);


    /**
     * \brief Do stuff that can not be done in the Evaluate call. All pairs call pre_evaluate before
     * Evaluate is called on one of them.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void pre_evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments) const {};

    /**
     * \brief Evaluate the geometry interaction of the line and the volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param segments (out) Vector with the segments of this line to volume pair.
     */
    virtual void evaluate(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<LineSegment<ScalarType>>& segments) const {};

    /**
     * \brief Return the pointer to the evaluation data of this pair.
     * @return Pointer to the evaluation data.
     */
    const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& GetEvaluationData() const
    {
      return line_to_3d_evaluation_data_;
    }

    /**
     * \brief Project a point in space to the volume element.
     * @param point (in) Point in space.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param xi (in/out) Parameter coordinates in the volume. The given values are the start values
     * for the Newton iteration.
     * @param projection_result (out) Flag for the result of the projection.
     * @return
     */
    void ProjectPointToOther(const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
        const ElementData<Volume, ScalarType>& element_data_volume,
        Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result) const;

    /**
     * \brief Get the intersection between the line and a surface in the volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param fixed_parameter (in) Index of parameter coordinate to be fixed on solid. In case
     * of tetraeder elements, fixed_parameter=3 represents the $r+s+t=1$ surface.
     * @param fixed_value (in) Value of fixed parameter.
     * @param eta (in/out) Parameter coordinate on the line. The given value is the start value
     * for the Newton iteration.
     * @param xi (in/out) Parameter coordinates in the volume. The given values are the start
     * values for the Newton iteration.
     * @param projection_result (out) Flag for the result of the intersection.
     */
    void intersect_line_with_surface(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        const unsigned int& fixed_parameter, const double& fixed_value, ScalarType& eta,
        Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result) const;

    /**
     * \brief Intersect a line with all surfaces of a volume.
     * @param element_data_line (in) Degrees of freedom for the line.
     * @param element_data_volume (in) Degrees of freedom for the volume.
     * @param intersection_points (out) vector with the found surface intersections.
     * @param eta_start (in) start value for parameter coordinate on line.
     * @param xi_start (in) start values for parameter coordinates in volume.
     */
    void intersect_line_with_other(const ElementData<Line, ScalarType>& element_data_line,
        const ElementData<Volume, ScalarType>& element_data_volume,
        std::vector<ProjectionPoint1DTo3D<ScalarType>>& intersection_points,
        const ScalarType& eta_start, const Core::LinAlg::Matrix<3, 1, ScalarType>& xi_start) const;

   protected:
    //! Link to the geometry evaluation container.
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> line_to_3d_evaluation_data_;
  };

  /**
   * \brief Project a point in space to the volume element.
   * @param point (in) Point in space.
   * @param element_data_volume (in) Degrees of freedom for the volume.
   * @param xi (in/out) Parameter coordinates in the volume. The given values are the start values
   * for the Newton iteration.
   * @param projection_result (out) Flag for the result of the projection.
   */
  template <typename ScalarType, typename Volume>
  void ProjectPointToVolume(const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
      const ElementData<Volume, ScalarType>& element_data_volume,
      Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result);
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
