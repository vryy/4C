/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_UTILITY_FUNCTIONS_HPP
#define FOUR_C_GEOMETRY_PAIR_UTILITY_FUNCTIONS_HPP


#include "baci_config.hpp"

#include "baci_geometry_pair_element.hpp"
#include "baci_geometry_pair_utility_classes.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_utils_fad.hpp"

#include <string>

BACI_NAMESPACE_OPEN

// Forward declarations.
namespace GEOMETRYPAIR
{
  enum class DiscretizationTypeGeometry;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Convert the enum DiscretizationTypeGeometry to a human readable string.
   * @param discretization_type (in)
   * @return Human readable string representation of the enum.
   */
  std::string DiscretizationTypeGeometryToString(
      const DiscretizationTypeGeometry discretization_type);

  /**
   * \brief Set a line-to-xxx segment from a segment with a different scalar type (all dependencies
   * of FAD variables will be deleted).
   *
   * @param vector_scalar_type_ptr (in) Pointer to a vector of arbitrary scalar type.
   * @param vector_double (out) Temp reference to double vector (this has to come from the outside
   * scope).
   * @return Pointer on the double vector.
   *
   * @tparam A Scalar type of segment_in.
   * @tparam B Scalar type of segment_out.
   */
  template <typename A, typename B>
  void CopySegment(const LineSegment<A>& segment_in, LineSegment<B>& segment_out)
  {
    // Add the start and end points.
    segment_out.GetStartPoint().SetFromOtherPointDouble(segment_in.GetStartPoint());
    segment_out.GetEndPoint().SetFromOtherPointDouble(segment_in.GetEndPoint());

    // Add the projection points.
    const auto n_points = segment_in.GetNumberOfProjectionPoints();
    const std::vector<ProjectionPoint1DTo3D<A>>& projection_points_in =
        segment_in.GetProjectionPoints();
    std::vector<ProjectionPoint1DTo3D<B>>& projection_points_out =
        segment_out.GetProjectionPoints();
    projection_points_out.resize(n_points);
    for (unsigned int i_point = 0; i_point < n_points; i_point++)
      projection_points_out[i_point].SetFromOtherPointDouble(projection_points_in[i_point]);
  }

  /**
   * @brief Add the current displacement state of a pair to an output stream
   */
  template <typename pair_type, typename scalar_type, typename... optional_type>
  void PrintPairInformation(std::ostream& out, const pair_type* pair,
      const ElementData<typename pair_type::t_line, scalar_type>& element_data_line,
      const ElementData<typename pair_type::t_other, scalar_type>& element_data_other)
  {
    using line = typename pair_type::t_line;
    using other = typename pair_type::t_other;

    constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
    out << std::setprecision(max_precision);

    const auto* face_element = dynamic_cast<const DRT::FaceElement*>(pair->Element2());
    if (face_element == nullptr)
    {
      out << "Pair consisting of the line with GID " << pair->Element1()->Id()
          << " and the other element GID " << pair->Element2()->Id() << "\n";
    }
    else
    {
      out << "Pair consisting of the line with GID " << pair->Element1()->Id()
          << " and the other element " << pair->Element2()->Id() << " with the parent element GID "
          << face_element->ParentElementId() << "\n";
    }
    out << "Line:";
    PrintElementData<line>::Print(element_data_line, out);
    out << "Other:";
    PrintElementData<other>::Print(element_data_other, out);
  }

}  // namespace GEOMETRYPAIR


BACI_NAMESPACE_CLOSE

#endif
