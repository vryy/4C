/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_UTILITY_FUNCTIONS_H
#define BACI_GEOMETRY_PAIR_UTILITY_FUNCTIONS_H


#include "baci_geometry_pair_utility_classes.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_utils_fad.H"

#include <string>


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
   * \brief Convert a pointer on a FAD vector to a pointer on a double vector.
   *
   * @param vector_scalar_type_ptr (in) Pointer to a vector of arbitrary scalar type.
   * @param vector_double (out) Temp reference to double vector (this has to come from the outside
   * scope).
   * @return Pointer on the double vector.
   *
   * @tparam scalar_type_normal Scalar type of the input vector.
   * @tparam n_dof Length of the in- and output vectors.
   */
  template <typename scalar_type_normal, unsigned int n_dof>
  CORE::LINALG::Matrix<n_dof, 1, double>* VectorPointerToVectorDouble(
      const CORE::LINALG::Matrix<n_dof, 1, scalar_type_normal>* vector_scalar_type_ptr,
      CORE::LINALG::Matrix<n_dof, 1, double>& vector_double)
  {
    if (vector_scalar_type_ptr == nullptr)
    {
      return nullptr;
    }
    else
    {
      vector_double = CORE::FADUTILS::CastToDouble(*vector_scalar_type_ptr);
      return &vector_double;
    }
  }

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
}  // namespace GEOMETRYPAIR


#endif
