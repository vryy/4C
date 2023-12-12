/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_UTILITY_FUNCTIONS_H
#define BACI_GEOMETRY_PAIR_UTILITY_FUNCTIONS_H


#include "baci_beam3_base.H"
#include "baci_geometry_pair_utility_classes.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_utils_fad.H"

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


  /**
   * @brief Print the current nodal normals of a surface element
   *
   * As the nodal normal are passed through the geometry pair functions via
   * "typename... optional_type" template arguments, this function will print something
   * in case nodal normals are provided.
   */
  template <typename other, typename scalar_type, typename... optional_type>
  void PrintNodalNormals(std::ostream& out,
      const CORE::LINALG::Matrix<3 * other::n_nodes_, 1, scalar_type>* nodal_normals,
      optional_type... optional_args)
  {
    if (nodal_normals == nullptr)
    {
      out << "No nodal normals are given\n";
    }
    else
    {
      out << "Averaged nodal normals: ";
      const auto& normal = *nodal_normals;
      CORE::FADUTILS::CastToDouble(normal).Print(out);
    }
  }

  /**
   * @brief Print the current nodal normals of a surface element
   *
   * As the nodal normal are passed through the geometry pair functions via
   * "typename... optional_type" template arguments, this function will print nothing if the other
   * element type does not provide nodal normals.
   */
  template <typename other, typename... optional_type>
  void PrintNodalNormals(std::ostream& out, optional_type... optional_args)
  {
  }

  /**
   * @brief Add the current displacement state of a pair to an output stream
   */
  template <typename pair_type, typename scalar_type, typename... optional_type>
  void PrintPairInformation(std::ostream& out, const pair_type* pair,
      const CORE::LINALG::Matrix<pair_type::t_line::n_dof_, 1, scalar_type>& q_line,
      const CORE::LINALG::Matrix<pair_type::t_other::n_dof_, 1, scalar_type>& q_other,
      optional_type... optional_args)
  {
    using other = typename pair_type::t_other;

    constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
    out << std::setprecision(max_precision);

    // Print the line information
    const auto* beam_element = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(pair->Element1());
    if (beam_element == nullptr)
      dserror(
          "The element pointer has to point to a valid beam element when evaluating the shape "
          "functions of a beam, as we need to get RefLength()!");

    out << "\nLine GID: " << pair->Element1()->Id();
    out << "\nLine ref length value: " << beam_element->RefLength();
    out << "\nLine position vector: ";
    CORE::FADUTILS::CastToDouble(q_line).Print(out);

    // Print the other information
    out << "Other GID: ";
    const auto* face_element = dynamic_cast<const DRT::FaceElement*>(pair->Element2());
    if (face_element == nullptr)
    {
      out << pair->Element2()->Id();
    }
    else
    {
      out << face_element->ParentElementId();
    }
    out << "\nOther position vector: ";
    CORE::FADUTILS::CastToDouble(q_other).Print(out);
    PrintNodalNormals<other>(out, optional_args...);
  }

}  // namespace GEOMETRYPAIR


BACI_NAMESPACE_CLOSE

#endif
