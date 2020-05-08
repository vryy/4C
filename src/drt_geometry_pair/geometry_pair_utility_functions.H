/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_UTILITY_FUNCTIONS_H
#define GEOMETRY_PAIR_UTILITY_FUNCTIONS_H


#include <string>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../headers/FAD_utils.H"


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
   * \brief Convert a pointer to a fixed size vector to a reference to a double vector.
   *
   * @param vector_scalar_type_ptr (in) Pointer to a vector of arbitrary scalar type.
   * @param vector_double (out) Reference to double vector.
   * @param valid_pointer (out) True if the given vector was a valid pointer, false otherwise.
   *
   * @tparam scalar_type_normal Scalar type of the input vector.
   * @tparam n_dof Length of the in- and output vectors.
   */
  template <typename scalar_type_normal, unsigned int n_dof>
  void VectorPointerToVectorDouble(
      const LINALG::Matrix<n_dof, 1, scalar_type_normal>* vector_scalar_type_ptr,
      LINALG::Matrix<n_dof, 1, double>& vector_double, bool& valid_pointer)
  {
    if (vector_scalar_type_ptr == nullptr)
    {
      valid_pointer = false;
    }
    else
    {
      valid_pointer = true;
      vector_double = FADUTILS::CastToDouble(*vector_scalar_type_ptr);
    }
  };
}  // namespace GEOMETRYPAIR


#endif
