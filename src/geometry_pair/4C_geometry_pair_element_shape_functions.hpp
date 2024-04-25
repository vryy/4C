/*----------------------------------------------------------------------*/
/*! \file

\brief Specify how to evaluate the shape functions for the elements defined in
baci_geometry_pair_element.H

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_ELEMENT_SHAPE_FUNCTIONS_HPP
#define FOUR_C_GEOMETRY_PAIR_ELEMENT_SHAPE_FUNCTIONS_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_geometry_pair_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  /**
   * \brief Struct to evaluate the shape functions for an element
   */
  template <typename element_type, typename enable = void>
  struct EvaluateShapeFunction
  {
    // Empty per default
  };

  /**
   * \brief Specialization for Lagrange elements and elements with averaged nodal normals (face
   * elements which require averaged nodal normals are currently always based on Lagrange elements
   * in BACI)
   */
  template <typename element_type>
  struct EvaluateShapeFunction<element_type,
      typename std::enable_if<IsLagrangeElement<element_type>::value_ ||
                              IsSurfaceAveragedNormalsElement<element_type>::value_>::type>
  {
    /**
     * \brief Evaluate the shape functions of the element at xi.
     *
     * The elements considered here, i.e., Lagrangian elements don't have specific shape function
     * data, so we allow this as an optional argument
     *
     * @param N (out) shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T, typename... not_needed_argument_type>
    static void Evaluate(V& N, const T& xi, const not_needed_argument_type&... not_needed_argument)
    {
      if constexpr (element_type::element_dim_ == 1)
      {
        CORE::FE::shape_function_1D(N, xi, element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 2)
      {
        CORE::FE::shape_function_2D(N, xi(0), xi(1), element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 3)
      {
        CORE::FE::shape_function_3D(N, xi(0), xi(1), xi(2), element_type::discretization_);
      }
      else
      {
        FOUR_C_THROW("Got unexpected element dimension %d", element_type::element_dim_);
      }
    }

    /**
     * \brief Evaluate the derivatives of the shape functions of the element at xi.
     *
     * The elements considered here, i.e., Lagrangian elements don't have specific shape function
     * data, so we allow this as an optional argument
     *
     * @param dN (out) derivatives of shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T, typename... not_needed_argument_type>
    static void EvaluateDeriv1(
        V& dN, const T& xi, const not_needed_argument_type&... not_needed_argument)
    {
      if constexpr (element_type::element_dim_ == 1)
      {
        CORE::FE::shape_function_1D_deriv1(dN, xi, element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 2)
      {
        CORE::FE::shape_function_2D_deriv1(dN, xi(0), xi(1), element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 3)
      {
        CORE::FE::shape_function_3D_deriv1(dN, xi(0), xi(1), xi(2), element_type::discretization_);
      }
      else
      {
        FOUR_C_THROW("Got unexpected element dimension %d", element_type::element_dim_);
      }
    }
  };

  /**
   * \brief Specialization for Hermite elements
   */
  template <>
  struct EvaluateShapeFunction<t_hermite>
  {
    /**
     * \brief Evaluate the shape functions of the element at xi.
     *
     * @param N (out) shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T>
    static void Evaluate(V& N, const T& xi, const ShapeFunctionData<t_hermite>& shape_function_data)
    {
      CORE::FE::shape_function_hermite_1D(
          N, xi, shape_function_data.ref_length_, t_hermite::discretization_);
    }

    /**
     * \brief Evaluate the derivatives of the shape functions of the element at xi.
     *
     * @param dN (out) derivatives of shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T>
    static void EvaluateDeriv1(
        V& dN, const T& xi, const ShapeFunctionData<t_hermite>& shape_function_data)
    {
      CORE::FE::shape_function_hermite_1D_deriv1(
          dN, xi, shape_function_data.ref_length_, t_hermite::discretization_);
    }
  };

  /**
   * \brief Specialization for NURBS elements
   */
  template <typename element_type>
  struct EvaluateShapeFunction<element_type,
      typename std::enable_if<IsNurbsElement<element_type>::value_>::type>
  {
    /**
     * \brief Evaluate the shape functions of the element at xi.
     *
     * @param N (out) shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T>
    static void Evaluate(
        V& N, const T& xi, const ShapeFunctionData<element_type>& shape_function_data)
    {
      if (shape_function_data.myknots_.size() == 0)
        FOUR_C_THROW(
            "Got shape function data with knot size 0 - did you forget to initialize the element "
            "data with GetElementData::Get?");

      if constexpr (element_type::element_dim_ == 2)
      {
        CORE::FE::NURBS::nurbs_get_2D_funct<typename V::scalar_type>(N, xi,
            shape_function_data.myknots_, shape_function_data.weights_,
            element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 3)
      {
        CORE::FE::NURBS::nurbs_get_3D_funct(N, xi, shape_function_data.myknots_,
            shape_function_data.weights_, element_type::discretization_);
      }
      else
      {
        FOUR_C_THROW("Got unexpected element dimension %d", element_type::element_dim_);
      }
    }

    /**
     * \brief Evaluate the derivatives of the shape functions of the element at xi.
     *
     * @param dN (out) derivatives of shape functions.
     * @param xi (in) Parameter coordinate on the element.
     * @param shape_function_data (in) Shape function data container.
     */
    template <typename V, typename T>
    static void EvaluateDeriv1(
        V& dN, const T& xi, const ShapeFunctionData<element_type>& shape_function_data)
    {
      if (shape_function_data.myknots_.size() == 0)
        FOUR_C_THROW(
            "Got shape function data with knot size 0 - did you forget to initialize the element "
            "data with GetElementData::Get?");

      using type_dummy = CORE::LINALG::Matrix<element_type::n_nodes_, 1, typename V::scalar_type>;
      type_dummy N_dummy;

      if constexpr (element_type::element_dim_ == 2)
      {
        CORE::FE::NURBS::nurbs_get_2D_funct_deriv<typename V::scalar_type>(N_dummy, dN, xi,
            shape_function_data.myknots_, shape_function_data.weights_,
            element_type::discretization_);
      }
      else if constexpr (element_type::element_dim_ == 3)
      {
        CORE::FE::NURBS::nurbs_get_3D_funct_deriv(N_dummy, dN, xi, shape_function_data.myknots_,
            shape_function_data.weights_, t_nurbs27::discretization_);
      }
      else
      {
        FOUR_C_THROW("Got unexpected element dimension %d", element_type::element_dim_);
      }
    }
  };

  /**
   * \brief Get the shape function matrix of the element at xi.
   *
   * Multiplying this shape function matrix with the element DOF vector results in the field
   * variable evaluated at xi.
   *
   * @param xi (in) Parameter coordinate on the element.
   * @param N (out) shape function matrix.
   * @param shape_function_data (in) Shape function data container.
   */
  template <typename element_type, typename T, typename scalar_type_result,
      typename... shape_function_data_type>
  inline void EvaluateShapeFunctionMatrix(
      CORE::LINALG::Matrix<element_type::spatial_dim_, element_type::n_dof_, scalar_type_result>& N,
      const T& xi, const shape_function_data_type&... shape_function_data)
  {
    // Matrix for shape function values.
    CORE::LINALG::Matrix<1, element_type::n_nodes_ * element_type::n_val_, scalar_type_result>
        N_flat(true);

    // Evaluate the shape function values.
    EvaluateShapeFunction<element_type>::Evaluate(N_flat, xi, shape_function_data...);

    // Fill up the full shape function matrix.
    N.Clear();
    for (unsigned int node = 0; node < element_type::n_nodes_; node++)
      for (unsigned int dim = 0; dim < element_type::spatial_dim_; dim++)
        for (unsigned int val = 0; val < element_type::n_val_; val++)
          N(dim, element_type::spatial_dim_ * element_type::n_val_ * node +
                     element_type::spatial_dim_ * val + dim) =
              N_flat(element_type::n_val_ * node + val);
  }

}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
