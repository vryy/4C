/*----------------------------------------------------------------------*/
/*! \file

\brief Functions to evaluate the templated elements

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_ELEMENT_EVALUATION_FUNCTIONS_HPP
#define FOUR_C_GEOMETRY_PAIR_ELEMENT_EVALUATION_FUNCTIONS_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_shape_functions.hpp"
#include "4C_geometry_pair_utility_functions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  /**
   * \brief Evaluate the field function in the element
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePosition(const T& xi, const ElementData<element_type, V>& element_data,
      Core::LinAlg::Matrix<element_type::spatial_dim_, 1, scalar_type>& r)
  {
    // Matrix for shape function values
    Core::LinAlg::Matrix<1, element_type::n_nodes_ * element_type::n_val_, scalar_type> N(true);

    // Evaluate the shape function values
    EvaluateShapeFunction<element_type>::evaluate(N, xi, element_data.shape_function_data_);

    // Calculate the field function
    r.clear();
    for (unsigned int node = 0; node < element_type::n_nodes_; node++)
      for (unsigned int dim = 0; dim < element_type::spatial_dim_; dim++)
        for (unsigned int val = 0; val < element_type::n_val_; val++)
          r(dim) += element_data.element_position_(
                        element_type::spatial_dim_ * element_type::n_val_ * node +
                        element_type::spatial_dim_ * val + dim) *
                    N(element_type::n_val_ * node + val);
  }

  /**
   * \brief Evaluate the derivative of the field function w.r.t xi in the element
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePositionDerivative1(const T& xi,
      const ElementData<element_type, V>& element_data,
      Core::LinAlg::Matrix<element_type::spatial_dim_, element_type::element_dim_, scalar_type>& dr)
  {
    // Matrix for shape function values
    Core::LinAlg::Matrix<element_type::element_dim_, element_type::n_nodes_ * element_type::n_val_,
        scalar_type>
        dN(true);

    // Evaluate the shape function values
    EvaluateShapeFunction<element_type>::EvaluateDeriv1(dN, xi, element_data.shape_function_data_);

    // Calculate the derivative of the field function
    dr.clear();
    for (unsigned int dim = 0; dim < element_type::spatial_dim_; dim++)
      for (unsigned int direction = 0; direction < element_type::element_dim_; direction++)
        for (unsigned int node = 0; node < element_type::n_nodes_; node++)
          for (unsigned int val = 0; val < element_type::n_val_; val++)
            dr(dim, direction) += element_data.element_position_(
                                      element_type::spatial_dim_ * element_type::n_val_ * node +
                                      element_type::spatial_dim_ * val + dim) *
                                  dN(direction, element_type::n_val_ * node + val);
  }

  /**
   * \brief Evaluate the geometrical normal (the one resulting from the FE approximation) on the
   * face element
   */
  template <typename surface, typename T, typename scalar_type_dof, typename scalar_type_result>
  void EvaluateFaceNormal(const T& xi,
      const ElementData<surface, scalar_type_dof>& element_data_surface,
      Core::LinAlg::Matrix<3, 1, scalar_type_result>& normal)
  {
    // Check at compile time if a surface (2D) element is given
    static_assert(
        surface::element_dim_ == 2, "EvaluateFaceNormal can only be called for 2D elements!");

    Core::LinAlg::Matrix<3, 2, scalar_type_result> dr;
    Core::LinAlg::Matrix<3, 1, scalar_type_result> dr_0;
    Core::LinAlg::Matrix<3, 1, scalar_type_result> dr_1;
    GEOMETRYPAIR::EvaluatePositionDerivative1<surface>(xi, element_data_surface, dr);
    for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
    {
      dr_0(i_dir) = dr(i_dir, 0);
      dr_1(i_dir) = dr(i_dir, 1);
    }
    normal.CrossProduct(dr_0, dr_1);
    normal.Scale(1.0 / Core::FADUtils::VectorNorm(normal));
    if constexpr (std::is_same<surface, t_nurbs9>::value)
    {
      // In the NURBS case some normals have to be flipped to point outward of the volume element
      normal.Scale(element_data_surface.shape_function_data_.surface_normal_factor_);
    }
  }

  /**
   * \brief Evaluate the normal at a position on the surface, depending on the surface type, either
   * the averaged nodal normal field or the kinematic normal field
   */
  template <typename surface, typename T, typename scalar_type_dof, typename scalar_type_result>
  void EvaluateSurfaceNormal(const T& xi,
      const ElementData<surface, scalar_type_dof>& element_data_surface,
      Core::LinAlg::Matrix<3, 1, scalar_type_result>& normal)
  {
    // Check at compile time if a surface (2D) element is given
    static_assert(
        surface::element_dim_ == 2, "EvaluateSurfaceNormal can only be called for 2D elements!");

    if constexpr (IsSurfaceAveragedNormalsElement<surface>::value_)
    {
      static_assert(surface::n_val_ == 1,
          "EvaluateSurfaceNormal with averaged nodal normals only makes sense for elements with "
          "one value per node and dimension");

      // Calculate the normal as a interpolation of nodal normals
      Core::LinAlg::Matrix<1, surface::n_nodes_, typename T::scalar_type> N(true);
      EvaluateShapeFunction<surface>::evaluate(N, xi, element_data_surface.shape_function_data_);
      normal.clear();
      for (unsigned int node = 0; node < surface::n_nodes_; node++)
        for (unsigned int dim = 0; dim < surface::spatial_dim_; dim++)
          normal(dim) +=
              element_data_surface.nodal_normals_(surface::spatial_dim_ * node + dim) * N(node);
      normal.Scale(1.0 / Core::FADUtils::VectorNorm(normal));
    }
    else
    {
      // No averaged normals, so evaluate the kinematic normals on the face element
      EvaluateFaceNormal(xi, element_data_surface, normal);
    }
  }

  /**
   * \brief Evaluate a position on the surface, where the third parameter coordinate is taken along
   * the surface normal
   */
  template <typename surface, typename scalar_type_xi, typename scalar_type_dof,
      typename scalar_type_result>
  void EvaluateSurfacePosition(const Core::LinAlg::Matrix<3, 1, scalar_type_xi>& xi,
      const ElementData<surface, scalar_type_dof>& element_data_surface,
      Core::LinAlg::Matrix<3, 1, scalar_type_result>& r)
  {
    // Check at compile time if a surface (2D) element is given
    static_assert(
        surface::element_dim_ == 2, "EvaluateSurfacePosition can only be called for 2D elements!");

    // Evaluate the normal
    Core::LinAlg::Matrix<3, 1, scalar_type_result> normal;
    EvaluateSurfaceNormal<surface>(xi, element_data_surface, normal);

    // Evaluate the position on the surface
    GEOMETRYPAIR::EvaluatePosition<surface>(xi, element_data_surface, r);

    // Add the normal part to the position
    normal.Scale(xi(2));
    r += normal;
  }

  /**
   * \brief Evaluate a triad on a plane curve
   *
   * The curve has to lie in the x-y plane. The first basis vector of the triad is the tangent along
   * the curve, the second basis vector is the tangent rotated clockwise around the z axis by pi/2.
   * The third basis vector is the Cartesian e_z basis vector.
   */
  template <typename line, typename scalar_type_xi, typename scalar_type_dof,
      typename scalar_type_result>
  void EvaluateTriadAtPlaneCurve(const scalar_type_xi xi,
      const ElementData<line, scalar_type_dof>& element_data_line,
      Core::LinAlg::Matrix<3, 3, scalar_type_result>& triad)
  {
    Core::LinAlg::Matrix<3, 1, scalar_type_result> tangent, cross_section_director_2,
        cross_section_director_3;
    EvaluatePositionDerivative1<line>(xi, element_data_line, tangent);

    if (std::abs(tangent(2)) > Constants::pos_tol)
      FOUR_C_THROW(
          "EvaluateTriadAtPlaneCurve: The tangent vector can not have a component in z direction! "
          "The component is %f!",
          Core::FADUtils::CastToDouble(tangent(2)));

    // Create the director vectors in the cross-section
    // Director 2 is the one in the y-axis (reference configuration)
    // Director 3 is the one in the z-axis (reference configuration)
    tangent.Scale(1. / Core::FADUtils::VectorNorm(tangent));
    cross_section_director_2.clear();
    cross_section_director_2(0) = -tangent(1);
    cross_section_director_2(1) = tangent(0);
    cross_section_director_3.clear();
    cross_section_director_3(2) = 1.;

    // Set the triad
    for (unsigned int dir = 0; dir < 3; dir++)
    {
      triad(dir, 0) = tangent(dir);
      triad(dir, 1) = cross_section_director_2(dir);
      triad(dir, 2) = cross_section_director_3(dir);
    }
  }

  /**
   * \brief Evaluate the Jacobi matrix for a volume element
   */
  template <typename volume, typename scalar_type>
  void EvaluateJacobian(const Core::LinAlg::Matrix<3, 1, scalar_type>& xi,
      const GEOMETRYPAIR::ElementData<volume, double>& X_volume,
      Core::LinAlg::Matrix<3, 3, scalar_type>& J)
  {
    // Check at compile time if a volume (3D) element is given
    static_assert(
        volume::element_dim_ == 3, "EvaluateJacobian can only be called for 3D elements!");

    // Get the derivatives of the reference position w.r.t the parameter coordinates. This is the
    // transposed Jacobi matrix.
    Core::LinAlg::Matrix<3, 3, scalar_type> dXdxi(true);
    EvaluatePositionDerivative1<volume>(xi, X_volume, dXdxi);
    J.clear();
    J.UpdateT(dXdxi);
  }

  /**
   * \brief Evaluate the 3D deformation gradient F in a volume element
   */
  template <typename volume, typename scalar_type_xi, typename scalar_type_dof,
      typename scalar_type_result>
  void EvaluateDeformationGradient(const Core::LinAlg::Matrix<3, 1, scalar_type_xi>& xi,
      const GEOMETRYPAIR::ElementData<volume, double>& X_volume,
      const GEOMETRYPAIR::ElementData<volume, scalar_type_dof>& q_volume,
      Core::LinAlg::Matrix<3, 3, scalar_type_result>& F)
  {
    // Check at compile time if a volume (3D) element is given
    static_assert(volume::element_dim_ == 3,
        "EvaluateDeformationGradient can only be called for 3D elements!");

    // Get the inverse of the Jacobian
    Core::LinAlg::Matrix<3, 3, scalar_type_xi> inv_J(true);
    EvaluateJacobian<volume>(xi, X_volume, inv_J);
    Core::LinAlg::Inverse(inv_J);

    // Get the derivatives of the shape functions w.r.t to the parameter coordinates
    Core::LinAlg::Matrix<volume::element_dim_, volume::n_nodes_ * volume::n_val_, scalar_type_xi>
        dNdxi(true);
    GEOMETRYPAIR::EvaluateShapeFunction<volume>::EvaluateDeriv1(
        dNdxi, xi, q_volume.shape_function_data_);

    // Transform to derivatives w.r.t physical coordinates
    Core::LinAlg::Matrix<volume::element_dim_, volume::n_nodes_ * volume::n_val_, scalar_type_xi>
        dNdX(true);
    for (unsigned int i_row = 0; i_row < 3; i_row++)
      for (unsigned int i_col = 0; i_col < volume::n_nodes_ * volume::n_val_; i_col++)
        for (unsigned int i_sum = 0; i_sum < 3; i_sum++)
          dNdX(i_row, i_col) += inv_J(i_row, i_sum) * dNdxi(i_sum, i_col);

    // Calculate F
    F.clear();
    for (unsigned int i_row = 0; i_row < 3; i_row++)
      for (unsigned int i_col = 0; i_col < 3; i_col++)
        for (unsigned int i_sum = 0; i_sum < volume::n_nodes_ * volume::n_val_; i_sum++)
          F(i_col, i_row) += dNdX(i_row, i_sum) * q_volume.element_position_(3 * i_sum + i_col);
  }

  /**
   * \brief Evaluate a position on the surface and its derivative w.r.t. xi.
   * @param element_data_surface (in) Degrees of freedom for the surface.
   * @param xi (in) Parameter coordinates on the surface (the first two are in the surface
   * parameter coordinates, the third one is in the normal direction).
   * @param r (out) Position on the surface.
   * @param dr (out) Derivative of the position on the surface, w.r.t xi.
   */
  template <typename scalar_type, typename surface>
  void EvaluateSurfacePositionAndDerivative(
      const ElementData<surface, scalar_type>& element_data_surface,
      const Core::LinAlg::Matrix<3, 1, scalar_type>& xi, Core::LinAlg::Matrix<3, 1, scalar_type>& r,
      Core::LinAlg::Matrix<3, 3, scalar_type>& dr)
  {
    // Create a nested FAD type.
    using FAD_outer = Sacado::ELRFad::SLFad<scalar_type, 3>;

    // Set up the AD variables.
    Core::LinAlg::Matrix<3, 1, FAD_outer> xi_AD;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++) xi_AD(i_dim) = FAD_outer(3, i_dim, xi(i_dim));
    Core::LinAlg::Matrix<3, 1, FAD_outer> r_AD;

    // Evaluate the position.
    EvaluateSurfacePosition<surface>(xi_AD, element_data_surface, r_AD);

    // Extract the return values from the AD types.
    for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
    {
      r(i_dir) = r_AD(i_dir).val();
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++) dr(i_dir, i_dim) = r_AD(i_dir).dx(i_dim);
    }
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 1D element
   */
  template <typename T>
  bool ValidParameter1D(const T& xi)
  {
    const double xi_limit = 1.0 + Constants::projection_xi_eta_tol;
    if (fabs(xi) < xi_limit)
      return true;
    else
      return false;
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 2D element
   */
  template <typename element_type, typename T>
  bool ValidParameter2D(const T& xi)
  {
    const double xi_limit = 1.0 + Constants::projection_xi_eta_tol;
    if (element_type::geometry_type_ == DiscretizationTypeGeometry::quad)
    {
      if (fabs(xi(0)) < xi_limit && fabs(xi(1)) < xi_limit) return true;
    }
    else if (element_type::geometry_type_ == DiscretizationTypeGeometry::triangle)
    {
      if (xi(0) > -Constants::projection_xi_eta_tol && xi(1) > -Constants::projection_xi_eta_tol &&
          xi(0) + xi(1) < 1.0 + Constants::projection_xi_eta_tol)
        return true;
    }
    else
    {
      FOUR_C_THROW("Wrong DiscretizationTypeGeometry given!");
    }

    // Default value
    return false;
  }

  /**
   * \brief Check if the parameter coordinate xi is in the valid range for a 3D element
   */
  template <typename element_type, typename T>
  bool ValidParameter3D(const T& xi)
  {
    const double xi_limit = 1.0 + Constants::projection_xi_eta_tol;
    if (element_type::geometry_type_ == DiscretizationTypeGeometry::hexahedron)
    {
      if (fabs(xi(0)) < xi_limit && fabs(xi(1)) < xi_limit && fabs(xi(2)) < xi_limit) return true;
    }
    else if (element_type::geometry_type_ == DiscretizationTypeGeometry::tetraeder)
    {
      if (xi(0) > -Constants::projection_xi_eta_tol && xi(1) > -Constants::projection_xi_eta_tol &&
          xi(2) > -Constants::projection_xi_eta_tol &&
          xi(0) + xi(1) + xi(2) < 1.0 + Constants::projection_xi_eta_tol)
        return true;
    }
    else
    {
      FOUR_C_THROW("Wrong DiscretizationTypeGeometry given!");
    }

    // Default value
    return false;
  }

  /**
   * \brief This struct sets parameter coordinates start values for local Newton iterations on an
   * element
   */
  template <DiscretizationTypeGeometry geometry_type>
  struct StartValues
  {
   public:
    /**
     * \brief Set the start values for a local Newton iterations on the element.
     *
     * This method will be specialized for each geometry type. By doing so the values can be set at
     * compile time.
     */
    template <typename T>
    static void set(T& xi)
    {
      std::string error_string =
          "This is the default implementation of StartValues::Set and should never be called, "
          "since we only want to call the templated versions. You are calling it with the "
          "DiscretizationTypeGeometry ";
      error_string += DiscretizationTypeGeometryToString(geometry_type);
      FOUR_C_THROW(error_string.c_str());
    }
  };

  /**
   * \brief Template specialization for line elements
   */
  template <>
  template <typename T>
  void StartValues<DiscretizationTypeGeometry::line>::set(T& xi)
  {
    xi = 0.0;
  }

  /**
   * \brief Template specialization for triangle surface elements
   */
  template <>
  template <typename T>
  void StartValues<DiscretizationTypeGeometry::triangle>::set(T& xi)
  {
    // We do not use xi.PutScalar(0.25) here, since this might be a surface element which has a
    // normal direction and we do not want to set an initial value in the normal direction
    xi.PutScalar(0.0);
    for (int i_dim = 0; i_dim < 2; i_dim++) xi(i_dim) = 0.25;
  }

  /**
   * \brief Template specialization for quad surface elements
   */
  template <>
  template <typename T>
  void StartValues<DiscretizationTypeGeometry::quad>::set(T& xi)
  {
    xi.PutScalar(0.0);
  }

  /**
   * \brief Template specialization for tetrahedra volume elements
   */
  template <>
  template <typename T>
  void StartValues<DiscretizationTypeGeometry::tetraeder>::set(T& xi)
  {
    xi.PutScalar(0.25);
  }

  /**
   * \brief Template specialization for hexahedron volume elements
   */
  template <>
  template <typename T>
  void StartValues<DiscretizationTypeGeometry::hexahedron>::set(T& xi)
  {
    xi.PutScalar(0.0);
  }

}  // namespace GEOMETRYPAIR


FOUR_C_NAMESPACE_CLOSE

#endif
