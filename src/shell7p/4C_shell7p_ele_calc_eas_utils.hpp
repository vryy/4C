/*! \file

\brief Helper functions for shell element EAS technology

\level 1
*/

#ifndef FOUR_C_SHELL7P_ELE_CALC_EAS_UTILS_HPP
#define FOUR_C_SHELL7P_ELE_CALC_EAS_UTILS_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_shell7p_ele_calc_lib.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS::SHELL::EAS
{
  /*!
   * @brief Evaluates the transformation matrix T0^{-1} which maps the M-matrix from centroid point
   * (parameter space) to the gaussian point (material configuration)
   * @tparam distype
   * @param akov (in) : Jacobian mapping evaluated at the gaussina point
   * @param akov0 (in) : Jacobian mapping evaluated at the element centroid
   *
   * @return CORE::LINALG::SerialDenseMatrix : Transformation matrix T0inv
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::SerialDenseMatrix EvaluateT0inv(
      const CORE::LINALG::Matrix<DRT::ELEMENTS::SHELL::DETAIL::num_dim,
          DRT::ELEMENTS::SHELL::DETAIL::num_dim>& akov,
      const CORE::LINALG::Matrix<DRT::ELEMENTS::SHELL::DETAIL::num_dim,
          DRT::ELEMENTS::SHELL::DETAIL::num_dim>& akon0)
  {
    CORE::LINALG::SerialDenseMatrix T0inv(DRT::ELEMENTS::SHELL::DETAIL::num_internal_variables,
        DRT::ELEMENTS::SHELL::DETAIL::num_internal_variables);
    double t11, t12, t13, t21, t22, t23, t31, t32, t33;
    // components of the transformation matrix T0^{-1}
    t11 = 0.0;
    t12 = 0.0;
    t13 = 0.0;
    t21 = 0.0;
    t22 = 0.0;
    t23 = 0.0;
    t31 = 0.0;
    t32 = 0.0;
    t33 = 1.0;
    for (int i = 0; i < DRT::ELEMENTS::SHELL::DETAIL::num_dim; ++i)
    {
      t11 += akov(0, i) * akon0(0, i);
      t12 += akov(0, i) * akon0(1, i);
      t21 += akov(1, i) * akon0(0, i);
      t22 += akov(1, i) * akon0(1, i);
    }

    T0inv(0, 0) = t11 * t11;
    T0inv(1, 0) = 2 * t11 * t21;
    T0inv(2, 0) = 2 * t11 * t31;
    T0inv(3, 0) = t21 * t21;
    T0inv(4, 0) = 2 * t21 * t31;
    T0inv(5, 0) = t31 * t31;

    T0inv(0, 1) = t11 * t12;
    T0inv(1, 1) = (t11 * t22 + t21 * t12);
    T0inv(2, 1) = (t11 * t32 + t31 * t12);
    T0inv(3, 1) = t21 * t22;
    T0inv(4, 1) = (t21 * t32 + t31 * t22);
    T0inv(5, 1) = t31 * t32;

    T0inv(0, 2) = t11 * t13;
    T0inv(1, 2) = (t11 * t23 + t21 * t13);
    T0inv(2, 2) = (t11 * t33 + t31 * t13);
    T0inv(3, 2) = t21 * t23;
    T0inv(4, 2) = (t21 * t33 + t31 * t23);
    T0inv(5, 2) = t31 * t33;

    T0inv(0, 3) = t12 * t12;
    T0inv(1, 3) = 2 * t12 * t22;
    T0inv(2, 3) = 2 * t12 * t32;
    T0inv(3, 3) = t22 * t22;
    T0inv(4, 3) = 2 * t22 * t32;
    T0inv(5, 3) = t32 * t32;

    T0inv(0, 4) = t12 * t13;
    T0inv(1, 4) = (t12 * t23 + t22 * t13);
    T0inv(2, 4) = (t12 * t33 + t32 * t13);
    T0inv(3, 4) = t22 * t23;
    T0inv(4, 4) = (t22 * t33 + t32 * t23);
    T0inv(5, 4) = t32 * t33;

    T0inv(0, 5) = t13 * t13;
    T0inv(1, 5) = 2 * t13 * t23;
    T0inv(2, 5) = 2 * t13 * t33;
    T0inv(3, 5) = t23 * t23;
    T0inv(4, 5) = 2 * t23 * t33;
    T0inv(5, 5) = t33 * t33;

    T0inv(6, 6) = t11 * t11;
    T0inv(7, 6) = 2 * t11 * t21;
    T0inv(8, 6) = 2 * t11 * t31;
    T0inv(9, 6) = t21 * t21;
    T0inv(10, 6) = 2 * t21 * t31;
    T0inv(11, 6) = t31 * t31;

    T0inv(6, 7) = t11 * t12;
    T0inv(7, 7) = (t11 * t22 + t21 * t12);
    T0inv(8, 7) = (t11 * t32 + t31 * t12);
    T0inv(9, 7) = t21 * t22;
    T0inv(10, 7) = (t21 * t32 + t31 * t22);
    T0inv(11, 7) = t31 * t32;

    T0inv(6, 8) = t11 * t13;
    T0inv(7, 8) = (t11 * t23 + t21 * t13);
    T0inv(8, 8) = (t11 * t33 + t31 * t13);
    T0inv(9, 8) = t21 * t23;
    T0inv(10, 8) = (t21 * t33 + t31 * t23);
    T0inv(11, 8) = t31 * t33;

    T0inv(6, 9) = t12 * t12;
    T0inv(7, 9) = 2 * t12 * t22;
    T0inv(8, 9) = 2 * t12 * t32;
    T0inv(9, 9) = t22 * t22;
    T0inv(10, 9) = 2 * t22 * t32;
    T0inv(11, 9) = t32 * t32;

    T0inv(6, 10) = t12 * t13;
    T0inv(7, 10) = (t12 * t23 + t22 * t13);
    T0inv(8, 10) = (t12 * t33 + t32 * t13);
    T0inv(9, 10) = t22 * t23;
    T0inv(10, 10) = (t22 * t33 + t32 * t23);
    T0inv(11, 10) = t32 * t33;

    T0inv(6, 11) = t13 * t13;
    T0inv(7, 11) = 2 * t13 * t23;
    T0inv(8, 11) = 2 * t13 * t33;
    T0inv(9, 11) = t23 * t23;
    T0inv(10, 11) = 2 * t23 * t33;
    T0inv(11, 11) = t33 * t33;

    for (int i = 0; i < DRT::ELEMENTS::SHELL::DETAIL::node_dof; ++i)
    {
      for (int j = 0; j < DRT::ELEMENTS::SHELL::DETAIL::node_dof; ++j)
      {
        T0inv(i, j + DRT::ELEMENTS::SHELL::DETAIL::node_dof) = 0.0;
        T0inv(i + DRT::ELEMENTS::SHELL::DETAIL::node_dof, j) = 0.0;
      }
    }
    return T0inv;
  }


  /*!
   * @brief Compute the matrix M that is the element-wise matrix of the shape functions for the
   * enhanced strains in the parameter space
   *
   * @tparam distype
   * @param xi_gp (in) : Coordinate in the parameter space
   * @param eas (in) : Number of eas parameters for each locking part
   * @param neas (in) : Total number of eas parameters
   * @return CORE::LINALG::SerialDenseMatrix : Enhanced strains shape function matrix in parameter
   * space M
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::SerialDenseMatrix EvaluateEASShapeFunctionsParameterSpace(
      const std::array<double, 2>& xi_gp, const STR::ELEMENTS::ShellLockingTypes& locking_types)
  {
    // evaluation of the shape function matrix to interpolate the enhanced strains alpha
    CORE::LINALG::SerialDenseMatrix M(
        DRT::ELEMENTS::SHELL::DETAIL::num_internal_variables, locking_types.total);

    const double xi = xi_gp[0];
    const double eta = xi_gp[1];
    int M_index = 0;

    const double xieta = xi * eta;
    const double xixi = xi * xi;
    const double etaeta = eta * eta;
    const double xixieta = xi * xieta;
    const double xietaeta = xieta * eta;
    const double xixietaeta = xieta * xieta;

    const int num_node = DRT::ELEMENTS::SHELL::DETAIL::num_node<distype>;
    if (num_node > 4)
    {
      // membrane locking: E_{11}, E_{12}, E_{22} constant
      switch (locking_types.membrane)
      {
        case 0:
          break;
        case 7:
          M(0, M_index) = eta - 3.0 * xixieta;
          M(0, M_index + 1) = etaeta - 3.0 * xixietaeta;
          M(3, M_index + 2) = xi - 3.0 * xietaeta;
          M(3, M_index + 3) = xixi - 3.0 * xixietaeta;
          M(1, M_index + 4) = eta - 3.0 * xixieta;
          M(1, M_index + 5) = xi - 3.0 * xietaeta;
          M(1, M_index + 6) = 1.0 - 3.0 * (xixi + etaeta) + 9.0 * xixietaeta;
          M_index += 7;
          break;
        case 9:
          M(0, M_index) = 1.0 - 3.0 * xixi;
          M(0, M_index + 1) = eta - 3.0 * xixieta;
          M(3, M_index + 2) = 1.0 - 3.0 * etaeta;
          M(3, M_index + 3) = xi - 3.0 * xietaeta;
          M(1, M_index + 4) = 1.0 - 3.0 * xixieta;
          M(1, M_index + 5) = 1.0 - 3.0 * xietaeta;
          M(1, M_index + 6) = eta - 3.0 * xixieta;
          M(1, M_index + 7) = xi - 3.0 * xietaeta;
          M(1, M_index + 8) = 1.0 - 3.0 * (xixi + etaeta) + 9.0 * xixietaeta;
          M_index += 9;
          break;
        case 11:
          M(0, M_index) = 1.0 - 3.0 * xixi;
          M(0, M_index + 1) = eta - 3.0 * xixieta;
          M(0, M_index + 2) = etaeta - 3.0 * xixietaeta;
          M(3, M_index + 3) = 1.0 - 3.0 * etaeta;
          M(3, M_index + 4) = xi - 3.0 * xietaeta;
          M(3, M_index + 5) = xixi - 3.0 * xixietaeta;
          M(1, M_index + 6) = 1.0 - 3.0 * xixi;
          M(1, M_index + 7) = 1.0 - 3.0 * etaeta;
          M(1, M_index + 8) = eta - 3.0 * xixieta;
          M(1, M_index + 9) = xi - 3.0 * xietaeta;
          M(1, M_index + 10) = 1.0 - 3.0 * (xixi + etaeta) + 9.0 * xixietaeta;
          M_index += 11;
          break;
        default:
          FOUR_C_THROW(
              "EAS Membrane locking: Only 0, 7, 9, 11 EAS modes are implemented. Given: "
              "%d",
              locking_types.membrane);
      }
      // bending locking: E_{11}, E_{12], E_{22} linear
      switch (locking_types.bending)
      {
        case 0:
          break;
        case 9:
          M(6, M_index) = 1.0 - 3.0 * xixi;
          M(6, M_index + 1) = eta - 3.0 * xixieta;
          M(9, M_index + 2) = 1.0 - 3.0 * etaeta;
          M(9, M_index + 3) = xi - 3.0 * xietaeta;
          M(7, M_index + 4) = 1.0 - 3.0 * xixi;
          M(7, M_index + 5) = 1.0 - 3.0 * etaeta;
          M(7, M_index + 6) = eta - 3.0 * xixieta;
          M(7, M_index + 7) = xi - 3.0 * xietaeta;
          M(7, M_index + 8) = 1.0 - 3.0 * (xixi + etaeta) + 9.0 * xixietaeta;
          M_index += 9;
          break;
        case 11:
          M(6, M_index) = 1.0 - 3.0 * xixi;
          M(6, M_index + 1) = eta - 3.0 * xixieta;
          M(6, M_index + 2) = etaeta - 3.0 * xixietaeta;
          M(9, M_index + 3) = 1.0 - 3.0 * etaeta;
          M(9, M_index + 4) = xi - 3.0 * xietaeta;
          M(9, M_index + 5) = xixi - 3.0 * xixietaeta;
          M(7, M_index + 6) = 1.0 - 3.0 * xixi;
          M(7, M_index + 7) = 1.0 - 3.0 * etaeta;
          M(7, M_index + 8) = eta - 3.0 * xixieta;
          M(7, M_index + 9) = xi - 3.0 * xietaeta;
          M(7, M_index + 10) = 1.0 - 3.0 * (xixi + etaeta) + 9.0 * xixietaeta;
          M_index += 11;
          break;
        default:
          FOUR_C_THROW("EAS bending part: Only 0, 9, 11 EAS modes are implemented. Given: %d",
              locking_types.bending);
      }
      // locking due to thickness changes E_{33} linear
      switch (locking_types.thickness)
      {
        case 0:
          break;
        case 1:
          M(11, M_index) = 1.0;
          M_index += 1;
          break;
        case 3:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M_index += 3;
          break;
        case 4:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M_index += 4;
          break;
        case 6:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = xixi;
          M(11, M_index + 5) = etaeta;
          M_index += 6;
          break;
        case 8:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = xixi;
          M(11, M_index + 5) = etaeta;
          M(11, M_index + 6) = xixieta;
          M(11, M_index + 7) = xietaeta;
          M_index += 8;
          break;
        case 9:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = 1.0 - 3.0 * xixi;
          M(11, M_index + 5) = 1.0 - 3.0 * etaeta;
          M(11, M_index + 6) = xixieta;
          M(11, M_index + 7) = xietaeta;
          M(11, M_index + 8) = 1.0 - 9.0 * xixietaeta;
          M_index += 9;
          break;
        default:
          FOUR_C_THROW(
              "EAS thickness locking: Only 0, 1, 3, 4, 5, 8, 9 EAS modes are implemented. Given: "
              "%d",
              locking_types.thickness);
      }
      // transverse shear strain locking: E_{13}, E_{23} constant
      switch (locking_types.transverse_shear_strain_const)
      {
        case 0:
          break;
        case 2:
          M(2, M_index) = eta - 3.0 * xixieta;
          M(4, M_index + 1) = xi - 3.0 * xietaeta;
          M_index += 2;
          break;
        case 4:
          M(2, M_index) = 1.0 - 3.0 * xixi;
          M(2, M_index + 1) = eta - 3.0 * xixieta;
          M(4, M_index + 2) = 1.0 - 3.0 * etaeta;
          M(4, M_index + 3) = xi - 3.0 * xietaeta;
          M_index += 4;
          break;
        case 6:
          M(2, M_index) = 1.0 - 3.0 * xixi;
          M(2, M_index + 1) = eta - 3.0 * xixieta;
          M(2, M_index + 2) = etaeta - 3.0 * xixietaeta;
          M(4, M_index + 3) = 1.0 - 3.0 * etaeta;
          M(4, M_index + 4) = xi - 3.0 * xietaeta;
          M(4, M_index + 5) = xixi - 3.0 * xixietaeta;
          M_index += 6;
          break;
        default:
          FOUR_C_THROW(
              "EAS transverse shear strain locking: Only 0, 2, 4, 6 EAS modes are implemented. "
              "Given: %d",
              locking_types.transverse_shear_strain_const);
      }
      // transverse shear strain locking: E_{13}, E_{23} linear
      switch (locking_types.transverse_shear_strain_lin)
      {
        case 0:
          break;
        case 2:
          M(8, M_index) = xixi;
          M(10, M_index + 1) = etaeta;
          M_index += 2;
          break;
        case 4:
          M(8, M_index) = xixi;
          M(8, M_index + 1) = xixietaeta;
          M(10, M_index + 2) = etaeta;
          M(10, M_index + 3) = xixietaeta;
          M_index += 4;
          break;
        case 6:
          M(8, M_index) = xixi;
          M(8, M_index + 1) = xixieta;
          M(8, M_index + 2) = xixietaeta;
          M(10, M_index + 3) = etaeta;
          M(10, M_index + 4) = xietaeta;
          M(10, M_index + 5) = xixietaeta;
          M_index += 6;
          break;
        default:
          FOUR_C_THROW(
              "EAS transverse shear strain locking: Only 0, 2, 4, 6 EAS modes are implemented. "
              "Given: %d",
              locking_types.transverse_shear_strain_lin);
      }
    }
    // four node element
    else if (num_node == 4)
    {
      // membrane locking: E_{11}, E_{22} constant
      switch (locking_types.membrane)
      {
        case 0:
          break;
        case 1:
          M(3, M_index) = eta;
          M_index += 1;
          break;
        case 2:
          M(1, M_index) = xi;
          M(1, M_index + 1) = eta;
          M_index += 2;
          break;
        case 3:
          M(1, M_index) = xi;
          M(1, M_index + 1) = eta;
          M(1, M_index + 2) = xieta;
          M_index += 3;
          break;
        case 4:
          M(0, M_index) = xi;
          M(3, M_index + 1) = eta;
          M(1, M_index + 2) = xi;
          M(1, M_index + 3) = eta;
          M_index += 4;
          break;
        case 5:
          M(0, M_index) = xi;
          M(3, M_index + 1) = eta;
          M(1, M_index + 2) = xi;
          M(1, M_index + 3) = eta;
          M(1, M_index + 4) = xieta;
          M_index += 5;
          break;
        case 7:
          M(0, M_index) = xi;
          M(3, M_index + 1) = eta;
          M(1, M_index + 2) = xi;
          M(1, M_index + 3) = eta;
          M(0, M_index + 4) = xieta;
          M(3, M_index + 5) = xieta;
          M(1, M_index + 6) = xieta;
          M_index += 7;
          break;
        default:
          FOUR_C_THROW(
              "EAS Membrane locking: Only 0, 1, 2, 3, 4, 5, 7 EAS modes are implemented. Given: %d",
              locking_types.membrane);
      }
      // bending part: E_{11}, E_{12], E_{22} linear
      switch (locking_types.bending)
      {
        case 0:
          break;
        case 4:
          M(6, M_index) = xi;
          M(9, M_index + 1) = eta;
          M(7, M_index + 2) = xi;
          M(7, M_index + 3) = eta;
          M_index += 4;
          break;
        case 5:
          M(6, M_index) = xi;
          M(9, M_index + 1) = eta;
          M(7, M_index + 2) = xi;
          M(7, M_index + 3) = eta;
          M(7, M_index + 4) = xieta;
          M_index += 5;
          break;
        case 6:
          M(6, M_index) = xixi;
          M(6, M_index + 1) = xixietaeta;
          M(9, M_index + 2) = etaeta;
          M(9, M_index + 3) = xixietaeta;
          M(7, M_index + 4) = xixi;
          M(7, M_index + 5) = etaeta;
          M_index += 6;
          break;
        case 7:
          M(6, M_index) = xi;
          M(9, M_index + 1) = eta;
          M(7, M_index + 2) = xi;
          M(7, M_index + 3) = eta;
          M(6, M_index + 4) = xieta;
          M(9, M_index + 5) = xieta;
          M(7, M_index + 6) = xieta;
          M_index += 7;
          break;
        default:
          FOUR_C_THROW("EAS bending part: Only 0, 4, 5, 7, 8 EAS modes are implemented. Given: %d",
              locking_types.bending);
      }
      // locking due to thickness changes E_{33} linear
      switch (locking_types.thickness)
      {
        case 0:
          break;
        case 1:
          M(11, M_index) = 1.0;
          M_index += 1;
          break;
        case 3:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M_index += 3;
          break;
        case 4:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M_index += 4;
          break;
        case 6:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = xixi;
          M(11, M_index + 5) = etaeta;
          M_index += 6;
          break;
        case 8:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = xixi;
          M(11, M_index + 5) = etaeta;
          M(11, M_index + 6) = xixieta;
          M(11, M_index + 7) = xietaeta;
          M_index += 8;
          break;
        case 9:
          M(11, M_index) = 1.0;
          M(11, M_index + 1) = xi;
          M(11, M_index + 2) = eta;
          M(11, M_index + 3) = xieta;
          M(11, M_index + 4) = 1.0 - 3.0 * xixi;
          M(11, M_index + 5) = 1.0 - 3.0 * etaeta;
          M(11, M_index + 6) = xixieta;
          M(11, M_index + 7) = xietaeta;
          M(11, M_index + 8) = 1.0 - 9.0 * xixietaeta;
          M_index += 9;
          break;
        default:
          FOUR_C_THROW(
              "EAS thickness locking: Only 0, 3, 4, 6, 8, 9 EAS modes are implemented. Given: %d",
              locking_types.thickness);
      }
      // transverse shear strain locking: E_{13}, E_{23} const
      switch (locking_types.transverse_shear_strain_const)
      {
        case 0:
          break;
        case 2:
          M(2, M_index) = xi;
          M(4, M_index + 1) = eta;
          M_index += 2;
          break;
        case 4:
          M(2, M_index) = xi;
          M(2, M_index + 1) = xieta;
          M(4, M_index + 2) = eta;
          M(4, M_index + 3) = xieta;
          M_index += 4;
          break;
        default:
          FOUR_C_THROW(
              "EAS transverse shear strain locking: Only 0, 2, 4 EAS modes are implemented. Given: "
              "%d",
              locking_types.transverse_shear_strain_const);
      }
      // transverse shear strain locking: E_{13}, E_{23} linear
      switch (locking_types.transverse_shear_strain_lin)
      {
        case 0:
          break;
        case 2:
          M(8, M_index) = xi;
          M(10, M_index + 1) = eta;
          M_index += 2;
          break;
        case 4:
          M(8, M_index) = xi;
          M(8, M_index + 1) = xieta;
          M(10, M_index + 2) = eta;
          M(10, M_index + 3) = xieta;
          M_index += 4;
          break;
        default:
          FOUR_C_THROW(
              "EAS transverse shear strain locking: Only 0, 2, 4 EAS modes are implemented. Given: "
              "%d",
              locking_types.transverse_shear_strain_lin);
      }
    }  // else if (iel==4)
    else
      FOUR_C_THROW("EAS only implemented for 4, 8 and 9 node elements");

    FOUR_C_ASSERT(M_index == locking_types.total,
        "Wrong total number of EAS parameters. Something went wrong.");

    return M;
  }

  /*!
   * @brief Update the current metric tensor due to EAS
   *
   * @tparam distype : discretization type
   * @param g (in/out) : An object holding the metrics and basis vectors
   * @param strain (in) : Strains
   * @param zeta (in) : Thickness coordinate of gaussian point (scaled via SDC)
   */
  template <CORE::FE::CellType distype>
  void UpdateCurrentMetricsEAS(DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& g,
      const CORE::LINALG::SerialDenseVector& strain, const double& zeta)
  {
    // update kovariant metric tensor
    // g_11 += 2 * (alpha_11 + zeta * beta_11)
    g.metric_kovariant_(0, 0) = g.metric_kovariant_(0, 0) + 2.0 * (strain(0) + zeta * strain(6));
    // g_12 += alpha_12 + zeta * beta_12
    g.metric_kovariant_(1, 0) = g.metric_kovariant_(1, 0) + strain(1) + zeta * strain(7);
    // g_13 += alpha_13 + zeta * beta_13
    g.metric_kovariant_(2, 0) = g.metric_kovariant_(2, 0) + strain(2) + zeta * strain(8);
    // g_22 += 2 * (alpha_22 + zeta * beta_22)
    g.metric_kovariant_(1, 1) = g.metric_kovariant_(1, 1) + 2.0 * (strain(3) + zeta * strain(9));
    // g_23 += alpha_23 + zeta * beta_23
    g.metric_kovariant_(2, 1) = g.metric_kovariant_(2, 1) + strain(4) + zeta * strain(10);
    // g_33 += 2 * (alpha_33 + zeta * beta_33)
    g.metric_kovariant_(2, 2) = g.metric_kovariant_(2, 2) + 2.0 * (strain(5) + zeta * strain(11));
    // g_21 = g_12
    g.metric_kovariant_(0, 1) = g.metric_kovariant_(1, 0);
    // g_31 = g_13
    g.metric_kovariant_(0, 2) = g.metric_kovariant_(2, 0);
    // g_32 = g_23
    g.metric_kovariant_(1, 2) = g.metric_kovariant_(2, 1);

    // re-evaluate kontravariant metric tensor
    g.metric_kontravariant_.Update(g.metric_kovariant_);

    double detJ = g.metric_kontravariant_.Invert();
    g.detJ_ = std::sqrt(detJ);
  }


  /*!
   * @brief Map the EAS shapefunction matrix M in the parameter space to M_gp in the global
   * configuration and return M_gp
   *
   *  M_gp = detJ0/detJ T0^{-T} M = J0^T M J0^{-T}
   *
   * @tparam distype : discretization type
   * @param neas (in) : Number of total EAS parameters
   * @param a_reference (in) : An object holding the current basis vectors and
   * metric tensors of the shell body
   * @param metrics_centroid_reference (in) : An object holding the current basis vectors and
   * metric tensors evaluated at the element center
   * @return CORE::LINALG::SerialDenseMatrix : EAS shapefunction matrix M
   */
  template <CORE::FE::CellType distype>
  CORE::LINALG::SerialDenseMatrix MapEasShapeFunctionsToGaussPoint(const int& neas,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& metrics_centroid_reference,
      CORE::LINALG::SerialDenseMatrix& M)
  {
    // get transformation matrix T0^-1 to map M from midpoint to gaussian point
    CORE::LINALG::SerialDenseMatrix T0inv =
        EvaluateT0inv<distype>(a_reference.kovariant_, metrics_centroid_reference.kontravariant_);
    // transform basis of M-matrix to gaussian point: M_gp = detJ0/ detJ * T0^-T * M
    CORE::LINALG::SerialDenseMatrix M_gp(
        DRT::ELEMENTS::SHELL::DETAIL::num_internal_variables, neas);
    M_gp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
        metrics_centroid_reference.detJ_ / a_reference.detJ_, T0inv, M, 0.0);
    return M_gp;
  }


  // Evaluates M_gp by setting up M and mapping M to M_gp via T0^{-T}
  template <CORE::FE::CellType distype>
  CORE::LINALG::SerialDenseMatrix EvaluateEasShapeFunctions(const std::array<double, 2>& xi_gp,
      const STR::ELEMENTS::ShellLockingTypes& locking_types,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& a_reference,
      const DRT::ELEMENTS::SHELL::BasisVectorsAndMetrics<distype>& metrics_centroid_reference)
  {
    CORE::LINALG::SerialDenseMatrix M =
        EAS::EvaluateEASShapeFunctionsParameterSpace<distype>(xi_gp, locking_types);
    CORE::LINALG::SerialDenseMatrix M_gp = MapEasShapeFunctionsToGaussPoint(
        locking_types.total, a_reference, metrics_centroid_reference, M);
    return M_gp;
  }


  // Evaluate enhanced EAS strains
  void EvaluateEasStrains(CORE::LINALG::SerialDenseVector& strain_enh,
      const CORE::LINALG::SerialDenseMatrix& alpha, CORE::LINALG::SerialDenseMatrix& M)
  {
    // evaluate enhanced strains = M * alpha to "unlock" element
    strain_enh.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, M, alpha, 0.0);
  }

  /*!
   * @brief Add EAS internal force contribution of one Gauss point
   *
   * @param LinvDTilde (in) : L*DTilde^{-1}
   * @param RTilde (in) : RTilde
   * @param force (in/out) : internal force vector where the local contribution is added to
   */
  void AddEASInternalForce(const CORE::LINALG::SerialDenseMatrix& LinvDTilde,
      const CORE::LINALG::SerialDenseMatrix& RTilde, CORE::LINALG::SerialDenseVector& force_vector)
  {
    // EAS internal force vector is: R - L * DTilde^-1 * RTilde (R is the usual unenhanced
    // internal force vector
    force_vector.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1., LinvDTilde, RTilde, 1.0);
  }

  /*!
   * @brief Add EAS stiffness matrix contribution of one Gauss point
   *
   * @param LinvDTilde(in) : L*DTilde^{-1}
   * @param transL (in) : L^T
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added
   * to
   */
  void AddEASStiffnessMatrix(const CORE::LINALG::SerialDenseMatrix& LinvDTilde,
      const CORE::LINALG::SerialDenseMatrix& transL,
      CORE::LINALG::SerialDenseMatrix& stiffness_matrix)
  {
    // EAS-stiffness matrix is: K- L * DTilde^-1 * L^T (K is the usual unenhanced stiffness
    // matrix)
    stiffness_matrix.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1., LinvDTilde, transL, 1.0);
  }
}  // namespace DRT::ELEMENTS::SHELL::EAS

FOUR_C_NAMESPACE_CLOSE

#endif
