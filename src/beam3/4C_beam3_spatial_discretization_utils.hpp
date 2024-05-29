/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief utility functions for spatial discretization of beam elements

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAM3_SPATIAL_DISCRETIZATION_UTILS_HPP
#define FOUR_C_BEAM3_SPATIAL_DISCRETIZATION_UTILS_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::UTILS::BEAM
{
  /** \brief evaluate shape functions at position \xi in element parameter space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionsAtXi(const T& xi, CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i,
      const CORE::FE::CellType& distype, double hermite_length_param = -1.0)
  {
    I_i.Clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape functions at xi
        CORE::FE::shape_function_1D(I_i, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape functions at xi: vpernode=2 means 3rd order, i.e. line2
        CORE::FE::shape_function_hermite_1D(
            I_i, xi, hermite_length_param, CORE::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate shape function derivatives at position \xi in element parameter space
   * [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionDerivsAtXi(const T& xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xi, const CORE::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xi.Clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        CORE::FE::shape_function_1D_deriv1(I_i_xi, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        CORE::FE::shape_function_hermite_1D_deriv1(
            I_i_xi, xi, hermite_length_param, CORE::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate second derivatives of shape function at position \xi in element parameter
   *         space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunction2ndDerivsAtXi(const T& xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xixi, const CORE::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xixi.Clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        CORE::FE::shape_function_1D_deriv2(I_i_xixi, xi, distype);
        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        CORE::FE::shape_function_hermite_1D_deriv2(
            I_i_xixi, xi, hermite_length_param, CORE::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate third derivatives of shape function at position \xi in element parameter
   *         space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunction3rdDerivsAtXi(const T& xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xixixi, const CORE::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    I_i_xixixi.Clear();

    switch (vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape function derivs at xi
        if (nnode >= 4)
          FOUR_C_THROW("Please implement 3rd derivatives of Lagrange shape functions!");

        break;
      }
      case 2:
      {
        FOUR_C_ASSERT(hermite_length_param != -1.0,
            "you must provide a length parameter in case of "
            "Hermite interpolation!");

        // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
        CORE::FE::shape_function_hermite_1D_deriv3(
            I_i_xixixi, xi, hermite_length_param, CORE::FE::CellType::line2);
        break;
      }
      default:
        FOUR_C_THROW("invalid value for vpernode (number of values per node) specified");
    }
  }

  /** \brief evaluate shape functions and its first derivatives at position \xi in element
   *         parameter space [-1,1]
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionsAndDerivsAtXi(const T& xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xi, const CORE::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    EvaluateShapeFunctionsAtXi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
    EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(xi, I_i_xi, distype, hermite_length_param);
  }

  /** \brief evaluate shape functions and its first and second derivatives at position \xi in
   *         element parameter space [-1,1]
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionsAndDerivsAnd2ndDerivsAtXi(const T& xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xi,
      CORE::LINALG::Matrix<1, vpernode * nnode, T>& I_i_xixi, const CORE::FE::CellType& distype,
      double hermite_length_param = -1.0)
  {
    EvaluateShapeFunctionsAtXi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
    EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(xi, I_i_xi, distype, hermite_length_param);
    EvaluateShapeFunction2ndDerivsAtXi<nnode, vpernode>(
        xi, I_i_xixi, distype, hermite_length_param);
  }

  /** \brief evaluate shape functions at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionsAllGPs(const CORE::FE::IntegrationPoints1D& gausspoints,
      std::vector<CORE::LINALG::Matrix<1, vpernode * nnode, T>>& I_i,
      const CORE::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    if (I_i.size() != (unsigned int)gausspoints.nquad)
      FOUR_C_THROW(
          "vector for individual shape functions to be evaluated at %d GPs has wrong size: %d",
          gausspoints.nquad, I_i.size());

    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const T xi_tilde = gausspoints.qxg[numgp][0];

      /* do a mapping into integration interval, i.e. coordinate transformation
       * note: this has no effect if integration interval is [-1;1] */
      const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                             (1.0 + xi_tilde) * integration_interval_upper_limit);

      EvaluateShapeFunctionsAtXi<nnode, vpernode>(xi, I_i[numgp], distype, hermite_length_param);
    }
  }

  /** \brief evaluate shape function derivatives at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionDerivsAllGPs(const CORE::FE::IntegrationPoints1D& gausspoints,
      std::vector<CORE::LINALG::Matrix<1, vpernode * nnode, T>>& I_i_xi,
      const CORE::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    if (I_i_xi.size() != (unsigned int)gausspoints.nquad)
      FOUR_C_THROW(
          "vector for individual shape function derivatives to be evaluated at %d GPs has "
          "wrong size: %d",
          gausspoints.nquad, I_i_xi.size());

    for (int numgp = 0; numgp < gausspoints.nquad; ++numgp)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const T xi_tilde = gausspoints.qxg[numgp][0];

      /* do a mapping into integration interval, i.e. coordinate transformation
       * note: this has no effect if integration interval is [-1;1] */
      const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                             (1.0 + xi_tilde) * integration_interval_upper_limit);

      EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(
          xi, I_i_xi[numgp], distype, hermite_length_param);
    }
  }

  /** \brief evaluate shape functions and derivatives at all specified Gauss points at once
   *
   */
  template <unsigned int nnode, unsigned int vpernode, typename T>
  void EvaluateShapeFunctionsAndDerivsAllGPs(const CORE::FE::IntegrationPoints1D& gausspoints,
      std::vector<CORE::LINALG::Matrix<1, vpernode * nnode, T>>& I_i,
      std::vector<CORE::LINALG::Matrix<1, vpernode * nnode, T>>& I_i_xi,
      const CORE::FE::CellType& distype, double hermite_length_param = -1.0,
      double integration_interval_lower_limit = -1.0, double integration_interval_upper_limit = 1.0)
  {
    EvaluateShapeFunctionsAllGPs<nnode, vpernode>(gausspoints, I_i, distype, hermite_length_param,
        integration_interval_lower_limit, integration_interval_upper_limit);
    EvaluateShapeFunctionDerivsAllGPs<nnode, vpernode>(gausspoints, I_i_xi, distype,
        hermite_length_param, integration_interval_lower_limit, integration_interval_upper_limit);
  }

  /** \brief assemble one shape function matrix, such that: r=N*d
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  void AssembleShapeFunctions(const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T>& N_i,
      CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, T>& N)
  {
    // assembly_N is just an array to help assemble the matrices of the shape functions
    // it determines, which shape function is used in which column of N
    std::array<std::array<unsigned int, 3 * numnodes * numnodalvalues>, 3> assembly_N{};

    /*
    Set number of shape functions for each 3*3 block:
    e.g. second order Lagrange polynomials (numnodes=3, numnodalvalues=1)
    int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
                             {0,1,0,0,2,0,0,3,0},
                             {0,0,1,0,0,2,0,0,3}};

    e.g. cubic Hermite polynomials (numnodes=2, numnodalvalues=2)
    int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
                             {0,1,0,0,2,0,0,3,0,0,4,0},
                             {0,0,1,0,0,2,0,0,3,0,0,4}};
    */

    for (unsigned int i = 0; i < numnodes * numnodalvalues; ++i)
    {
      assembly_N[0][3 * i] = i + 1;
      assembly_N[1][3 * i + 1] = i + 1;
      assembly_N[2][3 * i + 2] = i + 1;
    }

    // Assemble the matrices of the shape functions
    for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
      for (unsigned int j = 0; j < 3; ++j)
      {
        if (assembly_N[j][i] == 0)
          N(j, i) = 0.0;
        else
          N(j, i) = N_i(assembly_N[j][i] - 1);
      }
  }

  /** \brief assemble shape function matrices, such that: r=N*d, r_xi=N_xi*d, r_xixi=N_xixi*d
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  void AssembleShapeFunctionsAndDerivsAnd2ndDerivs(
      const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T>& N_i,
      const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T>& N_i_xi,
      const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, T>& N_i_xixi,
      CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, T>& N,
      CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_xi,
      CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_xixi)
  {
    AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i, N);
    AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i_xi, N_xi);
    AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i_xixi, N_xixi);
  }

  /** \brief interpolation of nodal DoFs based on given shape function values
   */
  template <unsigned int nnode, unsigned int vpernode, unsigned int ndim, typename T, typename T2>
  void CalcInterpolation(const CORE::LINALG::Matrix<ndim * vpernode * nnode, 1, T>& dof_vals,
      const CORE::LINALG::Matrix<1, vpernode * nnode, T2>& shapefcn_vals,
      CORE::LINALG::Matrix<ndim, 1, T>& result)
  {
    result.Clear();

    for (unsigned int idim = 0; idim < ndim; ++idim)
      for (unsigned int idofperdim = 0; idofperdim < vpernode * nnode; ++idofperdim)
        result(idim) += shapefcn_vals(idofperdim) * dof_vals(ndim * idofperdim + idim);
  }

  /** \brief evaluate length and derivative at all specified Gauss points at once
   */
  template <unsigned int nnode, unsigned int vpernode>
  std::tuple<double, double> IntegrateCenterlineArcLengthAndArcLengthDerivative(
      const CORE::FE::IntegrationPoints1D& gausspoints,
      const CORE::LINALG::Matrix<3 * vpernode * nnode, 1, double>& disp_centerline,
      const CORE::FE::CellType& distype, const double& reflength)
  {
    std::vector<CORE::LINALG::Matrix<1, nnode * vpernode, double>> H_i_xi(gausspoints.nquad);
    CORE::LINALG::Matrix<3, 1> r_xi;

    EvaluateShapeFunctionDerivsAllGPs<nnode, vpernode>(gausspoints, H_i_xi, distype, reflength);

    double int_length = 0.0, deriv_length = 0.0;

    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      double deriv_int = 0.0;

      CalcInterpolation<nnode, vpernode, 3, double>(disp_centerline, H_i_xi[numgp], r_xi);
      int_length += gausspoints.qwgt[numgp] * r_xi.Norm2();

      for (int dim = 0; dim < 3; dim++)
        deriv_int +=
            (disp_centerline(3 + dim) * H_i_xi[numgp](1) / reflength +
                disp_centerline(3 * vpernode * 1 + 3 + dim) * H_i_xi[numgp](3) / reflength) *
            r_xi(dim);

      deriv_length += gausspoints.qwgt[numgp] * deriv_int / r_xi.Norm2();
    }

    return {reflength - int_length, 1 - deriv_length};
  }
}  // namespace DRT::UTILS::BEAM

FOUR_C_NAMESPACE_CLOSE

#endif
