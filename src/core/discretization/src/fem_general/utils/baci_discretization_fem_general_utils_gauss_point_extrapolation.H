/*----------------------------------------------------------------------*/
/*! \file

\brief Provide functionality for extrapolation of any quantity from Gauss point to nodes via the
shape functions

\level 2

*----------------------------------------------------------------------*/
#ifndef BACI_DISCRETIZATION_FEM_GENERAL_UTILS_GAUSS_POINT_EXTRAPOLATION_H
#define BACI_DISCRETIZATION_FEM_GENERAL_UTILS_GAUSS_POINT_EXTRAPOLATION_H

#include "baci_lib_element.H"

namespace CORE::DRT::UTILS
{
  /*!
   * @brief Evaluates the linear operator that extrapolates quantities from the Gauss points to
   * the nodes.
   *
   * The evaluation happens only once per discretization type and integration rule. It will be
   * calculated only during the first execution and will then be stored internally for a faster
   * access. The procedure is as follows:
   *
   * (1) Get the base discretization type based on the given discretization type and the
   * integration rule. It is the discretization type of the same shape with the same or less
   * number of nodes compared to the number of Gauss points. Example: The base discretization type
   * of a HEX8-element that is integrated with 8 Gauss points is HEX8. The base discretization
   * type of a TET10-element with 4 Gauss points is TET4.
   *
   * (2) Evaluates the shape functions at the gauss points which is a matrix M=[N^i(xi_j)]
   *
   * (3) Either invert the matrix M (if it is a square matrix) or solve a least square algorithm
   *
   * (4) Extend the matrix such that the quantities on the additional nodes of the discretization
   * type compared to the base discretization type are interpolated using the shape functions of
   * the base discretization type.
   *
   * @note This method is not thread-safe and must not be called in a multi-threaded environment.
   *
   * @tparam distype Discretization type of the element
   * @param integration (in) : Gauss integration points
   * @return CORE::LINALG::SerialDenseMatrix
   */
  template <::DRT::Element::DiscretizationType distype, class GaussIntegration>
  CORE::LINALG::SerialDenseMatrix EvaluateGaussPointsToNodesExtrapolationMatrix(
      const GaussIntegration& integration);

  /*!
   * @brief Extrapolates Gauss point data to the nodes and assembles it into a global data vector
   *
   * @tparam distype discretization type
   * @tparam GaussIntegration type of the container that holds the gauss integration rule
   * @param ele (in) :  element
   * @param gp_data (in) : Data at the Gauss points
   * @param global_data (out) : global data
   * @param nodal_average (in) : Whether to compute nodal averages at the nodes which are shared
   * between elements. Only set this to true if every element contributes equally to the nodal
   * quantity.
   * @param integration (in) : Container that holds the integration points
   */
  template <::DRT::Element::DiscretizationType distype, class GaussIntegration>
  void ExtrapolateGPQuantityToNodesAndAssemble(const ::DRT::Element& ele,
      const CORE::LINALG::SerialDenseMatrix& gp_data, Epetra_MultiVector& global_data,
      bool nodal_average, const GaussIntegration& integration);
}  // namespace CORE::DRT::UTILS

#endif  // FEM_GENERAL_UTILS_GAUSS_POINT_EXTRAPOLATION_H
