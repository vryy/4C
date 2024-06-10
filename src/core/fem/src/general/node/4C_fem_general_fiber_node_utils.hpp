/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of utility functions for fiber interpolation

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_FIBER_NODE_UTILS_HPP
#define FOUR_C_FEM_GENERAL_FIBER_NODE_UTILS_HPP
#include "4C_config.hpp"

#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Nodes
{
  class NodalFiberHolder;


  /*!
   * \brief Projects the fibers and angles to the Gauss points with the shape functions
   *
   * Fibers are projected with the shape functions and then normalized to unit length. A
   * RAD-AXI-CIR coordinate system is orthogonalized with preserving the CIR direction over the
   * AXI direction.
   *
   * \tparam distype distype type of the element used
   * \param nodes list of nodes of the element
   * \param shapefcts list of shape functions evaluated at the gauss points
   * \param gpFiberHolder output of the projected fibers and angles
   */
  template <Core::FE::CellType distype>
  void ProjectFibersToGaussPoints(const Core::Nodes::Node* const* nodes,
      const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
      NodalFiberHolder& gpFiberHolder);

  /*!
   * @brief Projects an std::array<double, n> from the nodes with the shape functions into the
   * inner of the element
   *
   * @tparam distype type of the element used
   * @tparam dim dimension of the elemenet
   * @param quantity array of arrays to be projected with the shape functions. The first index
   * are the points, the second index the dimension of the quantity
   * @param shapefcts Shape functions
   * @param quantityProjected Output quantity projected with the shape functions
   */
  template <Core::FE::CellType distype, std::size_t dim>
  void ProjectQuantityWithShapeFunctions(
      const std::array<std::array<double, dim>, Core::FE::num_nodes<distype>> quantity,
      const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
      std::vector<Core::LinAlg::Matrix<dim, 1>>& quantityProjected);

  /*!
   * @brief Projectsd a double from the nodes with the shapefunctions into the innder of the
   * element
   *
   * @tparam distype type of the element used
   * @param quantity array of doubles to be projected with the shape functions.
   * @param shapefcts Shape functions
   * @param quantityProjected Output quantity projected with the shape functions
   */
  template <Core::FE::CellType distype>
  void ProjectQuantityWithShapeFunctions(
      const std::array<double, Core::FE::num_nodes<distype>> quantity,
      const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
      std::vector<double>& quantityProjected);

  /*!
   * \brief Check whether all nodes of the element have fibers.
   *
   * \tparam distype discretization type
   * \param nodes Pointer to the nodes of the element
   * \return true All nodes have fibers
   * \return false At least one one does not have a fiber
   */
  template <Core::FE::CellType distype>
  bool HaveNodalFibers(const Core::Nodes::Node* const* nodes);
}  // namespace Core::Nodes

FOUR_C_NAMESPACE_CLOSE

#endif
