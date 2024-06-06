/*! \file
\brief Helpers for solid element nullspace computation
\level 0
*/

#ifndef FOUR_C_SO3_NULLSPACE_HPP
#define FOUR_C_SO3_NULLSPACE_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Nodes
{
  class Node;
}

namespace Discret::ELEMENTS
{
  /*!
    \brief Helper function for the nodal nullspace of solid elements in 3D

  \param node (in):    node to calculate the nullspace on
     \param x0 (in):      center of discretization
                      */
  Core::LinAlg::SerialDenseMatrix ComputeSolid3DNullSpace(
      const Core::Nodes::Node& node, const double* x0);

  /*!
   \brief Helper function for the nodal nullspace of solid elements in 2D

    \param node (in):    node to calculate the nullspace on
    \param x0 (in):      center of discretization
  */
  Core::LinAlg::SerialDenseMatrix ComputeSolid2DNullSpace(
      const Core::Nodes::Node& node, const double* x0);
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
