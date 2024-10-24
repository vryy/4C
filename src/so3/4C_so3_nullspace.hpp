// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_NULLSPACE_HPP
#define FOUR_C_SO3_NULLSPACE_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Nodes
{
  class Node;
}

namespace Discret::Elements
{
  /*!
    \brief Helper function for the nodal nullspace of solid elements in 3D

  \param node (in):    node to calculate the nullspace on
     \param x0 (in):      center of discretization
                      */
  Core::LinAlg::SerialDenseMatrix compute_solid_3d_null_space(
      const Core::Nodes::Node& node, const double* x0);

  /*!
   \brief Helper function for the nodal nullspace of solid elements in 2D

    \param node (in):    node to calculate the nullspace on
    \param x0 (in):      center of discretization
  */
  Core::LinAlg::SerialDenseMatrix compute_solid_2d_null_space(
      const Core::Nodes::Node& node, const double* x0);
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
