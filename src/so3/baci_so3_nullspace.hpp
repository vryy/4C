/*! \file
\brief Helpers for solid element nullspace computation
\level 0
*/

#ifndef FOUR_C_SO3_NULLSPACE_HPP
#define FOUR_C_SO3_NULLSPACE_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Node;
}

namespace DRT::ELEMENTS
{
  /*!
    \brief Helper function for the nodal nullspace of solid elements in 3D

  \param node (in):    node to calculate the nullspace on
     \param x0 (in):      center of discretization
                      */
  CORE::LINALG::SerialDenseMatrix ComputeSolid3DNullSpace(const DRT::Node& node, const double* x0);

  /*!
   \brief Helper function for the nodal nullspace of solid elements in 2D

    \param node (in):    node to calculate the nullspace on
    \param x0 (in):      center of discretization
  */
  CORE::LINALG::SerialDenseMatrix ComputeSolid2DNullSpace(const DRT::Node& node, const double* x0);
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
