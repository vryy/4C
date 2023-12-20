/*! \file
\brief Helpers for solid element nullspace computation
\level 0
*/

#ifndef BACI_SO3_NULLSPACE_H
#define BACI_SO3_NULLSPACE_H

#include "baci_config.H"

#include "baci_linalg_serialdensematrix.H"

BACI_NAMESPACE_OPEN

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

BACI_NAMESPACE_CLOSE

#endif
