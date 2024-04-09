/*! \file
\brief Helpers for fluid element nullspace computation
\level 0
*/

#ifndef FOUR_C_FLUID_ELE_NULLSPACE_HPP
#define FOUR_C_FLUID_ELE_NULLSPACE_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Node;
}

namespace FLD
{
  /*!
  \brief Helper function for the nodal nullspace of fluid elements

  \param node (in):    node to calculate the nullspace on
    \param numdof (in):  number of degrees of freedom
    \param dimnsp (in):  dimension of the nullspace
                         */
  CORE::LINALG::SerialDenseMatrix ComputeFluidNullSpace(
      const DRT::Node& node, const int numdof, const int dimnsp);
}  // namespace FLD

BACI_NAMESPACE_CLOSE

#endif
