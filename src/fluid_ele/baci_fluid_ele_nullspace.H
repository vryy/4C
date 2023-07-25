/*! \file
\brief Helpers for fluid element nullspace computation
\level 0
*/

#ifndef BACI_FLUID_ELE_NULLSPACE_H
#define BACI_FLUID_ELE_NULLSPACE_H

#include <Teuchos_SerialDenseMatrix.hpp>

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
  Teuchos::SerialDenseMatrix<int, double> ComputeFluidNullSpace(
      const DRT::Node& node, const int numdof, const int dimnsp);
}  // namespace FLD

#endif
