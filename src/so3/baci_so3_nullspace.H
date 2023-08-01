/*! \file
\brief Helpers for solid element nullspace computation
\level 0
*/

#ifndef BACI_SO3_NULLSPACE_H
#define BACI_SO3_NULLSPACE_H

#include <Teuchos_SerialDenseMatrix.hpp>

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
  Teuchos::SerialDenseMatrix<int, double> ComputeSolid3DNullSpace(
      const DRT::Node& node, const double* x0);

  /*!
   \brief Helper function for the nodal nullspace of solid elements in 2D

    \param node (in):    node to calculate the nullspace on
    \param x0 (in):      center of discretization
  */
  Teuchos::SerialDenseMatrix<int, double> ComputeSolid2DNullSpace(
      const DRT::Node& node, const double* x0);
}  // namespace DRT::ELEMENTS

#endif
