/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the nullspace calculation at
       node level

\level 0

*/
/*---------------------------------------------------------------------*/
#ifndef LINALG_UTILS_NULLSPACE_H
#define LINALG_UTILS_NULLSPACE_H

#include <Teuchos_SerialDenseMatrix.hpp>

#include "discret.H"

namespace LINALG
{
  // forward declarations
  class Discretization;

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


  /*!
    \brief Helper function for the nodal nullspace of fluid elements

    \param node (in):    node to calculate the nullspace on
    \param numdof (in):  number of degrees of freedom
    \param dimnsp (in):  dimension of the nullspace
  */
  Teuchos::SerialDenseMatrix<int, double> ComputeFluidNullSpace(
      const DRT::Node& node, const int numdof, const int dimnsp);

}  // namespace LINALG

#endif
