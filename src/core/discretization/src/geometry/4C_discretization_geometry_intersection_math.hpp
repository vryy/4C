/*----------------------------------------------------------------------*/
/*! \file

\brief collection of math tools for the interface computation of two
       curved meshes
       !WARNING: Except Tolerances not used at the moment
       (remove this comment and change level as soon as this functionality is tested again!)

\level 3

*----------------------------------------------------------------------*/


#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_INTERSECTION_MATH_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_INTERSECTION_MATH_HPP


#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN


namespace Core::Geo
{
  //! named tolerance for easy search/grep
  const double TOL14 = 1e-14;
  //! named tolerance for easy search/grep
  const double TOL13 = 1e-13;
  //! named tolerance for easy search/grep
  const double TOL12 = 1e-12;
  //! named tolerance for easy search/grep
  const double TOL7 = 1e-7;
  //! named tolerance for easy search/grep
  const double TOL6 = 1e-6;
  //! named tolerance for easy search/grep
  const double TOL2 = 1e-2;
  //! named tolerance for easy search/grep
  const double TOLPLUS8 = 1e8;
  //! large number to start computations of nearest distance in tree
  const double LARGENUMBER = 1e30;
}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
