/*----------------------------------------------------------------------*/
/*! \file

\brief collection of math tools for the interface computation of two
       curved meshes
       !WARNING: Except Tolerances not used at the moment
       (remove this comment and change level as soon as this functionality is tested again!)

\level 3

*----------------------------------------------------------------------*/


#ifndef BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_MATH_H
#define BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_MATH_H


#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensevector.H"

#include <Teuchos_SerialDenseSolver.hpp>


namespace CORE::GEO
{
  //! tolerance used for residuals in Newton-methods
  const double TOL14 = 1e-14;

  //! tolerance used for residuals in Newton-methods
  const double TOL13 = 1e-13;
  const double TOL12 = 1e-12;

  //! tolerance used for residuals in Newton-methods
  const double TOL10 = 1e-10;

  //! named tolerance for easy search/grep
  const double TOL7 = 1e-7;
  //! named tolerance for easy search/grep
  const double TOL6 = 1e-6;
  //! named tolerance for easy search/grep
  const double TOL5 = 1e-5;
  //! named tolerance for easy search/grep
  const double TOL4 = 1e-4;
  //! named tolerance for easy search/grep
  const double TOL3 = 1e-3;
  //! named tolerance for easy search/grep
  const double TOL2 = 1e-2;
  //! named tolerance for easy search/grep
  const double TOLPLUS8 = 1e8;
  //! named tolerance for easy search/grep
  const double TOLPLUS15 = 1e15;
  //! large number to start computations of nearest distance in tree
  const double LARGENUMBER = 1e30;


  /*!
  \brief sum up the absolute value of all entries of a vector
  \param v (in) : CORE::LINALG::Matrix<3,1>
  \return sum
  */
  inline double SumOfFabsEntries(const CORE::LINALG::Matrix<3, 1>& v)
  {
    return fabs(v(0)) + fabs(v(1)) + fabs(v(2));
  }


  /*!
  \brief computes the cross product of two CORE::LINALG::Matrix<3,1> a x b
  \param a (in) : arbitrary vector
  \param b (in) : arbitrary vector
  \return cross product of two CORE::LINALG::Matrix<3,1> a x b
  */
  inline CORE::LINALG::Matrix<3, 1> computeCrossProduct(
      const CORE::LINALG::Matrix<3, 1>& a, const CORE::LINALG::Matrix<3, 1>& b)
  {
    CORE::LINALG::Matrix<3, 1> c;

    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);

    return c;
  }

  /*!
  \brief computes the cross product of two CORE::LINALG::Matrix<3,1,T> a x b
  \param a (in) : arbitrary vector
  \param b (in) : arbitrary vector
  \return cross product of two CORE::LINALG::Matrix<3,1,T> a x b
  */
  template <typename T>
  CORE::LINALG::Matrix<3, 1, T> computeCrossProductT(
      const CORE::LINALG::Matrix<3, 1, T>& a, const CORE::LINALG::Matrix<3, 1, T>& b)
  {
    CORE::LINALG::Matrix<3, 1, T> c;

    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);

    return c;
  }
}  // namespace CORE::GEO


#endif  // DISCRETIZATION_GEOMETRY_INTERSECTION_MATH_H
