/*----------------------------------------------------------------------*/
/*! \file

\brief collection of math tools for the interface computation of two
       curved meshes
       !WARNING: Except Tolerances not used at the moment
       (remove this comment and change level as soon as this functionality is tested again!)

\level 3

*----------------------------------------------------------------------*/


#ifndef INTERSECTION_MATH_H_
#define INTERSECTION_MATH_H_


#include "../linalg/linalg_fixedsizematrix.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseVector.h"


namespace GEO
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



  inline double XSIGN(const double a, const double b) { return (b >= 0.0 ? fabs(a) : -fabs(a)); }



  /*!
  \brief calculate a*a
  \param a
  \return a^2
   */
  template <typename T>
  inline T sqr(const T a)
  {
    return a * a;
  }



  /*!
  \brief  computes the Theorem of Pythagoras
          (a^2 + b^2)^(1/2)
          (modified from NUMERICAL RECIPES)

  \param a    (in)        : a
  \param b    (in)        : b
  \return result of the Theorem of Pythagoras
  */
  inline double pythagoras(const double a, const double b)
  {
    // note: sqrt(0.0) is properly defined and returns 0.0 as expected
    //      no check for zero needed
    return sqrt(sqr(a) + sqr(b));
  }


  /*!
  \brief sum up the absolute value of all entries of a vector
  \param v (in) : LINALG::Matrix<3,1>
  \return sum
  */
  inline double SumOfFabsEntries(const LINALG::Matrix<3, 1>& v)
  {
    return fabs(v(0)) + fabs(v(1)) + fabs(v(2));
  }



  /*!
  \brief computes the cross product of two LINALG::Matrix<3,1> a x b
  \param a (in) : arbitrary vector
  \param b (in) : arbitrary vector
  \return cross product of two LINALG::Matrix<3,1> a x b
  */
  LINALG::Matrix<3, 1> computeCrossProduct(
      const LINALG::Matrix<3, 1>& a, const LINALG::Matrix<3, 1>& b);


  /*!
  \brief computes the cross product of two LINALG::Matrix<3,1,T> a x b
  \param a (in) : arbitrary vector
  \param b (in) : arbitrary vector
  \return cross product of two LINALG::Matrix<3,1,T> a x b
  */
  template <typename T>
  LINALG::Matrix<3, 1, T> computeCrossProductT(
      const LINALG::Matrix<3, 1, T>& a, const LINALG::Matrix<3, 1, T>& b)
  {
    LINALG::Matrix<3, 1, T> c;

    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);

    return c;
  }



  /*-------------------------------------------------------------------*
  |  VON MIR GESCHRIEBEN                                               |
  *--------------------------------------------------------------------*/

  /*!
      \brief  Given a matrix A[1..m][1..n], this method computes its
              singular value decomposition, A =U ·W ·V T .
              The matrix U replaces a on output. The diagonal matrix of
              singular values W is output as a vector W[1..n].
              The matrix V (not the transpose V T ) is output as V[1..n][1..n].
              (modified from NUMERICAL RECIPES)

      \param A    (in/out)        : system matrix / matrix U
      \param W    (out)           : diagonal matrix stored in a vector
      \param V    (out)           : V matrix (not its transpose)
      */

  // Selbe Funktion wie svdcmp nur SerialDenseMatrizen als Eingabeparameter

  void svdcmpSerialDense(
      Epetra_SerialDenseMatrix& A, Epetra_SerialDenseMatrix& W, Epetra_SerialDenseMatrix& V);


  /*!
  \brief  Tests the singular value decomposition

  \param A    (in)        : system matrix
  \param U    (in)        : decomposed matrix U
  \param W    (in)        : diagonal matrix W stored in  vector
  \param V    (in)        : decomposed matrix V
  \param dim  (in)        : dimension
  */
  void test_svdcmp(LINALG::Matrix<3, 3>& A, LINALG::Matrix<3, 3>& U, LINALG::Matrix<3, 1>& W,
      LINALG::Matrix<3, 3>& V, int dim);


}  // namespace GEO


#endif  // INTERSECTION_MATH_H_
