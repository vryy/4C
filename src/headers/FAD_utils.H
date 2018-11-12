/*----------------------------------------------------------------------*/
/*!
\file FAD_utils.H

\brief A set of utility functions for Forward Automatic Differentiation

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/


#ifndef FAD_UTILS_H
#define FAD_UTILS_H

#include "../linalg/linalg_fixedsizematrix.H"
#include <Sacado.hpp>

namespace FADUTILS
{
  /*!
  \brief Cast of FAD to double
  */
  template <typename scalar_type>
  inline double CastToDouble(scalar_type a)
  {
    return a.val();
  }

  /*!
  \brief Template specialization for double
  */
  template <>
  inline double CastToDouble(double a)
  {
    return a;
  }

  /*!
  \brief Cast of a FAD matrix to a double matrix
  */
  template <typename type, unsigned int dim1, unsigned int dim2>
  LINALG::TMatrix<double, dim1, dim2> CastToDouble(LINALG::TMatrix<type, dim1, dim2> a)
  {
    LINALG::TMatrix<double, dim1, dim2> b(true);

    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim2; j++)
      {
        b(i, j) = CastToDouble(a(i, j));
      }
    }
    return b;
  }

  /*!
  \brief Calculate signum function of FAD or double quantity
  */
  template <typename type>
  double Signum(type a)
  {
    if (a >= 0.0)
      return 1.0;
    else
      return -1.0;
  }

  /*!
  \brief Calculate square root of a scalar FAD quantity. If a compiler Error is thrown then calling
  this function, check if the template argument is explicitly stated in the function call, i.e.
  FADUTILS::sqrt<my_AD_type>(...)
  */
  template <typename scalar_type>
  inline scalar_type sqrt(scalar_type a)
  {
    // avoid non-differentiable point of square root function by using conditional
    /* Todo should we use a tolerance/threshold here? */
    if (a == 0.0)
      return 0.0;
    else
      return std::sqrt(a);
  }

  /*!
  \brief Template specialization for double
  */
  template <>
  inline double sqrt(double a)
  {
    return std::sqrt(a);
  }

  /*!
  \brief Calculate Norm of a scalar FAD quantity. If a compiler Error is thrown then calling this
  function, check if the template argument is explicitly stated in the function call, i.e.
  FADUTILS::Norm<my_AD_type>(...)
  */
  template <typename scalar_type>
  inline scalar_type Norm(scalar_type a)
  {
    return FADUTILS::sqrt<scalar_type>(a * a);
  }

  /*!
  \brief Calculate Norm of a FAD vector
  */
  template <typename scalar_type, unsigned int length>
  scalar_type VectorNorm(LINALG::TMatrix<scalar_type, length, 1> v)
  {
    scalar_type norm_squared = 0.0;
    for (unsigned int i = 0; i < length; i++)
    {
      norm_squared += v(i) * v(i);
    }

    return FADUTILS::sqrt<scalar_type>(norm_squared);
  }

  /*!
  \brief Template specialization for double
  */
  template <unsigned int length>
  double VectorNorm(LINALG::TMatrix<double, length, 1> v)
  {
    return v.Norm2();
  }

  //! Calculates the Norm of a FAD vector, since .Norm2() is not available for FAD vectors
  // Todo this function is obsolete
  template <typename T>
  T Norm(LINALG::TMatrix<T, 3, 1> v)
  {
    T norm_squared = 0.0;
    for (int i = 0; i < 3; i++)
    {
      norm_squared += v(i) * v(i);
    }

    return std::sqrt(norm_squared);
  }

  /*!
  \brief Calculate inner product of two FAD or double vectors
  */
  // Todo this function is obsolete, use Dot of LINALG::TMatrix instead
  template <typename type>
  type ScalarProduct(LINALG::TMatrix<type, 3, 1> a, LINALG::TMatrix<type, 3, 1> b)
  {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  /*!
  \brief Calculate difference of two FAD or double vectors
  */
  // Todo this function is obsolete, use Update of LINALG::TMatrix instead
  template <typename type>
  LINALG::TMatrix<type, 3, 1> DiffVector(
      LINALG::TMatrix<type, 3, 1> a, LINALG::TMatrix<type, 3, 1> b)
  {
    LINALG::TMatrix<type, 3, 1> c(true);
    for (int i = 0; i < 3; i++) c(i) = a(i) - b(i);

    return c;
  }

  /*!
  \brief Calculate vector product of two FAD or double vectors
  */
  // Todo this function is obsolete, use CrossProduct of LINALG::TMatrix instead
  template <typename type>
  LINALG::TMatrix<type, 3, 1> VectorProduct(
      LINALG::TMatrix<type, 3, 1> first_vector, LINALG::TMatrix<type, 3, 1> second_vector)
  {
    LINALG::TMatrix<type, 3, 1> result_vector;
    result_vector.Clear();
    LINALG::TMatrix<type, 3, 3> S_first_vector;
    S_first_vector.Clear();

    S_first_vector(0, 0) = 0.0;
    S_first_vector(0, 1) = -first_vector(2);
    S_first_vector(0, 2) = first_vector(1);
    S_first_vector(1, 0) = first_vector(2);
    S_first_vector(1, 1) = 0.0;
    S_first_vector(1, 2) = -first_vector(0);
    S_first_vector(2, 0) = -first_vector(1);
    S_first_vector(2, 1) = first_vector(0);
    S_first_vector(2, 2) = 0.0;

    result_vector.Multiply(S_first_vector, second_vector);

    return result_vector;
  }

}  // namespace FADUTILS

#endif
