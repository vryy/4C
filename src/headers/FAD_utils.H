/*----------------------------------------------------------------------*/
/*! \file

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
  /**
   * \brief Cast of FAD to double.
   *
   * This function will be used for FAD types, since the double version is specialized (see next
   * function). We pass the value of the FAD variable to this function again. By doing so this
   * function can also handle nested FAD types, i.e. higher order derivatives.
   */
  template <typename scalar_type>
  inline double CastToDouble(scalar_type a)
  {
    return CastToDouble(a.val());
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
  LINALG::Matrix<dim1, dim2, double> CastToDouble(LINALG::Matrix<dim1, dim2, type> a)
  {
    LINALG::Matrix<dim1, dim2, double> b(true);

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
  \brief Calculate square root of a scalar FAD quantity. If a compiler Error is thrown when calling
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
  \brief Calculate Norm of a scalar FAD quantity. If a compiler Error is thrown when calling this
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
  scalar_type VectorNorm(LINALG::Matrix<length, 1, scalar_type> v)
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
  double VectorNorm(LINALG::Matrix<length, 1, double> v)
  {
    return v.Norm2();
  }

  //! Calculates the Norm of a FAD vector, since .Norm2() is not available for FAD vectors
  // Todo this function is obsolete
  template <typename T>
  T Norm(LINALG::Matrix<3, 1, T> v)
  {
    T norm_squared = 0.0;
    for (int i = 0; i < 3; i++)
    {
      norm_squared += v(i) * v(i);
    }

    return FADUTILS::sqrt(norm_squared);
  }

  /*!
  \brief Calculate inner product of two FAD or double vectors
  */
  // Todo this function is obsolete, use Dot of LINALG::Matrix instead
  template <typename type>
  type ScalarProduct(LINALG::Matrix<3, 1, type> a, LINALG::Matrix<3, 1, type> b)
  {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  /*!
  \brief Calculate difference of two FAD or double vectors
  */
  // Todo this function is obsolete, use Update of LINALG::Matrix instead
  template <typename type>
  LINALG::Matrix<3, 1, type> DiffVector(LINALG::Matrix<3, 1, type> a, LINALG::Matrix<3, 1, type> b)
  {
    LINALG::Matrix<3, 1, type> c(true);
    for (int i = 0; i < 3; i++) c(i) = a(i) - b(i);

    return c;
  }

  /*!
  \brief Calculate vector product of two FAD or double vectors
  */
  // Todo this function is obsolete, use CrossProduct of LINALG::Matrix instead
  template <typename type>
  LINALG::Matrix<3, 1, type> VectorProduct(
      LINALG::Matrix<3, 1, type> first_vector, LINALG::Matrix<3, 1, type> second_vector)
  {
    LINALG::Matrix<3, 1, type> result_vector;
    result_vector.Clear();
    LINALG::Matrix<3, 3, type> S_first_vector;
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


  /**
   * \brief Build nested Fad type for computing derivatives up to order N.
   *
   * This function is a slightly modified version of MakeFad from
   * trilinos-source/packages/sacado/example/high_order_example.cpp
   *
   * @tparam N Number of nested derivatives, i.e. order of derivatives to compute.
   * @tparam BaseFadType Basic Fad type of this nested type.
   */
  template <int N, typename BaseFadType>
  struct HigherOrderFadType
  {
    //! Nested type of this Fad type.
    typedef typename HigherOrderFadType<N - 1, BaseFadType>::type nested_type;

    //! Fad type of this object.
    typedef typename Sacado::mpl::apply<BaseFadType, nested_type>::type type;
  };

  /**
   * \brief Specialization for the last nested type, i.e. first order derivatives.
   */
  template <typename BaseFadType>
  struct HigherOrderFadType<1, BaseFadType>
  {
    typedef BaseFadType type;
  };


  /**
   * \brief Create a variable of a nested FAD type.
   *
   * This function is a slightly modified version of MakeFad from
   * trilinos-source/packages/sacado/example/high_order_example.cpp
   *
   * @tparam FadType Fad type of the variable to create.
   */
  template <typename FadType>
  struct HigherOrderFadValue
  {
    //! Type of nested value.
    typedef typename FadType::value_type nested_type;

    /**
     * \brief Set a value of this Fad type.
     * @param n (in) Number of dependent variables.
     * @param i (in) Index of the current dependent variable.
     * @param x (in) Value of the current dependent variable.
     * @return Created Fad type variable.
     */
    static FadType apply(const int n, const int i, const double x)
    {
      return FadType(n, i, HigherOrderFadValue<nested_type>::apply(n, i, x));
    }
  };

  /**
   * \brief Specialization for the last nested type, i.e. double.
   */
  template <>
  struct HigherOrderFadValue<double>
  {
    static double apply(const int n, const int i, const double x) { return x; }
  };

}  // namespace FADUTILS

#endif
