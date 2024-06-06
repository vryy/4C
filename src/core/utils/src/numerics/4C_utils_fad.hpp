/*----------------------------------------------------------------------*/
/*! \file

\brief A set of utility functions for Forward Automatic Differentiation

\level 3


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_UTILS_FAD_HPP
#define FOUR_C_UTILS_FAD_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FADUtils
{
  namespace Details
  {
    template <typename T, typename AlwaysVoid = void>
    constexpr bool is_double_convertible = false;

    template <typename T>
    constexpr bool is_double_convertible<T,
        std::void_t<decltype(static_cast<double>(std::declval<std::remove_cv_t<T>>()))>> = true;
  }  // namespace Details

  /*!
   * \brief Overload of CastToDouble() for any type that is convertible to double.
   */
  template <typename ScalarType,
      std::enable_if_t<Details::is_double_convertible<ScalarType>, bool> = true>
  inline double CastToDouble(ScalarType a)
  {
    return static_cast<double>(a);
  }

  /**
   * \brief Cast of FAD to double.
   *
   * This function will be used exclusively for FAD or FAD expression types and can also handle
   * nested FAD types, i.e., higher order derivatives.
   *
   * @note We recursively call this function with the value of the FAD variable @p a until we reach
   * the innermost double value.
   */
  template <typename FADType,
      typename =
          std::enable_if_t<Sacado::IsFad<FADType>::value || Sacado::IsExpr<FADType>::value, void*>>
  inline double CastToDouble(FADType a)
  {
    return CastToDouble(a.val());
  }


  /*!
  \brief Cast of a FAD matrix to a double matrix
  */
  template <typename type, unsigned int dim1, unsigned int dim2>
  Core::LinAlg::Matrix<dim1, dim2, double> CastToDouble(Core::LinAlg::Matrix<dim1, dim2, type> a)
  {
    Core::LinAlg::Matrix<dim1, dim2, double> b(true);

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
  Core::FADUtils::sqrt<my_AD_type>(...)
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
  \brief Overload for double.
  */
  inline double sqrt(double a) { return std::sqrt(a); }

  /*!
  \brief Calculate Norm of a scalar FAD quantity. If a compiler Error is thrown when calling this
  function, check if the template argument is explicitly stated in the function call, i.e.
  Core::FADUtils::Norm<my_AD_type>(...)
  */
  template <typename scalar_type>
  inline scalar_type Norm(scalar_type a)
  {
    return Core::FADUtils::sqrt<scalar_type>(a * a);
  }

  /*!
  \brief Calculate Norm of a FAD vector
  */
  template <typename scalar_type, unsigned int length>
  scalar_type VectorNorm(Core::LinAlg::Matrix<length, 1, scalar_type> v)
  {
    scalar_type norm_squared = 0.0;
    for (unsigned int i = 0; i < length; i++)
    {
      norm_squared += v(i) * v(i);
    }

    return Core::FADUtils::sqrt<scalar_type>(norm_squared);
  }

  /*!
  \brief Template specialization for double
  */
  template <unsigned int length>
  double VectorNorm(Core::LinAlg::Matrix<length, 1, double> v)
  {
    return v.Norm2();
  }

  //! Calculates the Norm of a FAD vector, since .Norm2() is not available for FAD vectors
  // Todo this function is obsolete
  template <typename T>
  T Norm(Core::LinAlg::Matrix<3, 1, T> v)
  {
    T norm_squared = 0.0;
    for (int i = 0; i < 3; i++)
    {
      norm_squared += v(i) * v(i);
    }

    return Core::FADUtils::sqrt(norm_squared);
  }

  /*!
  \brief Calculate inner product of two FAD or double vectors
  */
  // Todo this function is obsolete, use Dot of Core::LinAlg::Matrix instead
  template <typename type>
  type ScalarProduct(Core::LinAlg::Matrix<3, 1, type> a, Core::LinAlg::Matrix<3, 1, type> b)
  {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  /*!
  \brief Calculate difference of two FAD or double vectors
  */
  // Todo this function is obsolete, use Update of Core::LinAlg::Matrix instead
  template <typename type>
  Core::LinAlg::Matrix<3, 1, type> DiffVector(
      Core::LinAlg::Matrix<3, 1, type> a, Core::LinAlg::Matrix<3, 1, type> b)
  {
    Core::LinAlg::Matrix<3, 1, type> c(true);
    for (int i = 0; i < 3; i++) c(i) = a(i) - b(i);

    return c;
  }

  /*!
  \brief Calculate vector product of two FAD or double vectors
  */
  // Todo this function is obsolete, use CrossProduct of Core::LinAlg::Matrix instead
  template <typename type>
  Core::LinAlg::Matrix<3, 1, type> VectorProduct(
      Core::LinAlg::Matrix<3, 1, type> first_vector, Core::LinAlg::Matrix<3, 1, type> second_vector)
  {
    Core::LinAlg::Matrix<3, 1, type> result_vector;
    result_vector.Clear();
    Core::LinAlg::Matrix<3, 3, type> S_first_vector;
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

}  // namespace Core::FADUtils

FOUR_C_NAMESPACE_CLOSE

#endif
