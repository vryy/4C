/*---------------------------------------------------------------------------*/
/*! \file
\brief utils for particle interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_UTILS_H
#define PARTICLE_INTERACTION_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "lib_dserror.H"

#include <cmath>

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  namespace UTILS
  {
    /**
     *  \brief provide an efficient method to determine the power with integer exponents
     */
    template <class T, int N>
    struct helper
    {
      static_assert(N >= 0, "The exponent must be positive!");
      static constexpr T Pow(const T x)
      {
        return ((N % 2) == 0 ? 1 : x) * helper<T, (N / 2)>::Pow(x * x);
      }
    };

    template <class T>
    struct helper<T, 0>
    {
      static constexpr T Pow(const T x) { return 1; }
    };

    /**
     *  \brief helper function
     *
     *  when you use this helper function there will be no need to explicitly insert the class type
     */
    template <int N, class T>
    T constexpr Pow(T const x)
    {
      return helper<T, N>::Pow(x);
    }

    //! @name collection of three dimensional vector operations
    //@{

    /**
     *  \brief clear vector c
     */
    template <class T>
    inline void VecClear(T* c)
    {
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = 0.0;
    }

    /**
     *  \brief set vector a to vector c
     */
    template <class T>
    inline void VecSet(T* c, const T* a)
    {
      c[0] = a[0];
      c[1] = a[1];
      c[2] = a[2];
    }

    /**
     *  \brief add vector a to vector c
     */
    template <class T>
    inline void VecAdd(T* c, const T* a)
    {
      c[0] += a[0];
      c[1] += a[1];
      c[2] += a[2];
    }

    /**
     *  \brief subtract vector a from vector c
     */
    template <class T>
    inline void VecSub(T* c, const T* a)
    {
      c[0] -= a[0];
      c[1] -= a[1];
      c[2] -= a[2];
    }

    /**
     *  \brief scale vector c
     */
    template <class T>
    inline void VecScale(T* c, const T fac)
    {
      c[0] *= fac;
      c[1] *= fac;
      c[2] *= fac;
    }

    /**
     *  \brief scale vector a and set to vector c
     */
    template <class T>
    inline void VecSetScale(T* c, const T fac, const T* a)
    {
      c[0] = fac * a[0];
      c[1] = fac * a[1];
      c[2] = fac * a[2];
    }

    /**
     *  \brief scale vector a and add to vector c
     */
    template <class T>
    inline void VecAddScale(T* c, const T fac, const T* a)
    {
      c[0] += fac * a[0];
      c[1] += fac * a[1];
      c[2] += fac * a[2];
    }

    /**
     *  \brief set cross product of vector a and vector b to vector c
     */
    template <class T>
    inline void VecSetCross(T* c, const T* a, const T* b)
    {
      c[0] = a[1] * b[2] - a[2] * b[1];
      c[1] = a[2] * b[0] - a[0] * b[2];
      c[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     *  \brief add cross product of vector a and vector b to vector c
     */
    template <class T>
    inline void VecAddCross(T* c, const T* a, const T* b)
    {
      c[0] += a[1] * b[2] - a[2] * b[1];
      c[1] += a[2] * b[0] - a[0] * b[2];
      c[2] += a[0] * b[1] - a[1] * b[0];
    }

    /**
     *  \brief return scalar product of vector a and vector b
     */
    template <class T>
    inline T VecDot(const T* a, const T* b)
    {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /**
     *  \brief return 2-norm of vector a
     */
    template <class T>
    inline T VecNormTwo(const T* a)
    {
      return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    //@}

    //! @name methods for construction of three dimensional vector space
    //@{

    /**
     *  \brief construct orthogonal unit surface tangent vectors from given unit surface normal
     */
    template <class T>
    inline void UnitSurfaceTangents(const T* n, T* t1, T* t2)
    {
#ifdef DEBUG
      if (std::abs(1.0 - VecNormTwo(n)) > 1.0e-14)
        dserror("given unit surface normal not normalized!");
#endif

      if ((std::abs(n[0]) <= std::abs(n[1])) and (std::abs(n[0]) <= std::abs(n[2])))
      {
        t1[0] = 0.0;
        t1[1] = -n[2];
        t1[2] = n[1];
      }
      else if (std::abs(n[1]) <= std::abs(n[2]))
      {
        t1[0] = -n[2];
        t1[1] = 0.0;
        t1[2] = n[0];
      }
      else
      {
        t1[0] = -n[1];
        t1[1] = n[0];
        t1[2] = 0.0;
      }

      VecScale(t1, 1.0 / VecNormTwo(t1));

      VecSetCross(t2, n, t1);
    }

    //@}

    //! @name methods for linear transition in a given interval
    //@{

    /**
     *  \brief linear transition function
     */
    inline double LinTrans(const double x, const double x1, const double x2)
    {
#ifdef DEBUG
      if (not(std::abs(x2 - x1) > 1.0e-14)) dserror("danger of division by zero!");
#endif

      if (x < x1) return 0.0;
      if (x > x2) return 1.0;
      return (x - x1) / (x2 - x1);
    }

    /**
     *  \brief complementary linear transition function
     */
    inline double CompLinTrans(const double x, const double x1, const double x2)
    {
      return 1.0 - LinTrans(x, x1, x2);
    }

    //@}

  }  // namespace UTILS

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
