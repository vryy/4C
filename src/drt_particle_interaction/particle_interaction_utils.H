/*---------------------------------------------------------------------------*/
/*!
\brief utils for particle interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_UTILS_H
#define PARTICLE_INTERACTION_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 01/2019 |
 *---------------------------------------------------------------------------*/
#include <cmath>

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  namespace UTILS
  {
    /**
     *  \brief providing a more efficient method to determine the power with integer exponents
     */
    template <class T, int N>
    struct helper
    {
      static_assert(N >= 0, "The exponent must be positive!");
      static constexpr T pow(const T x)
      {
        return ((N % 2) == 0 ? 1 : x) * helper<T, (N / 2)>::pow(x * x);
      }
    };

    template <class T>
    struct helper<T, 0>
    {
      static constexpr T pow(const T x) { return 1; }
    };

    /**
     *  \brief helper function
     *
     *  when you use this helper function there will be no need to explicitly insert the class type
     */
    template <int N, class T>
    T constexpr pow(T const x)
    {
      return helper<T, N>::pow(x);
    }

    //! @name collection of three dimensional vector operations
    //@{

    /**
     *  \brief clear vector c
     */
    template <class T>
    inline void vec_clear(T* c)
    {
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = 0.0;
    }

    /**
     *  \brief set vector a to vector c
     */
    template <class T>
    inline void vec_set(T* c, const T* a)
    {
      c[0] = a[0];
      c[1] = a[1];
      c[2] = a[2];
    }

    /**
     *  \brief add vector a to vector c
     */
    template <class T>
    inline void vec_add(T* c, const T* a)
    {
      c[0] += a[0];
      c[1] += a[1];
      c[2] += a[2];
    }

    /**
     *  \brief subtract vector a from vector c
     */
    template <class T>
    inline void vec_sub(T* c, const T* a)
    {
      c[0] -= a[0];
      c[1] -= a[1];
      c[2] -= a[2];
    }

    /**
     *  \brief scale vector c
     */
    template <class T>
    inline void vec_scale(T* c, const T fac)
    {
      c[0] *= fac;
      c[1] *= fac;
      c[2] *= fac;
    }

    /**
     *  \brief scale vector a and set to vector c
     */
    template <class T>
    inline void vec_setscale(T* c, const T fac, const T* a)
    {
      c[0] = fac * a[0];
      c[1] = fac * a[1];
      c[2] = fac * a[2];
    }

    /**
     *  \brief scale vector a and add to vector c
     */
    template <class T>
    inline void vec_addscale(T* c, const T fac, const T* a)
    {
      c[0] += fac * a[0];
      c[1] += fac * a[1];
      c[2] += fac * a[2];
    }

    /**
     *  \brief set cross product of vector a and vector b to vector c
     */
    template <class T>
    inline void vec_setcross(T* c, const T* a, const T* b)
    {
      c[0] = a[1] * b[2] - a[2] * b[1];
      c[1] = a[2] * b[0] - a[0] * b[2];
      c[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     *  \brief add cross product of vector a and vector b to vector c
     */
    template <class T>
    inline void vec_addcross(T* c, const T* a, const T* b)
    {
      c[0] += a[1] * b[2] - a[2] * b[1];
      c[1] += a[2] * b[0] - a[0] * b[2];
      c[2] += a[0] * b[1] - a[1] * b[0];
    }

    /**
     *  \brief return scalar product of vector a and vector b
     */
    template <class T>
    inline T vec_dot(const T* a, const T* b)
    {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /**
     *  \brief return 2-norm of vector a
     */
    template <class T>
    inline T vec_norm2(const T* a)
    {
      return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    //@}

  }  // namespace UTILS

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
