/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_utils.H

\brief utils for particle interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_UTILS_H
#define PARTICLE_INTERACTION_UTILS_H

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
      static constexpr T pow(const T x) { return helper<T, N - 1>::pow(x) * x; }
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
  }  // namespace UTILS

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
