
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_UTILS_H
#define FSI_UTILS_H

#include <set>
#include <string>
#include <functional>

#include "../drt_lib/drt_discret.H"


namespace FSI
{
  namespace Utils
  {

    /// return the first slot of a pair
    /*!
      To be used with stl algorithms.

      This should be part of stl but is not (yet). So we define our
      own version.
    */
    template <typename pair_type>
    struct select1st: public std::unary_function<const pair_type &,
                                                 const typename pair_type::first_type &>
    {
      const typename pair_type::first_type &operator()(const pair_type &v) const
        {
          return v.first;
        }
    };

    /// return the second slot of a pair
    /*!
      To be used with stl algorithms.

      This should be part of stl but is not (yet). So we define our
      own version.
    */
    template <typename pair_type>
    struct select2nd: public std::unary_function<const pair_type &,
                                                 const typename pair_type::second_type &>
    {
      const typename pair_type::second_type &operator()(const pair_type &v) const
        {
          return v.second;
        }
    };

 }
}

#endif
#endif
#endif
