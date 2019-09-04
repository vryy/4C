/*----------------------------------------------------------------------*/
/*! \file

\brief General functions for the particle dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_UTILS_H_
#define PARTICLE_UTILS_H_
/*----------------------------------------------------------------------*/
/* headers */
#include <math.h>


namespace PARTICLE
{
  class Utils
  {
   public:
    //! Radius -> volume
    static double Radius2Volume(double radius)
    {
      return (4.0 / 3.0) * M_PI * radius * radius * radius;
    }

    //! Volume -> radius
    static double Volume2Radius(double volume)
    {
      return std::pow((3.0 / 4.0) * M_1_PI * volume, 1.0 / 3.0);
    }
  };
}  // namespace PARTICLE

#endif /* PARTICLE_UTILS_H_ */
