/*---------------------------------------------------------------------*/
/*!
\file base_vol.cpp

\brief used in boundary cell integration

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "base_vol.H"

/*---------------------------------------------------------------------------------------------------------*
 *   Returns the actual base function to be integrated over the volume to form the moment fitting
 *matrix    *
 *----------------------------------------------------------------------------------------------------------*/
double GEO::CUT::base_function(std::vector<double> coordi, int base_num)
{
  if (base_num == 1)  // f(x,y,z) = 1.0
    return 1.0;
  if (base_num == 2)  // f(x,y,z) = x
    return coordi[0];
  if (base_num == 3)  // f(x,y,z) = y
    return coordi[1];
  if (base_num == 4)  // f(x,y,z) = z
    return coordi[2];

  if (base_num == 5)  // f(x,y,z) = x^2
    return coordi[0] * coordi[0];
  if (base_num == 6)  // f(x,y,z) = xy
    return coordi[0] * coordi[1];
  if (base_num == 7)  // f(x,y,z) = xz
    return coordi[0] * coordi[2];
  if (base_num == 8)  // f(x,y,z) = y^2
    return coordi[1] * coordi[1];
  if (base_num == 9)  // f(x,y,z) = yz
    return coordi[1] * coordi[2];
  if (base_num == 10)  // f(x,y,z) = z^2
    return coordi[2] * coordi[2];

  if (base_num == 11)  // f(x,y,x) = x^3
    return coordi[0] * coordi[0] * coordi[0];
  if (base_num == 12)  // f(x,y,x) = x^2*y
    return coordi[0] * coordi[0] * coordi[1];
  if (base_num == 13)  // f(x,y,x) = x^2*z
    return coordi[0] * coordi[0] * coordi[2];
  if (base_num == 14)  // f(x,y,x) = xy^2
    return coordi[0] * coordi[1] * coordi[1];
  if (base_num == 15)  // f(x,y,x) = xyz
    return coordi[0] * coordi[1] * coordi[2];
  if (base_num == 16)  // f(x,y,x) = xz^2
    return coordi[0] * coordi[2] * coordi[2];
  if (base_num == 17)  // f(x,y,x) = y^3
    return coordi[1] * coordi[1] * coordi[1];
  if (base_num == 18)  // f(x,y,x) = y^2*z
    return coordi[1] * coordi[1] * coordi[2];
  if (base_num == 19)  // f(x,y,x) = yz^2
    return coordi[1] * coordi[2] * coordi[2];
  if (base_num == 20)  // f(x,y,x) = z^3
    return coordi[2] * coordi[2] * coordi[2];

  if (base_num == 21)  // f(x,y,z) = x^4
    return coordi[0] * coordi[0] * coordi[0] * coordi[0];
  if (base_num == 22)  // f(x,y,z) = x^3*y
    return coordi[0] * coordi[0] * coordi[0] * coordi[1];
  if (base_num == 23)  // f(x,y,z) = x^3*z
    return coordi[0] * coordi[0] * coordi[0] * coordi[2];
  if (base_num == 24)  // f(x,y,z) = x^2*y^2
    return coordi[0] * coordi[0] * coordi[1] * coordi[1];
  if (base_num == 25)  // f(x,y,z) = x^2*yz
    return coordi[0] * coordi[0] * coordi[1] * coordi[2];
  if (base_num == 26)  // f(x,y,z) = x^2*z^2
    return coordi[0] * coordi[0] * coordi[2] * coordi[2];
  if (base_num == 27)  // f(x,y,z) = xy^3
    return coordi[0] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 28)  // f(x,y,z) = xy^2*z
    return coordi[0] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 29)  // f(x,y,z) = xyz^2
    return coordi[0] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 30)  // f(x,y,z) = xz^3
    return coordi[0] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 31)  // f(x,y,z) = y^4
    return coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 32)  // f(x,y,z) = y^3*z
    return coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 33)  // f(x,y,z) = y^2*z^2
    return coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 34)  // f(x,y,z) = yz^3
    return coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 35)  // f(x,y,z) = z^4
    return coordi[2] * coordi[2] * coordi[2] * coordi[2];

  if (base_num == 36)  // f(x,y,z) = x^5
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[0];
  if (base_num == 37)  // f(x,y,z) = x^4*y
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[1];
  if (base_num == 38)  // f(x,y,z) = x^4*z
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[2];
  if (base_num == 39)  // f(x,y,z) = x^3*y^2
    return coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[1];
  if (base_num == 40)  // f(x,y,z) = x^3*yz
    return coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[2];
  if (base_num == 41)  // f(x,y,z) = x^3*z^2
    return coordi[0] * coordi[0] * coordi[0] * coordi[2] * coordi[2];
  if (base_num == 42)  // f(x,y,z) = x^2*y^3
    return coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 43)  // f(x,y,z) = x^2*y^2*z
    return coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 44)  // f(x,y,z) = x^2*yz^2
    return coordi[0] * coordi[0] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 45)  // f(x,y,z) = x^2*z^3
    return coordi[0] * coordi[0] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 46)  // f(x,y,z) = x*y^4
    return coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 47)  // f(x,y,z) = x*y^3*z
    return coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 48)  // f(x,y,z) = x*y^2*z^2
    return coordi[0] * coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 49)  // f(x,y,z) = x*y*z^3
    return coordi[0] * coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 50)  // f(x,y,z) = x*z^4
    return coordi[0] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 51)  // f(x,y,z) = y^5
    return coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 52)  // f(x,y,z) = y^4*z
    return coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 53)  // f(x,y,z) = y^3*z^2
    return coordi[1] * coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 54)  // f(x,y,z) = y^2*z^3
    return coordi[1] * coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 55)  // f(x,y,z) = y*z^4
    return coordi[1] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 56)  // f(x,y,z) = z^5
    return coordi[2] * coordi[2] * coordi[2] * coordi[2] * coordi[2];

  if (base_num == 57)  // f(x,y,z) = x^6
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[0];
  if (base_num == 58)  // f(x,y,z) = x^5*y
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[1];
  if (base_num == 59)  // f(x,y,z) = x^5*z
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[2];
  if (base_num == 60)  // f(x,y,z) = x^4*y^2
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[1];
  if (base_num == 61)  // f(x,y,z) = x^4*yz
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[2];
  if (base_num == 62)  // f(x,y,z) = x^4*z^2
    return coordi[0] * coordi[0] * coordi[0] * coordi[0] * coordi[2] * coordi[2];
  if (base_num == 63)  // f(x,y,z) = x^3*y^3
    return coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 64)  // f(x,y,z) = x^3*y^2*z
    return coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 65)  // f(x,y,z) = x^3*y*z^2
    return coordi[0] * coordi[0] * coordi[0] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 66)  // f(x,y,z) = x^3*z^3
    return coordi[0] * coordi[0] * coordi[0] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 67)  // f(x,y,z) = x^2*y^4
    return coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 68)  // f(x,y,z) = x^2*y^3*z
    return coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 69)  // f(x,y,z) = x^2*y^2*z^2
    return coordi[0] * coordi[0] * coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 70)  // f(x,y,z) = x^2*yz^3
    return coordi[0] * coordi[0] * coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 71)  // f(x,y,z) = x^2*z^4
    return coordi[0] * coordi[0] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 72)  // f(x,y,z) = xy^5
    return coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 73)  // f(x,y,z) = xy^4*z
    return coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 74)  // f(x,y,z) = xy^3*z^2
    return coordi[0] * coordi[1] * coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 75)  // f(x,y,z) = xy^2*z^3
    return coordi[0] * coordi[1] * coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 76)  // f(x,y,z) = xy*z^4
    return coordi[0] * coordi[1] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 77)  // f(x,y,z) = xz^5
    return coordi[0] * coordi[2] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 78)  // f(x,y,z) = y^6
    return coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[1];
  if (base_num == 79)  // f(x,y,z) = y^5*z
    return coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[2];
  if (base_num == 80)  // f(x,y,z) = y^4*z^2
    return coordi[1] * coordi[1] * coordi[1] * coordi[1] * coordi[2] * coordi[2];
  if (base_num == 81)  // f(x,y,z) = y^3*z^3
    return coordi[1] * coordi[1] * coordi[1] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 82)  // f(x,y,z) = y^2*z^4
    return coordi[1] * coordi[1] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 83)  // f(x,y,z) = y*z^5
    return coordi[1] * coordi[2] * coordi[2] * coordi[2] * coordi[2] * coordi[2];
  if (base_num == 84)  // f(x,y,z) = z^6
    return coordi[2] * coordi[2] * coordi[2] * coordi[2] * coordi[2] * coordi[2];

  //  std::dserror("the base function required is not defined");
  return 0.0;
}
