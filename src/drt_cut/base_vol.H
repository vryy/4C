/*---------------------------------------------------------------------*/
/*!

\brief used in boundary cell integration

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#ifndef BASE_VOL_H
#define BASE_VOL_H

#include <vector>

namespace GEO
{
  namespace CUT
  {
    /*!
    \brief Returns the actual base function to be integrated over the volume to form the moment
    fitting matrix
    */
    double base_function(std::vector<double> coordi, int base_num);
  }  // namespace CUT
}  // namespace GEO
#endif
