/*---------------------------------------------------------------------*/
/*! \file

\brief used in boundary cell integration

\level 3


*----------------------------------------------------------------------*/

#ifndef BACI_CUT_BASE_VOL_H
#define BACI_CUT_BASE_VOL_H

#include "baci_config.H"

#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    /*!
    \brief Returns the actual base function to be integrated over the volume to form the moment
    fitting matrix
    */
    double base_function(std::vector<double> coordi, int base_num);
  }  // namespace CUT
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif
