/*---------------------------------------------------------------------*/
/*! \file

\brief used in boundary cell integration

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_BASE_VOL_HPP
#define FOUR_C_CUT_BASE_VOL_HPP

#include "4C_config.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE

#endif
