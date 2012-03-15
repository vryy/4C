#ifndef BASE_VOL_H
#define BASE_VOL_H

#include<vector>

namespace GEO
{
  namespace CUT
  {
  /*!
  \brief Returns the actual base function to be integrated over the volume to form the moment fitting matrix
  */
    double base_function(std::vector<double> coordi,int base_num);
  }
}
#endif

