#ifndef BASE_VOL_H
#define BASE_VOL_H

#include<vector>

namespace GEO
{
  namespace CUT
  {
    //returns function which is used to form the matrix of moment fitting equations
    double base_function(std::vector<double> coordi,int base_num);
  }
}
#endif

