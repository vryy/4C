#ifndef DRT_PAROBJECTREGISTER_H
#define DRT_PAROBJECTREGISTER_H

#ifdef __cplusplus /* stuff which is c++ specific */

#include <string>

namespace DRT
{
  std::string ParObjectList();
}

#endif

/*! print list of known parobjects */
void PrintParObjectList();


#endif
