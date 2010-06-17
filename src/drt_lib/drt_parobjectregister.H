#ifndef DRT_PAROBJECTREGISTER_H
#define DRT_PAROBJECTREGISTER_H

#ifdef __cplusplus /* stuff which is c++ specific */

#include <string>

namespace DRT
{
  std::string ParObjectList();
}

#endif

/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

  /*! print list of known parobjects */
  void PrintParObjectList();

#ifdef __cplusplus
}
#endif

#endif
