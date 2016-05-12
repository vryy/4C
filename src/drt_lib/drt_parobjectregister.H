/*!----------------------------------------------------------------------
\file drt_parobjectregister.H

\brief registration functionality for parallel objects

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DRT_PAROBJECTREGISTER_H
#define DRT_PAROBJECTREGISTER_H


#include <string>

namespace DRT
{
  std::string ParObjectList();
}

/*! print list of known parobjects */
void PrintParObjectList();


#endif
