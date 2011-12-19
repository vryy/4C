/*!----------------------------------------------------------------------
\file beam3contact_utils.cpp
\brief Utility functions for beam contact

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "beam3contact_utils.H"


int sgn(double skalar)
{
  if (skalar < 0) return -1;
  else return 1;
}

#endif /*CCADISCRET*/
