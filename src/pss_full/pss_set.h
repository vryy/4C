/*!
\file
\brief A very simple set implementation.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is really the most basic implementation possible. The aim is to
keep the code clear. As soon as we see the need for something more
elaborated we can easily change this code.

*/

#ifndef PSS_SET_H
#define PSS_SET_H

#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*/
/*!
  \brief Sorted set of integers.

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef struct _INTSET
{
  INT size;			/* number of entries allocated */
  INT count;			/* number of entries used */
  INT* value;			/* entry array */
} INTSET;


void intset_init(INTSET* s, INT size);
void intset_destroy(INTSET* s);
void intset_add(INTSET* s, INT v);
INT intset_contains(INTSET* s, INT v);
void intset_clear(INTSET* s);
INT intset_get(INTSET* s, INT i);

#endif
