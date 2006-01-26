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

#include "pss_set.h"

void intset_init(INTSET* s, INT size)
{
#ifdef DEBUG
  dstrc_enter("intset_init");
#endif

  s->size = size;
  s->count = 0;
  s->value = (INT*)CCAMALLOC(size*sizeof(INT));

#ifdef DEBUG
  dstrc_exit();
#endif
}


void intset_destroy(INTSET* s)
{
#ifdef DEBUG
  dstrc_enter("intset_destroy");
#endif
  CCAFREE(s->value);
#ifdef DEBUG
  dstrc_exit();
#endif
}


void intset_add(INTSET* s, INT v)
{
#ifdef DEBUG
  dstrc_enter("intset_add");
#endif
  if (s->count == s->size)
  {
    s->size *= 2;
    s->value = CCAREALLOC(s->value, s->size*sizeof(INT));
  }
  s->value[s->count++] = v;
  qsort(s->value, s->count, sizeof(INT), cmp_int);
#ifdef DEBUG
  dstrc_exit();
#endif
}


INT intset_contains(INTSET* s, INT v)
{
  return bsearch(&v, s->value, s->count, sizeof(INT), cmp_int) != NULL;
}


void intset_clear(INTSET* s)
{
  s->count = 0;
}


INT intset_get(INTSET* s, INT i)
{
#ifdef DEBUG
  dstrc_enter("intset_get");
#endif

  dsassert((i >= 0) && (i < s->count), "index out of range");

#ifdef DEBUG
  dstrc_exit();
#endif

  return s->value[i];
}
