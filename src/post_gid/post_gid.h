/*!
\file
\brief Postprocessing utility that takes ccarat output and produces a
GiD input file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

The general idea is that we cannot load the whole result data into
memory at once.

\author u.kue
\date 09/04

*/

#ifndef POST_GID_H
#define POST_GID_H

#include "../post_common/post_common.h"

#ifdef D_SHELL9

/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a shell9 discretization.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void write_shell9_field(PROBLEM_DATA* problem, INT num);

#endif


#endif
