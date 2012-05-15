/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/


#ifndef MATH_PROTOTYPES_H
#define MATH_PROTOTYPES_H

/*----------------------------------------------------------------------*
 |  math1.c                                               m.gee 11/01   |
 *----------------------------------------------------------------------*/
void math_matvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor);


#endif
