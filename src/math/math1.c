/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |  r(I) = A(I,K)*b(K)*factor -----  r = A*b*factor          m.gee 6/01 |
 |  or                                                                  |
 |  r(I) += A(I,K)*b(K)*factor                                          |
 *----------------------------------------------------------------------*/
void math_matvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor)
{
INT i,k;
DOUBLE sum;
/*----------------------------------------------------------------------*/
if (init==0)
{
   for (i=0; i<ni; i++) r[i]=0.0;
}
for (i=0; i<ni; i++)
{
   sum=0.0;
   for (k=0; k<nk; k++) sum += A[i][k]*b[k];
   r[i] += sum*factor;
}
/*----------------------------------------------------------------------*/
return;
} /* end of math_matvecdense */

