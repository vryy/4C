/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_disd' which calclate displacement
       derivatives for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief compute displacement derivatives

<pre>                                                              al 06/02
This routine calcuates the displacement derivatives for an 3D-hex-element.

</pre>
\param  bop   DOUBLE**   (i)   b-operator matrix
\param edis   DOUBLE*    (i)   element displacements
\param disd   DOUBLE*    (o)   displacement derivatives
\param  iel      INT     (i)   number of element nodes

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_disd(DOUBLE    **bop, /* b-operator matrix                      */
             DOUBLE    *edis, /* element displacements                  */
             DOUBLE    *disd, /* displacement derivatives               */
             INT         iel) /* number of element nodes                */
{
/*----------------------------------------------------------------------*/
INT    i, cc;
DOUBLE h1, h2, h3, ed1, ed2, ed3;
#ifdef DEBUG
dstrc_enter("c1_disd");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<9; i++) disd[i] = 0.0;
  cc=0;
  for (i=0; i<iel; i++)
  {
    h1 = bop[0][cc+0];
    h2 = bop[1][cc+1];
    h3 = bop[2][cc+2];

    ed1=edis[cc+0];
    ed2=edis[cc+1];
    ed3=edis[cc+2];

    cc+=3;

    disd[0] += h1*ed1;
    disd[1] += h2*ed2;
    disd[2] += h3*ed3;
    disd[3] += h2*ed1;
    disd[4] += h1*ed2;
    disd[5] += h3*ed2;
    disd[6] += h2*ed3;
    disd[7] += h3*ed1;
    disd[8] += h1*ed3;
  } /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_disd */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
