/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_bope' which calculates B-operator
       at gausian points for nonlocal equiv. strains of a gradient
       enhanced wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief   calculates B-operator at gausian points for nonlocal equiv. strains
         for the gradient enhanced wall element

\param **bope           DOUBLE      O: B-operator
\param **derive         DOUBLE      I: Derivatives of Ansatz-functions
\param **xjm            DOUBLE      I: Jacobian Matrix
\param det              DOUBLE      I: Determinant of Jacobian

\return void
\sa calling:   nothing;
    called by: wgestatic_ke();

*----------------------------------------------------------------------*/
void wge_bope(DOUBLE    **bope,
              DOUBLE    **derive,
              DOUBLE    **xjm,
              DOUBLE      det)
{
#ifdef D_WALLGE

INT inode;
DOUBLE dum;
DOUBLE xji[2][2];

#ifdef DEBUG
dstrc_enter("wge_bope");
#endif
/*----------------------------------------------------------------------*/

/*---------------------------------------------- inverse of jacobian ---*/
dum = 1.0/det;
xji[0][0] = xjm[1][1]* dum;
xji[0][1] =-xjm[0][1]* dum;
xji[1][0] =-xjm[1][0]* dum;
xji[1][1] = xjm[0][0]* dum;
/*----------------------------- get operator b of global derivatives ---*/
for (inode=0; inode<4; inode++)
{
  bope[0][inode] += xji[0][0] * derive[0][inode]
                 +  xji[0][1] * derive[1][inode];
  bope[1][inode] += xji[1][0] * derive[0][inode]
                 +  xji[1][1] * derive[1][inode];
} /* end of loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_bope */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
