/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_eps' which evaluates the linear/nonlinear
       strains from displacement derivatives for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | evaluates linear/nonlinear strains from               al    9/01     |
 | displacement derivatives                                             |
 *----------------------------------------------------------------------*/
void w1_eps( DOUBLE   *disd,
             WALL_TYPE wtype,
             DOUBLE   *eps)
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE u11, u12, u21, u22, u33;
#ifdef DEBUG
dstrc_enter("w1_eps");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++) eps[i] = 0.0;

  u11 = disd[0];
  u22 = disd[1];
  u12 = disd[2];
  u21 = disd[3];
  eps[0] = u11;
  eps[1] = u22;
  eps[2] = u12 + u21;

  if(wtype==plane_strain)
  {
    eps[3] = 0.;
  }

  if(wtype==rotat_symmet)
  {
    u33 = disd[4];
    eps[3] = u33;
  }
  else
  {
    u33 = 0.0;
    eps[3] = u33;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_eps */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
