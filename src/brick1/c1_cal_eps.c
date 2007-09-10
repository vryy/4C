/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_eps' which evaluates linear/nonlinear
                                     strains from displacement derivatives

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief compute evaluates linear/nonlinear strains

<pre>                                                              al 06/02
This routine evaluates linear/nonlinear strains from displacement derivatives
for an 3D-hex-element.

</pre>
\param disd     DOUBLE*  (i)   displacement derivatives
\param eps      DOUBLE*  (o)   strain vector
\param iform       INT   (i)   index for nonlinear formulation

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint() and routines for plasticity

*----------------------------------------------------------------------*/
void c1_eps( DOUBLE   *disd,                /* displacement derivatives */
             DOUBLE   *eps,                 /* strain vector            */
             INT      iform)         /* index for nonlinear formulation */
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE fac, half;
DOUBLE u11, u12, u13, u21, u22, u23, u31, u32, u33;
DOUBLE dn[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_eps");
#endif
/*----------------------------------------------------------------------*/
  half = 1./2.;
/*---------------------------------------------- linear strain terms ---*/
  u11 = disd[0];
  u22 = disd[1];
  u33 = disd[2];
  u12 = disd[3];
  u21 = disd[4];
  u23 = disd[5];
  u32 = disd[6];
  u13 = disd[7];
  u31 = disd[8];

  eps[0] = u11;
  eps[1] = u22;
  eps[2] = u33;
  eps[3] = u12 + u21;
  eps[4] = u23 + u32;
  eps[5] = u13 + u31;
/*--------------------------------------- add nonlinear strain terms ---*/
  if (iform>1)
  {
    dn[0] = half * (u11*u11 + u21*u21 + u31*u31);
    dn[1] = half * (u12*u12 + u22*u22 + u32*u32);
    dn[2] = half * (u13*u13 + u23*u23 + u33*u33);
    dn[3] = u11*u12 + u21*u22 + u31*u32;
    dn[4] = u12*u13 + u22*u23 + u32*u33;
    dn[5] = u11*u13 + u21*u23 + u31*u33;
/*----------------------------------------------------------------------*
 |   green-lagrange strains --> add      nonlinear strain terms         |
 |   almansi strains        --> subtract nonlinear strain terms         |
 *----------------------------------------------------------------------*/
      fac = 1.;
      if (iform==3) fac = -1.;
      for (i=0; i<6; i++) eps[i] += fac*dn[i];
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_eps */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/

#endif
