/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1consig' which condenses the stress vector
       plane strain --> plane stress
 contains the routine 'w1concep' which condeses the constitutive tensor
       plane strain --> plane stress
 contains the routine 'w1de33' which evaluates the incremental strain
       in thickness direction for wall element (planes stress)

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------|
|      topic: condensed stress vector                                    |
|             plane strain --> plane stress                              |
|-----------------------------------------------------------------------*/
void w1consig (DOUBLE **d,/* current material matrix components d14-d44 */
               DOUBLE  *sigma,  /* current stresses (local)             */
               DOUBLE  *sigmac) /* condensed stress vector              */
{
/*----------------------------------------------------------------------*/
    static INT i;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_enter("w1consig");
#endif
/*---------------------------- condensed stress- or backstressvector ---*/
    for (i=0; i<4; i++) sigmac[i] = sigma[i] - sigma[3] * d[i][3]/d[3][3];
    sigmac[3] = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1consig */
/*-----------------------------------------------------------------------|
|      topic : condense the constitutive tensor                          |
|              plane strain --> plane stress                             |
|-----------------------------------------------------------------------*/
void w1concep (DOUBLE **d)   /* material matrix to be calculated        */
{
/*----------------------------------------------------------------------*/
  static INT    i, j;
  static DOUBLE fac, dm[3][3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_enter("w1concep");
#endif
/*---------------------------------- condensed constitutive tensor d ---*/
    fac = 1. / d[3][3];
    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
      {
          dm[i][j] = fac * d[3][i] * d[j][3];
    }}

    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
      {
          d[i][j] -= dm[i][j];
    }}
/*----------------------------------------------------------------------*/
    for (j=0; j<4; j++)
    {
      d[3][j] = 0.;
      d[j][3] = 0.;
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1concep */
/*-----------------------------------------------------------------------|
|    topic: evaluate the incremental strain in thickness direction       |
|           for wall element (plane stress)                              |
|-----------------------------------------------------------------------*/
void w1de33(DOUBLE *sigi,  /* stresses from last iteration step         */
            DOUBLE *epsi,  /* strains from last iteration step          */
            DOUBLE *di,    /* components d41,d42,d43,d44 of the         */
                           /* const. tensor from last iteration step    */
            DOUBLE *strain)/* current strains (local)                   */
{
/*----------------------------------------------------------------------*/
    static DOUBLE dezz;
    static INT i;
    static DOUBLE de33;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_enter("w1de33");
#endif
/*------------------------------------------ incremental strain de44 ---*/
    de33 = 0.;
    for (i = 0; i < 3; i++)
    {
	de33 += di[i] * (strain[i] - epsi[i]);
    }
    dezz = -(de33 + sigi[3]) / di[3];
/*-------------------------------------------------- total strain[3] ---*/
    strain[3] = dezz + epsi[3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
    dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1de33 */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
#endif
