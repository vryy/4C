/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9_tvbo: which calculates the B-Oparator for a shell9 element


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief B-Operator for compatible strains

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the B-Operator for compatible strains
</pre>
\param  DOUBLE   **bop     (o)  B-Operator matrix
\param  DOUBLE    *funct   (i)  shape functions at GP
\param  DOUBLE   **deriv   (i)  shape function derivatives at GP
\param  INT        iel     (i)  number of nodes to this element
\param  INT        numdf   (i)  number of dofs to one node
\param  DOUBLE  ***akov    (i)  kovariant basis vectors at reference layer of each kinematic layer
\param  DOUBLE  ***a3kvp   (i)  partial derivatives of a3_L for each kinematic layer
\param  INT        num_klay(i)  number of kin layers to this element
\param  INT        klay    (i)  actual kinematic layer
\param  INT        nsansq  (i)  number of collocation points (ANS)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_tvbo(DOUBLE    **bop,
             DOUBLE     *funct,
             DOUBLE    **deriv,
             INT         iel,
             INT         numdf,
             DOUBLE   ***akov,
             DOUBLE   ***a3kvp,
             INT         num_klay,    /* number of kin layers to this element */
             INT         klay,        /* actual kin layer */
             DOUBLE      condfac,
             INT         nsansq)
{
INT    jlay,inode,node_start;
DOUBLE a1x,a1y,a1z,a2x,a2y,a2z;
DOUBLE a3xi,a3yi,a3zi;
DOUBLE a3xj,a3yj,a3zj;
DOUBLE a31xj,a31yj,a31zj,a32xj,a32yj,a32zj;
DOUBLE pk,pk1,pk2;
INT    idof, jdof;
DOUBLE fac, fac1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_tvbo");
#endif
/*-------- Initialize bob[12][NUMDOF_SHELL9*MAXNOD_SHELL9] to 0.0 ------------------*/
/*for (i=0; i<12; i++)
{
   for (j=0; j<numdof_shell9*MAXNOD_SHELL9; j++) bop[i][j]= 0.0;
}*/
/*----------------------------------------------------------------------*/
  idof = 3 + (klay) * 3;
  a3xi=akov[0][2][klay];
  a3yi=akov[1][2][klay];
  a3zi=akov[2][2][klay];

/*---------- loop over all kinematic layers ---------------------------*/
for (jlay=0; jlay<num_klay; jlay++)
{
  /* is jlay on trajectory of reference layer to klay (=ilay) ? */
  fac1 = s9notr(num_klay,klay,jlay);
  if (fac1 == 0.0) continue;

  /* independent part of xsi */
  fac = s9ksi(num_klay,klay,jlay,condfac);
  jdof = 3 + (jlay) * 3;

  a1x=akov[0][0][jlay];
  a1y=akov[1][0][jlay];
  a1z=akov[2][0][jlay];
  a2x=akov[0][1][jlay];
  a2y=akov[1][1][jlay];
  a2z=akov[2][1][jlay];
  a3xj=akov[0][2][jlay];
  a3yj=akov[1][2][jlay];
  a3zj=akov[2][2][jlay];
  a31xj=a3kvp[0][0][jlay];
  a31yj=a3kvp[1][0][jlay];
  a31zj=a3kvp[2][0][jlay];
  a32xj=a3kvp[0][1][jlay];
  a32yj=a3kvp[1][1][jlay];
  a32zj=a3kvp[2][1][jlay];

/*----- loop over all nodes -----------------------------------------------*/
   for (inode=0; inode<iel; inode++)
   {
      pk  = funct[inode];
      pk1 = deriv[0][inode];
      pk2 = deriv[1][inode];

      node_start = inode*numdf;

      if (klay == jlay)
      {
          /* e11 (const) */
          bop[0][node_start+0] += pk1*a1x;
          bop[0][node_start+1] += pk1*a1y;
          bop[0][node_start+2] += pk1*a1z;

          /* e12 (const) */
          bop[1][node_start+0] += pk2*a1x + pk1*a2x;
          bop[1][node_start+1] += pk2*a1y + pk1*a2y;
          bop[1][node_start+2] += pk2*a1z + pk1*a2z;

          /* e22 (const) */
          bop[2][node_start+0] += pk2*a2x;
          bop[2][node_start+1] += pk2*a2y;
          bop[2][node_start+2] += pk2*a2z;

          if (nsansq == 0) /*do this only if there is no Querschub-ANS*/
          {
            /* e13 (const) */
            bop[3][node_start+0] += pk1*a3xi;
            bop[3][node_start+1] += pk1*a3yi;
            bop[3][node_start+2] += pk1*a3zi;

            bop[3][node_start+idof+0] += pk *a1x;
            bop[3][node_start+idof+1] += pk *a1y;
            bop[3][node_start+idof+2] += pk *a1z;

            /* e23 (const) */
            bop[4][node_start+0] += pk2*a3xi;
            bop[4][node_start+1] += pk2*a3yi;
            bop[4][node_start+2] += pk2*a3zi;

            bop[4][node_start+idof+0] += pk *a2x;
            bop[4][node_start+idof+1] += pk *a2y;
            bop[4][node_start+idof+2] += pk *a2z;
          }

          /* e33 (const) */
          bop[5][node_start+idof+0] += pk *a3xi;
          bop[5][node_start+idof+1] += pk *a3yi;
          bop[5][node_start+idof+2] += pk *a3zi;

          /* e11 (lin) */
          bop[6][node_start+0]= pk1*a31xj;
          bop[6][node_start+1]= pk1*a31yj;
          bop[6][node_start+2]= pk1*a31zj;

          bop[6][node_start+jdof+0]= pk1*a1x;
          bop[6][node_start+jdof+1]= pk1*a1y;
          bop[6][node_start+jdof+2]= pk1*a1z;

          /* e12 (lin) */
          bop[7][node_start+0]= pk1*a32xj + pk2*a31xj;
          bop[7][node_start+1]= pk1*a32yj + pk2*a31yj;
          bop[7][node_start+2]= pk1*a32zj + pk2*a31zj;

          bop[7][node_start+jdof+0]= pk1* a2x + pk2* a1x;
          bop[7][node_start+jdof+1]= pk1* a2y + pk2* a1y;
          bop[7][node_start+jdof+2]= pk1* a2z + pk2* a1z;

          /* e22 (lin) */
          bop[8][node_start+0]= pk2*a32xj;
          bop[8][node_start+1]= pk2*a32yj;
          bop[8][node_start+2]= pk2*a32zj;

          bop[8][node_start+jdof+0]= pk2*a2x;
          bop[8][node_start+jdof+1]= pk2*a2y;
          bop[8][node_start+jdof+2]= pk2*a2z;

          /* e13 (lin) */
          bop[9][node_start+jdof+0]= pk1*a3xi;
          bop[9][node_start+jdof+1]= pk1*a3yi;
          bop[9][node_start+jdof+2]= pk1*a3zi;

          bop[9][node_start+idof+0] += pk *a31xj;
          bop[9][node_start+idof+1] += pk *a31yj;
          bop[9][node_start+idof+2] += pk *a31zj;

          /* e23 (lin) */
          bop[10][node_start+jdof+0]= pk2*a3xi;
          bop[10][node_start+jdof+1]= pk2*a3yi;
          bop[10][node_start+jdof+2]= pk2*a3zi;

          bop[10][node_start+idof+0] += pk *a32xj;
          bop[10][node_start+idof+1] += pk *a32yj;
          bop[10][node_start+idof+2] += pk *a32zj;

          /* e33 (lin) -> 0.0 */
      }
      else if (klay != jlay)
      {
          /* e11 (const) */
          bop[0][node_start+0] += fac*pk1*a31xj;
          bop[0][node_start+1] += fac*pk1*a31yj;
          bop[0][node_start+2] += fac*pk1*a31zj;

          bop[0][node_start+jdof+0] += fac*pk1*a1x;
          bop[0][node_start+jdof+1] += fac*pk1*a1y;
          bop[0][node_start+jdof+2] += fac*pk1*a1z;

          /* e12 (const) */
          bop[1][node_start+0] += fac * (pk2*a31xj + pk1*a32xj);
          bop[1][node_start+1] += fac * (pk2*a31yj + pk1*a32yj);
          bop[1][node_start+2] += fac * (pk2*a31zj + pk1*a32zj);

          bop[1][node_start+jdof+0] += fac * (pk1*a2x + pk2*a1x);
          bop[1][node_start+jdof+1] += fac * (pk1*a2y + pk2*a1y);
          bop[1][node_start+jdof+2] += fac * (pk1*a2z + pk2*a1z);

          /* e22 (const) */
          bop[2][node_start+0] += fac * pk2*a32xj;
          bop[2][node_start+1] += fac * pk2*a32yj;
          bop[2][node_start+2] += fac * pk2*a32zj;

          bop[2][node_start+jdof+0] += fac * pk2*a2x;
          bop[2][node_start+jdof+1] += fac * pk2*a2y;
          bop[2][node_start+jdof+2] += fac * pk2*a2z;

          if (nsansq == 0) /*do this only if there is no Querschub-ANS*/
          {
            /* e13 (const) */
            bop[3][node_start+jdof+0] += fac * pk1*a3xi;
            bop[3][node_start+jdof+1] += fac * pk1*a3yi;
            bop[3][node_start+jdof+2] += fac * pk1*a3zi;

            bop[3][node_start+idof+0] += fac * pk*a31xj;
            bop[3][node_start+idof+1] += fac * pk*a31yj;
            bop[3][node_start+idof+2] += fac * pk*a31zj;

            /* e23 (const) */
            bop[4][node_start+jdof+0] += fac * pk2*a3xi;
            bop[4][node_start+jdof+1] += fac * pk2*a3yi;
            bop[4][node_start+jdof+2] += fac * pk2*a3zi;

            bop[4][node_start+idof+0] += fac * pk*a32xj;
            bop[4][node_start+idof+1] += fac * pk*a32yj;
            bop[4][node_start+idof+2] += fac * pk*a32zj;
          }
      }

   } /* end of loop over nodes */

} /* end of loop over all kinematic layers*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_tvbo */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
#endif
