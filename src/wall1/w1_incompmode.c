/*!----------------------------------------------------------------------
\file
\brief

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

/*----------------------------------------------------------------------*
 | stiffness matrix  knc        (8*4)                           ah 9/02 |
 *----------------------------------------------------------------------*/
void w1_knc(DOUBLE  **knc,             /* stiffness knc= BT C G         */
            DOUBLE  **bop,             /* operator matrix               */
            DOUBLE   *gop,             /* additional opperator matrix   */
            DOUBLE  **d,               /* constitutive matrix           */
            DOUBLE    fac)             /* integration factor            */
{
INT            i, j, k, l, m;
DOUBLE         dum;
DOUBLE         cg[4];
DOUBLE         g[3][4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_knc");
#endif
/*------------------------------------------- vector gop -> matrix g ---*/
for(i=0;i<2;i++)
{
g[0][2*i]   = gop[2*i];
g[0][2*i+1] = 0;
g[1][2*i]   = 0;
g[1][2*i+1] = gop[2*i+1];
g[2][2*i]   = gop[2*i+1];
g[2][2*i+1] = gop[2*i];
}
/*----------------------------------------------------------------------*/
for (j=0; j<4; j++)
{
  for (k=0; k<3; k++)
  {
   cg[k] = 0.0 ;
   for (l=0; l<3; l++)
   {
    cg[k] += d[k][l] * g[l][j]*fac ;
   }
  }
  for (i=0; i<8; i++)
  {
    dum = 0.0 ;
    for (m=0; m<3; m++)
    {
     dum += bop[m][i] * cg[m] ;
    }
     knc[i][j] += dum ;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_knc */

/*----------------------------------------------------------------------*
 | stiffness matrix  knn      (4*4)                             ah 9/02 |
 *----------------------------------------------------------------------*/
void w1_knn(DOUBLE  **knn,             /* stiffness knn= GT C G         */
            DOUBLE   *gop,             /* additional opperator matrix   */
            DOUBLE  **d,               /* constitutive matrix           */
            DOUBLE    fac)             /* integration factor            */
{
INT            i, j, k, l, m;
DOUBLE         dum;
DOUBLE         cg[4];
DOUBLE         g[3][4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_knn");
#endif
/*------------------------------------------- vector gop -> matrix g ---*/
for(i=0;i<2;i++)
{
g[0][2*i]   = gop[2*i];
g[0][2*i+1] = 0;
g[1][2*i]   = 0;
g[1][2*i+1] = gop[2*i+1];
g[2][2*i]   = gop[2*i+1];
g[2][2*i+1] = gop[2*i];
}
/*----------------------------------------------------------------------*/
for (j=0; j<4; j++)
{
  for (k=0; k<3; k++)
  {
   cg[k] = 0.0 ;
   for (l=0; l<3; l++)
   {
    cg[k] += d[k][l] * g[l][j]*fac ;
   }
  }
  for (i=0; i<4; i++)
  {
    dum = 0.0 ;
    for (m=0; m<3; m++)
    {
     dum += g[m][i] * cg[m] ;
    }
     knn[i][j] += dum ;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_knn */

/*----------------------------------------------------------------------*
 | evaluate internal element forces due to incomp. modes    ah 9/02     |
 *----------------------------------------------------------------------*/
void w1_fintn(DOUBLE  *F,              /* stress                        */
              DOUBLE   fac,            /* integration factor            */
              DOUBLE  *gop,            /* additional opperator matrix   */
              DOUBLE  *fintn)          /* INT forces due to inc modes   */
{
/*----------------------------------------------------------------------*/
INT j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_fintn");
#endif
/*----------------------------------------------------------------------*/
  for (j=0; j<2; j++)
  {
    fintn[2*j]  += fac * ( gop[2*j] * F[0] + gop[2*j+1] * F[2] );
    fintn[2*j+1]+= fac * ( gop[2*j] * F[2] + gop[2*j+1] * F[1] );
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_fintn */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | static kondensation: Kele = K- knc*inverse(knn)*kcn          ah 9/02 |
 *----------------------------------------------------------------------*/
void  w1_stat_cond(DOUBLE **knninv,/*I:stiffness inverse(knn) (4x4)     */
                   DOUBLE **knc,   /*I: knc = BT C G (8x4)              */
                   DOUBLE **deltak,/*O: knc*inverse(knn)*kcn            */
                   DOUBLE  *fintn, /*I: fintn=GT*singma                 */
                   DOUBLE  *deltaf,/*O: knc*inverse(knn)*fintn          */
                   ELEMENT *ele)   /*I: actual element                  */
{
INT            i, j, k, l, m;
DOUBLE         dumK,dumF;
DOUBLE         help[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_stat_cond");
#endif
/*--------------------------------------- calculate Knc*inverse (Knn) ---*/
for (i=0; i<8; i++)
{
  for (k=0; k<4; k++)
  {
   help[k] = 0.0 ;
   for (l=0; l<4; l++)
   {
    help[k] +=  knc[i][l] * knninv[l][k];
   }
  }
/*------------------------------------------------- calculate delta K ---*/
  for (j=0; j<8; j++)
  {
    dumK = 0.0 ;
    for (m=0; m<4; m++)
    {
     dumK += help[m] * knc[j][m];
    }
     deltak[i][j] += dumK ;
   }
/*---------------------------------------------- calculate delta Fint ---*/
   if (fintn)
   {
     dumF = 0.0;
     for (m=0; m<4; m++)
     {
       dumF += help[m] * fintn[m];
     }
     deltaf[i] += dumF;
   }
}
/*---- store inv(Knn), Knc,fintn at element working array for update ---*/

 for(i=0;i<4;i++)
 {
    if(fintn)
    {
    ele->e.w1->elewa[0].imodewa[0].fintn[i]=fintn[i];
    }
    for(j=0;j<4;j++)
    ele->e.w1->elewa[0].imodewa[0].knninv[i][j]=knninv[i][j];
    for(j=0;j<8;j++)
    ele->e.w1->elewa[0].imodewa[0].knc[j][i]=knc[j][i];
 }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_stat_cond */

/*----------------------------------------------------------------------*
 | update of internal dof alpha                                  ah 9/02 |
 *----------------------------------------------------------------------*/
void  w1_updalpha(DOUBLE  *alpha,  /*O: internal dof for incomp modes        */
                  ELEMENT *ele,    /*I:actual element                        */
                  DOUBLE **knc,    /*I: mixed stiffness                      */
                  DOUBLE **knninv, /*I: inverse of incomp. stiffness         */
                  DOUBLE  *fintn,  /*I: INT. forces of inc. modes            */
                  INT      istore) /*I: flag:update after loadstep->istore=1 */
{
INT            i, j, k, l;
DOUBLE         dum,yip0,yip1,yip2,yip3;
DOUBLE         help[4];
DOUBLE         deltad[8];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_updalpha");
#endif
/*---------------------- check wether alpha is to be evaluated or not --*/
yip0    = ele->e.w1->elewa[0].ipwa[0].yip;
yip1    = ele->e.w1->elewa[0].ipwa[0].yip;
yip2    = ele->e.w1->elewa[0].ipwa[0].yip;
yip3    = ele->e.w1->elewa[0].ipwa[0].yip;
if(yip0*yip1 < 0. ||  yip0*yip2 < 0. || yip0*yip3 < 0.)
{
 dserror("Widerspruch ob globaler Praed oder Korrektor");
}
/*--------  --------------------------------------------------------------*/
if(yip0>0 || istore == 1)
{
 for (i=0;i<4;i++)
 alpha[i] = ele->e.w1->elewa[0].imodewa[0].alpha[i];
}
else
{
 for (i=0; i<4; i++)
 {
/*----------------------------------- read the necessary information ---*/
/*--- sol_residual.a.da = Deltad im Praed. und = delta d im Korrektor --*/
  deltad[2*i]   = ele->node[i]->sol_residual.a.da[0][0];
  deltad[2*i+1] = ele->node[i]->sol_residual.a.da[0][1];
  fintn[i]=ele->e.w1->elewa[0].imodewa[0].fintn[i];
  for(j=0;j<4;j++)
   knninv[i][j]=ele->e.w1->elewa[0].imodewa[0].knninv[i][j];
  for(j=0;j<8;j++)
   knc[j][i]=ele->e.w1->elewa[0].imodewa[0].knc[j][i];
 }
/*-------------------------------------------- evaluate delta alpha ---*/
 for (i=0; i<4; i++)
 {
  dum=0.0;
  for (k=0; k<4; k++)
  {
   help[k] = 0.0;
   for (l=0; l<8; l++)
   {
    help[k] += knc[l][k] * deltad[l];
   }
   help[k] += fintn[k];
   dum += knninv[i][k] * help[k];
  }
  alpha[i] -= dum ;
/*---------------------------------- evaluate and store actual alpha ---*/
  alpha[i] += ele->e.w1->elewa[0].imodewa[0].alpha[i];
  ele->e.w1->elewa[0].imodewa[0].alpha[i] = alpha[i];
 }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_updalpha */


#endif



























