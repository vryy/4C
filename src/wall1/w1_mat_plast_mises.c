/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_mat_plast_mises' which calculates the
       constitutive matrix - forces - linear elastic - von Mises - 2D
       (plane stress, plane strain)

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
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
 | constitutive matrix - forces - linear elastic- von Mises - 2D al 9/01|
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_plast_mises(DOUBLE ym,
                        DOUBLE pv,
                        DOUBLE sigy,
                        DOUBLE hard,
                        DOUBLE gf,
                        DOUBLE betah,
                        ELEMENT   *ele,
                        WALL_TYPE wtype,
                        DOUBLE **bop,
                        DOUBLE  *gop,
                        DOUBLE  *alpha,
                        INT ip,
                        DOUBLE *stress,
                        DOUBLE **d,
                        INT istore,/* controls storing of new stresses to wa */
                        INT newval)/* controls evaluation of new stresses    */
{
/*----------------------------------------------------------------------*/
INT i,j;
INT yip;
INT isoft;
INT iupd;
DOUBLE sum, epstn, ft;
DOUBLE disd[5];
DOUBLE sig[4];
DOUBLE eps[4];
DOUBLE strain[4];
DOUBLE delsig[4];
DOUBLE deleps[4];
DOUBLE tau[4];
DOUBLE qn[4];
DOUBLE tol = 1.0E-10;
DOUBLE dlam;
DOUBLE dia;

DOUBLE tauc[4];
#ifdef DEBUG
dstrc_enter("w1_mat_plast_mises");
#endif
/*----------------------------------------------------------------------*/
  dia = ele->e.w1->elewa[0].dia;
  iupd=0;
/*---------------------------------------------- hardening/softening ---*/
 if(fabs(gf)>tol)
 {
   isoft= 1;
   hard =gf;
 }
 else isoft=0;
/*------------ original global elastic matrix for current point -> D ---*/
  w1_mat_linel(ym, pv, wtype, d);
/*--------------------------------- compute displacement derivatives ---*/
  w1_disd (ele,bop,gop,alpha,wtype,disd) ;
/*------------------------------------- get actual strains -> strain ---*/
  w1_eps (disd,wtype,strain);
/*----------------------------- get old values -> sig, eps,epstn,yip ---*/
  for (i=0; i<4; i++)
  {
    sig[i] = ele->e.w1->elewa[0].ipwa[ip].sig[i];
    eps[i] = ele->e.w1->elewa[0].ipwa[ip].eps[i];
    qn[i]  = ele->e.w1->elewa[0].ipwa[ip].qn[ i];
  }
  epstn  = ele->e.w1->elewa[0].ipwa[ip].epstn;
  yip    = ele->e.w1->elewa[0].ipwa[ip].yip;
/*----------------------------------------------------------------------*/
  if(newval==1)
  {
    for (i=0; i<4; i++)  stress[i] = sig[i];
    goto end;
  }
/*----------------------------------------------------------------------*/
for (i=0; i<4; i++)
{
  stress[i] = 0.0;
  tau[i]    = 0.0;
}
/*-----------------------------------------------------------------------|
|  YIP > 0  STRESSES ARE AVAILABLE FROM LAST UPDATE  (GLOBAL PRAEDIKTOR) |
|      = 1  E L A S T I C                                                |
|      = 2  P L A S T I C                                                |
|     UPDATE FLAG MUST SET TO STORE CHANGE OF PARAMETER YIP              |
|     NO CHANGES HAVE BEEN MADE ON STRESS STATE                          |
|-----------------------------------------------------------------------*/
      if(yip>0)
      {
        for (i=0; i<4; i++)
        {
          stress[i] = sig[i];
          tau[i]    = sig[i];
        }
        if(yip==1)
        {
          yip=-yip;
        }
        else
        {
          dlam=0.;
          w1mapl(ym, hard, betah, sigy, pv, dia, tau, isoft, &epstn, &dlam, d, wtype);
          yip=-yip;
        }
        iupd=1;
        goto end;
      }
/*-----------------------------------------------------------------------|
|   1. CALCULATE INCREMENTAL STRAINS     DELEPS                          |
|   2. CALCULATE STRESS INCREMENT ASSUMING ELASTIC BEHAVIOUR             |
|   3. CALCULATE TOTAL STRESS                                            |
|   4. CHECK STRESS DEVIATOR AGAINST CURRENT YIELD SURFACE               |
|-----------------------------------------------------------------------*/

  for (i=0; i<4; i++) deleps[i] = strain[i] - eps[i];

  for (i=0; i<4; i++) delsig[i] = 0.0;
  for (i=0; i<4; i++) for (j=0; j<4; j++) delsig[i] += d[i][j]*deleps[j];

  for (i=0; i<4; i++) tau[i] = sig[i] + delsig[i] - qn[i];

  sum = 0.0;
  for (i=0; i<4; i++) sum += sig[i] * delsig[i];

  if(sum<0.0)
  {
    for (i=0; i<4; i++) stress[i] = tau[i] + qn[i];
    yip=1;
    goto end;
  }
/*-------------------------------------------------- yield condition ---*/
  w1yilcr(ym, hard, betah, sigy, epstn, isoft, dia, tau, &ft);
/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<tol)
  {
    yip = 1;
    for (i=0; i<4; i++) stress[i] = tau[i] + qn[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else
  {
    yip = 2;

    for (i=0; i<4; i++) tauc[i]=tau[i]; /*praedictor stresses -> tauc*/

    w1radi(ym, hard, betah, sigy, pv, dia, tau, qn, isoft, &epstn, &dlam, wtype);

    switch(wtype)
    {
    case plane_stress:/*new stresses*/
       w1mapl(ym, hard, betah, sigy, pv, dia, tau, isoft, &epstn, &dlam, d, wtype);
    break;
    case plane_strain:/*praedictor stresses*/
       w1mapl(ym, hard, betah, sigy, pv, dia, tauc, isoft, &epstn, &dlam, d, wtype);
    break;
   /*============================================== rotational symmetry ===*/
    default:
      dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
    break;
    }

    for (i=0; i<4; i++) stress[i] = tau[i] + qn[i];
  }
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
end:
/*----------------------------------------------------------------------*/
  if(istore==1 || iupd==1)
  {
    for (i=0; i<4; i++)
    {
      ele->e.w1->elewa[0].ipwa[ip].sig[i] = stress[i];
      ele->e.w1->elewa[0].ipwa[ip].eps[i] = strain[i];
      ele->e.w1->elewa[0].ipwa[ip].qn[ i] = qn[i] ;
    }
    ele->e.w1->elewa[0].ipwa[ip].epstn = epstn;
    ele->e.w1->elewa[0].ipwa[ip].yip   = yip  ;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_mises */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
