#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - Drucker Prager - 2D     al    9/01    |
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_plast_dp(   double ym,
                        double pv,
                        double ALFAT,
                        double sigy,
                        double hard,
                        double phi,
                        ELEMENT   *ele,
                        WALL_TYPE wtype,
                        double **bop,
                        int ip,
                        double *stress,       
                        double **d,
                        int istore,/* controls storing of new stresses to wa */
                        int newval)/* controls evaluation of new stresses    */
{
/*----------------------------------------------------------------------*/
int i,j,k;
int yip;
int iupd;
double e1, e2, e3, a1, b1, c1, sum, epstn, ft;
double disd[5];
double sig[4];
double eps[4];
double strain[4];
double delsig[4];
double deleps[4];
double tau[4];
double qn[4];
double tol = 1.0E-10;
double dlam;
double sigym;
double rad;
double dia = 0.;
int    isoft = 0;
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_dp");
#endif
/*----------------------------------------------------------------------*/
  iupd=0;
/*----------------------------------------------------------------------*/
  rad   = atan(1.)/45.;
  phi   = phi * rad;          
/*------------ original global elastic matrix for current point -> D ---*/
  w1_mat_linel(ym, pv, wtype, d);
/*--------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,wtype,disd) ;                  
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
|     YIP > 0  STRESSES ARE AVAILABLE FROM LAST UPDATE                   |
|         = 1  E L A S T I C                                             |
|         = 2  P L A S T I C                                             |
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
          w1mapl(ym, hard, sigy, pv, dia, tau, isoft, &epstn, &dlam, d, wtype);
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
  w1yilcr_dp(ym, hard, phi, sigy, &sigym, epstn, tau, &ft);
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
    w1radi_dp(ym, hard, phi, sigy, pv, tau, qn, &epstn, &sigym, &dlam, wtype);
    w1mapl(ym, hard, sigy, pv, dia, tau, isoft, &epstn, &dlam, d, wtype);
    
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
} /* end of w1_mat_plast_dp */
