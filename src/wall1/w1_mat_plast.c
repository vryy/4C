#include "../headers/standardtypes.h"
#include "wall1.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- von Mises - 2D al 9/01|
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_plast_mises(double ym,
                        double pv,
                        double ALFAT,
                        double sigy,
                        double hard,
                        double gf,
                        ELEMENT   *ele,
                        WALL_TYPE wtype,
                        double **bop,
                        int ip,
                        double *stress,       
                        double **d,
                        int istore)
{
/*----------------------------------------------------------------------*/
int i,j,k;
int yip;
int isoft;
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
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_mises");
#endif
/*----------------------------------------------------------------------*/
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
          w1mapl(ym, hard, sigy, pv, tau, isoft, &epstn, &dlam, d, wtype);
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
  w1yilcr(ym, hard, sigy, epstn, isoft, tau, &ft);
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
    w1radi(ym, hard, sigy, pv, tau, qn, isoft, &epstn, &dlam, wtype);
    w1mapl(ym, hard, sigy, pv, tau, isoft, &epstn, &dlam, d, wtype);
    
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
/*--------------------------------------- consistent material matrix ---*/
/*      CALL MXCR8 (C,NUMEPS,NUMEPS,D) */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_mises */
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
                        int istore)
{
/*----------------------------------------------------------------------*/
int i,j,k;
int yip;
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
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_dp");
#endif
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
/*      IF (NINT(YIP).GT.0) THEN
        DO 20 I=1,4
        STRESS(I) = SIG(I)
   20   TAU(I) = SIG(I)   
C
        IF (NINT(YIP).EQ.1) THEN
          YIP=-YIP
        ELSE
          CALL W1LSS (TAU,G,GI,0)
          IF (INDMAT.GT.0) THEN
            DLAM = ZERO
            CALL W1MAPL (TAU,D,PROPM(1,MAT),DLAM,NUMCM,ITYPE,DIA,EPSTN)
            CALL W1GLD  (D,G)
          ENDIF
          YIP=-YIP
        ENDIF
        IUPD=1
        GO TO 900
      ENDIF
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
    w1mapl(ym, hard, sigy, pv, tau, &epstn, &dlam, d, wtype);
    
    for (i=0; i<4; i++) stress[i] = tau[i] + qn[i];
  }
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
  if(istore==1)
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
/*--------------------------------------- consistent material matrix ---*/
/*      CALL MXCR8 (C,NUMEPS,NUMEPS,D) */
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_dp */
