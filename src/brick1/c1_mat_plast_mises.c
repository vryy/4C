#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- von Mises - 3D al 9/01|
 *----------------------------------------------------------------------*/
void c1_mat_plast_mises(double ym,      /* young's modulus              */
                        double pv,      /* poisson's ratio              */
                        double alfat,   /* temperature expansion factor */
                        double uniax,   /* yield stresse                */
                        double fhard,   /* hardening modulus            */
                        double gf,      /* fracture energy              */
                        ELEMENT   *ele, /* actual element               */
                        double **bop,   /* derivative operator          */
                        int ip,         /* integration point Id         */
                        double *stress, /*ele stress (-resultant) vector*/      
                        double **d,     /* material matrix              */
                        double  *disd,  /* displacement derivatives     */
                        double g[6][6], /* transformation matrix        */
                        double gi[6][6],/* inverse of g                 */
                        int istore,     /* controls storing of stresses */
                        int newval)     /* controls eval. of stresses   */
{
/*----------------------------------------------------------------------*/
int i,j,k;
int yip;
int isoft;
int iupd;
double e1, e2, e3, a1, b1, c1, sum, epstn, ft;
double sig[6];
double eps[6];
double strain[6];
double delsig[6];
double deleps[6];
double tau[6];
double tol = 1.0E-10;
double dlam;


double yld, sm, sx, sy, sz, sxy, sxz, syz, sig2, hard;
double expo  = 0.;
double alpha = 0.;
#ifdef DEBUG 
dstrc_enter("c1_mat_plast_mises");
#endif
/*----------------------------------------------------------------------*/
  iupd=0;
/*------------ original global elastic matrix for current point -> D ---*/
  c1_mat_linel(ym, pv, d);
/*------------------------ transform local to global material matrix ---*/
  c1gld (d,g);
/*------------------------------------- get actual strains -> strain ---*/
  c1_eps (disd,strain,1);
/*----------------------------- get old values -> sig, eps,epstn,yip ---*/
  for (i=0; i<6; i++)
  {
    sig[i] = ele->e.c1->elewa[0].ipwa[ip].sig[i];
    eps[i] = ele->e.c1->elewa[0].ipwa[ip].eps[i];
  }
  epstn  = ele->e.c1->elewa[0].ipwa[ip].epstn;
  yip    = ele->e.c1->elewa[0].ipwa[ip].yip;
/*----------------------------------------------------------------------*/
  if(newval==1)
  {
    for (i=0; i<6; i++)  stress[i] = sig[i];
    goto end;
  }
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++)
  {
    stress[i] = 0.0;
    tau[i]    = 0.0;
  }
/*-----------------------------------------------------------------------|
|     yip > 0  stresses are available from last update                   |
|         = 1  e l a s t i c                                             |
|         = 2  p l a s t i c                                             |
|     update flag must set to store change of parameter yip              |
|     no changes have been made on stress state                          |
|-----------------------------------------------------------------------*/
      if(yip>0)
      {
        for (i=0; i<6; i++)
        {
          stress[i] = sig[i];
          strain[i] = eps[i];
          tau[i]    = sig[i];
        }
        if(yip==1)
        {
          yip=-yip;
        }
        else
        {
          dlam=0.;
          /*transform stresses to local coordinate system */
          c1trss2local (tau,gi);

          sm = (tau[0]+tau[1]+tau[2])/3.;
          sx = tau[0] - sm;
          sy = tau[1] - sm;
          sz = tau[2] - sm;
          sxy = tau[3];
          sxz = tau[4];
          syz = tau[5];
          sig2 = ( sx*sx + sy*sy + sz*sz)/2. + sxy*sxy + sxz*sxz + syz*syz;
          sig2*=2.;
          sig2=sqrt(sig2);
 
          c1matp1 (ym, fhard, uniax, pv, sig2,tau,epstn,dlam,d);
          c1gld (d,g);/* transform local to global material matrix */

          yip=-yip;
        }
        iupd=1;
        goto end;
      }
/*-------------------------------------------------- yield condition ---*/
/*------ determine uniaxial yield stress (with non-linear hardening) ---*/
  yld = uniax + fhard*epstn + alpha*(1.-exp(-expo*epstn));


/*-----------------------------------------------------------------------|
|   1. calculate incremental strains     deleps                          |
|   2. calculate stress increment assuming elastic behaviour             |
|   3. calculate total stress                                            |
|   4. check stress deviator against current yield surface               |
|-----------------------------------------------------------------------*/
  for (i=0; i<6; i++) deleps[i] = strain[i] - eps[i];
  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += d[i][j]*deleps[j];
  
  for (i=0; i<6; i++) tau[i] = sig[i] + delsig[i];
/*------------------------------------------------ verify plasticity ---*/
  sm  = (tau[0]+tau[1]+tau[2])/3.;
  sx  = tau[0] - sm;
  sy  = tau[1] - sm;
  sz  = tau[2] - sm;
  sxy = tau[3];
  sxz = tau[4];
  syz = tau[5];
/*-------------------------------------------------- yield condition ---*/
  sig2 = ( sx*sx + sy*sy + sz*sz)/2. + sxy*sxy + sxz*sxz + syz*syz;
  ft = sig2 - yld*yld/3.; 
/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<=tol) 
  {
    yip = 1;
    for (i=0; i<6; i++) stress[i] = tau[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    yip = 2;
    /*transform stresses to local coordinate system */
    c1trss2local (tau,gi);
    /* evaluate new stresses with stress projection */
    for (i=0; i<6; i++) stress[i] = tau[i];
    sig2*=2.;
    sig2=sqrt(sig2);
    c1rad1 (ym, fhard, uniax, pv, sig2, tau, &epstn, &dlam);
    
    /* evaluate new material matrix if requested    */
    
    c1matp1 (ym, fhard, uniax, pv, sig2,stress,epstn,dlam,d);
    c1gld (d,g);/* transform local to global material matrix */

    tau[0] += sm;
    tau[1] += sm;
    tau[2] += sm;
    /*transform stresses to global coordinate system */
    c1trss2global (tau,g);
    
    for (i=0; i<6; i++) stress[i] = tau[i];
  }
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
end:
/*----------------------------------------------------------------------*/
  if(istore==1 || iupd==1)
  {
    for (i=0; i<6; i++)
    {
      ele->e.c1->elewa[0].ipwa[ip].sig[i] = stress[i];
      ele->e.c1->elewa[0].ipwa[ip].eps[i] = strain[i];
    }
    ele->e.c1->elewa[0].ipwa[ip].epstn = epstn;
    ele->e.c1->elewa[0].ipwa[ip].yip   = yip  ;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_mat_plast_mises */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void c1rad1(double e,        /* young's modulus                         */
            double fhard,    /* hardening modulus                       */
            double uniax,    /* yield stresse                           */
            double vnu,      /* poisson's ratio                         */
            double sig2,
            double *dev,     /* elastic predicor projected onto yield   */
            double *epstn,   /* equivalent uniaxial plastic strain      */
            double *dlam)    /* increment of plastic multiplier         */
{
/*----------------------------------------------------------------------*/
int i;
int isoft1 = 0;
int nsoft  = 1;
double ro23, ro32, bmu, sm, epst, esig, f, fobj, dfdl, dum, dhard;
double rnorm[6];
double alpha = 0.;
double expo  = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1rad1");
#endif
/*----------------------------------------------------------------------*/
  ro32 = sqrt(3./2.);
  ro23 = 1./ro32;
  bmu = e/(2.+2.*vnu);
/*-------------------------------------------- initialized variables ---*/
  *dlam = 0.;
  fobj  = 1.;
/*----------------------------------------------------------------------*/
  sm = (dev[0]+dev[1]+dev[2])/3.;
  dev[0] = dev[0]-sm;
  dev[1] = dev[1]-sm;
  dev[2] = dev[2]-sm;
  dev[3] = dev[3];
  dev[4] = dev[4];
  dev[5] = dev[5];

  for (i=0;i<6;i++) rnorm[i]=dev[i]/sig2;
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
  i=0;
/* for non-linear hardening rule */
L500:
  ++i;
/* new plastic uniaxial deformation */

  epst = *epstn + ro23* *dlam;

/* new hardening rate */
                        
  dhard = fhard + alpha*(expo*exp(-expo*epst));

/*  new uniaxial yield stress */

  esig = ro23*(uniax + fhard*epst + alpha*(1.-exp(-expo*epst)));

/*  apply von mises yield criteria */

  f= sig2 - esig -2.*bmu* *dlam*fobj;

/*  derivative of the yield criteria with respect to plastic increment */

  dfdl = -2.*bmu*(fobj+dhard/(3.*bmu));

/*  new plastic increment */

  *dlam = *dlam - f/dfdl;

/* check convergence */

  if (esig == 0.) {
      if (fabs(f) > 1e-5) {
          if (i > 30) {
              dserror("local iteration exceeds limit");
          }
          goto L500;
      }
  } else {
      if ((dum = f / esig, fabs(dum)) > 1e-5) {
          if (i > 30) {
              dserror("local iteration exceeds limit");
          }
          goto L500;
      }
  }
/*----------------------------------------------------------------------*/
/* update plastic strain */

  *epstn = epst;

/* new stress state */  
  for (i=0;i<6;i++) dev[i] -= 2.*bmu*(*dlam)*rnorm[i]*fobj;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1rad1 */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void c1matp1(double e,       /* young's modulus                         */
             double fhard,   /* hardening modulus                       */
             double uniax,   /* yield stresse                           */
             double vnu,     /* poisson's ratio                         */
             double sig2,
             double *tau,    /* current stresses (local)                */
             double epstn,   /* equivalent uniaxial plastic strain      */ 
             double dlam,    /* increment of plastic multiplier         */
             double **cc)    /* material matrix to be calculated        */
{
/*----------------------------------------------------------------------*/
int i, j, k, l;
double alpha, expo, fobj, sq23, sm, sqsig, b, fac1, fac2, fac3;
double fkh, g, hard;
double rn[3][3];
double c[3][3][3][3];
double gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1matp1");
#endif
/*----------------------------------------------------------------------*/
  alpha = 0.;
  expo  = 0.;
  fobj  = 1.;
  sq23 = sqrt(2./3.);
/*----------------------------------------------------------------------*/
  g     = e/(2. + 2. * vnu);
  fkh   = e/(3. *(1. - 2.*vnu));
/*----------------------------------------------------------------------*/
/* derivative from yld with respect to epstn */
/* ramberg-osgood hardening law or another nonlinear hardening law */
  hard =fhard + alpha*(expo*exp(-expo*epstn));
/*----------------------------------------------------------------------*/
  sm = (tau[0]+tau[1]+tau[2])/3.;
  rn[0][0] = tau[0]-sm;
  rn[1][1] = tau[1]-sm;
  rn[2][2] = tau[2]-sm;
  rn[0][1] = tau[3];
  rn[1][0] = tau[3];
  rn[1][2] = tau[4];
  rn[2][1] = tau[4];
  rn[0][2] = tau[5];
  rn[2][0] = tau[5];
/*----------------------------------------------------------------------*/
  sqsig   = 0.;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      sqsig+=rn[i][j]*rn[i][j];
    }
  }
  sqsig=sqrt(sqsig);

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      rn[i][j]=rn[i][j]/sqsig;
    }
  }
/*----------------------------------------------------------------------*/
  b  = 1.-2.*g*dlam/sig2*fobj;
  fac1 = fkh -2.*g*b/3.;
  fac2 = g*b; 
  fac3 = 2.*g*(fobj*1./(fobj+hard/(3.*g)) - 1. + b );

  /*
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        for (l=0; l<3; l++)
        {
          c[i][j][k][l]= 
                   fac1 * gk[i][j] * gk[k][l] +
                   fac2 *(gk[i][k] * gk[j][l] + gk[i][l] * gk[j][k]) -
                   fac3 * rn[i][j] * rn[k][l];
  }}}}

  cc[0][0]=c[0][0][0][0];                                                     
  cc[0][1]=c[1][1][0][0];                                                     
  cc[0][2]=c[2][2][0][0];                                                     
  cc[0][3]=c[0][1][0][0];                                                     
  cc[0][4]=c[1][2][0][0];                                                     
  cc[0][5]=c[0][2][0][0];                                                     
                                                                         
  cc[1][0]=c[1][1][0][0];                                                     
  cc[1][1]=c[1][1][1][1];                                                     
  cc[1][2]=c[2][2][1][1];                                                     
  cc[1][3]=c[0][1][1][1];                                                     
  cc[1][4]=c[1][2][1][1];                                                     
  cc[1][5]=c[0][2][1][1];                                                     
                                                                         
  cc[2][0]=c[2][2][0][0];                                                     
  cc[2][1]=c[2][2][1][1];                                                     
  cc[2][2]=c[2][2][2][2];                                                     
  cc[2][3]=c[0][1][2][2];                                                     
  cc[2][4]=c[1][2][2][2];                                                     
  cc[2][5]=c[0][2][2][2];                                                     
                                                                         
  cc[3][0]=c[0][1][0][0];                                                     
  cc[3][1]=c[0][1][1][1];                                                     
  cc[3][2]=c[0][1][2][2];                                                     
  cc[3][3]=c[0][1][0][1];                                                     
  cc[3][4]=c[1][2][0][1];                                                     
  cc[3][5]=c[0][2][0][1];                                                     
                                                                         
  cc[4][0]=c[1][2][0][0];                                                     
  cc[4][1]=c[1][2][1][1];                                                     
  cc[4][2]=c[1][2][2][2];                                                     
  cc[4][3]=c[1][2][0][1];                                                     
  cc[4][4]=c[1][2][1][2];                                                     
  cc[4][5]=c[1][2][0][2];                                                     
                                                                         
  cc[5][0]=c[0][2][0][0];                                                     
  cc[5][1]=c[0][2][1][1];                                                     
  cc[5][2]=c[0][2][2][2];                                                     
  cc[5][3]=c[0][2][0][1];                                                     
  cc[5][4]=c[0][2][1][2];                                                     
  cc[5][5]=c[0][2][0][2];                                              
  */
  for (k=0; k<3; k++)
  {
    for (l=0; l<3; l++)
    {
      c[0][0][k][l]= fac1*gk[0][0]*gk[k][l] +
                     fac2*(gk[0][k]*gk[0][l]+gk[0][l]*gk[0][k]) -
                     fac3*rn[0][0]*rn[k][l];
      
      c[1][0][k][l]= fac1*gk[1][0]*gk[k][l] +
                     fac2*(gk[1][k]*gk[0][l]+gk[1][l]*gk[0][k]) -
                     fac3*rn[1][0]*rn[k][l];
      
      c[1][1][k][l]= fac1*gk[1][1]*gk[k][l] +
                     fac2*(gk[1][k]*gk[1][l]+gk[1][l]*gk[1][k]) -
                     fac3*rn[1][1]*rn[k][l];
      
      c[2][0][k][l]= fac1*gk[2][0]*gk[k][l] +
                     fac2*(gk[2][k]*gk[0][l]+gk[2][l]*gk[0][k]) -
                     fac3*rn[2][0]*rn[k][l];
      
      c[2][1][k][l]= fac1*gk[2][1]*gk[k][l] +
                     fac2*(gk[2][k]*gk[1][l]+gk[2][l]*gk[1][k]) -
                     fac3*rn[2][1]*rn[k][l];
      
      c[2][2][k][l]= fac1*gk[2][2]*gk[k][l] +
                     fac2*(gk[2][k]*gk[2][l]+gk[2][l]*gk[2][k]) -
                     fac3*rn[2][2]*rn[k][l];
  }}
/**/
 cc[0][0]=c[0][0][0][0];                                                   
 cc[0][1]=c[0][0][1][1];                                                   
 cc[0][2]=c[0][0][2][2];                                                  
 cc[0][3]=c[0][0][1][0];                                                  
 cc[0][4]=c[0][0][2][1];                                                  
 cc[0][5]=c[0][0][2][0];                                                  
                                              
 cc[1][0]=c[1][1][0][0];                                                    
 cc[1][1]=c[1][1][1][1];                                                   
 cc[1][2]=c[1][1][2][2];                                                  
 cc[1][3]=c[1][1][1][0];                                                  
 cc[1][4]=c[1][1][2][1];                                                  
 cc[1][5]=c[1][1][2][0];                                                   
                                
 cc[2][0]=c[2][2][0][0];                                                     
 cc[2][1]=c[2][2][1][1];                                                   
 cc[2][2]=c[2][2][2][2];                                                  
 cc[2][3]=c[2][2][1][0];                                                  
 cc[2][4]=c[2][2][2][1];                                                 
 cc[2][5]=c[2][2][2][0];                                                    
                                
 cc[3][0]=c[1][0][0][0];                                                    
 cc[3][1]=c[1][0][1][1];                                                   
 cc[3][2]=c[1][0][2][2];                                                   
 cc[3][3]=c[1][0][1][0];                                                  
 cc[3][4]=c[1][0][2][1];                                                  
 cc[3][5]=c[1][0][2][0];                                                    
                                
 cc[4][0]=c[2][1][0][0];                                                    
 cc[4][1]=c[2][1][1][1];                                                   
 cc[4][2]=c[2][1][2][2];                                                   
 cc[4][3]=c[2][1][1][0];                                                  
 cc[4][4]=c[2][1][2][1];                                                  
 cc[4][5]=c[2][1][2][0];                                                   
                                
 cc[5][0]=c[2][0][0][0];                                                    
 cc[5][1]=c[2][0][1][1];                                                   
 cc[5][2]=c[2][0][2][2];                                                   
 cc[5][3]=c[2][0][1][0];                                                 
 cc[5][4]=c[2][0][2][1];                                                  
 cc[5][5]=c[2][0][2][0];                                                     
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif  
return; 
} /* end of c1matp1 */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
