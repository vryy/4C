/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_mat_plast_mises' to establish local material law
       stress-strain law for mises material for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1
#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief establish local material

<pre>                                                              al 06/02
This routine to establish local material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param        ym  DOUBLE    (i) young's modulus
\param        pv  DOUBLE    (i) poisson's ratio
\param     alfat  DOUBLE    (i) temperature expansion factor
\param     uniax  DOUBLE    (i) yield stresse
\param     fhard  DOUBLE    (i) hardening modulus
\param        gf  DOUBLE    (i) fracture energy
\param       ele  ELEMENT*  (i) actual element
\param       bop  DOUBLE**  (i) derivative operator
\param        ip  INT       (i) integration point Id
\param    stress  DOUBLE*   (o) ele stress (-resultant) vector
\param         d  DOUBLE**  (o) material matrix
\param      disd  DOUBLE*   (i) displacement derivatives
\param   g[6][6]  DOUBLE    (i) transformation matrix
\param  gi[6][6]  DOUBLE    (i) inverse of g
\param    istore  INT       (i) controls storing of stresses
\param    newval  INT       (i) controls eval. of stresses

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_mat_plast_mises(DOUBLE ym,      /* young's modulus              */
                        DOUBLE pv,      /* poisson's ratio              */
                        DOUBLE uniax,   /* yield stresse                */
                        DOUBLE fhard,   /* hardening modulus            */
                        ELEMENT   *ele, /* actual element               */
                        INT ip,         /* integration point Id         */
                        DOUBLE *stress, /*ele stress (-resultant) vector*/
                        DOUBLE **d,     /* material matrix              */
                        DOUBLE  *disd,  /* displacement derivatives     */
                        DOUBLE g[6][6], /* transformation matrix        */
                        DOUBLE gi[6][6],/* inverse of g                 */
                        INT istore,     /* controls storing of stresses */
                        INT newval)     /* controls eval. of stresses   */
{
/*----------------------------------------------------------------------*/
INT i,j;
INT yip;
INT iupd;
DOUBLE epstn, ft;
DOUBLE sig[6];
DOUBLE eps[6];
DOUBLE strain[6];
DOUBLE delsig[6];
DOUBLE deleps[6];
DOUBLE tau[6];
DOUBLE tol = 1.0E-10;
DOUBLE dlam;


DOUBLE yld, sm, sx, sy, sz, sxy, sxz, syz, sig2;
DOUBLE expo  = 0.;
DOUBLE alpha = 0.;
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

          c1matp1 (ym, fhard, pv, sig2,tau,epstn,dlam,d);
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

    c1matp1 (ym, fhard, pv, sig2,stress,epstn,dlam,d);
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
/*!----------------------------------------------------------------------
\brief radial return for elements with mises material model

<pre>                                                              al 06/02
This routine performs a
radial return for elements with mises material model.
Nonlinear Ramberg-Osgood -Hardening law.

</pre>
\param e     DOUBLE  (i) young's modulus
\param fhard DOUBLE  (i) hardening modulus
\param uniax DOUBLE  (i) yield stresse
\param vnu   DOUBLE  (i) poisson's ratio
\param sig2  DOUBLE  (i) equivalent stress
\param dev   DOUBLE* (i) elastic predicor projected onto yield
\param epstn DOUBLE* (i) equivalent uniaxial plastic strain
\param dlam  DOUBLE* (i) increment of plastic multiplier

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1rad1(DOUBLE e,        /* young's modulus                         */
            DOUBLE fhard,    /* hardening modulus                       */
            DOUBLE uniax,    /* yield stresse                           */
            DOUBLE vnu,      /* poisson's ratio                         */
            DOUBLE sig2,     /* equivalent stress                       */
            DOUBLE *dev,     /* elastic predicor projected onto yield   */
            DOUBLE *epstn,   /* equivalent uniaxial plastic strain      */
            DOUBLE *dlam)    /* increment of plastic multiplier         */
{
/*----------------------------------------------------------------------*/
INT i;
INT isoft1 = 0;
INT nsoft  = 1;
DOUBLE ro23, ro32, bmu, sm, epst, esig, f, fobj, dfdl, dum, dhard;
DOUBLE rnorm[6];
DOUBLE alpha = 0.;
DOUBLE expo  = 0.;
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
/*!----------------------------------------------------------------------
\brief forms the elasto-plastic consistent tangent material tenso

<pre>                                                              al 06/02
This routine forms the elasto-plastic consistent tangent material tensor.

</pre>
\param e     DOUBLE  (i) young's modulus
\param fhard DOUBLE  (i) hardening modulus
\param uniax DOUBLE  (i) yield stresse
\param vnu   DOUBLE  (i) poisson's ratio
\param sig2  DOUBLE  (i) equivalent stress
\param tau   DOUBLE* (i) current stresses (local)
\param epstn DOUBLE* (i) equivalent uniaxial plastic strain
\param dlam  DOUBLE* (i) increment of plastic multiplier
\param cc    DOUBLE**(o) material matrix to be calculated

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1matp1(DOUBLE e,       /* young's modulus                         */
             DOUBLE fhard,   /* hardening modulus                       */
             DOUBLE vnu,     /* poisson's ratio                         */
             DOUBLE sig2,
             DOUBLE *tau,    /* current stresses (local)                */
             DOUBLE epstn,   /* equivalent uniaxial plastic strain      */
             DOUBLE dlam,    /* increment of plastic multiplier         */
             DOUBLE **cc)    /* material matrix to be calculated        */
{
/*----------------------------------------------------------------------*/
INT i, j, k, l;
DOUBLE alpha, expo, fobj, sq23, sm, sqsig, b, fac1, fac2, fac3;
DOUBLE fkh, g, hard;
DOUBLE rn[3][3];
DOUBLE c[3][3][3][3];
DOUBLE gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
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
#endif
/*! @} (documentation module close)*/
