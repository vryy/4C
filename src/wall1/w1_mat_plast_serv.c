/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1yilcr' which calculates the yield criterion
       for material model - Plasticity: Mises
 contains the routine 'w1radi' - radial return for element with mises
       material model
 contains the routine 'w1mapl' which forms the elasto-plastic consistent
       tangent material tensor for comined, linear isotropic and kinematic
       hardening for mises material model
 contains the routine 'w1yilcr_dp' which calculates the yield criterion
       for material model - Plasticity: Drucker-Prager
 contains the routine 'w1radi_dp' - radial return for element with 
       Drucker-Prager material model
 contains the routine 'w1maplg' which forms the elasto-plastic consistent
       tangent material tensor for Drucker-Prager material model
       (generalized for the multisurface problem: including the apex)
 contains the routine 'w1radcap' - radial return for element with 
       elastoplastic material model (concrete) - Cap region
 contains the routine 'w1mplcap' which forms the elasto-plastic consistent
       tangent material tensor for elastoplastic material model (concrete)
       - Cap region
 contains the routine 'w1cradi' - radial return for element with elasto-
       plastic material model (concrete)
 contains the routine 'w1mapl2' which forms the elasto-plastic consistent
       tangent material tensor 
 contains the routine 'w1cradms' - radial return for element with elasto-
       plastic material model (concrete) - multisurface problem
 contains the routine 'w1_mat_rebar' which calculates the constitutive
       matrix - forces - plastic - rebar
 contains the routine 'w1cdia' which calculates the element diameter
       (equivalent length)
 contains the routine 'w1cpreva' which calculates some prevalues
 contains the routine 'w1conver' which converts some arrays
 contains the routine 'w1arcs' which computes the average crack spacing
 contains the routine 'w1cpreds' which calculates the gradients of the 
       elastic predictor deviatoric/hydrostatic stresses
 and contains some mathematical routines
      

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   TOPIC : WALL1 - YIELD CRITERION PLANE_STRESS                       |
 |           -   -   -- -  --                                           |
 |           YIELD - CRITERION:  PHI > 0 DRUCKER - PRAGER               |
 |                               PHI = 0 MISES                          |
 |           FOR MATERIAL MODEL 3 - PLASTICITY                          |
 *----------------------------------------------------------------------*/
void w1yilcr(DOUBLE E, 
             DOUBLE Eh,
             DOUBLE betah,
             DOUBLE sigy,
             DOUBLE epstn,
             INT    isoft,
             DOUBLE dia,
             DOUBLE *tau,
             DOUBLE *ft)
{
/*----------------------------------------------------------------------*/
DOUBLE J2, sx, sy, sz, sxy, sigym, hards, epstmax;
/*DOUBLE betah = 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1yilcr");
#endif
/*----------------------------------------------------------------------*/
  sx  = tau[0];
  sy  = tau[1];
  sxy = tau[2];
  sz  = tau[3];
/*----------------------------------------------------------------------*/
  if(isoft==0)
  {
    hards = E * Eh / (E-Eh);
  }
  else
  {
    epstmax = (2. * Eh)/(sigy * dia);
    hards   = -(sigy*sigy * dia)/(2. * Eh); 
  }
/*----------------------------------------------------------------------*/
/*    Beruecksichtigung des Verfestigungsverhaltens */
  if(isoft==0)
  {
    sigym = sigy + betah * hards * epstn;
  }
  else
  {
    if(epstn<epstmax)  sigym = sigy + betah * hards * epstn;
    else sigym = 0.01*sigy;    
  }
/*-------------------------------------------------------- von mises ---*/
  J2 = 1./3. * (sx*sx + sy*sy + sz*sz -sx*sy - sx*sz - sy*sz) + sxy*sxy;
  *ft = sqrt(3.* J2) - sigym;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1yilcr */

/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                              modified by sh  8/02    |
 *----------------------------------------------------------------------*/
void w1radi(DOUBLE e, 
            DOUBLE eh,
            DOUBLE betah,
            DOUBLE sigy,
            DOUBLE vnu,
            DOUBLE dia,
            DOUBLE *sigma,
            DOUBLE *qn,
            INT    isoft,
            DOUBLE *epstn,
            DOUBLE *dlam,
            WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
INT i;
INT isoft1 = 0;
INT nsoft  = 1;
DOUBLE half, ro23, q13, g, c1, c2, xsi1, xsi2, xsi3, hard2, hards;
DOUBLE f, f1, f2, f3, f4, dfdl, esig, epst, fr, det, epstmax, df, dfi; 
DOUBLE J2, xsi4, stps, dlf, dlfi, dum;
DOUBLE hd11, hd21, hd12, hd22, hd33, hd44;
DOUBLE hm11, hm21, hm12, hm22, hm33;
DOUBLE dm11, dm21, dm12, dm22, dm33;
/*DOUBLE betah = 1.0;
/**/
DOUBLE fac;
DOUBLE x1,x2,x3,x4;
DOUBLE ft;
/**/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1radi");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(2./3.);
    q13  = 1./3.;
/*----------------------------------------------------------------------*/
  if(isoft==0)
  {
    hards = e * eh / (e - eh);
    hard2 = hards;
  }
  else
  {
    epstmax = (2. * eh)/(sigy * dia);
    epstmax = fabs(epstmax);
    if(*epstn<epstmax)
    {
      hards   = -(sigy*sigy * dia)/(2. * eh);
    }
    else
    {
      hards   = 0.0001;
    }
    hard2 = -(sigy*sigy * dia)/(2. * eh);
  }
 /*----------------------------------------------------------------------*/
   g = e / (vnu * 2. + 2.);

   fac= ro23 * ( 1. - betah ) * hards;
/*--------------------------------------------- initialize variables ---*/
   i = 0;
   *dlam = 0.;

/*===================================================== plane stress ===*/
  switch(wtype)
  {
  case plane_stress:
/*----------------------------------------------------------------------*/
    c1 = e / (3. - vnu * 3.) + (1. - betah) * 2. * hards / 3.;
    c2 = g * 2. + (1. - betah) * 2. * hards / 3.;

    xsi1 = sigma[0] + sigma[1];
    xsi2 = sigma[1] - sigma[0];
    xsi3 = sigma[2];
/*-----------------------------------------------------------------------|
|       iterate for new plastic multiplier                               |
|       max. number of iteration is set by imax                          |
|-----------------------------------------------------------------------*/
    f2 = xsi2 * xsi2 / 2. + xsi3 * xsi3 * 2.;
/*------------------- projection on yield surface - newton iteration ---*/

L500:
	++i;

	f1 = xsi1 / (c1 * *dlam + 1.);
	f3 = c2 * *dlam + 1.;
	f = f1 * f1 / 6. + f2 / (f3 * f3);
	fr = sqrt(f);

	epst = *epstn + *dlam * fr * ro23;

/*------------ isotropic hardening for Mises yield criteria ------------*/
/*--------------- HARD2 (Ausgangswert) statt HARDS (evtl. upgedatet) ---*/
	if (isoft1 == 0) {
	    if (nsoft == 1) {
		esig = sigy + betah * hard2 * epst;
	    } else {
		esig = sigy * exp(-epst / epstmax);
	    }
	} else {
	    if (nsoft == 1) {
		esig = sigy - sigy * (epst / epstmax);
	    } else {
		esig = sigy * exp(-epst / epstmax);
	    }
	}
/*------------------------------- minimale Werte fuer ESIG und HARDS ---*/
	if (esig < sigy * .01) {
	    esig = sigy * .01;
	    hards = 1e-4;
	}
/*-------------------------------------------------- yield criteria ----*/
	f4 = esig * esig / 3.;
	f = f / 2. - f4;
/*------------------------------------ derivative of yield criteria ----*/
	dfdl = -c1 * f1 * f1 / ((c1 * *dlam + 1.) * 6.) - c2 * f2 / 
              (f3 * f3 * f3);
	f4 = esig * (betah * hards) * ro23 * 2. / 3.;
	dfdl = (fr - f4 * *dlam) * dfdl / fr - f4 * fr;

	*dlam -= f / dfdl;
/*------------------------------------------------- Abbruchkriterium ---*/
	if ((dum = f / esig, fabs(dum)) > 1e-5) {
	    if (i > 100) dserror("i>100!");         
	    goto L500;
	} else {
	    *epstn = epst;
	    df = *dlam / ((1. - betah) * 2.*hards * *dlam / 3.+ 1.) / 3.;
	}
/*-------- modified constitutive tensor (compliance) HM = H + DF * P ---*/
	hm11 = 1. / e + df * 2.;
	hm21 = -vnu / e - df;
	hm12 = -vnu / e - df;
	hm22 = 1. / e + df * 2.;
	hm33 = 1. / g + df * 6.;
/*------------------------------ inverse (material stiffness tensor) ---*/
	det = hm11 * hm22 - hm12 * hm21;

	dm11 = hm22 / det;
	dm21 = -hm21 / det;
	dm12 = -hm12 / det;
	dm22 = hm11 / det;
	dm33 = 1. / hm33;
/*-------------------------- new stress state STRES = DM * H * SIGMA ---*/
	dfi = 1. / ((1. - betah) * 2. * hards * *dlam / 3. + 1.);

	hm11 = 1. / e;
	hm21 = -vnu / e;
	hm12 = -vnu / e;
	hm22 = 1. / e;
	hm33 = 1. / g;

	f1 = hm11 * sigma[0] + hm12 * sigma[1];
	f2 = hm21 * sigma[0] + hm22 * sigma[1];
	f3 = hm33 * sigma[2];

	sigma[0] = dfi * (dm11 * f1 + dm12 * f2);
	sigma[1] = dfi * (dm21 * f1 + dm22 * f2);
	sigma[2] = dfi * (dm33 * f3);

	for (i = 0; i < 3; i++) {
	    qn[i] += sigma[i] * 2. * (1. - betah) * hards * *dlam / 3.;
	}
/*----------------------------------------------------------------------*/
  break;
  case plane_strain:

/*===================================================== plane strain ===*/
	xsi1 = sigma[0];
	xsi2 = sigma[1];
	xsi3 = sigma[3]; /*!!!*/
	xsi4 = sigma[2];
/*---------------------------------- increment of plastic multiplier ---*/
   J2 = 1./3. * (xsi1*xsi1 + xsi2*xsi2 + xsi3*xsi3 - 
                 xsi1*xsi2 - xsi1*xsi3 - xsi2*xsi3) + xsi4*xsi4;
   stps = sqrt(2.* J2);    /*norm of devaiatoric praedictor stresses*/     

/*-------- direct evaluation of dlam  see Simo & Hughes p.120f ---------*/
/*------------- only possible for linear hardening ---------------------*/

   esig = ro23 * (sigy + betah * hards * *epstn);
   ft = sqrt (2.*J2) - esig;
   *dlam = ft/(2.*g + (2./3.) * hards);

/*-------- update of uniaxial plastic strain, stresses & backstress ----*/
   *epstn = *epstn + ro23 * *dlam;

/*-------- components of gradient n ------------------------------------*/
   x1= (2. * xsi1 - xsi2 - xsi3 ) / 3. / stps;
   x2= (2. * xsi2 - xsi1 - xsi3 ) / 3. / stps;
   x3= (2. * xsi3 - xsi1 - xsi2 ) / 3. / stps;
   x4= xsi4 / stps;

/*-- see equation (3.3.3) in Simo & Hughes on p.121 --------------------*/
   sigma[0] = sigma[0] - (2.*g + ro23*fac)* *dlam*x1;   
   sigma[1] = sigma[1] - (2.*g + ro23*fac)* *dlam*x2;   
   sigma[3] = sigma[3] - (2.*g + ro23*fac)* *dlam*x3;   
   sigma[2] = sigma[2] - (2.*g + ro23*fac)* *dlam*x4;   
         
/*-- backstress vector: qn = qn_trial + fac2 * dlam * grad (3.3.1) -----*/
   qn[0]  = qn[0] + ro23*fac * *dlam * x1;
   qn[1]  = qn[1] + ro23*fac * *dlam * x2;
   qn[3]  = qn[3] + ro23*fac * *dlam * x3;
   qn[2]  = qn[2] + ro23*fac * *dlam * x4;
   
/*----------------------------------------------------------------------*/
  break;
/*============================================== rotational symmetry ===*/
  default:
    dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1radi */

/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 |   for combined, linear isotropic and kinematic hardening             |
 |                                                modified sh 8/02      |
 |   see Simo & Huges: plane_strain: S.124                              |
 |                     plane_stress: S.130                              |
 *----------------------------------------------------------------------*/
void w1mapl(DOUBLE e, 
            DOUBLE eh,
            DOUBLE betah,
            DOUBLE sigy,
            DOUBLE vnu,
            DOUBLE dia,
            DOUBLE *tau,       /* !! predictor stresses if plane_strain */
            INT    isoft,
            DOUBLE *epstn,
            DOUBLE *dlam,
            DOUBLE **d,
            WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
INT i, j;
INT nsoft  = 1;
DOUBLE gamma1, gamma2, beta,x1,x2,x3,x4,abeta,fact,fact1;
DOUBLE vect[4];
DOUBLE dm[4][4];
DOUBLE xsi1, xsi2, xsi3, xsi4, hards, g, epstmax, df, det, dum;
DOUBLE d11, d21, d12, d22, d33;
DOUBLE hm11, hm21, hm12, hm22, hm33, hm44;
DOUBLE dm11, dm21, dm12, dm22, dm33;
/*DOUBLE betah = 1.;*/

DOUBLE k,fac,fac1,fac2;
DOUBLE stps,J2;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1mapl");
#endif
/*----------------------------------------------------------------------*/
    if (isoft == 0) {
	hards = e * eh / (e - eh);
    } else {
	if (nsoft == 1) {
	    epstmax = eh * 2. / (sigy * dia);
	    epstmax = fabs(epstmax);
	    if (*epstn < epstmax) {
                /* Computing 2nd power */
		dum = sigy;
		hards = -(dum * dum * dia) / (eh * 2.);
	    } else {
		hards = 1e-4;
	    }
	} else if (nsoft == 2) {
	    epstmax = eh * 1. / (sigy * dia);
	    epstmax = abs(epstmax);
	    hards = -sigy * exp(-(*epstn) / epstmax) / epstmax;
	} else {
            dserror("WRONG NUMBER OF NSOFT");
	}
    }
    g = e / (vnu * 2. + 2.);           /*Schubmodul*/
    k = e / ((1. - vnu * 2.) * 3.);    /*Kompressionsmodul*/

/*===================================================== plane stress ===*/
  switch(wtype)
  {
  case plane_stress:
/*----------------------------------------------------------------------*/
	gamma1 = (1. - betah) * 2. * hards * *dlam / 3. + 1.;
	gamma2 = 1. - betah * 2. * hards * *dlam / 3.;

	xsi1 = tau[0];
	xsi2 = tau[1];
	xsi3 = tau[2];

	beta = (betah * hards * gamma1 + (1. - betah) * hards * gamma2) *
               2. * gamma1 / (gamma2 * 3.);
	beta = beta * 2. * (xsi1 * xsi1 + xsi2 * xsi2 - xsi1 * 
               xsi2 + xsi3 * xsi3 * 3.) / 3.;

	df = *dlam / ((1. - betah) * 2. * hards * *dlam / 3. + 1.) / 3.;
/*-----------------------------------------------------------------------|
|        inverse modified constitutive tensor (compliance)               |
|                 D  =  D(INV) + P * DLAM                                |
|-----------------------------------------------------------------------*/
	d11 = 1. / e + df * 2.;
	d21 = -vnu / e - df;
	d12 = -vnu / e - df;
	d22 = 1. / e + df * 2.;
	d33 = 1. / g + df * 6.;
/*-----------------------------------------------------------------------|
|        inverse material stiffness tensor     DM = INV (D(INV) + P*DLA  |
|-----------------------------------------------------------------------*/
	det = d11 * d22 - d12 * d21;

	dm11 =  d22 / det;
	dm21 = -d21 / det;
	dm12 = -d12 / det;
	dm22 =  d11 / det;
	dm33 =   1. / d33;
/*-----------------------------------------------------------------------|
|        product  P * TAU                                                |
|-----------------------------------------------------------------------*/
	x1 = (xsi1 * 2. - xsi2) / 3.;
	x2 = (-xsi1 + xsi2 * 2.) / 3.;
	x3 = xsi3 * 2.;
/*-----------------------------------------------------------------------|
|        product  DM * P * TAU                                           |
|-----------------------------------------------------------------------*/
	vect[0] = dm11 * x1 + dm12 * x2;
	vect[1] = dm21 * x1 + dm22 * x2;
	vect[2] = dm33 * x3;
/*-----------------------------------------------------------------------|
|        scalar product  TAU * P * DM * P * TAU                          |
|-----------------------------------------------------------------------*/
	abeta = x1 * vect[0] + x2 * vect[1] + x3 * vect[2];
	abeta = 1. / (abeta + beta);
/*-----------------------------------------------------------------------|
|        form consistent 4x4 material matrix                             |
|-----------------------------------------------------------------------*/
        d[0][0]=dm11;
        d[0][1]=dm12;
        d[0][2]=0.0;
        d[1][0]=dm21;
        d[1][1]=dm22;
        d[1][2]=0.0;
        d[2][0]=0.0;
        d[2][1]=0.0;
        d[2][2]=dm33;

	for (i = 0; i < 3; ++i) {
	    for (j = i; j <3; ++j) {
               d[i][j] -= abeta * vect[i] * vect[j];
		if (j > i)  d[j][i] = d[i][j];
	    }
	}
/*----------------------------------------------------------------------*/
  break;
  case plane_strain:

/*===================================================== plane strain ===*/
   xsi1 = tau[0];
   xsi2 = tau[1];
   xsi3 = tau[3];
   xsi4 = tau[2];
/*------- Grad n  -> Normiert  -----------------------------------------*/
   J2 = 1./3. * (xsi1*xsi1 + xsi2*xsi2 + xsi3*xsi3 - 
                 xsi1*xsi2 - xsi1*xsi3 - xsi2*xsi3) + xsi4*xsi4;
   stps = sqrt(2.* J2); /*norm of devaiatoric predictor stresses*/     

   x1= (2. * xsi1 - xsi2 - xsi3 ) / 3. / stps;
   x2= (2. * xsi2 - xsi1 - xsi3 ) / 3. / stps;
   x3= (2. * xsi2 - xsi1 - xsi2 ) / 3. / stps;
   x4= xsi4 / stps;
/*------- Vorfaktoren  -------------------------------------------------*/
   fac  = (2.*g* *dlam)/ (stps);
   fac1 = 1. - fac;
   fac2 = 1./(1.+ hards/(3.*g)) - fac ;      
/*-------- Cep_alg nach Simo & Hughes  ---------------------------------*/

   dm[0][0] = k + 2.*g*fac1*( 2./3.) - 2.*g*fac2*x1*x1;
   dm[1][1] = k + 2.*g*fac1*( 2./3.) - 2.*g*fac2*x2*x2;
   dm[2][2] = k + 2.*g*fac1*( 2./3.) - 2.*g*fac2*x3*x3;
   dm[3][3] =     2.*g*fac1*( 1./2 ) - 2.*g*fac2*x4*x4;
           
   dm[0][1] = k + 2.*g*fac1*(-1./3.) - 2.*g*fac2*x1*x2;
   dm[1][0] = dm[0][1];
   dm[0][2] = k + 2.*g*fac1*(-1./3.) - 2.*g*fac2*x1*x3;
   dm[2][0] = dm[0][2];
   dm[1][2] = k + 2.*g*fac1*(-1./3.) - 2.*g*fac2*x2*x3;
   dm[2][1] = dm[1][2];
     
   dm[0][3] =                        - 2.*g*fac2*x1*x4;
   dm[3][0] = dm[0][3];
   dm[1][3] =                        - 2.*g*fac2*x2*x4;
   dm[3][1] = dm[1][3];
   dm[2][3] =                        - 2.*g*fac2*x3*x4;
   dm[3][2] = dm[2][3];
   
   /*Sorting back rows/columns 2 <-> 3 */
   d[0][0] = dm[0][0];
   d[1][0] = dm[1][0]; 
   d[2][0] = dm[3][0];
   d[3][0] = dm[2][0];
   d[0][1] = dm[0][1];
   d[1][1] = dm[1][1];
   d[2][1] = dm[3][1];
   d[3][1] = dm[2][1];
   d[0][2] = dm[0][3];
   d[1][2] = dm[1][3];
   d[2][2] = dm[3][3]; 
   d[3][2] = dm[2][3];
   d[0][3] = dm[0][2];
   d[1][3] = dm[1][2];
   d[2][3] = dm[3][2];
   d[3][3] = dm[2][2];              
/*----------------------------------------------------------------------*/
  break;
/*============================================== rotational symmetry ===*/
  default:
    dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1mapl */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   TOPIC : WALL1 - YIELD CRITERION PLANE_STRESS                       |
 |           -   -   -- -  --                                           |
 |           YIELD - CRITERION:  PHI > 0 DRUCKER - PRAGER               |
 |                               PHI = 0 MISES                          |
 *----------------------------------------------------------------------*/
void w1yilcr_dp(DOUBLE E, 
                DOUBLE Eh,
                DOUBLE phi,
                DOUBLE sigy,
                DOUBLE *sigym,
                DOUBLE epstn,
                DOUBLE *tau,
                DOUBLE *ft)
{
/*----------------------------------------------------------------------*/
DOUBLE coh, sx, sy, sxy, hards, alpha, yld;
DOUBLE betah = 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1yilcr_dp");
#endif
/*----------------------------------------------------------------------*/
  sx  = tau[0];
  sy  = tau[1];
  sxy = tau[2];
/*----------------------------------------------------------------------*/
  coh = sigy / 2.;
/*----------------------------------------------------------------------*/
  hards = E * Eh / (E-Eh);
/*----------------------------------------------------------------------*/
/*     Vergleichsspannung Drucker - Prager */
  yld = 6. * coh * cos(phi) / (3. - sin(phi));
/*----------------------------------------------------------------------*/
/*     Beruecksichtigung des Verfestigungsverhaltens */
  *sigym = yld + betah * hards * epstn;
/*----------------------------------------------------------------------*/
/*     Vorwert bei I1 */
  alpha = 1. / 3. * 6. * sin(phi) / (3. - sin(phi));
/*----------------------------------------------------------------------*/
/*     Fliessbedingung  Drucker - Prager */
  *ft = sqrt(sx*sx-sx*sy+sy*sy+3.0*(sxy*sxy)) + alpha*(sx+sy)- *sigym;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1yilcr_dp */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1radi_dp(DOUBLE e, 
               DOUBLE eh,
               DOUBLE phi,
               DOUBLE sigy,
               DOUBLE vnu,
               DOUBLE *sigma,
               DOUBLE *qn,
               DOUBLE *epstn,
               DOUBLE *sigym,
               DOUBLE *dlam,
               WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
INT i;
INT isoft1 = 0;
INT nsoft  = 1;
DOUBLE half, ro23, q13, g, c1, c2, g1, xsi1, xsi2, xsi3, hard2, hards;
DOUBLE f, f1, f2, f3, f4, dfdl, esig, epst, fr, det, epstmax, df, dfi; 
DOUBLE dum, xsi4, stps, dlf, dlfi, coh, y, alpha, devinv, alph, fac;
DOUBLE hd11, hd21, hd12, hd22, hd33, hd44;
DOUBLE hm11, hm21, hm12, hm22, hm33;
DOUBLE dm11, dm21, dm12, dm22, dm33;
DOUBLE betah = 1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1radi_dp");
#endif
/*----------------------------------------------------------------------*/
  half = .5;
  ro23 = sqrt(.66666666666666663);
  q13 = .33333333333333331;
/*----------------------------------------------------------------------*/
  hards = e * eh / (e - eh);
  hard2 = hards;
  g     = e / (vnu * 2. + 2.);
  coh   = sigy/2.;  
/*----------------------------------------------------------------------*/
/*     MAXIMALE HYDROSTATISCHER SPANNUNGSZUSTAND */
  if(phi!=0.) y=coh/tan(phi); 
/*===================================================== plane stress ===*/
  switch(wtype)
  {
  case plane_stress:
/*----------------------------------------------------------------------*/
    c1 = e / (3. - vnu * 3.) + (1. - betah) * 2. * hards / 3.;
    c2 = g * 2. + (1. - betah) * 2. * hards / 3.;
    g1 = g + (1. - betah) * hards / 3.;

    xsi1 = sigma[0] + sigma[1];
    xsi2 = sigma[1] - sigma[0];
    xsi3 = sigma[2];
/*-----------------------------------------------------------------------|
|       iterate for new plastic multiplier                               |
|       max. number of iteration is set by imax                          |
|-----------------------------------------------------------------------*/
    f2 = xsi2 * xsi2 / 2. + xsi3 * xsi3 * 2.;
/*--------------------------------------------- initialize variables ---*/
    i = 0;
    *dlam = 0.;
/*------------------- projection on yield surface - newton iteration ---*/

L500:
	++i;

	f1 = xsi1 / (c1 * *dlam + 1.);
	f3 = c2 * *dlam + 1.;
/*----------------------------------------------------------------------*/
/*       Vorwert bei I1  (1. Invariante) */
        alpha = 1. / 3. * 6. * sin(phi) / (3. - sin(phi)); 
/*----------------------------------------------------------------------*/
/*       Wurzel 3J2  (Deviatorinvariante) */
        devinv =  sqrt(3./2.*(f1*f1/6. + f2/(f3*f3))); 
/*----------------------------------------------------------------------*/
/*       Vergleichsspannung */
        f  =  devinv + alpha * f1 ;  
        fr = sqrt(f);
/*----------------------------------------------------------------------*/
/*       Beiwert fuer plastische Verzerrung */
          alph  =  6. * sin(phi) / (3. - sin(phi));
/*----------------------------------------------------------------------*/
/*       Vorfaktor der plastischen Verzerrung */
          fac = (alph + 1. / sqrt(3.)) / sqrt(3. * alph * alph + .5);
/*----------------------------------------------------------------------*/
	  epst = *epstn + *dlam * fr * fac;
/*----------------------------------------------------------------------*/
/*       isotrope Materialverfestigung DRUCKER PRAGER */    
          esig = *sigym  + betah*hards*epst; 
/*----------------------------------------------------------------------*/
/*       Fliessbedingung */
          f4 = esig;
          f  = f - f4;
/*----------------------------------------------------------------------*/
/*       Ableitung der Fliessbedingung */
	dfdl = 3./2.*(-c1*f1*f1/(6.*(1.+c1* *dlam)) - c2*f2/(f3*f3*f3)) /
                 devinv - alpha*c1*xsi1/(1+c1* *dlam) * (1+c1* *dlam);

	f4 = esig * (betah * hards) * ro23 * 2. / 3.;
	dfdl = (fr - f4 * *dlam) * dfdl / fr - f4 * fr;

	*dlam -= f / dfdl;
/*------------------------------------------------- Abbruchkriterium ---*/
/*       Abbruchkriterium */
          if (esig>(0.001* *sigym)) 
          {
            if (fabs(f/esig)>0.0001 ) 
            {
	      if (i > 30) dserror("i>30!");         
	      goto L500;
            }
            else
            {
              *epstn = epst;
	      df = *dlam /(1.+2.*(1. - betah) * hards* *dlam/3)/3.;
            }  
          }
          else
          {
            sigy = 0.001 * *sigym; 
            *epstn = epst;
            df = *dlam /(1.+2.*(1. - betah) * hards* *dlam/3)/3.;  
          }
/*-------- modified constitutive tensor (compliance) HM = H + DF * P ---*/
	hm11 = 1. / e + df * 2.;
	hm21 = -vnu / e - df;
	hm12 = -vnu / e - df;
	hm22 = 1. / e + df * 2.;
	hm33 = 1. / g + df * 6.;
/*------------------------------ inverse (material stiffness tensor) ---*/
	det = hm11 * hm22 - hm12 * hm21;

	dm11 = hm22 / det;
	dm21 = -hm21 / det;
	dm12 = -hm12 / det;
	dm22 = hm11 / det;
	dm33 = 1. / hm33;
/*-------------------------- new stress state STRES = DM * H * SIGMA ---*/
	hm11 = 1. / e;
	hm21 = -vnu / e;
	hm12 = -vnu / e;
	hm22 = 1. / e;
	hm33 = 1. / g;

	f1 = hm11 * sigma[0] + hm12 * sigma[1];
	f2 = hm21 * sigma[0] + hm22 * sigma[1];
	f3 = hm33 * sigma[2];

	sigma[0] = dm11 * f1 + dm12 * f2;
	sigma[1] = dm21 * f1 + dm22 * f2;
	sigma[2] = dm33 * f3;

	for (i = 0; i < 3; i++) {
	    qn[i] += sigma[i] * 2. * (1. - betah) * hards * *dlam / 3.;
	}
/*----------------------------------------------------------------------*/
  break;
  case plane_strain:
/*===================================================== plane strain ===*/
	xsi1 = sigma[0];
	xsi2 = sigma[1];
	xsi3 = sigma[3]; /*!!!*/
	xsi4 = sigma[2];
/*---------------------------------- increment of plastic multiplier ---*/
	dum = q13 * (xsi1 * (xsi1 * 2. - xsi2 - xsi3) + xsi2 * (xsi2 * 2. - 
		xsi1 - xsi3) + xsi3 * (xsi3 * 2. - xsi1 - xsi2)) + xsi4 * 2. *
		 xsi4;
	stps = pow(dum, half);
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
	dlf = *dlam * 2. * (g + (1. - betah) * hards / 3.) + 1.;
/*--------------------------------- new plastic uniaxial deformation ---*/
	epst = *epstn + *dlam / dlf * ro23 * stps;
	esig = sigy + betah * hards * epst;
/*----------------------------------- apply von mises yield criteria ---*/
	f = 1. / (dlf * 2. * dlf) * stps * stps - q13 * esig * esig;
/*--------- derivative of the yield criteria with respect to plastic ---*/
/*---------                                                increment ---*/
	dfdl = -stps / (dlf * dlf) * (stps / dlf * 2. * (g + (1.-betah) * 
		hards / 3.)+pow(ro23, 3.) * esig * betah * hards);
/*-------------------------------------------- new plastic increment ---*/
	*dlam -= f / dfdl;
/*------------------------------------------------ check convergence ---*/
	if (esig == 0.) {
	    if (fabs(f) > 1e-5) {
		if (i > 30) {
		    dserror("local iteration exceeds limit");
		}
		goto L600;
	    }
	} else {
	    if ((dum = f / esig, fabs(dum)) > 1e-5) {
		if (i > 30) {
		    dserror("local iteration exceeds limit");
		}
		goto L600;
	    }
	}
	*epstn = epst;
/*-----------------------------------------------------------------------|
|        NEW STRESS STATE  SIGMA = DM * H * SIGMA(PRED)                  |  
|        WITH MODIFIED CONSTITUTIVE TENSOR (COMPLIANCE)                  |
|        HM = H + DLAM * P                                               |
|        HD = INV[H]*INV[D]                                              |
|        IN 3D-CONTINUUM: HD11=HD22=HD33; HD12=HD13=HD21=HD23=HD31=HD32  |
|-----------------------------------------------------------------------*/
	dlf = *dlam * 2. * (g + (1. - betah) * hards / 3.) + 1.;
	dlfi = 1. / dlf;

	hd11 = q13 * (dlfi * 2. + 1.);
	hd12 = q13 * (1. - dlfi);
	hd44 = dlfi;

	sigma[0] = hd11 * xsi1 + hd12 * xsi2 + hd12 * xsi3;
	sigma[1] = hd12 * xsi1 + hd11 * xsi2 + hd12 * xsi3;
	sigma[3] = hd12 * xsi1 + hd12 * xsi2 + hd11 * xsi3;
	sigma[2] = hd44 * xsi4;

	for (i = 0; i < 3; ++i) {
	    qn[i] += sigma[i] * 2. * (1. - betah) * hards * *dlam / 3.;
	}
/*----------------------------------------------------------------------*/
  break;
/*============================================== rotational symmetry ===*/
  default:
    dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1radi_dp */
/*----------------------------------------------------------------------*
 |                                                       al    9/01     |
 |    calculate element diameter (equivalent length) for one element    |	     
 *----------------------------------------------------------------------*/
void w1cdia(ELEMENT   *ele, 
            W1_DATA   *data,
            DOUBLE    *funct_h,
            DOUBLE   **deriv_h,
            DOUBLE   **xjm_h)
{
INT                 i,j,k;            /* some loopers */
INT                 nir,nis;          /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
const INT           numdf  = 2;
const INT           numeps = 3;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE              det, exp, dia;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1cdia");
#endif
/*------------------------------------------- integration parameters ---*/
w1intg(ele,data,1);
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
iel     = ele->numnp;
nd      = numdf * iel;
/*================================================ integration loops ===*/
fac = 0.;
for (lr=0; lr<nir; lr++)
{
   /*=============================== gaussian point and weight at it ===*/
   e1   = data->xgrr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      /*============================ gaussian point and weight at it ===*/
      e2   = data->xgss[ls];
      facs = data->wgts[ls];
      /*------------------------- shape functions and their derivatives */
      w1_funct_deriv(funct_h,deriv_h,e1,e2,ele->distyp,1);
      /*------------------------------------ compute jacobian matrix ---*/       
      w1_jaco (funct_h,deriv_h,xjm_h,&det,ele,iel);                         
      fac += facr * facs * det; 
      /*----------------------------------------------------------------*/
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */

/*------------------------------- mofification factor exp = alpha**2 ---*/
  if(iel==4)       exp = 2.;
  else if(iel==8)  exp = 1.;
  else if(iel==9)  exp = 1.;
  else             exp = 2./3.;
  
  dia = sqrt(exp * fac);
  
  ele->e.w1->elewa[0].dia = dia;
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1cdia */
/*-----------------------------------------------------------------------|
|     topic : wall1 - concrete material prevalues           al  2/02     |
|             -   -   -                 -----                            |
|             yield - criterion:  drucker-prager / spherical cap         |
|-----------------------------------------------------------------------*/
void w1cpreva (DOUBLE *epst,     /* equivalent uniaxial plastic strain  */
               DOUBLE *sigym,    /* [0] uniaxial tension yield stress   */
                                 /* [1] yield stress "inverted cone"    */
                                 /* [2] uniaxial compr. yield stress    */
               DOUBLE *alpha,    /* factor for the first invariants     */
	       DOUBLE *hards,    /* hardening modulus                   */
               DOUBLE *e,        /* young modulus                       */
               DOUBLE *g,        /* shear modulus                       */
               DOUBLE *vnu,      /* poisson's ratio                     */
               DOUBLE *com,      /* bilk modulus                        */
               DOUBLE *fcm,      /* compressive strenght                */
               DOUBLE *gc,       /* compression fracture energy         */
               DOUBLE *ftm,      /* tensile strenght                    */
               DOUBLE *gt,       /* tensile fracture energy             */
               DOUBLE *gamma1,   /* fitting factor yield function 1     */
               DOUBLE *gamma2,   /* fitting factor yield function 2     */
               DOUBLE *dfac,     /* damage factor                       */
               DOUBLE *dia,      /* equivalent element length           */
	       DOUBLE *acrs,     /* average crack spacing               */
               DOUBLE *cappaet,  /* max. elastic tension strain         */
               DOUBLE *cappaut,  /* tensile fracture strain             */
               DOUBLE *cappae,   /* max. elastic compression strain     */
               DOUBLE *cappauc,  /* compressive fracture strain         */
               DOUBLE *sig3,     /* equivalent compressive stress       */
               DOUBLE *fbd,      /* tension stiffening stress           */
               DOUBLE **d,       /* elastic material matrix             */
               DOUBLE *sig,      /* stresses from last update           */
               DOUBLE *eps)      /* strains from last update            */
{                                
/*-------------------------------------------------- local variables ---*/
    static INT i, j;
    static DOUBLE dum;
    static DOUBLE vect[4], dfac1, dfac2, epsn1, epsn2;
    static DOUBLE sigym1, sigym2, q23, fac, dam, ro23, ro54, dia2, sig1;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("w1cpreva");
  #endif
/*----------------------------------------------------------------------*/
    ro23 = sqrt(.66666666666666663);
    ro54 = sqrt(54.);
    q23 = .66666666666666663;
/*------------------------------------------------------- initialize ---*/
    sigym1 = 0.;
    sigym2 = 0.;
    *dfac = 1.;

    for (i = 0; i < 4; i++) 
    {
	vect[i ] = 0.;
	sigym[i] = 0.;
	alpha[i] = 0.;
	hards[i] = 0.;
    }
/*----------------------------------------------------- shear/compr. ---*/
    *g = *e / (*vnu * 2. + 2.);
    *com = *e / ((1. - *vnu * 2.) * 3.);
/*----------------------------------------------------------------------*/
/* Objektivitaet bei grosser equivalenter Elementlaenge gewaehrleisten  */ 
    dia2 = *gt * *e / pow(*ftm,2.);
    if (*dia > dia2) 
    {
	*ftm = sqrt(*gt * *e / *dia);
    }

    dia2 = *gc * 3. * *e / (pow(*fcm, 2.) * 4.);
    if (*dia > dia2) 
    {
	*fcm = sqrt(*gc * 3. * *e / (*dia * 4.));
    }
/*----------------------------------------------------------------------*/
/* Vorwerte fuer: Zug-/Zug-Druck-Bereiche */
    alpha[0] = ro23 * (*gamma1 * *fcm - *ftm) / (*gamma1 * *fcm + *ftm);
    *cappaet = *ftm / *e;
    *cappaut = *gt / (*acrs * *ftm);
/*----------------------------------------------------------------------*/
/* FBD = 1.05 * FTM */
    *fbd = *ftm * .9;
    sig1 = *gamma1 * 2. * *fcm * *ftm / (*gamma1 * *fcm + *ftm);
    sigym[0] = sig1 * exp(-epst[0] / *cappaut);

    if (sigym[0] > sig1 * .01) 
    {
	hards[0] = -sigym[0] / *cappaut;
    }
    else 
    {
	sigym[0] = sig1 * .01;
	hards[0] = 1e-4;
    }
/*----------------------------------------------------------------------*/
/* Vorwerte fuer: APEX */
    alpha[1] = -1. / (alpha[0] * 3.);
    sigym[1] = -sigym[0] / (alpha[0] * alpha[0] * 3.);
    hards[1] = -hards[0] / (alpha[0] * alpha[0] * 3.);
/*----------------------------------------------------------------------*/
/* Vorwerte fuer: Druck-Bereich */
    alpha[2] = ro23 * (*gamma2 - 1.) / (*gamma2 * 2. - 1.);
    if (alpha[2] == 0.) 
    {
	alpha[2] = 1e-4;
    }
    *cappae = *fcm * 4. / (*e * 3.);
    *cappauc = *cappae + *gc * 3. / (*dia * 2. * *fcm);
    *sig3 = *fcm * *gamma2 / (*gamma2 * 2. - 1.);
    if (epst[1] < *cappae) 
    {
	fac = epst[1] / *cappae;
	*sig3 /= 3.;
	sigym[2] = *sig3 * (fac * 4. + 1. - pow(fac, 2.) * 2.);
	hards[2] = *sig3 * 4. * (1. - fac) / *cappae;
    } 
    else if (epst[1] >= *cappae && epst[1] < *cappauc) 
    {
	fac = epst[1] - *cappae;
	fac /= *cappauc - *cappae;
	sigym[2] = *sig3 * (1. - pow(fac, 2.));
	dum = *cappauc - *cappae;
	hards[2] = *sig3 * -2. * (epst[1] - *cappae) / (dum * dum);
    } 
    else if (epst[1] >= *cappauc) 
    {
	sigym[2] = 1e-4;
	hards[2] = 1e-4;
    }

    if (sigym[2] < *sig3 * .03) 
    {
	sigym[2] = *sig3 * .03;
	hards[2] = 1e-4;
    }

/*----------------------------------------------------------------------*/
/* Vorwerte fuer: CAP-Bereich */
    alpha[3] = -(ro54 * alpha[2] + 2.) * *gamma2;
    sigym[3] = sqrt(q23 + alpha[2] * alpha[2] * 6.) * *gamma2 * sigym[2];
    hards[3] = hards[2] * sigym[3] / sigym[2];
/*----------------------------------------------------------------------*/
/* Schaedigungsfaktor (ingenieurmaessiger Ansatz) */
    if (dam > 0. && dam <= 1. && (epst[0] > 0. || epst[1] > 0.)) 
    {
	dum = ro23 + alpha[0];
	epsn1 = pow(dum, 2.) * epst[0] / ro23;
	dum = ro23 - alpha[2];
	epsn2 = pow(dum, 2.) * epst[1] / ro23;
	dfac1 = sigym[0] / (sigym[0] + dam * *e * epsn1);
	dfac2 = sigym[2] / (sigym[2] + dam * *e * epsn2);
	if (epst[0] > 0. && epst[1] <= 0.) 
        {
	    *dfac = dfac1;
	} 
        else if (epst[0] <= 0. && epst[1] > 0.) 
        {
	    *dfac = dfac2;
	} 
        else if (epst[0] > 0. && epst[1] > 0.) 
        {
	    *dfac = (dfac1 * epsn1 + dfac2 * epsn2) / (epsn1 + epsn2);
	} 
        else {
	    *dfac = 1.;
	}
/*----------------------------------------------------------------------*/
/* Schaedigungsfaktor (mit "Freier Energie nach Helmholtz") */
/*                SIG : EPS    */
/*     DFAC =  --------------- */
/* 	       EPS : D : EPS   */
    } 
    else if (dam > 1. && dam <= 2. && (epst[0] > 0. || epst[1] > 0.)) 
    {
	dfac1 = 0.;
	dfac2 = 0.;
/*----------------------------------------------------------------------*/
/* VECT  = DM : DEPS */
	for (j = 0; j < 4; j++) 
        {
	    for (i = 0; i < 4; i++) 
            {
		vect[j] += d[j][i] * eps[i];
	    }
	}
/*----------------------------------------------------------------------*/
/* DFAC1 = SIG  : EPS */
/* DFAC2 = VECT : EPS */
	eps[2] *= 2.;
	for (i = 0; i < 4; i++) 
        {
	    dfac1 += sig[i] * eps[i];
	    dfac2 += vect[i] * eps[i];
	}
	if (dfac2 != 0.) 
        {
	    *dfac = dfac1 / dfac2;
	    *dfac = (*dfac - 1.) * (dam - 1.) + 1.;
	}
	eps[2] /= 2.;
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ;
} /* end of w1cpreva */
/*----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------|
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)                                               |
|-----------------------------------------------------------------------*/
void w1cradi (DOUBLE *sigma,  /* elastic predictor projected onto yield surface */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain             */ 
              DOUBLE *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              INT    yip,     /* stress state   1 =elastic   >=2 =plastic       */ 
              DOUBLE *alpha,  /* factor for the first invariants                */ 
              DOUBLE *ft,     /* yield condition                                */ 
              DOUBLE *e,      /* young modulus                                  */ 
              DOUBLE *g,      /* shear modulus                                  */ 
	      DOUBLE *com,    /* bulk modulus                                   */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress                */ 
              DOUBLE *hards,  /* plastic modulus                                */ 
              DOUBLE *sigy,   /* actual uniaxial yield stress                   */ 
              DOUBLE *dn,     /* gradient components in deviatoric direction    */ 
              DOUBLE *dcom,   /* gradient components in hdrostatic direction    */ 
              DOUBLE *grad,   /* total gradient                                 */ 
              DOUBLE *devsig, /* deviatoric predictor stresses                  */ 
              DOUBLE *sm,     /* hydrostatic predictor stresses                 */ 
              DOUBLE *fcm,    /* compressive strenght                           */ 
              DOUBLE *gc,     /* compression fracture energy                    */ 
              DOUBLE *ftm,    /* tensile strenght                               */ 
	      DOUBLE *gt,     /* tensile fracture energy                        */ 
              DOUBLE *gamma1, /* fitting factor yield function 1                */ 
              DOUBLE *gamma2, /* fitting factor yield function 2                */ 
              DOUBLE *dia,    /* equivalent element length                      */ 
              DOUBLE *acrs)   /* average crack spacing                          */ 
{                             
/*----------------------------------------------------------------------*/
    DOUBLE dum;               
/*----------------------------------------------------------------------*/
    static DOUBLE dfdl, half, epst;
    static DOUBLE f;
    static INT i;
    static DOUBLE q13, dalpha, cappae, fac, ro23, smi[4], cappauc, 
	    devsigi[4], cappaut, sig1, sig3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1cradi");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
/*   Beiwert fuer plastische Verzerrung */
    dalpha = *alpha * *alpha;
/*--------------------------------------------- initialize variables ---*/
    i = 0;
    *dlam = 0.;

/*------------------------------------- plane strain / plane stress  ---*/
    if (wtype==plane_stress || wtype==plane_strain) 
    {
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
        epst = *epstn + *dlam;
/*------------------------------------------------- new yield stress ---*/
/*     Zug-/Zug-Druck-Bereiche */
	if (yip == 2) {
	    cappaut = *gt / (*acrs * *ftm);
	    sig1 = *gamma1 * 2. * *fcm * *ftm / (*gamma1 * *fcm + *ftm);
	    *sigy = sig1 * exp(-epst / cappaut);

	    if (*sigy > sig1 * .01) {
		*hards = -(*sigy) / cappaut;
	    } else {
		*sigy = sig1 * .01;
		*hards = 1e-4;
	    }
/*     Druck-Bereich */
	} else if (yip == 4) {
	    cappae = *fcm * 4. / (*e * 3.);
	    cappauc = cappae + *gc * 3. / (*dia * 2. * *fcm);
	    sig3 = *fcm * *gamma2 / (*gamma2 * 2. - 1.);

	    if (epst < cappae) {
		fac = epst / cappae;
		sig3 /= 3.;
		*sigy = sig3 * (fac * 4. + 1. - pow(fac, 2.) * 2.);
		*hards = sig3 * 4. * (1. - fac) / cappae;
	    } else if (epst >= cappae && epst < cappauc) {
		fac = epst - cappae;
		fac /= cappauc - cappae;
		*sigy = sig3 * (1. - pow(fac, 2.));
		if (*sigy < 1e-4) {
		    *sigy = 1e-4;
		}
		dum = cappauc - cappae;
		*hards = sig3 * -2. * (epst - cappae) / (dum * dum);
	    } else if (epst >= cappauc) {
		*sigy = sig3 * 1e-4;
		*hards = sig3 * 1e-4;
	    }

	    if (*sigy < sig3 * .03) {
		*sigy = sig3 * .03;
		*hards = 1e-4;
	    }

	}
/*------------------------------ apply drucker prager yield criteria ---*/
	f = *ft - (*g * 2. + *com * 9 * dalpha) * *dlam + ro23 * (*sigym - *
		sigy);
/*---- derivative of the yield criteria with respect to plastic --------*/
/*---- increment                                                --------*/
	dfdl = -(*g * 2. + *com * 9 * dalpha) - ro23 * *hards;
	if (dfdl >= 0.)   dserror("CHECK THE SOFTENING PARAMETER");        
/*------------------------------------------- new plastic increment  ---*/
	*dlam -= f / dfdl;
/*------------------------------------------------ check convergence ---*/
	if (*sigy == 0.) {
	    if (fabs(f) > 1e-8) {
		if (i > 60) {
	           dserror("CONVERGENCE CHECK FAILED");        
		}
		goto L600;
	    }
	} else {
	    if ((dum = f / *sigy, fabs(dum)) > 1e-8) {
		if (i > 60) {
	           dserror("CONVERGENCE CHECK FAILED");        
		}
		goto L600;
	    }
	}
	*epstn = epst;
/*-------- new stress state  sigma = dm * h * sigma(pred) --------------*/
/*-------- (considering the back stress vector)           --------------*/
	for (i = 0; i < 4; i++) {
	    devsigi[i] = devsig[i] - *g * 2. * dn[i] * *dlam;
	    smi[i] = sm[i] - *com * 3. * dcom[i] * *dlam;
	}
	sigma[0] = devsigi[0] + smi[0];
	sigma[1] = devsigi[1] + smi[1];
	sigma[3] = devsigi[2] + smi[2];
	sigma[2] = devsigi[3] + smi[3];
/*---------------------------------------------- rotational symmetry ---*/
    } else  {

	dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");        
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1cradi */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void w1mapl2(DOUBLE *tau,      /* current stresses (local)              */
             DOUBLE **d,       /* material matrix to be calculated      */
             DOUBLE *dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* type of problem                       */
             DOUBLE *alpha,    /* neigungswinkel der fliessflaechen     */
             DOUBLE *emod,     /* elastizitaetsmodul                    */
             DOUBLE *g,        /* schubmodul                            */
             DOUBLE *com,      /* kompressionsmodul                     */
             DOUBLE *betah,    /* factor for isotrop/kinemat. hardening */
             DOUBLE *hards,    /* plastic hardeningmodulus              */
             DOUBLE *dn,       /* gradient components in dev. direction */
             DOUBLE *grad,     /* total gradient                        */
             DOUBLE *dev)      /* norm of the dev. predictor stresses   */
{
/*----------------------------------------------------------------------*/
static DOUBLE half, fact, vect[4], stps, fact1,dum;
static INT i, j;
static DOUBLE x[4];
static DOUBLE stpsi, dm[4][4], hm[4][4], q13, q23; 
static DOUBLE dalpha, add, fac, fkg, fki, ro23; 
static DOUBLE xsi1, xsi2, xsi3, xsi4;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1mapl2");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
    q23 = .66666666666666663;

    dalpha = *alpha * *alpha;
    /* concrete */
    fact = ro23 * *hards;
/*----------------------------------------------------------------------*/
  if(wtype==plane_stress || wtype==plane_strain)
  {
/*----------------------------------------------------------------------*/
    xsi1 = tau[0];
    xsi2 = tau[1];
    xsi3 = tau[3];
    xsi4 = tau[2];

    dum = q13 * (xsi1 * (xsi1 * 2. - xsi2 - xsi3) 
               + xsi2 * (xsi2 * 2. - xsi1 - xsi3) 
               + xsi3 * (xsi3 * 2. - xsi1 - xsi2)) 
               + xsi4 * xsi4 * 2.;
    stps = sqrt(dum);

    if (stps < 9.9999999999999998e-13) 
    {
        fkg = *g * 2.;
        fki = 1.;
        stpsi = 0.;
    } else {
        fkg = *g * 2. / stps;
        fki = *g * 2. * *dlam;
        stpsi = stps - *g * 2. * *dlam;
    }
/*----------------- product:   GRAD = P * TAU / |S|  +  ALPHA * KRON ---*/
    x[0] = grad[0];
    x[1] = grad[1];
    x[2] = grad[2];
    x[3] = grad[3] * 2.;
/*-----------------------------------------------------------------------|
|        inverse modified constitutive tensor (compliance)               |
|                 D  =  D(INV) + P * DLAM                                |
|        inverse material stiffness tensor     DM = INV (D(INV) + P*DLAM |
|-----------------------------------------------------------------------*/
    for (i=0; i<4; i++)for (j=0; j<4; j++) hm[i][j] = 0.;
    hm[0][0] =  *com + fkg*(fki*dn[0]*dn[0] + stpsi * q23);
    hm[0][1] =  *com + fkg*(fki*dn[0]*dn[1] - stpsi * q13);
    hm[0][2] =  *com + fkg*(fki*dn[0]*dn[2] - stpsi * q13);
    hm[0][3] =         fkg*(fki*dn[0]*dn[3]); 

    hm[1][1] =  *com + fkg*(fki*dn[1]*dn[1] + stpsi * q23);
    hm[1][2] =  *com + fkg*(fki*dn[1]*dn[2] - stpsi * q13);
    hm[1][3] =         fkg*(fki*dn[1]*dn[3]); 

    hm[2][2] =  *com + fkg*(fki*dn[2]*dn[2] + stpsi * q23);
    hm[2][3] =         fkg*(fki*dn[2]*dn[3]); 

    hm[3][3] =         fkg*(fki*dn[3]*dn[3] + stpsi * half);  

    for (i=1; i<4; i++)
      for (j=0; j<i; j++)
            hm[i][j] = hm[j][i];
/*--------------------------------------------- product:   HM * GRAD ---*/
    for (i=0; i<4; i++)
    {
      add = 0.;
      for (j=0; j<4; j++)
      {
        add = add + hm[i][j]*x[j];
      }
      vect[i] = add;
    }
/*------------------------------- scalar product:   GRAD * HM * GRAD ---*/
    fact1 = 0.;
    for (i=0; i<4; i++) fact1 = fact1 + x[i]*vect[i];
    fact1 = 1. / (fact1 + fact);
/*------------------------------ form consistent 4x4 material matrix ---*/
    for (i=0; i<4; i++)for (j=0; j<4; j++) dm[i][j] = 0.;
    for (i=0; i<4; i++)for (j=0; j<4; j++) d[ i][j] = 0.;
    
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
            dm[i][j] = hm[i][j] - fact1 * vect[i] * vect[j];

  
    d[0][0] = dm[0][0];
    d[1][0] = dm[1][0]; 
    d[2][0] = dm[3][0];
    d[3][0] = dm[2][0];

    d[0][1] = dm[0][1];
    d[1][1] = dm[1][1];
    d[2][1] = dm[3][1];
    d[3][1] = dm[2][1];

    d[0][2] = dm[0][3];
    d[1][2] = dm[1][3];
    d[2][2] = dm[3][3]; 
    d[3][2] = dm[2][3];

    d[0][3] = dm[0][2];
    d[1][3] = dm[1][2];
    d[2][3] = dm[3][2];
    d[3][3] = dm[2][2];
/*----------------------------------------------------------------------*/
  }else{
         dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1mapl2 */

/*-----------------------------------------------------------------------|
|       topic : radial return for multisurface problems                  |
|               (material model: concrete)                               |
|-----------------------------------------------------------------------*/
void w1cradms(DOUBLE *sigma,  /* elastic predictor projected onto yield surface */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain             */ 
              DOUBLE *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              INT     yip,    /* stress state   1 =elastic   >=2 =plastic       */ 
              DOUBLE *alpha,  /* factor for the first invariants                */ 
              DOUBLE *ft,     /* yield condition                                */ 
              DOUBLE  e,      /* young modulus                                  */ 
              DOUBLE  g,      /* shear modulus                                  */ 
	      DOUBLE  com,    /* bulk modulus                                   */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress                */ 
              DOUBLE *hards,  /* plastic modulus                                */ 
              DOUBLE  dn[3][4],  /* gradient components in deviatoric direction */ 
              DOUBLE  dcom[3][4],/* gradient components in hdrostatic direction */ 
              DOUBLE *devsig, /* deviatoric predictor stresses                  */ 
              DOUBLE *sm,     /* hydrostatic predictor stresses                 */ 
              DOUBLE  fcm,    /* compressive strenght                           */ 
              DOUBLE  gc,     /* compression fracture energy                    */ 
              DOUBLE  ftm,    /* tensile strenght                               */ 
	      DOUBLE  gt,     /* tensile fracture energy                        */ 
              DOUBLE  gamma1, /* fitting factor yield function 1                */ 
              DOUBLE  gamma2, /* fitting factor yield function 2                */ 
              DOUBLE  dia,    /* equivalent element length                      */ 
              DOUBLE  acrs)   /* average crack spacing                          */ 
{                             
/*----------------------------------------------------------------------*/
    DOUBLE dum;               
/*----------------------------------------------------------------------*/
    DOUBLE dum1, dum2;

    static DOUBLE half, soli[4]	, epst[2], sigy[3];
    static DOUBLE f;
    static INT i;
    static DOUBLE dhard[2];
    static INT kflag;
    static DOUBLE c1, c2, alpha1, alpha2, df[4];
    static DOUBLE q13, cappae, fac, dlambda[2], det, dev;
    static DOUBLE cappaut, sig1, sig2, sig3;
    static DOUBLE fti[2], hyd, ro23, cappauc, sol[4];
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
/*--------------------------------------------- initialize variables ---*/
    i = 0;
    dlambda[0] = 0.;
    dlambda[1] = 0.;
    c1 = 1.;
    c2 = 1.;

    if (yip == 2) {
	alpha1 = alpha[0];
	alpha2 = alpha[1];
    } else if (yip == 3) {
	alpha1 = alpha[0];
	alpha2 = alpha[2];
    }

    df[0] = g * 2. + com * 9 * alpha1 * alpha1;
    df[2] = g * 2. + com * 9 * alpha1 * alpha2;
    df[1] = df[2];
    df[3] = g * 2. + com * 9 * alpha2 * alpha2;

/*------------------------------------- plane strain / plane stress  ---*/
    if (wtype==plane_stress || wtype==plane_strain) 
    {
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
	if (yip == 2) {
	    epst[0] = epstn[0] + dlambda[0];
	    epst[1] = epstn[1];
	} else if (yip == 3) {
	    epst[0] = epstn[0] + dlambda[0];
	    epst[1] = epstn[1] + dlambda[1];
	}
/*------------------------------------------------- new yield stress ---*/
/*   Zug-/Zug-Druck-Bereiche */

	cappaut = gt / (acrs * ftm);
	sig1 = gamma1 * 2. * fcm * ftm / (gamma1 * fcm + ftm);
	sigy[0] = sig1 * exp(-epst[0] / cappaut);

	if (epst[0] == 0.) {
	    hards[0] = 0.;
	} else if (sigy[0] > sig1 * .01) {
	    hards[0] = -sigy[0] / cappaut;
	} else {
	    sigy[0] = sig1 * .01;
	    hards[0] = 1e-4;
	}

/*----------------------------------------------------------------------*/
/*   Vorwerte fuer: APEX */

	if (yip == 2) {
	    dum1 = alpha[0];
	    sigy[1] = -sigy[0] / (dum1 * dum1 * 3.);
	    dum1 = alpha[0];
	    hards[1] = -hards[0] / (dum1 * dum1 * 3.);
	    sigy[2] = 1e-4;

/*----------------------------------------------------------------------*/
/*   Vorwerte fuer: Druck-Bereich */

	} else if (yip == 3) {
	    cappae = fcm * 4. / (e * 3.);
	    cappauc = cappae + gc * 3. / (dia * 2. * fcm);
	    sig3 = fcm * gamma2 / (gamma2 * 2. - 1.);
	    sigy[1] = 1e-4;

	    if (epst[1] == 0.) {
		sig3 /= 3.;
		sigy[2] = sig3;
		hards[2] = 0.;
	    } else if (epst[1] < cappae && epst[1] > 0.) {
		fac = epst[1] / cappae;
		sig3 /= 3.;
		sigy[2] = sig3 * (fac * 4. + 1. - pow(fac, 2.) * 2.);
		hards[2] = sig3 * 4. * (1. - fac) / cappae;
	    } else if (epst[1] >= cappae && epst[1] < cappauc) {
		fac = epst[1] - cappae;
		fac /= cappauc - cappae;
		sigy[2] = sig3 * (1. - pow(fac, 2.));
		if (sigy[2] < 1e-4) {
		    sigy[2] = 1e-4;
		}

		dum1 = cappauc - cappae;
		hards[2] = sig3 * -2. * (epst[1] - cappae) / (dum1 * dum1);
	    } else if (epst[1] >= cappauc) {
		sigy[2] = 1e-4;
		hards[2] = 1e-4;
	    }

	    if (sigy[2] < sig3 * .03) {
		sigy[2] = sig3 * .03;
		hards[2] = 1e-4;
	    }
	}
/*------------------------------ apply drucker prager yield criteria ---*/
	if (yip == 2) {
	    fti[0] = ft[0] + ro23 * (sigym[0] - sigy[0]) - (df[0] * dlambda[0]
		     + df[2] * dlambda[1]);
	    fti[1] = ft[1] + ro23 * (sigym[1] - sigy[1]) - (df[1] * dlambda[0]
		     + df[3] * dlambda[1]);
	    dhard[0] = hards[0];
	    dhard[1] = hards[1];
	    sig1 = sigy[0];
	    sig2 = sigy[1];
	} else if (yip == 3) {
	    fti[0] = ft[0] + ro23 * (sigym[0] - sigy[0]) - (df[0] * dlambda[0]
		     + df[2] * dlambda[1]);
	    fti[1] = ft[2] + ro23 * (sigym[2] - sigy[2]) - (df[1] * dlambda[0]
		     + df[3] * dlambda[1]);
	    dhard[0] = hards[0];
	    dhard[1] = hards[2];
	    sig1 = sigy[0];
	    sig2 = sigy[2];
	}
/*---- derivative of the yield criteria with respect to plastic --------*/
/*---- increment                                                --------*/
	sol[0] = -df[0] - ro23 * dhard[0];
	sol[2] = -df[2];
	if (yip == 2) {
	    sol[1] = -df[1] - ro23 * dhard[1];
	    sol[3] = -df[3];
	} else if (yip == 3) {
	    sol[1] = -df[1];
	    sol[3] = -df[3] - ro23 * dhard[1];
	}

	if (sol[0] >= 0.) {
	           dserror("CHECK THE SOFTENING PARAMETER");        
	}
/*--------------------------------- modify the multisurface problem  ---*/
L700:
	kflag = 0;
	fti[0] = c1 * fti[0] + (1. - c1) * dlambda[0];
	fti[1] = c2 * fti[1] + (1. - c2) * dlambda[1];

	sol[0] = c1 * sol[0] + (1. - c1);
	sol[2] = c1 * sol[2];
	sol[1] = c2 * sol[1];
	sol[3] = c2 * sol[3] + (1. - c2);
/*---------------------------------- solution of the equation system ---*/
	det = sol[0] * sol[3] - sol[1] * sol[2];
	soli[0] = sol[3] / det;
	soli[1] = -sol[1] / det;
	soli[2] = -sol[2] / det;
	soli[3] = sol[0] / det;

	dlambda[0] -= soli[0] * fti[0] + soli[2] * fti[1];
	dlambda[1] -= soli[1] * fti[0] + soli[3] * fti[1];

	dum1 = fti[0] / sig1;

	dum2 = fti[1] / sig2;
	f = sqrt(dum1 * dum1 + dum2 * dum2);
/*----------------------------------- check the multisurface problem ---*/
	if (yip == 2) {
	    if (dlambda[0] < 0.) {
		c1 = 0.;
		kflag = 1;
	    } else if (dlambda[1] > 0.) {
		c2 = 0.;
		kflag = 1;
	    }
	} else if (yip == 3) {
	    if (dlambda[0] < 0.) {
		c1 = 0.;
		kflag = 1;
	    } else if (dlambda[1] < 0.) {
		c2 = 0.;
		kflag = 1;
	    }
	}

	if (kflag == 1) {
	    dlambda[0] = 0.;
	    dlambda[1] = 0.;
	    fti[0] = ft[0];
	    if (yip == 2) {
		fti[1] = ft[1];
	    } else if (yip == 3) {
		fti[1] = ft[2];
	    }
	    goto L700;
	}
/*------------------------------------------------ check convergence ---*/
	if (f > 1e-8) {
	    if (i > 60) {
		 dserror("CHECK CONVERGENCE"); 
	    }
	    goto L600;
	}

	epstn[0] = epst[0];
	epstn[1] = epst[1];

	dlam[0] = dlambda[0];
	dlam[1] = dlambda[1];

	if (yip == 2) {
	    dlam[0] = dlambda[0] + dlambda[1];
	    dlam[1] = 0.;
	} else if (yip == 3) {
	    dlam[0] = dlambda[0];
	    dlam[1] = dlambda[1];
	}
/*-------- new stress state  sigma = dm * h * sigma(pred) --------------*/
/*-------- (considering the back stress vector)           --------------*/
	for (i = 0; i < 4; ++i) {
	    if (yip == 2) {
		devsig[i] = devsig[i] - g * 2. * dn[0][i] * dlambda[0] - g 
			* 2. * dn[1][i] * dlambda[1];
		sm[i] = sm[i] - com * 3. * dcom[0][i] * dlambda[0] - com * 
			3. * dcom[1][i] * dlambda[1];
	    } else if (yip == 3) {
		devsig[i] = devsig[i] - g * 2. * dn[0][i] * dlambda[0] - g 
			* 2. * dn[2][i] * dlambda[1];
		sm[i] = sm[i] - com * 3. * dcom[0][i] * dlambda[0] - com * 
			3. * dcom[2][i] * dlambda[1];
	    }
	}

	sigma[0] = devsig[0] + sm[0];
	sigma[1] = devsig[1] + sm[1];
	sigma[3] = devsig[2] + sm[2];
	sigma[2] = devsig[3] + sm[3];

	w1pres(&sigma[0], &devsig[0], &sm[0], &dev, &hyd);
	w1yicsr(dev, hyd, sigy[0], alpha[0], &ft[0]);
	w1yicsr(dev, hyd, sigy[1], alpha[1], &ft[1]);
	w1yicsr(dev, hyd, sigy[2], alpha[2], &ft[2]);
/*----------------------------------------------------------------------*/
  }else{
         dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1cradms */
/*-----------------------------------------------------------------------|
|     topic: convert some arrays                                         |
|-----------------------------------------------------------------------*/
      void w1conver(DOUBLE *dlam,       /*  incr. of plastic multiplier */             
                    DOUBLE *alpha,      /*  factor for the first inv.   */             
                    DOUBLE *hards,      /*  plastischer hardeningmodul  */             
                    DOUBLE dn[3][4],    /*  gradient comp. in dev. dir. */  
                    DOUBLE grad[3][4],  /*  total gradient              */             
                    DOUBLE *dlamc,      /*  -|                          */             
                    DOUBLE *alphac,     /*   |                          */             
                    DOUBLE *hardsc,     /*   |-  converted arrays       */             
                    DOUBLE dnc[2][4],   /*   |                          */             
                    DOUBLE gradc[2][4]) /*  -|                          */             
{
/*----------------------------------------------------------------------*/
    static INT i;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_enter("w1conver");
    #endif
/*----------------------------------------------------------------------*/
      *dlamc = dlam[0] + dlam[1];
/*----------------------------------------------------------------------*/
      alphac[0] = alpha[0];
      alphac[1] = alpha[2];
/*----------------------------------------------------------------------*/
      hardsc[0] = hards[0];
      hardsc[1] = hards[2];
/*----------------------------------------------------------------------*/
    for (i=0; i<4; i++)
    {
        dnc[0][i]   = dn[0][i];
        dnc[1][i]   = dn[2][i];
        gradc[0][i] = grad[0][i];
        gradc[1][i] = grad[2][i];
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1conver */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------|
 |   forms the elasto-plastic consistent tangent material tensor        |
 |   for wall element   (drucker - prager)                              |
 |   (generalized for the multisurface problem: including the apex)     |
 *----------------------------------------------------------------------*/
void w1maplg(DOUBLE *tau,      /* current stresses (local)              */
             DOUBLE **d,       /* material matrix to be calculated      */
             DOUBLE  dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* type of problem                       */
             INT     yip,      /* stress state  1=elastic 2=plastic     */
             DOUBLE  emod,     /* elastizitaetsmodul                    */
             DOUBLE  g,        /* schubmodul                            */
             DOUBLE  com,      /* kompressionsmodul                     */
             DOUBLE *hards,    /* plastic hardeningmodulus              */
             DOUBLE  dn[2][4], /* gradient components in dev. direction */
             DOUBLE  grad[2][4])/*norm of the dev. predictor stresses   */
{
/*----------------------------------------------------------------------*/
    DOUBLE dum;

    static INT i, j, k;
    static DOUBLE fact[2][2], ghme[4][2], stps;
    static DOUBLE e[2][2];
    static DOUBLE e2[2][2], stpsi, dm[4][4];
    static DOUBLE  hm[4][4], q13, q23, dalpha, facmin;
    static DOUBLE fak, fkg, fki, ghm[4][2], det, ro23; 
    static DOUBLE tra[4][4], xsi1, xsi2, xsi3, xsi4;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1maplg");
#endif
/*----------------------------------------------------------------------*/
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
    q23 = .66666666666666663;
/*----------------------------------------------------------------------*/
    fact[0][0] = ro23 * hards[0];
    fact[0][1] = 0.;
    if (yip == 2) {
        fact[1][0] = ro23 * hards[1];
        fact[1][1] = 0.;
    } else if (yip == 3) {
        fact[1][0] = 0.;
        fact[1][1] = ro23 * hards[1];
    }
/*----------------------------------------------------------------------*/
  if(wtype==plane_stress || wtype==plane_strain)
  {

	xsi1 = tau[0];
	xsi2 = tau[1];
	xsi3 = tau[3];
	xsi4 = tau[2];

	dum = q13 * (xsi1 * (xsi1 * 2. - xsi2 - xsi3)
                   + xsi2 * (xsi2 * 2. - xsi1 - xsi3)
                   + xsi3 * (xsi3 * 2. - xsi1 - xsi2))
                   + xsi4 * xsi4 * 2.;
	stps = sqrt(dum);

	if (stps < 1e-5) {
	    fkg = g * 2.;
	    fki = 1.;
	    stpsi = 0.;
	} else {
	    fkg = g * 2. / stps;
	    fki = g * 2. * dlam;
	    stpsi = stps - g * 2. * dlam;
	}
/*------------ product   GRAD = P * TAU / |S|  +  ALPHA(1) * KRON ------*/
/*------------                =   GRAD(DEV)    +     GRAD(COM)    ------*/
	for (j = 0; j < 2; j++) {
	    grad[j][0] = grad[j][0];
	    grad[j][1] = grad[j][1];
	    grad[j][2] = grad[j][2];
	    grad[j][3] *= 2.;
	}
/*-----------------------------------------------------------------------|
|        inverse modified constitutive tensor (compliance)               |
|                 D  =  D(INV) + P * DLAM                                |
|        inverse material stiffness tensor     DM = INV (D(INV) + P*DLAM |
|-----------------------------------------------------------------------*/
    for (i=0; i<4; i++)for (j=0; j<4; j++) hm[i][j] = 0.;
    hm[0][0] =  com + fkg*(fki*dn[0][0]*dn[0][0] + stpsi * q23);
    hm[0][1] =  com + fkg*(fki*dn[0][0]*dn[0][1] - stpsi * q13);
    hm[0][2] =  com + fkg*(fki*dn[0][0]*dn[0][2] - stpsi * q13);
    hm[0][3] =        fkg*(fki*dn[0][0]*dn[0][3]); 

    hm[1][1] =  com + fkg*(fki*dn[0][1]*dn[0][1] + stpsi * q23);
    hm[1][2] =  com + fkg*(fki*dn[0][1]*dn[0][2] - stpsi * q13);
    hm[1][3] =        fkg*(fki*dn[0][1]*dn[0][3]); 

    hm[2][2] =  com + fkg*(fki*dn[0][2]*dn[0][2] + stpsi * q23);
    hm[2][3] =        fkg*(fki*dn[0][2]*dn[0][3]); 

    hm[3][3] =        fkg*(fki*dn[0][3]*dn[0][3] + stpsi * .5);  

    for (i=1; i<4; i++)
      for (j=0; j<i; j++)
            hm[i][j] = hm[j][i];
        for (i=0; i<4; i++)for (j=0; j<4; j++) hm[i][j] = 0.;
/*-------------------------------------- HARDENING MATRIX E (2X2) ------*/
    for (i=0; i<2; i++)for (j=0; j<2; j++) e[i][j] = 0.;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
            e[i][j] = fact[i][j];
/*-----------------------------------------------------------------------|
|        DENOMINATOR:                                                    |
|        E2 (2X2) = (GRADT * HM) * GRAD  +  E                            |
|                 =     GHMT     * GRAD  +  E                            |
|                 =             E2       +  E                            |
|        AND ITS INVERSE INV(E2)                                         |
|-----------------------------------------------------------------------*/
        for (i=0; i<4; i++)for (j=0; j<2; j++) ghm[i][j] = 0.;
        for (i=0; i<2; i++)for (j=0; j<2; j++) e2[i][j] = 0.;
        
        for (i=0; i<4; i++)
          for (j=0; j<4; j++) 
            for (k=0; k<2; k++) 
              ghm[i][j] = hm[i][k]* grad[k][j];

        for (i=0; i<4; i++)
          for (j=0; j<2; j++) 
            for (k=0; k<2; k++) 
              ghm[i][j] = grad[k][i]* e2[k][j];
	
        for (i=0; i<2; i++)for (j=0; j<2; j++) e[i][j] += e2[i][j];

	det = e[0][0] * e[1][1] - e[1][0] * e[0][1];
	e2[0][0] =  e[1][1] / det;
	e2[1][0] = -e[1][0] / det;
	e2[0][1] = -e[0][1] / det;
	e2[1][1] =  e[0][0] / det;
/*-----------------------------------------------------------------------|
|       numerator:                                                       |
|       TRA(4,4) =  GHM * INV(E2) * GHMT                                 |
|-----------------------------------------------------------------------*/
        for (i=0; i<4; i++)for (j=0; j<2; j++) ghme[i][j] = 0.;
        for (i=0; i<4; i++)for (j=0; j<4; j++)  tra[i][j] = 0.;

        for (i=0; i<4; i++)
          for (j=0; j<2; j++) 
            for (k=0; k<2; k++) 
              ghme[i][j] = ghm[i][k]* e2[k][j];

        for (i=0; i<4; i++)
          for (j=0; j<2; j++) 
            for (k=0; k<4; k++) 
	    tra[i][j] += ghme[i][j] * ghm[j][k];

/*------------------------------ form consistent 4x4 material matrix ---*/
    for (i=0; i<4; i++)for (j=0; j<4; j++) dm[i][j] = 0.;
    for (i=0; i<4; i++)for (j=0; j<4; j++) d[ i][j] = 0.;

    fak = -1.;
    for (i=0; i<4; i++)for (j=0; j<4; j++)  tra[i][j] *= fak;
    for (i=0; i<2; i++)for (j=0; j<2; j++) dm[i][j] = hm[i][j] + tra[i][j];
    
    /* initialization at the apex */
    facmin = 1e-5;
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
        if(fabs(dm[i][j])>facmin) facmin = dm[i][j];
    
    if(facmin < emod * 1e-5)
    {
      for (i=0; i<4; i++)for (j=0; j<4; j++) dm[i][j] = 0.;
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            dm[i][j] = emod * 1e-5;
    }   

    d[0][0] = dm[0][0];
    d[1][0] = dm[1][0]; 
    d[2][0] = dm[3][0];
    d[3][0] = dm[2][0];

    d[0][1] = dm[0][1];
    d[1][1] = dm[1][1];
    d[2][1] = dm[3][1];
    d[3][1] = dm[2][1];

    d[0][2] = dm[0][3];
    d[1][2] = dm[1][3];
    d[2][2] = dm[3][3]; 
    d[3][2] = dm[2][3];

    d[0][3] = dm[0][2];
    d[1][3] = dm[1][2];
    d[2][3] = dm[3][2];
    d[3][3] = dm[2][2];
/*----------------------------------------------------------------------*/
  }else{
         dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1maplg */
/*----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------|
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)      CAP-REGION                               |
|-----------------------------------------------------------------------*/
void w1radcap(DOUBLE *sigma,  /* elastic predictor projected onto yield surface */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain             */ 
              DOUBLE *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              DOUBLE *alpha,  /* factor for the first invariants                */ 
              DOUBLE *e,      /* young modulus                                  */ 
              DOUBLE *g,      /* shear modulus                                  */ 
	      DOUBLE *com,    /* bulk modulus                                   */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress                */ 
              DOUBLE *hards,  /* plastic modulus                                */ 
              DOUBLE *grad,   /* total gradient                                 */ 
              DOUBLE *devsig, /* deviatoric predictor stresses                  */ 
              DOUBLE *sm,     /* hydrostatic predictor stresses                 */ 
              DOUBLE *dev,     /* norm of the deviatoric predictor stresses     */ 
              DOUBLE *hyd,     /* 1st invariant  of the predictor stresses      */ 
              DOUBLE *hydn,    /* 1st invariant  of the new stresses            */ 
              DOUBLE *fcm,    /* compressive strenght                           */ 
              DOUBLE *gc,     /* compression fracture energy                    */ 
              DOUBLE *ftm,    /* tensile strenght                               */ 
	      DOUBLE *gt,     /* tensile fracture energy                        */ 
              DOUBLE *gamma2, /* fitting factor yield function 2                */ 
              DOUBLE *dia)    /* equivalent element length                      */ 
{
/*----------------------------------------------------------------------*/
    DOUBLE dum1, dum2;               
/*----------------------------------------------------------------------*/
    static INT i;
    
    static DOUBLE dadk, dadl, dbdl, dbdk, half, fact[4], dcom[4]; 
    static DOUBLE soli[2][2], epst, root;
    static DOUBLE a, b, f;
    static DOUBLE hardi, hardr, rdevc, rcomc, rdevn;
    static DOUBLE depst, rcomn, dn[4], q13, q23, cappae, q29, fac;
    static DOUBLE det, fti[2], ro23, smi[4], ro54, cappauc, sol[2][2];
    static DOUBLE devsigi[4], fac1, fac2, sig3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1radcap");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    ro54 = sqrt(54.);
    q13 = .33333333333333331;
    q23 = .66666666666666663;
    q29 = .22222222222222221;
/*--------------------------------------------- initialize variables ---*/
    i = 0;
    *dlam = 0.;
    depst = 0.;

    fact[0] = 1.;
    fact[1] = 1.;
    fact[2] = 1.;
    fact[3] = 0.;
/*------------------------------------- plane strain / plane stress  ---*/
    if (wtype==plane_stress || wtype==plane_strain) 
    {
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L600:
	++i;
/*--------------------------------- new plastic uniaxial deformation ---*/
	epst = *epstn + depst;
/*------------------------------------------------- new yield stress ---*/
/*     Zug-/Zug-Druck-Bereiche */
	cappae = *fcm * 4. / (*e * 3.);
	cappauc = cappae + *gc * 3. / (*dia * 2. * *fcm);
	sig3 = *fcm * *gamma2 / (*gamma2 * 2. - 1.);

	if (epst < cappae) {
	    fac = epst / cappae;
	    sig3 /= 3.;
	    sigym[2] = sig3 * (fac * 4. + 1. - pow(fac, 2.) * 2.);
	    hards[2] = sig3 * 4. * (1. - fac) / cappae;
	} else if (epst >= cappae && epst < cappauc) {
	    fac = epst - cappae;
	    fac /= cappauc - cappae;
	    sigym[2] = sig3 * (1. - pow(fac, 2.));
	    if (sigym[2] < 1e-6) {
		sigym[2] = 1e-6;
	    }
	    dum1 = cappauc - cappae;
	    hards[2] = sig3 * -2. * (epst - cappae) / (dum1 * dum1);
	} else if (epst >= cappauc) {
	    sigym[2] = 1e-6;
	    hards[1] = 1e-6;
	}

	if (sigym[2] < sig3 * .03) {
	    sigym[2] = sig3 * .03;
	    hards[2] = 1e-6;
	}

/*     CAP-Bereich */

	alpha[3] = -(ro54 * alpha[2] + 2.) * *gamma2;
	dum1 = alpha[2];
	sigym[3] = sqrt(q23 + dum1 * dum1 * 6.) * *gamma2 * sigym[2];
	hards[3] = hards[2] * sigym[3] / sigym[2];

	hardr = hards[3];
	hardi = hards[2] * alpha[3];
	fac1 = hardi / hardr;
	fac2 = 1.;
/*----------- apply cap yield criteria and associated hardening rule ---*/
/*----------- for a local newton-iteration                           ---*/
	rdevc = *dev;
	rdevn = sigym[3] + *g * 2. * *dlam;
	rcomc = *hyd - alpha[3] * sigym[2];
	rcomn = sigym[3] + *com * *dlam;

	a = rdevc / rdevn;
	b = rcomc / rcomn;
	root = sqrt(pow(a, 2.) + pow(b, 2.) / 9);

	fti[0] = pow(a, 2.) + pow(b, 2.) / 9 - 1.;
	fti[1] = depst - *dlam * (fac1 * b / (root * 9) + fac2);
/*------------- derivative of both equations with respect to plastic ---*/
/*------------- increment and the internal strain increment          ---*/
	dadl = -rdevc * 2. * *g / pow(rdevn, 2.);
	dadk = -rdevc * hardr / pow(rdevn, 2.);
	dbdl = -rcomc * *com / pow(rcomn, 2.);

	dum1 = rcomn;
	dbdk = -(rcomn * hardi + rcomc * hardr) / (dum1 * dum1);

	sol[0][0] = pow(rdevc, 2.) * -2. * 2. * *g / pow(rdevn, 3.) 
                  - q29 * rcomc * (rcomc * *com) / pow(rcomn, 3.);
	sol[0][1] = pow(rdevc, 2.) * -2. * hardr / pow(rdevn, 3.) 
		- q29 * rcomc * (hardi * rcomn + rcomc * hardr) / pow(
		rcomn, 3.);

	sol[1][0] = -(fac1 * b / (root * 9) + fac2) - (root * dbdl - b * (a * 
		dadl + b * dbdl / 9) / root) / (pow(root, 2.) * 9) * *
		dlam * fac1;
	sol[1][1] = 1. - (root * dbdk - b * (a * dadk + b * dbdk / 9) / root) / (
		pow(root, 2.) * 9) * *dlam * fac1;

	if (sol[0][0] >= 0.) dserror("CHECK THE SOFTENING PARAMETER");        
/*----------------------------------- solution of the equation system---*/
	det = sol[0][0] * sol[1][1] - sol[1][0] * sol[0][1];
	soli[0][0] =  sol[1][1] / det;
	soli[1][0] = -sol[1][0] / det;
	soli[0][1] = -sol[0][1] / det;
	soli[1][1] =  sol[0][0] / det;

	*dlam -= soli[0][0] * fti[0] + soli[0][1] * fti[1];
	depst -= soli[1][0] * fti[0] + soli[1][1] * fti[1];

	dum1 = fti[0] / 1e-6;
	dum2 = fti[1] / 1e-6;
	f = sqrt(dum1 * dum1 + dum2 * dum2);
/*------------------------------------------------ check convergence ---*/
	if (fabs(f) > 1e-6) {
	    if (i > 60)  dserror("CONVERGENCE CHECK FAILED");
	    goto L600;
	} else {
	    *epstn = epst;
	}
/*-------- new stress state  sigma                        --------------*/
/*-------- (considering the back stress vector)           --------------*/
	for (i = 0; i < 4; ++i) {
	    devsigi[i] = devsig[i] * sigym[3] / rdevn;
	    smi[i] = *hyd * sigym[3] + alpha[3] * sigym[2] * *com * *dlam;
	    smi[i] = fact[i] * smi[i] / (rcomn * 3.);
	}

	sigma[0] = devsigi[0] + smi[0];
	sigma[1] = devsigi[1] + smi[1];
	sigma[3] = devsigi[2] + smi[2];
	sigma[2] = devsigi[3] + smi[3];

	*hydn = smi[0] + smi[1] + smi[2];
/*--------------------------------------- components ot the gradient ---*/
	for (i = 0; i < 4; ++i) {
	    dn[i] = devsigi[i] / sigym[3];
	    dcom[i] = fact[i] * (*hydn - alpha[3] * sigym[2]) / (
		    sigym[3] * 9);
	    grad[i] = dn[i] + dcom[i];
	}
/*---------------------------------------------- rotational symmetry ---*/
    } else  {

	dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");        
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1radcap */
/*----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------|
|       topic: forms the  e l a s t o - p l a s t i c                    |
|              c o n s i s t e n t  tangent material tensor              |
|              for wall element   (drucker - prager)                     |
|-----------------------------------------------------------------------*/
void w1mplcap(DOUBLE *tau,      /* current stresses (local)              */
              DOUBLE **d,       /* material matrix to be calculated      */
              DOUBLE *dlam,     /* increment of plastic multiplier       */
              WALL_TYPE wtype,  /* type of problem                       */
              DOUBLE *alpha,    /* factor for the first invariants       */
              DOUBLE *emod,     /* young modulus                         */
              DOUBLE *vnu,      /* poisson's ratio                       */
              DOUBLE *hards,    /* plastic hardeningmodulus              */
              DOUBLE *sigym,    /* [0] uniaxial tension yield stress     */
              DOUBLE *grad,     /* total gradient                        */
              DOUBLE *cappae,   /* max. elastic compression strain       */
              DOUBLE *cappauc,  /* compressive fracture strain           */
              DOUBLE *epst,     /* equivalent uniaxial plastic strain    */
              DOUBLE *sig3)     /* equivalent compressive stress         */
{
/*----------------------------------------------------------------------*/
    DOUBLE dum;               
/*----------------------------------------------------------------------*/
    static INT i, j ,k;
    
    static DOUBLE diag[4][4], half, hard, dkkf, hinv[4][4];
    static DOUBLE hydr[4], fack1, quot, a[2][2];
    static DOUBLE q[4][4], facth, hardi, gradk[4];
    static DOUBLE dkapp, hardr, fkapp;
    static DOUBLE theta[4][4];
    static DOUBLE dconst, const1, const2, const3, ai[2][2];
    static DOUBLE dm[4][4], q13, q23, qn[4][4], dhards;
    static DOUBLE dn1[4], dn2[4], dkf, grd[4], det, sig, ro12;
    static DOUBLE ro13, hyd, ro16, fac1, fac2, fac3, val3, xsi1;
    static DOUBLE xsi2, xsi3, xsi4;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1mplcap");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro12 = sqrt(.5);
    ro13 = sqrt(.33333333333333331);
    ro16 = sqrt(.16666666666666666);
    q13 = .33333333333333331;
    q23 = .66666666666666663;

    if (*dlam < 1e-10) { 
	*dlam = 1e-10;
    }
/*------------------------------------- plane strain / plane stress  ---*/
    if (wtype==plane_stress || wtype==plane_strain) 
    {
	xsi1 = tau[0];
	xsi2 = tau[1];
	xsi3 = tau[3];
	xsi4 = tau[2];

	hyd = xsi1 + xsi2 + xsi3;

	fack1 = alpha[3] * sigym[2] / sigym[3];
	quot = (hyd - alpha[3] * sigym[2]) / sigym[3];
	facth = fack1 * quot / 9;
	fkapp = facth + 1.;
	dkapp = fkapp * *dlam;

	hardi = alpha[3] * hards[2];
	hardr = hards[3];

	sig = *sig3 * 3.;
	if (*epst < *cappae) {
	    dhards = sig * -4. / (*cappae* *cappae);
	} else if (*epst >= *cappae && *epst < *cappauc) {
	    dum = *cappauc - *cappae;
	    dhards = sig * -2. / (dum * dum);
	} else if (*epst >= *cappauc) {
	    dhards = 0.;
	    hardr = 1e-4;
	}
	hard = dhards * dkapp + hardr;
/*---------------------------------------------- gradient components ---*/
	grad[0] = grad[0];
	grad[1] = grad[1];
	grad[2] = grad[2];
	grad[3] *= 2.;

	for (i=0; i<4; i++) hydr[i] = 1.;
	hydr[3] = 0.;
/*-----------------------------------------------------------------------|
|        inverse modified constitutive tensor (compliance)               |
|-----------------------------------------------------------------------*/
	w1_init44(diag, 0.);
	w1_init44(hinv, 0.);
	w1_init44(q   , 0.);
	w1_init44(qn  , 0.);
	w1_init44(dm  , 0.);

/* 	DEFINE MODAL MATRIX Q */

	q[0][0] = ro16;
	q[1][0] = ro16;
	q[2][0] = ro16 * -2.;
	q[0][1] =-ro12;
	q[1][1] = ro12;
	q[0][2] = ro13;
	q[1][2] = ro13;
	q[2][2] = ro13;
	q[3][3] = 1.;
/*-----------------------------------------------------------------------|
|        define inverse diagonal matrix diag:                            |
|        DIAG = INV( LAMC + LAMP*DLAM/SIGYM(4) + LAMD*DLAM/SIGYM(4)/9.0) |
|-----------------------------------------------------------------------*/
	fac1 = (*vnu + 1.) / *emod;
	val3 = (1. - *vnu * 2.) / (*vnu + 1.);
	fac2 = *dlam / sigym[3];
	fac3 = fac2 / 9;

	dum = 1. /(fac1 + fac2);
	diag[0][0] = dum;
	diag[1][1] = dum;
	diag[2][2] = 1. /(fac1 * val3 + fac3 * 3.);
	diag[3][3] = 1. /((fac1 + fac2) * 2.);
/*-----------------------------------------------------------------------|
|        compute matrix: H(INV) = Q * DIAG * QT                          |
|-----------------------------------------------------------------------*/
	w1_AxB_444(q, diag, qn);
	w1_AxBT_444(qn, q, hinv);
/*-----------------------------------------------------------------------|
|        form consistent modified stiffness matrix THETA(4X4):           |
|                                                                        |
|                         (H(INV) : GRAD) x (H(INV) : GRAD)              |
|         THETA = H(INV) + --------------------------------              |
|                         INV(FAC2) - GRAD : H(INV) : GRAD               |
|-----------------------------------------------------------------------*/
	w1_AxB_441(hinv, grad, grd);

	dconst = 0.;
	for (i = 0; i < 4; ++i) {
	    dconst += grad[i] * grd[i];
	}

	dum   =  1. / fac2 - dconst;
	dconst = 1. / dum;

	w1_init44(qn, 0.);
	w1_AxBT_414(grd, grd, qn);
	w1_AxFAC_44(qn, qn, dconst);
	w1_AaddB_44(hinv, qn, theta);
/*-----------------------------------------------------------------------|
|    form components of matrix A(2X2):                                   |
|        |-                 -|                                           |
|        |  A(1,1)   A(1,2)  |                                           |
|    A = |                   |   ----> INV(A)                            |
|        |  A(2,1)   A(2,2)  |                                           |
|        |-                 -|                                           |
|                                                                        |
|    where:                                                              |
|    A(1,1) = GRAD  : THETA : GRAD                                       |
|    A(1,2) = GRAD  : THETA : GRADK - DKF / DLAM                         |
|    A(2,1) = GRADK : THETA : GRAD  - DKF / DLAM = A(1,2)                |
|    A(2,2) = GRADK : THETA : GRADK - DKKF/ DLAM - HARD / DLAM**2        |
|-----------------------------------------------------------------------*/

	dkf = -fkapp * hardr;
	dkkf = -(hardi*hardi) * (1. - (quot*quot) / 9) / (
		sigym[3] * 9) - fkapp * dhards;

	for (i = 0; i < 4; i++) gradk[i] = 0.;
	for (i = 0; i < 4; i++) {
	    gradk[i] = -(hydr[i] - grad[i] * quot) * hardi / (sigym[3] * 9);
	}

	w1_AxB_441(theta, grad , dn1);
	w1_AxB_441(theta, gradk, dn2);

	const1 = 0.;
	for (i = 0; i < 4; i++) {
	    const1 += grad[i] * dn1[i];
	}
	a[0][0] = const1;

	const2 = 0.;
	for (i = 0; i < 4; i++) {
	    const2 += grad[i] * dn2[i];
	}
	a[0][1] = const2 - dkf / *dlam;
	a[1][0] = a[0][1];

	const3 = 0.;
	for (i = 0; i < 4; i++) {
	    const3 += gradk[i] * dn2[i];
	}
	a[1][1] = const3 - dkkf / *dlam - hard / (*dlam * *dlam);

	det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	ai[0][0] =  a[1][1] / det;
	ai[1][0] = -a[1][0] / det;
	ai[0][1] = -a[0][1] / det;
	ai[1][1] =  a[0][0] / det;
/*-----------------------------------------------------------------------|
|        form consistent 4x4 material matrix                             |
|-----------------------------------------------------------------------*/
        for (i=0; i<4; i++)for (j=0; j<4; j++) d[i][j] = 0.;
	w1_init44(dm, 0.);

	w1_init44(qn, 0.);
	w1_AxBT_414(dn1, dn1, qn);
	dum = -ai[0][0];
	w1_AxFAC_44(qn, qn, dum);
	w1_AaddB_44(dm, qn, dm);

	w1_init44(qn, 0.);
	w1_AxBT_414(dn1, dn2, qn);
	dum = -ai[0][1];
	w1_AxFAC_44(qn, qn, dum);
	w1_AaddB_44(dm, qn, dm);

	w1_init44(qn, 0.);
	w1_AxBT_414(dn2, dn1, qn);
	dum = -ai[1][0];
	w1_AxFAC_44(qn, qn, dum);
	w1_AaddB_44(dm, qn, dm);

	w1_init44(qn, 0.);
	w1_AxBT_414(dn2, dn2, qn);
	dum = -ai[1][1];
	w1_AxFAC_44(qn, qn, dum);
	w1_AaddB_44(dm, qn, dm);

	w1_AaddB_44(dm, theta, dm);


        d[0][0] = dm[0][0];
        d[1][0] = dm[1][0]; 
        d[2][0] = dm[3][0];
        d[3][0] = dm[2][0];

        d[0][1] = dm[0][1];
        d[1][1] = dm[1][1];
        d[2][1] = dm[3][1];
        d[3][1] = dm[2][1];

        d[0][2] = dm[0][3];
        d[1][2] = dm[1][3];
        d[2][2] = dm[3][3]; 
        d[3][2] = dm[2][3];

        d[0][3] = dm[0][2];
        d[1][3] = dm[1][2];
        d[2][3] = dm[3][2];
        d[3][3] = dm[2][2];
/*---------------------------------------------- rotational symmetry ---*/
    } else  {

	dserror("ROTATIONAL SYMMETRY FOR PLASTICITY NOT IMPLEMENTED");        
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1mplcap */
/*-----------------------------------------------------------------------|
| compute average crack spacing for element wall1     (model: concrete)  |
|-----------------------------------------------------------------------*/
void w1acrs(DOUBLE  hte,    /* section properties                       */
               INT  maxreb, /* max. number of rebars                    */
            DOUBLE *stress, /* element stresses                         */
            DOUBLE *angle,  /* angle of the principal concrete stress   */
            DOUBLE *reb_area,
            DOUBLE *reb_ang,
            DOUBLE *reb_so,
            DOUBLE *reb_ds,
            DOUBLE *reb_rgamma,
            DOUBLE *thick,  /* thickness ot the structure               */
            DOUBLE dia,     /* equivalent element length                */
            DOUBLE **xjm,
            DOUBLE *acrs)   /* average crack spacing                    */
{
/*----------------------------------------------------------------------*/
    static INT i, k, id;
    static DOUBLE dum, rad, fps, sps, aps; 
    static DOUBLE sig[8], dl, alphas, fac, area, alphar, so, ds, rho, rls;
    static DOUBLE dalpha, dls[2], dgamma;              
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_enter("w1acrs");
    #endif
/*-------------------------------------------- initialized variables ---*/
  for (i=0; i<4; i++)
  {
    sig[i]   = stress[i];
    sig[i+4] = 0.;
  }

  *thick = hte;
  *angle  = 0.;
  dl     = dia;

  rad   = atan(1.)/45.;
/*----------------------------------------------------------------------*/
  dum = sig[0]+sig[1]+sig[2];
  if (fabs(dum)>0.000001) 
  {/*w1acrs01*/
/*----------------------------------------------------------------------*/
    w1_mami(sig, &fps,&sps,&aps);
    sig[4] = fps;
    sig[5] = sps;
    sig[6] = aps;
    sig[7] = 0.;
    *angle = sig[6];
/*------------------------------------ compute average crack spacing ---*/      
    for (i=0; i<2; i++)
    {/*800*/
      alphas = sig[6] + 90. * (DOUBLE)(i);
      if(alphas<0.) alphas = alphas + 180.;
      fac = 0.;

      id = 0;           
      for (k=0; k<maxreb; k++)
      {
	area   = reb_area[k];
	if (area>0.) id++;
      }
	    
      if (id>0) 
      {
        if(id>1) id = 2;
        for (k=0; k<id; k++)
        {
          area   = reb_area[k];
          alphar = reb_ang[k];
          so     = reb_so[k];
          ds     = reb_ds[k];
          dgamma = reb_rgamma[k];
       
          rho    = area / *thick;
          rls    = 2.*(2.*so + ds/(dgamma*rho))/3.;
          dalpha = fabs(alphas - alphar) * rad;
       
          if (id==1) fac = 1. / rls;
          else if (id>1) fac += fabs(cos(dalpha))/rls;
          
        }
      }
       
  	  if (fabs(fac)<0.000001) fac = 1. / dia;
          dls[i] = 1. / fac;
  	  if (dls[i]>dia) dls[i] = dia;
    }/*800*/	
    dl = (dls[0] + dls[1]) * .5;
/*----------------------------------------------------------------------*/
  }/*w1acrs01*/
/*------------------------------- modified equivalent element length ---*/
   if(dia<dl) *acrs = dia;
   else       *acrs = dl;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1acrs */
/*-----------------------------------------------------------------------|
|    topic: calculate the gradients of the elastic predictor             |              
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1cpreds(DOUBLE *sigym, /* uniaxial predictor yield stress         */
              DOUBLE *alpha, /* factor for the first invariants         */
              DOUBLE *devsig,/* deviatoric predictor stresses           */
              DOUBLE  hyd,   /* 1. invariant of the predictor stresses  */
              DOUBLE *grad)  /* total gradient                          */
{
/*----------------------------------------------------------------------*/
    static INT i;
    static DOUBLE dn[4], fact[4], dcom[4];
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_enter("w1cpreds");
    #endif
/*-------------------------------------------- initialized variables ---*/
    for (i = 0; i < 3; i++) fact[i]  = 1.;
    fact[3]  = 0.;
/*--------------------------------------- components ot the gradient ---*/
    for (i=0;i<4;i++)
    {
      dn[i]   = devsig[i] / sigym[3];
      dcom[i] = fact[i] * (hyd-alpha[3]*sigym[2])/(sigym[3]*9.) ;
      grad[i] = dn[i] + dcom[i];
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1cpreds */
/*----------------------------------------------------------------------*/






/*----------------------------------------------------------------------*
 |   WALL1 : array = value                                              |
 *----------------------------------------------------------------------*/
void w1_init44(DOUBLE a[4][4],DOUBLE val)
{
/*----------------------------------------------------------------------*/
INT i,j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_init44");
#endif
/*----------------------------------------------------------------------*/
      for (i=0; i<4; i++)for (j=0; j<4; j++) a[i][j] = val;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_init44 */
/*----------------------------------------------------------------------*
 |   WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                        |
 *----------------------------------------------------------------------*/
void w1_AxB_444(DOUBLE a[4][4],DOUBLE b[4][4],DOUBLE r[4][4])
{
/*----------------------------------------------------------------------*/
INT i,j,k;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AxB_444");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)for (j=0; j<4; j++)  r[i][j] = 0.;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) 
        r[i][j] += a[i][k]* b[k][j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AxB_444 */
/*----------------------------------------------------------------------*
 |   WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                        |
 *----------------------------------------------------------------------*/
void w1_AxBT_444(DOUBLE a[4][4],DOUBLE b[4][4],DOUBLE r[4][4])
{
/*----------------------------------------------------------------------*/
INT i,j,k;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AxBT_444");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)for (j=0; j<4; j++)  r[i][j] = 0.;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) 
        r[i][j] += a[i][k]* b[j][k];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AxBT_444 */
/*----------------------------------------------------------------------*
 |   WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                        |
 *----------------------------------------------------------------------*/
void w1_AxB_441(DOUBLE a[4][4],DOUBLE b[4],DOUBLE r[4])
{
/*----------------------------------------------------------------------*/
INT i,k;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AxB_441");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)  r[i] = 0.;

  for (i=0; i<4; i++)
      for (k=0; k<4; k++) 
        r[i] += a[i][k]* b[k];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AxB_441 */
/*----------------------------------------------------------------------*
 |   WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                        |
 *----------------------------------------------------------------------*/
void w1_AxBT_414(DOUBLE a[4],DOUBLE b[4],DOUBLE r[4][4])
{
/*----------------------------------------------------------------------*/
INT i,j,k;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AxBT_414");
#endif
/*----------------------------------------------------------------------*/

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 

        r[i][j] = a[i]* b[j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AxBT_414 */
/*----------------------------------------------------------------------*
 |   WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                        |
 *----------------------------------------------------------------------*/
void w1_AxFAC_44(DOUBLE a[4][4], DOUBLE r[4][4], DOUBLE fac)
{
/*----------------------------------------------------------------------*/
INT i,j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AxBT_414");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
        r[i][j] = a[i][j]* fac;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AxBT_414 */
/*----------------------------------------------------------------------*
 |   WALL1 : ADDITION   R = A+B   (FULL MATRICES)                       |
 *----------------------------------------------------------------------*/
void w1_AaddB_44(DOUBLE a[4][4], DOUBLE b[4][4], DOUBLE r[4][4])
{
/*----------------------------------------------------------------------*/
INT i,j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_AaddB_44");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
        r[i][j] = a[i][j] + b[i][j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_AaddB_44 */
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - plastic - rebar                al 9/01|
 *----------------------------------------------------------------------*/
void w1_mat_rebar(ELEMENT   *ele,
                  MATERIAL  *mat, 
                  DOUBLE   **bop, /*derivative operator */
                  DOUBLE    *gop, /*derivative operator */
                  DOUBLE    *alpha, /*derivative operator */
                  DOUBLE   **xjm, /* jacobian matrix    */
                  DOUBLE *stress,       
                  DOUBLE     **d,   
                  INT         ip,
                  INT       lanz,
                  INT     istore)
{
/*----------------------------------------------------------------------*/
INT    yip, ryip, iupd, nstiff;
DOUBLE sig, eps, epstn, rad, alfr, areb, den, prod, c, ac, arad;
DOUBLE x1r, x2r, x1s, x2s;
DOUBLE straino, strainrb, stresso, stressrb, dstrain;
DOUBLE ca, sa, epste, eh, sigy, stiff1; 
DOUBLE rsig, reps, repstn, rarea, alfrr, alfpri, ca2, sa2, emod, emod1, sigyv,
factor;  
DOUBLE disd[5];
DOUBLE strain[4];
DOUBLE cappaet, cappaut, angle, thick, fbd, epsy, dalpha, ecappaet;
DOUBLE eccappaet, dcappaet, roh, ecappaut, dcappaut;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_mat_rebar");
#endif
/*--------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,gop,alpha,ele->e.w1->wtype,disd) ;                  
/*------------------------------------- get actual strains -> strain ---*/
  w1_eps (disd,ele->e.w1->wtype,strain);
/*----------------------------- get old values -> stressrb, strainrb ---*/
  stressrb    = ele->e.w1->elewa[0].ipwa[ip].rsig[lanz-1];
  strainrb    = ele->e.w1->elewa[0].ipwa[ip].reps[lanz-1];
  epstn  = ele->e.w1->elewa[0].ipwa[ip].repstn[lanz-1];
  yip    = ele->e.w1->elewa[0].ipwa[ip].ryip[lanz-1];
/*---------------------------------- evaluate new stresses -> stress ---*/

/*-------------------------------------------------- set rebar angle ---*/
  rad  = atan(1.)/45.;
  alfr = mat->m.pl_epc->reb_ang[lanz-1];  /* input angle of rebar in degrees*/
  areb = mat->m.pl_epc->reb_area[lanz-1]; /* input area of rebar in degrees */

  if(alfr>=360.0) 
  {
    x1r = xjm[0][0];
    x2r = xjm[0][1];
    x1s = xjm[1][0];
    x2s = xjm[1][1];
 
    prod=x1r*x1s + x2r*x2s; 
    den=(x1r*x1r+x2r*x2r)*(x1s*x1s+x2s*x2s);
    den=sqrt(den);
    c=prod/den;
    ac=acos(c);
    arad=45./atan(1.);
    alfrr=ac;
    alfpri=alfrr*arad; /* angle to be printed in degrees */
  }
  else
  {
    alfpri=alfr;
    if(alfr>=90.0) alfr = alfr - 180.;
    rad=atan(1.)/45.;
    alfrr=alfr*rad;
  }
/*--------------------- define stress and strain from last load step ---*/
  straino  = strainrb;
  stresso  = stressrb;
/*------------------------------ determine strain in rebar direction ---*/
  ca       = cos(alfrr);                                       
  sa       = sin(alfrr);
  ca2      = ca*ca;
  sa2      = sa*sa;
  strainrb = ca2*strain[0] + sa2*strain[1] + ca*sa*strain[2];
  dstrain  = strainrb - straino;
  emod  = mat->m.pl_epc->reb_emod[lanz-1];
  emod1 = emod;
  eh    = mat->m.pl_epc->reb_hard[lanz-1];
  sigy  = mat->m.pl_epc->reb_sigy[lanz-1];
/*----------------------------------- considering tension stiffening ---*/
  stiff1 = 0.;
  epste  = sigy/emod;
  nstiff=mat->m.pl_epc->nstiff;
  if (fabs(strainrb)<=0.000000001) nstiff=0;

   
  if (nstiff==1) 
  {
    cappaet = ele->e.w1->elewa[0].ipwa[ip].stcappae ;
    cappaut = ele->e.w1->elewa[0].ipwa[ip].stcappaut;
    angle   = ele->e.w1->elewa[0].ipwa[ip].stalpha;  
    thick   = ele->e.w1->elewa[0].ipwa[ip].stthick;  
    fbd     = ele->e.w1->elewa[0].ipwa[ip].stfbd;    
   
    dalpha   = (angle - alfr) * rad;
    epsy     = strainrb; 
    ecappaet = 0.01 * cappaet;
    ecappaut = 0.01 * cappaut;
    dcappaet = (cos(dalpha));
    dcappaet *= dcappaet;
    dcappaet *= cappaet;

    dcappaut  = cos(dalpha);
    dcappaut *= dcappaut;
    dcappaut *= cappaut;
    if(ecappaet>dcappaet) cappaet  = ecappaet;
    else                  cappaet  = dcappaet;
    if(ecappaut>dcappaut) cappaut  = ecappaut;
    else                  cappaut  = dcappaut;
 
    roh      = areb / thick;
   
    if (epsy>cappaet && epsy<(3.*cappaut)) 
    {
      stiff1  = (3. * cappaut - cappaet);
      stiff1 *= stiff1;
      stiff1 *=roh;
      stiff1  = fbd / stiff1;
      stiff1 *= (6. * cappaut - 2. * epsy);
      if (stiff1>(2.*emod1)) stiff1 = 2.*emod1;
    }
  }

/*--------------------------------------------------------praedictor ---*/  
  iupd = 0;
  if (yip>0) 
  { 
    stressrb = stresso; 
    if (yip==2||yip==-2) emod = eh;
    yip    = -yip;
    iupd   = 1;
    goto local_stresses;
  }
/*---------------------------------------------------------------------*/   
/* 1. calculate praedictor stress from incremental strains using emod  */
  emod     = emod1 + stiff1;
  stressrb = dstrain * emod + stresso; 
/* 2. yield criterion */
  if (yip==1||yip==-1) sigyv = sigy;
  else                 sigyv = sigy+(fabs(epstn)-epste)*eh;
/*-------------------------------------------------- yield condition ---*/  	  
  if (sigyv>fabs(stressrb)) 
  { 
/*---------------------------------------------------------- elastic ---*/  
    yip = 1;
  }
/*---------------------------------------------------------- plastic ---*/  
  else
  {
    yip   = 2;
    epstn = strainrb;       

    if (stressrb>=0.) factor = 1.;
    else factor =-1.;

    stressrb = (sigy+(fabs(strainrb)-epste)*eh)*factor;
    emod     = eh;
  }          
/*----------------------------------------------------------------------*/
local_stresses:
/*--------------------------------------------------- local stresses ---*/
/*-------------------------------------------------- global stresses ---*/
      stress[0]=ca2*stressrb;
      stress[1]=sa2*stressrb;
      stress[3]=0.;
      stress[2]=ca*sa*stressrb;     
/*-------------------------------------------- local material matrix ---*/ 
/*------------------------------------------- global material matrix ---*/
      d[0][0]=ca2*ca2*emod;
      d[0][1]=ca2*sa2*emod;
      d[0][2]=ca2*ca*sa*emod;       
      d[0][3]=0.;
      
      d[1][0]=d[0][1];
      d[1][1]=sa2*sa2*emod;
      d[1][2]=ca*sa2*sa*emod;  
      d[1][3]=0.;
      
      d[2][0]=d[0][2];
      d[2][1]=d[1][2];
      d[2][2]=d[0][1];
      d[2][3]=0.;
      
      d[3][0]=d[0][3];
      d[3][1]=d[1][3];
      d[3][2]=d[2][3]; 
      d[3][3]=0.;
/*------------------------------------------------- store new values ---*/
  if(istore==1||iupd==1)
  {
    ele->e.w1->elewa[0].ipwa[ip].rsig[lanz-1]    = stressrb;
    ele->e.w1->elewa[0].ipwa[ip].reps[lanz-1]    = strainrb;
    ele->e.w1->elewa[0].ipwa[ip].repstn[lanz-1]  = epstn   ;
    ele->e.w1->elewa[0].ipwa[ip].ryip[lanz-1]    = yip     ;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_rebar */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/










