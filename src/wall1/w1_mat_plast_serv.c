#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   TOPIC : WALL1 - YIELD CRITERION PLANE_STRESS                       |
 |           -   -   -- -  --                                           |
 |           YIELD - CRITERION:  PHI > 0 DRUCKER - PRAGER               |
 |                               PHI = 0 MISES                          |
 |           FOR MATERIAL MODEL 3 - PLASTICITY                          |
 *----------------------------------------------------------------------*/
void w1yilcr(double E, 
             double Eh,
             double sigy,
             double epstn,
             int    isoft,
             double dia,
             double *tau,
             double *ft)
{
/*----------------------------------------------------------------------*/
double sx, sy, sxy, sigym, hards, epstmax;
double betah = 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1yilcr");
#endif
/*----------------------------------------------------------------------*/
  sx  = tau[0];
  sy  = tau[1];
  sxy = tau[2];
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
  *ft = sqrt(sx*sx - sx*sy + sy*sy+ 3.0 * (sxy*sxy)) - sigym;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1yilcr */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   radial return for elements with mises material model               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1radi(double e, 
            double eh,
            double sigy,
            double vnu,
            double dia,
            double *sigma,
            double *qn,
            int    isoft,
            double *epstn,
            double *dlam,
            WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
int i;
int isoft1 = 0;
int nsoft  = 1;
double half, ro23, q13, g, c1, c2, g1, xsi1, xsi2, xsi3, hard2, hards;
double f, f1, f2, f3, f4, dfdl, esig, epst, fr, det, epstmax, df, dfi; 
double dum, xsi4, stps, dlf, dlfi;
double hd11, hd21, hd12, hd22, hd33, hd44;
double hm11, hm21, hm12, hm22, hm33;
double dm11, dm21, dm12, dm22, dm33;
double betah = 1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1radi");
#endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
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
	f = f1 * f1 / 6. + f2 / (f3 * f3);
	fr = sqrt(f);

	epst = *epstn + *dlam * fr * ro23;

/*------------ isotrope Materialverfestigung bei von Mises Bedingung ---*/
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
/*-------------------------------------------------- Fliessbedingung ---*/
	f4 = esig * esig / 3.;
	f = f / 2. - f4;
/*------------------------------------ Ableitung der Fliessbedingung ---*/
	dfdl = -c1 * f1 * f1 / ((c1 * *dlam + 1.) * 6.) - c2 * f2 / 
              (f3 * f3 * f3);
	f4 = esig * (betah * hards) * ro23 * 2. / 3.;
	dfdl = (fr - f4 * *dlam) * dfdl / fr - f4 * fr;

	*dlam -= f / dfdl;
/*------------------------------------------------- Abbruchkriterium ---*/
	if ((dum = f / esig, fabs(dum)) > 1e-5) {
	    if (i > 30) dserror("i>30!");         
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
} /* end of w1radi */
/*----------------------------------------------------------------------*
 |                                                        al    9/01    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 *----------------------------------------------------------------------*/
void w1mapl(double e, 
            double eh,
            double sigy,
            double vnu,
            double dia,
            double *tau,
            int    isoft,
            double *epstn,
            double *dlam,
            double **d,
            WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
int i, j;
int nsoft  = 1;
double fkh, fkg, gamma1, gamma2, beta,x1,x2,x3,abeta,fact,x4,fact1;
double vect[4];
double dm[4][4];
double xsi1, xsi2, xsi3, hards, g, epstmax, df, det, dum, xsi4;
double d11, d21, d12, d22, d33;
double hm11, hm21, hm12, hm22, hm33, hm44;
double dm11, dm21, dm12, dm22, dm33;
double betah = 1.;
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
    g = e / (vnu * 2. + 2.);
    fkh = e / ((1. - vnu * 2.) * 3.);
    fkg = g * 2. / (((g + (1. - betah) * hards / 3.) * 2. * *dlam + 1.) * 3.);
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
	fact = 1. - betah * 2. * hards * *dlam / 3.;
	fact = betah * 2. * hards / 3. / fact;

	xsi1 = tau[0];
	xsi2 = tau[1];
	xsi3 = tau[3];
	xsi4 = tau[2];

	fact = fact * 2. * (xsi1 * xsi1 + xsi2 * xsi2 + xsi3 * xsi3 - xsi1 * 
		xsi2 - xsi1 * xsi3 - xsi2 * xsi3 + xsi4 * xsi4 * 3.) / 3.;

/*-----------------------------------------------------------------------|
|        inverse modified constitutive tensor (compliance)               |
|                 D  =  D(INV) + P * DLAM                                |
|        inverse material stiffness tensor     DM = INV (D(INV) + P*DLAM |
|           (HM11=HM22=HM33,HM12=HM23=HM21=HM32=HM13,H44)                |
|-----------------------------------------------------------------------*/
	hm11 = fkh + fkg * 2.;
	hm12 = fkh - fkg;
	hm44 = fkg * 1.5;
/*-----------------------------------------------------------------------|
|        product  P * TAU                                                |
|-----------------------------------------------------------------------*/
	x1 = (xsi1 * 2. - xsi2 - xsi3) / 3.;
	x2 = (-xsi1 + xsi2 * 2. - xsi3) / 3.;
	x3 = (-xsi1 - xsi2 + xsi3 * 2.) / 3.;
	x4 = xsi4 * 2.;
/*-----------------------------------------------------------------------|
|        product  HM * P * TAU                                           |
|-----------------------------------------------------------------------*/
	vect[0] = hm11 * x1 + hm12 * x2 + hm12 * x3;
	vect[1] = hm12 * x1 + hm11 * x2 + hm12 * x3;
	vect[2] = hm12 * x1 + hm12 * x2 + hm11 * x3;
	vect[3] = hm44 * x4;
/*-----------------------------------------------------------------------|
|        scalar product  TAU * P * HM * P * TAU                          |
|-----------------------------------------------------------------------*/
	fact1 = x1 * vect[0] + x2 * vect[1] + x3 * vect[2] + x4 * vect[3];
	fact1 = 1. / (fact1 + fact);
/*-----------------------------------------------------------------------|
|        form consistent 4x4 material matrix                             |
|-----------------------------------------------------------------------*/
	for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    dm[i][j] = 0.;
	     d[i][j] = 0.;
	}}
        dm[0][0] = hm11;
        dm[0][1] = hm12;
        dm[0][2] = hm12;
        dm[1][0] = hm12;
        dm[1][1] = hm11;
        dm[1][2] = hm12;
        dm[2][0] = hm12;
        dm[2][1] = hm12;
        dm[2][2] = hm11;
        dm[3][3] = hm44;

	for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    dm[i][j] -= fact1 * vect[i] * vect[j];
	}}

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
void w1yilcr_dp(double E, 
                double Eh,
                double phi,
                double sigy,
                double *sigym,
                double epstn,
                double *tau,
                double *ft)
{
/*----------------------------------------------------------------------*/
double coh, sx, sy, sxy, hards, alpha, yld;
double betah = 1.;
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
void w1radi_dp(double e, 
               double eh,
               double phi,
               double sigy,
               double vnu,
               double *sigma,
               double *qn,
               double *epstn,
               double *sigym,
               double *dlam,
               WALL_TYPE wtype)
{
/*----------------------------------------------------------------------*/
int i;
int isoft1 = 0;
int nsoft  = 1;
double half, ro23, q13, g, c1, c2, g1, xsi1, xsi2, xsi3, hard2, hards;
double f, f1, f2, f3, f4, dfdl, esig, epst, fr, det, epstmax, df, dfi; 
double dum, xsi4, stps, dlf, dlfi, coh, y, alpha, devinv, alph, fac;
double hd11, hd21, hd12, hd22, hd33, hd44;
double hm11, hm21, hm12, hm22, hm33;
double dm11, dm21, dm12, dm22, dm33;
double betah = 1.0;
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
            double    *funct_h,
            double   **deriv_h,
            double   **xjm_h)
{
int                 i,j,k;            /* some loopers */
int                 nir,nis;          /* num GP in r/s/t direction */
int                 lr, ls;           /* loopers over GP */
int                 iel;              /* numnp to this element */
int                 nd;
const int           numdf  = 2;
const int           numeps = 3;

double              fac;
double              e1,e2,e3;         /*GP-coords*/
double              facr,facs,fact;   /* weights at GP */
double              det, exp, dia;
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
void w1cpreva (double *epst,     /* equivalent uniaxial plastic strain  */
               double *sigym,    /* [0] uniaxial tension yield stress   */
                                 /* [1] yield stress "inverted cone"    */
                                 /* [2] uniaxial compr. yield stress    */
               double *alpha,    /* factor for the first invariants     */
	       double *hards,    /* hardening modulus                   */
               double *e,        /* young modulus                       */
               double *g,        /* shear modulus                       */
               double *vnu,      /* poisson's ratio                     */
               double *com,      /* bilk modulus                        */
               double *fcm,      /* compressive strenght                */
               double *gc,       /* compression fracture energy         */
               double *ftm,      /* tensile strenght                    */
               double *gt,       /* tensile fracture energy             */
               double *gamma1,   /* fitting factor yield function 1     */
               double *gamma2,   /* fitting factor yield function 2     */
               double *dfac,     /* damage factor                       */
               double *dia,      /* equivalent element length           */
	       double *acrs,     /* average crack spacing               */
               double *cappaet,  /* max. elastic tension strain         */
               double *cappaut,  /* tensile fracture strain             */
               double *cappae,   /* max. elastic compression strain     */
               double *cappauc,  /* compressive fracture strain         */
               double *sig3,     /* equivalent compressive stress       */
               double *fbd,      /* tension stiffening stress           */
               double **d,       /* elastic material matrix             */
               double *sig,      /* stresses from last update           */
               double *eps)      /* strains from last update            */
{                                
/*-------------------------------------------------- local variables ---*/
    static int i, j;
    static double dum;
    static double vect[4], dfac1, dfac2, epsn1, epsn2;
    static double sigym1, sigym2, q23, fac, dam, ro23, ro54, dia2, sig1;
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
void w1cradi (double *sigma,  /* elastic predictor projected onto yield surface */ 
              double *epstn,  /* equivalent uniaxial plastic strain             */ 
              double *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              int    yip,     /* stress state   1 =elastic   >=2 =plastic       */ 
              double *alpha,  /* factor for the first invariants                */ 
              double *ft,     /* yield condition                                */ 
              double *e,      /* young modulus                                  */ 
              double *g,      /* shear modulus                                  */ 
	      double *com,    /* bulk modulus                                   */ 
              double *sigym,  /* uniaxial predictor yield stress                */ 
              double *hards,  /* plastic modulus                                */ 
              double *sigy,   /* actual uniaxial yield stress                   */ 
              double *dn,     /* gradient components in deviatoric direction    */ 
              double *dcom,   /* gradient components in hdrostatic direction    */ 
              double *grad,   /* total gradient                                 */ 
              double *devsig, /* deviatoric predictor stresses                  */ 
              double *sm,     /* hydrostatic predictor stresses                 */ 
              double *fcm,    /* compressive strenght                           */ 
              double *gc,     /* compression fracture energy                    */ 
              double *ftm,    /* tensile strenght                               */ 
	      double *gt,     /* tensile fracture energy                        */ 
              double *gamma1, /* fitting factor yield function 1                */ 
              double *gamma2, /* fitting factor yield function 2                */ 
              double *dia,    /* equivalent element length                      */ 
              double *acrs)   /* average crack spacing                          */ 
{                             
/*----------------------------------------------------------------------*/
    double dum;               
/*----------------------------------------------------------------------*/
    static double dfdl, half, epst;
    static double f;
    static int i;
    static double q13, dalpha, cappae, fac, ro23, smi[4], cappauc, 
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
void w1mapl2(double *tau,      /* current stresses (local)              */
             double **d,       /* material matrix to be calculated      */
             double *dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* type of problem                       */
             double *alpha,    /* neigungswinkel der fliessflaechen     */
             double *emod,     /* elastizitaetsmodul                    */
             double *g,        /* schubmodul                            */
             double *com,      /* kompressionsmodul                     */
             double *betah,    /* factor for isotrop/kinemat. hardening */
             double *hards,    /* plastic hardeningmodulus              */
             double *dn,       /* gradient components in dev. direction */
             double *grad,     /* total gradient                        */
             double *dev)      /* norm of the dev. predictor stresses   */
{
/*----------------------------------------------------------------------*/
static double half, fact, vect[4], stps, fact1,dum;
static int i, j;
static double x[4];
static double stpsi, dm[4][4], hm[4][4], q13, q23; 
static double dalpha, add, fac, fkg, fki, ro23; 
static double xsi1, xsi2, xsi3, xsi4;
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
void w1cradms(double *sigma,  /* elastic predictor projected onto yield surface */ 
              double *epstn,  /* equivalent uniaxial plastic strain             */ 
              double *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              int     yip,    /* stress state   1 =elastic   >=2 =plastic       */ 
              double *alpha,  /* factor for the first invariants                */ 
              double *ft,     /* yield condition                                */ 
              double  e,      /* young modulus                                  */ 
              double  g,      /* shear modulus                                  */ 
	      double  com,    /* bulk modulus                                   */ 
              double *sigym,  /* uniaxial predictor yield stress                */ 
              double *hards,  /* plastic modulus                                */ 
              double  dn[3][4],  /* gradient components in deviatoric direction */ 
              double  dcom[3][4],/* gradient components in hdrostatic direction */ 
              double *devsig, /* deviatoric predictor stresses                  */ 
              double *sm,     /* hydrostatic predictor stresses                 */ 
              double  fcm,    /* compressive strenght                           */ 
              double  gc,     /* compression fracture energy                    */ 
              double  ftm,    /* tensile strenght                               */ 
	      double  gt,     /* tensile fracture energy                        */ 
              double  gamma1, /* fitting factor yield function 1                */ 
              double  gamma2, /* fitting factor yield function 2                */ 
              double  dia,    /* equivalent element length                      */ 
              double  acrs)   /* average crack spacing                          */ 
{                             
/*----------------------------------------------------------------------*/
    double dum;               
/*----------------------------------------------------------------------*/
    double dum1, dum2;

    static double half, soli[4]	, epst[2], sigy[3];
    static double f;
    static int i;
    static double dhard[2];
    static int kflag;
    static double c1, c2, alpha1, alpha2, df[4];
    static double q13, cappae, fac, dlambda[2], det, dev;
    static double cappaut, sig1, sig2, sig3;
    static double fti[2], hyd, ro23, cappauc, sol[4];
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
      void w1conver(double *dlam,       /*  incr. of plastic multiplier */             
                    double *alpha,      /*  factor for the first inv.   */             
                    double *hards,      /*  plastischer hardeningmodul  */             
                    double dn[3][4],    /*  gradient comp. in dev. dir. */  
                    double grad[3][4],  /*  total gradient              */             
                    double *dlamc,      /*  -|                          */             
                    double *alphac,     /*   |                          */             
                    double *hardsc,     /*   |-  converted arrays       */             
                    double dnc[2][4],   /*   |                          */             
                    double gradc[2][4]) /*  -|                          */             
{
/*----------------------------------------------------------------------*/
    static int i;
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
void w1maplg(double *tau,      /* current stresses (local)              */
             double **d,       /* material matrix to be calculated      */
             double  dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* type of problem                       */
             int     yip,      /* stress state  1=elastic 2=plastic     */
             double  emod,     /* elastizitaetsmodul                    */
             double  g,        /* schubmodul                            */
             double  com,      /* kompressionsmodul                     */
             double *hards,    /* plastic hardeningmodulus              */
             double  dn[2][4], /* gradient components in dev. direction */
             double  grad[2][4])/*norm of the dev. predictor stresses   */
{
/*----------------------------------------------------------------------*/
    double dum;

    static int i, j, k;
    static double fact[2][2], ghme[4][2], stps;
    static double e[2][2];
    static double e2[2][2], stpsi, dm[4][4];
    static double  hm[4][4], q13, q23, dalpha, facmin;
    static double fak, fkg, fki, ghm[4][2], det, ro23; 
    static double tra[4][4], xsi1, xsi2, xsi3, xsi4;
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
void w1radcap(double *sigma,  /* elastic predictor projected onto yield surface */ 
              double *epstn,  /* equivalent uniaxial plastic strain             */ 
              double *dlam,   /* incremental plastic multiplier                 */ 
              WALL_TYPE wtype,/* type of problem                                */ 
              double *alpha,  /* factor for the first invariants                */ 
              double *e,      /* young modulus                                  */ 
              double *g,      /* shear modulus                                  */ 
	      double *com,    /* bulk modulus                                   */ 
              double *sigym,  /* uniaxial predictor yield stress                */ 
              double *hards,  /* plastic modulus                                */ 
              double *grad,   /* total gradient                                 */ 
              double *devsig, /* deviatoric predictor stresses                  */ 
              double *sm,     /* hydrostatic predictor stresses                 */ 
              double *dev,     /* norm of the deviatoric predictor stresses     */ 
              double *hyd,     /* 1st invariant  of the predictor stresses      */ 
              double *hydn,    /* 1st invariant  of the new stresses            */ 
              double *fcm,    /* compressive strenght                           */ 
              double *gc,     /* compression fracture energy                    */ 
              double *ftm,    /* tensile strenght                               */ 
	      double *gt,     /* tensile fracture energy                        */ 
              double *gamma2, /* fitting factor yield function 2                */ 
              double *dia)    /* equivalent element length                      */ 
{
/*----------------------------------------------------------------------*/
    double dum1, dum2;               
/*----------------------------------------------------------------------*/
    static int i;
    
    static double dadk, dadl, dbdl, dbdk, half, fact[4], dcom[4]; 
    static double soli[2][2], epst, root;
    static double a, b, f;
    static double hardi, hardr, rdevc, rcomc, rdevn;
    static double depst, rcomn, dn[4], q13, q23, cappae, q29, fac;
    static double det, fti[2], ro23, smi[4], ro54, cappauc, sol[2][2];
    static double devsigi[4], fac1, fac2, sig3;
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
void w1mplcap(double *tau,      /* current stresses (local)              */
              double **d,       /* material matrix to be calculated      */
              double *dlam,     /* increment of plastic multiplier       */
              WALL_TYPE wtype,  /* type of problem                       */
              double *alpha,    /* factor for the first invariants       */
              double *emod,     /* young modulus                         */
              double *vnu,      /* poisson's ratio                       */
              double *hards,    /* plastic hardeningmodulus              */
              double *sigym,    /* [0] uniaxial tension yield stress     */
              double *grad,     /* total gradient                        */
              double *cappae,   /* max. elastic compression strain       */
              double *cappauc,  /* compressive fracture strain           */
              double *epst,     /* equivalent uniaxial plastic strain    */
              double *sig3)     /* equivalent compressive stress         */
{
/*----------------------------------------------------------------------*/
    double dum;               
/*----------------------------------------------------------------------*/
    static int i, j ,k;
    
    static double diag[4][4], half, hard, dkkf, hinv[4][4];
    static double hydr[4], fack1, quot, a[2][2];
    static double q[4][4], facth, hardi, gradk[4];
    static double dkapp, hardr, fkapp;
    static double theta[4][4];
    static double dconst, const1, const2, const3, ai[2][2];
    static double dm[4][4], q13, q23, qn[4][4], dhards;
    static double dn1[4], dn2[4], dkf, grd[4], det, sig, ro12;
    static double ro13, hyd, ro16, fac1, fac2, fac3, val3, xsi1;
    static double xsi2, xsi3, xsi4;
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
void w1acrs(double  hte,    /* section properties                       */
               int  maxreb, /* max. number of rebars                    */
            double *stress, /* element stresses                         */
            double *angle,  /* angle of the principal concrete stress   */
            double *reb_area,
            double *reb_ang,
            double *reb_so,
            double *reb_ds,
            double *reb_rgamma,
            double *thick,  /* thickness ot the structure               */
            double dia,     /* equivalent element length                */
            double **xjm,
            double *acrs)   /* average crack spacing                    */
{
/*----------------------------------------------------------------------*/
    static int i, k, id;
    static double dum, rad, fps, sps, aps; 
    static double sig[8], dl, alphas, fac, area, alphar, so, ds, rho, rls;
    static double dalpha, dls[2], dgamma;              
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
      alphas = sig[6] + 90. * (double)(i);
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
void w1cpreds(double *sigym, /* uniaxial predictor yield stress         */
              double *alpha, /* factor for the first invariants         */
              double *devsig,/* deviatoric predictor stresses           */
              double  hyd,   /* 1. invariant of the predictor stresses  */
              double *grad)  /* total gradient                          */
{
/*----------------------------------------------------------------------*/
    static int i;
    static double dn[4], fact[4], dcom[4];
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
void w1_init44(double a[4][4],double val)
{
/*----------------------------------------------------------------------*/
int i,j;
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
void w1_AxB_444(double a[4][4],double b[4][4],double r[4][4])
{
/*----------------------------------------------------------------------*/
int i,j,k;
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
void w1_AxBT_444(double a[4][4],double b[4][4],double r[4][4])
{
/*----------------------------------------------------------------------*/
int i,j,k;
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
void w1_AxB_441(double a[4][4],double b[4],double r[4])
{
/*----------------------------------------------------------------------*/
int i,k;
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
void w1_AxBT_414(double a[4],double b[4],double r[4][4])
{
/*----------------------------------------------------------------------*/
int i,j,k;
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
void w1_AxFAC_44(double a[4][4], double r[4][4], double fac)
{
/*----------------------------------------------------------------------*/
int i,j;
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
void w1_AaddB_44(double a[4][4], double b[4][4], double r[4][4])
{
/*----------------------------------------------------------------------*/
int i,j;
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
                  double   **bop, /*derivative operator */
                  double   **xjm, /* jacobian matrix    */
                  double *stress,       
                  double     **d,   
                  int         ip,
                  int       lanz,
                  int     istore)
{
/*----------------------------------------------------------------------*/
int    yip, ryip, iupd, nstiff;
double sig, eps, epstn, rad, alfr, areb, den, prod, c, ac, arad;
double x1r, x2r, x1s, x2s;
double straino, strainrb, stresso, stressrb, dstrain;
double ca, sa, epste, eh, sigy, stiff1; 
double rsig, reps, repstn, rarea, alfrr, alfpri, ca2, sa2, emod, emod1, sigyv,
factor;  
double disd[5];
double strain[4];
double cappaet, cappaut, angle, thick, fbd, epsy, dalpha, ecappaet;
double eccappaet, dcappaet, roh, ecappaut, dcappaut;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_mat_rebar");
#endif
/*--------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,ele->e.w1->wtype,disd) ;                  
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












