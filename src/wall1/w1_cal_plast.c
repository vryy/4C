#include "../headers/standardtypes.h"
#include "wall1.h"
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
             double *tau,
             double *ft)
{
/*----------------------------------------------------------------------*/
double sx, sy, sxy, sigym, hards, epstmax;
double betah = 1.;
double dia = 1.414213562373095;
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
double dia = 1.414213562373095;
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
  case plain_stress:
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
  case plain_strain:
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
double dia = 1.414213562373095;
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
  case plain_stress:
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
  case plain_strain:
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
  case plain_stress:
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
  case plain_strain:
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
