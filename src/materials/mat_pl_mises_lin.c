/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - mat_pl_mises_lin_main: which calculates the constitutive matrix for a
                          plasticity model based on the von Mises yield
                          criterion, with a combined, isotropic and kinematic
                          linear hardening law (see Simo & Hughes).
                          This routine is formulated in cartesian coordinate
                          system, general 3D with the sorting [11,22,33,12,23,13]
 - mat_pl_mises_lin_yilcr: which calculates the von Mises yield criterion
 - mat_pl_mises_lin_radi: which performs the radial return algorithm
 - mat_pl_mises_lin_mapl: which calculates the elasto-plastic consistent
                          material matrix
*----------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

/*! 
\addtogroup MAT 
*//*! @{ (documentation module open)*/ 


/*!----------------------------------------------------------------------
\brief consitutive matrix for von Mises-Plasticity Model with linear 
       hardening                                  

<pre>                                                            sh 03/03
This routine calculates the constitutive matrix and forces for a plasticity
model, based on the von Mises yield criterion with a combined isotropic
and kinematic, linear hardening law. (see Simo & Hughes)
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    ym      (i)  young's modulus 
\param  DOUBLE    pv      (i)  poisson's ration               
\param  DOUBLE    sigy    (i)  uniaxial yield stress         
\param  DOUBLE    hard    (i)  hardening modulus                 
\param  DOUBLE    gf      (i)  fracture energy
\param  DOUBLE    betah   (i)  controls the isotropic/kinematic hardening
\param  DOUBLE   *stress  (o)  vector of stresses [11,22,33,12,23,13]
\param  DOUBLE   *strain  (i)  actual strains from displacements  [11,22,33,12,23,13]
\param  DOUBLE  **d       (o)  constitutive matrix          
\param  INT      *iupd   (i/o) controls update of new stresses to wa         
\param  INT      *yip    (i/o) flag if global predictor step an if last step was plastic/elastic                 
\param  DOUBLE   *epstn  (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE    sig[6]  (i)  stresses from WA  [11,22,33,12,23,13]
\param  DOUBLE    eps[6]  (i)  strains from WA  [11,22,33,12,23,13]
\param  DOUBLE    qn[6]   (i)  backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE    dia     (i)  internal length parameter from WA

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_mat_plast_mises()     [s9_mat_plast_mises.c]

*----------------------------------------------------------------------*/
void mat_pl_mises_lin_main(
             DOUBLE   ym,       /*!< young's modulus */
             DOUBLE   pv,       /*!< poisson's ration */
             DOUBLE   sigy,     /*!< uniaxial yield stress */
             DOUBLE   hard,     /*!< hardening modulus */
             DOUBLE   gf,       /*!< fracture energy */
             DOUBLE   betah,    /*!< controls the iso/kin hardening */
             DOUBLE  *stress,   /*!< ele stress (-resultant) vector */      
             DOUBLE   strain[6],/*!< actual strains from displacements */
             DOUBLE **d,        /*!< material matrix 3D */
             INT     *iupd,     /*!< controls update of new stresses to wa */
             INT     *yip,      /*!< from WA */
             DOUBLE  *epstn,    /*!< from WA */
             DOUBLE   sig[6],   /*!< stresses from WA */
             DOUBLE   eps[6],   /*!< strains from WA */
             DOUBLE   qn[6],    /*!< backstress vector from WA */
             DOUBLE   dia)      /*!< internal length parameter from WA */
{
/*----------------------------------------------------------------------*/
INT    i,j;
INT    isoft;
DOUBLE sum, ft;
DOUBLE delsig[6];
DOUBLE deleps[6];
DOUBLE tau[6];
DOUBLE dlam;
DOUBLE taupre[6];

#ifdef DEBUG 
dstrc_enter("mat_pl_mises_lin_main");
#endif
/*----------------------------------------------------------------------*/
  *iupd=0;
/*---------------------------------------------- hardening/softening ---*/
 if(fabs(gf)>EPS10)
 {
   isoft= 1;
   hard =gf;
 }
 else isoft=0;

/*----- get linear elastic isotropic material tensor -------------------*/
   mat_el_iso(ym, pv, d); 

/*-----------------------------------------------------------------------|
|     YIP > 0  STRESSES ARE AVAILABLE FROM LAST UPDATE                   |
|         = 1  E L A S T I C                                             |
|         = 2  P L A S T I C                                             |
|     UPDATE FLAG MUST SET TO STORE CHANGE OF PARAMETER YIP              |
|     NO CHANGES HAVE BEEN MADE ON STRESS STATE                          |
|-----------------------------------------------------------------------*/
      if(*yip>0)
      {
        for (i=0; i<6; i++)
        {
          stress[i] = sig[i];
          tau[i]    = sig[i];
        }

        if(*yip==1)           /*elastic*/
        {
          *yip  = - *yip;
        }
        else                  /*plastic*/
        {
          dlam=0.;
          mat_pl_mises_lin_mapl(ym, hard, sigy, pv, dia, tau, 
                                isoft, epstn, dlam, d);
          *yip= - *yip;
        }
        *iupd=1;
        goto end;
      }
/*-----------------------------------------------------------------------|
|   1. CALCULATE INCREMENTAL STRAINS     DELEPS                          |
|   2. CALCULATE STRESS INCREMENT ASSUMING ELASTIC BEHAVIOUR             |
|   3. CALCULATE TOTAL STRESS                                            |
|   4. CHECK STRESS DEVIATOR AGAINST CURRENT YIELD SURFACE               |
|-----------------------------------------------------------------------*/
  
  for (i=0; i<6; i++) deleps[i] = strain[i] - eps[i]; 

  deleps[3] = 2. * deleps[3];
  deleps[4] = 2. * deleps[4];
  deleps[5] = 2. * deleps[5];

  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += d[i][j]*deleps[j];
  
  for (i=0; i<6; i++) tau[i] = sig[i] + delsig[i] - qn[i];
  /*store predictor stresses in taupre -> needed in mapl_mises_lin*/
  for (i=0; i<6; i++) taupre[i] = tau[i];
  
  sum = 0.0;
  for (i=0; i<6; i++) sum += sig[i] * delsig[i];
  if(sum<0.0)
  {
    for (i=0; i<6; i++) stress[i] = tau[i] + qn[i];
    *yip=1;
    goto end;
  }

/*---------yield condition - Mises 3D - linear hardening/softening   ---*/
  mat_pl_mises_lin_yilcr(ym, hard, betah, sigy, *epstn, isoft, dia, tau, &ft);
  
/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<EPS10) 
  {
    *yip = 1;
    for (i=0; i<6; i++) stress[i] = tau[i] + qn[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    *yip = 2;
    /* return -> new stresses, dlam, epstn */
    mat_pl_mises_lin_radi(ym, hard, betah, sigy, pv, dia, tau, 
                          qn, isoft, epstn, &dlam);
    /* algorithmic elasto plastic tangent material matrix */
    mat_pl_mises_lin_mapl(ym, hard, sigy, pv, dia, taupre, 
                          isoft, epstn, dlam, d);
    
    for (i=0; i<6; i++) stress[i] = tau[i] + qn[i];
  }
/*----------------------------------------------------------------------*/
end:
/*Store values into Working-Array --> outside of this routine*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_mises_lin_main */





/*!----------------------------------------------------------------------
\brief von Mises yield criterion

<pre>                                                            sh 07/02
This routine calculates the von mises yield criterion, formulated with
a combined linear hardening law
</pre>
\param  DOUBLE    E       (i)  young's modulus  
\param  DOUBLE    Eh      (i)  hardening modulus             
\param  DOUBLE    betah   (i)  controls the isotropic/kinematic hardening
\param  DOUBLE    sigy    (i)  uniaxial yield stress 
\param  DOUBLE   *ft      (o)  yield criterion 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_mises_lin_main()  [mat_pl_mises_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_mises_lin_yilcr(DOUBLE  E, 
                            DOUBLE  Eh,
                            DOUBLE  betah,
                            DOUBLE  sigy,
                            DOUBLE  epstn,
                            INT     isoft,
                            DOUBLE  dia,
                            DOUBLE *tau,
                            DOUBLE *ft)
{
/*----------------------------------------------------------------------*/
DOUBLE sx, sy, sz, sxy, sxz, syz;
DOUBLE sigym, hards, epstmax;
/*DOUBLE betah = 1.;*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_mises_lin_yilcr");
#endif
/*----------------------------------------------------------------------*/
  sx  = tau[0];
  sy  = tau[1];
  sz  = tau[2];
  sxy = tau[3];
  syz = tau[4];
  sxz = tau[5];
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
  *ft = sqrt(sx*sx + sy*sy + sz*sz - sx*sy - sy*sz - sx*sz
           + 3.0 * (sxy*sxy) + 3.0 * (sxz*sxz) + 3.0 * (syz*syz)) 
           - sigym;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_mises_lin_yilcr */



/*!----------------------------------------------------------------------
\brief return algorithm for mises material model

<pre>                                                            sh 08/02
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for the mises material model with combined 
linear iso/kin. hardenig law -> see Simo & Hughes p.120f 
(plane strain -> quasi 3D)
</pre>
\param  DOUBLE    e        (i)  young's modulus  
\param  DOUBLE    eh       (i)  hardening modulus             
\param  DOUBLE    betah    (i)  controls the isotropic/kinematic hardening
\param  DOUBLE    sigy     (i)  uniaxial yield stress 
\param  DOUBLE    vnu      (i)  poisson's ration  
\param  DOUBLE    dia      (i)  internal length parameter from WA
\param  DOUBLE   *sigma   (i/o) trial stresses to be projected  
\param  DOUBLE   *qn       (i)  backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE   *epstn   (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE   *dlam     (o)   plastic multiplier 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_mises_lin_main()  [mat_pl_mises_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_mises_lin_radi(DOUBLE  e, 
                           DOUBLE  eh,
                           DOUBLE  betah,
                           DOUBLE  sigy,
                           DOUBLE  vnu,
                           DOUBLE  dia,
                           DOUBLE *sigma,
                           DOUBLE *qn,
                           INT     isoft,
                           DOUBLE *epstn,
                           DOUBLE *dlam)
{
/*----------------------------------------------------------------------*/
INT    i;
DOUBLE half, ro23, q13, xsi1, xsi2, xsi3, xsi4, xsi5, xsi6, hard2, hards;
DOUBLE esig, epstmax; 
DOUBLE J2, stps, g;
/**/
DOUBLE fac;
DOUBLE x[6];
DOUBLE ft;
/**/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_mises_lin_radi");
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
/*--------------------------------------------- initialize dlam --------*/
   *dlam = 0.;
/*----------------------------------------------------------------------*/
	xsi1 = sigma[0];
	xsi2 = sigma[1];
	xsi3 = sigma[2];
	xsi4 = sigma[3];
      xsi5 = sigma[4];
      xsi6 = sigma[5];
/*---------------------------------- increment of plastic multiplier ---*/
   J2 = 1./3. * (xsi1*xsi1 + xsi2*xsi2 + xsi3*xsi3 - 
                 xsi1*xsi2 - xsi1*xsi3 - xsi2*xsi3) 
                 + xsi4*xsi4 + xsi5*xsi5 + xsi6*xsi6;
   stps = sqrt(2.* J2);    /*norm of devaiatoric praedictor stresses*/     

/*-------- direct evaluation of dlam  see Simo & Hughes p.120f ---------*/
/*------------- only possible for linear hardening ---------------------*/

   esig = ro23 * (sigy + betah * hards * *epstn);
   ft = sqrt (2.*J2) - esig;
   *dlam = ft/(2.*g + (2./3.) * hards);

/*-------- update of uniaxial plastic strain, stresses & backstress ----*/
   *epstn = *epstn + ro23 * *dlam;

/*-------- components of gradient n ------------------------------------*/
   x[0]= (2. * xsi1 - xsi2 - xsi3 ) / 3. / stps;
   x[1]= (2. * xsi2 - xsi1 - xsi3 ) / 3. / stps;
   x[2]= (2. * xsi3 - xsi1 - xsi2 ) / 3. / stps;
   x[3]= xsi4 / stps;
   x[4]= xsi5 / stps;
   x[5]= xsi6 / stps;

/*-- update of stresses and backstress   -------------------------------*/
/*-- see equation (3.3.3) & (3.3.1) in Simo & Hughes on p.121 ----------*/
   for(i=0; i<6; i++) 
   {
     sigma[i] = sigma[i] - (2.*g + ro23*fac)* *dlam*x[i];
     qn[i]    = qn[i]    + ro23*fac * *dlam * x[i];
   }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_mises_lin_radi */

/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 von Mises plasticity                                   

<pre>                                                            sh 08/02
This routine calculates the elasto-plastic consistent material tangent
for the von Mises plasticity with combined, linear iso/kin. hardening
law. -> see Simo & Hughes: plane_strain: p.124 
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    e        (i)  young's modulus  
\param  DOUBLE    eh       (i)  hardening modulus             
\param  DOUBLE    sigy     (i)  uniaxial yield stress 
\param  DOUBLE    vnu      (i)  poisson's ration  
\param  DOUBLE    dia      (i)  internal length parameter from WA
\param  DOUBLE   *tau      (i)  trial stresses   
\param  DOUBLE   *qn       (i)  backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE   *epstn   (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE   *dlam     (o)  plastic multiplier 
\param  DOUBLE  **d       (i/o) material matrix to be calculated             

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_mises_lin_mapl(DOUBLE   e, 
                           DOUBLE   eh,
                           DOUBLE   sigy,
                           DOUBLE   vnu,
                           DOUBLE   dia,
                           DOUBLE  *tau,   /*Praediktorspannungen*/
                           INT      isoft,
                           DOUBLE  *epstn,
                           DOUBLE   dlam,
                           DOUBLE **d)
{
/*----------------------------------------------------------------------*/
INT    i, j;
INT    nsoft  = 1;
DOUBLE x[6];
DOUBLE xsi1, xsi2, xsi3, xsi4, xsi5, xsi6, hards, g, epstmax, dum;
DOUBLE k,fac,fac1,fac2;
DOUBLE stps,J2;

DOUBLE A[6][6] = {{1.,1.,1.,0.,0.,0.},  
                  {1.,1.,1.,0.,0.,0.},  
                  {1.,1.,1.,0.,0.,0.},  
                  {0.,0.,0.,0.,0.,0.},  
                  {0.,0.,0.,0.,0.,0.},  
                  {0.,0.,0.,0.,0.,0.}}; 

DOUBLE B[6][6] = {{ 2./3.,-1./3.,-1./3., 0.  , 0.  , 0.  },
                  {-1./3., 2./3.,-1./3., 0.  , 0.  , 0.  },
                  {-1./3.,-1./3., 2./3., 0.  , 0.  , 0.  },
                  { 0.   , 0.   , 0.   ,1./2., 0.  , 0.  },
                  { 0.   , 0.   , 0.   , 0.  ,1./2., 0.  },
                  { 0.   , 0.   , 0.   , 0.  , 0.  ,1./2.}};
                  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_mises_lin_mapl");
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
/*===================================================== plane strain ===*/
   xsi1 = tau[0];
   xsi2 = tau[1];
   xsi3 = tau[2];
   xsi4 = tau[3];
   xsi5 = tau[4];
   xsi6 = tau[5];
/*------- Grad n  -> Normiert  -----------------------------------------*/
   J2 = 1./3. * (xsi1*xsi1 + xsi2*xsi2 + xsi3*xsi3 - 
                 xsi1*xsi2 - xsi1*xsi3 - xsi2*xsi3)
                 + xsi4*xsi4 + xsi5*xsi5 + xsi6*xsi6;
   stps = sqrt(2.* J2); /*norm of devaiatoric praedictor stresses*/     

   x[0]= (2. * xsi1 - xsi2 - xsi3 ) / 3. / stps;
   x[1]= (2. * xsi2 - xsi1 - xsi3 ) / 3. / stps;
   x[2]= (2. * xsi3 - xsi1 - xsi2 ) / 3. / stps;
   x[3]= xsi4 / stps;
   x[4]= xsi5 / stps;
   x[5]= xsi6 / stps;
/*------- Vorfaktoren  -------------------------------------------------*/
   fac  = (2.*g* dlam)/ (stps);
   fac1 = 1. - fac;
   fac2 = 1./(1.+ hards/(3.*g)) - fac ;      
/*-------- Cep_alg nach Simo & Hughes  S.124 ---------------------------*/

   for (i=0; i<6; i++) for (j=0; j<6; j++) 
   {
     d[i][j] = k*A[i][j]+ 2.*g*fac1*B[i][j] - 2.*g*fac2*x[i]*x[j];
   }
            
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_mises_lin_mapl */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/

