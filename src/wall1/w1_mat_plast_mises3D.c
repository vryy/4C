/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_mat_plast_mises_3D' which calculates 
 the constitutive matrix and forces using a von Mises-Plasticity 
 with a combined linear isotropic/kinematic hardening law
 the material formulation is 3D, so that the 2D-conditions from a plane
 calculation (Wall1-Element) have to be blown up ('w1mat_trans_up before) 
 calling 'mat_plast_mises_3D', which is a Element independent 3D-Routine 
 for calculating the 3D constitutive matrix. The calculated 3D-based 
 values have to be condesed back to either plane_stress or plane_strain
 conditions ('w1mat_trans_down').              
 (rotational symmetry is not implemented)   
 contains the routine 'mat_plast_mises_3D' which calcualtes the
 constitutive matrix - forces - von Mises - 3D which is element type 
 independent  
 contains the routine 'mat_linel3D' which calculates the local material
 law for 3D isotropic material   
 contains the routine 'yilcr_mises_lin' which calculates the yield 
 criterion for Mises Plasticity - 3D (combined linear hardening law)   
 contains the routine 'radi_mises_lin' which makes the radial return
 for von Mises Plasticity - 3D                  
 contains the routine 'mapl_mises_lin' which forms the elasto-plastic
 consistent tangent material tensor - 3D
 
*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"


/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculates the constitutive matrix for a von Mises-Plasticity

<pre>                                                             sh 08/02 
This routine calculates the constitutive matrix and forces using a 3D
von Mises-Plasticity with a combined linear isotropic/kinematic hardening 
law. Needed routines: w1mat_trans_up (2D->3D), w1mat_trans_down (3D->2D),
mat_plast_mises_3D (constitutive matrix, general 3D, element independent).
Works for plane_strain & plane_stress. Rotational symmetrie is not 
implemented yet.

</pre>
\param *d     DOUBLE     (o)   the constitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: w1_call_mat()

*----------------------------------------------------------------------*/
void w1_mat_plast_mises_3D(
                 DOUBLE     ym,      /*!< young's modulus              */
                 DOUBLE     pv,      /*!< poisson's ratio              */
                 DOUBLE     ALFAT,   /*!< temperature expansion factor */
                 DOUBLE     sigy,    /*!< yield stresse                */
                 DOUBLE     hard,    /*!< hardening modulus            */
                 DOUBLE     gf,      /*!< fracture energy              */
                 DOUBLE     betah,   /*!< controls the iso/kin hard.   */
                 ELEMENT   *ele,     /*!< actual element               */
                 WALL_TYPE  wtype,   /*!< plane stress/strain...       */
                 DOUBLE   **bop,     /*!< derivative operator          */
                 DOUBLE    *gop,
                 DOUBLE    *alpha,
                 INT        ip,      /*!< integration point Id         */
                 DOUBLE    *stress,  /*!< vector of stresses condensed */
                 DOUBLE   **d,       /*!< constitutive matrix          */
                 INT        istore,  /*!< controls storing of stresses */
                 INT        newval)  /*!< controls eval. of stresses   */
{
INT i;
DOUBLE stress3D[6];  /*actual stresses*/
DOUBLE strain3D[6];  /*actual strains from displacements*/
DOUBLE sig3D[6];     /*stresses from last update -> WA [4]->[6]*/
DOUBLE eps3D[6];     /*strains from last updat -> WA [4]->[6]*/
DOUBLE qn3D[6];      /*backstress vector from last update -> WA [4]->[6]*/
DOUBLE sig[4];       /*stresses from last update -> WA*/
DOUBLE eps[4];       /*strains from last update -> WA*/
DOUBLE sigi[4];      /*stress from last iteration step (not condensed*/
DOUBLE epsi[4];      /*strains from last iteration step (with e33)*/
DOUBLE di[4];        /*components d41,d42,d43,d44 of th const. tensor from*/
INT iupd=0;
INT yip,yipc;
DOUBLE epstn;
DOUBLE dia;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_mises_3D");
#endif
/*----------------------------------------------------------------------*/
    w1mat_trans_up(ym,
                   pv,
                   ele,
                   wtype,
                   bop,
                   gop,
                   alpha,
                   ip,
                   stress,   /*stress to be calculated [4]*/
                   stress3D, /*stress to be calculated [6]*/
                   strain3D,
                   sig3D,    /*11,22,33,12,23,13*/
                   eps3D,    /*11,22,33,12,23,13*/
                   qn3D,     /*11,22,33,12,23,13*/
                   newval);
    
    /*Werte aus Elementinformationen -> Material soll unabhaengig von Ele sein*/
    dia   = ele->e.w1->elewa[0].dia;
    yip   = ele->e.w1->elewa[0].ipwa[ip].yip;
    yipc  = yip;       /*make a copy of yip for correct update */
    epstn = ele->e.w1->elewa[0].ipwa[ip].epstn;
    /**/
   if(newval!=1)  /*Check of yield criteria with possible return*/
   {
    /*Aufruf der Materialroutine, allg. 3D*/
    mat_plast_mises_3D(ym,
                       pv,
                       ALFAT,
                       sigy,
                       hard,
                       gf,
                       betah,
                       stress3D, /*stress*/
                       d,        /*Material-Matrix to be calculated 3D*/
                       &iupd,    /*to be modified*/
/*zusaetzliche Parameter*/                       
                       &yip,     /*to be modified*/ 
                       &epstn,   /*to be modified*/
                       strain3D,
                       sig3D,
                       eps3D,
                       qn3D,
                       dia);
   
   }
    w1mat_trans_down(d,         /*Material-Matrix to be condensed 3D->2D*/
                     ele,
                     wtype,
                     ip,
                     yipc,
                     stress,     /*condensed stresses*/
                     sig3D,      /*11,22,33,12,23,13*/
                     eps3D,      /*11,22,33,12,23,13*/
                     stress3D,   /*to be condensed [6] -> [4]*/
                     strain3D,   /*to be condensed [6] -> [4]*/
                     qn3D);      /*to be condensed [6] -> [4]*/
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
   if(istore==1 || iupd==1)
   {
     for (i=0; i<4; i++)
     {
       ele->e.w1->elewa[0].ipwa[ip].sig[i] = stress3D[i];  
       ele->e.w1->elewa[0].ipwa[ip].eps[i] = strain3D[i];
       ele->e.w1->elewa[0].ipwa[ip].qn[ i] = qn3D[i] ;
     }
     ele->e.w1->elewa[0].ipwa[ip].epstn = epstn;
     ele->e.w1->elewa[0].ipwa[ip].yip   = yip  ;
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_mises_3D */



/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- von Mises - 3D sh 7/02|
 |                Test for 3D Materialmodell in Wall-Element            |
 *----------------------------------------------------------------------*/
void mat_plast_mises_3D(DOUBLE ym,      /*young's modulus*/
                        DOUBLE pv,      /*poisson's ration*/
                        DOUBLE ALFAT,   /*temperature expansion factor*/
                        DOUBLE sigy,    /*uniaxial yield stress*/
                        DOUBLE hard,    /*hardening modulus*/
                        DOUBLE gf,      /*fracture energy*/
                        DOUBLE betah,   
                        DOUBLE *stress, /*ele stress (-resultant) vector*/      
                        DOUBLE **d3D,   /*material matrix 3D*/
                        INT    *iupd,   /*controls update of new stresses to wa */
/*zusaetzliche Uebergabeparameter*/
                        INT    *yip,     /*from WA*/
                        DOUBLE *epstn,    /*from WA*/
                        DOUBLE strain3D[6],/*actual strains from displacements*/
                        DOUBLE sig3D[6],   /*stresses from WA*/
                        DOUBLE eps3D[6],   /*strains from WA*/
                        DOUBLE qn3D[6],    /*backstress vector from WA*/
                        DOUBLE dia)      /*internal length parameter from WA*/
{
/*----------------------------------------------------------------------*/
INT i,j;
INT isoft;
DOUBLE sum, ft;
DOUBLE delsig[6];
DOUBLE deleps[6];
DOUBLE tau[6];
DOUBLE dlam;
DOUBLE taupre[6];

#ifdef DEBUG 
dstrc_enter("mat_plast_mises_3D");
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
/*--------- 3D original global elastic matrix for current point -> D ---*/
   mat_linel3D(ym, pv, d3D); 
/*----------------------------------------------------------------------*/
/*Initialize*/
for (i=0; i<6; i++)
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
      if(*yip>0)
      {
        for (i=0; i<6; i++)
        {
          stress[i] = sig3D[i];
          tau[i]    = sig3D[i];
        }

        if(*yip==1)           /*elastic*/
        {
          *yip  = - *yip;
        }
        else                  /*plastic*/
        {
          dlam=0.;
          mapl_mises_lin(ym, hard, betah, sigy, pv, dia, tau, isoft, epstn, 
                         dlam, d3D);
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
  
  for (i=0; i<6; i++) deleps[i] = strain3D[i] - eps3D[i]; 
  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += d3D[i][j]*deleps[j];
  
  for (i=0; i<6; i++) tau[i] = sig3D[i] + delsig[i] - qn3D[i];
  /*store predictor stresses in taupre -> needed in mapl_mises_lin*/
  for (i=0; i<6; i++) taupre[i] = tau[i];
  
  sum = 0.0;
  for (i=0; i<6; i++) sum += sig3D[i] * delsig[i];
  if(sum<0.0)
  {
    for (i=0; i<6; i++) stress[i] = tau[i] + qn3D[i];
    *yip=1;
    goto end;
  }

/*---------yield condition - Mises 3D - linear hardening/softening   ---*/
  yilcr_mises_lin(ym, hard, betah, sigy, *epstn, isoft, dia, tau, &ft);
  
/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<EPS10) 
  {
    *yip = 1;
    for (i=0; i<6; i++) stress[i] = tau[i] + qn3D[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    *yip = 2;
    /* return -> new stresses, dlam, epstn */
    radi_mises_lin(ym, hard, betah, sigy, pv, dia, tau, qn3D, isoft, epstn, &dlam);
    /* algorithmic elasto plastic tangent material matrix */
    mapl_mises_lin(ym, hard, betah, sigy, pv, dia, taupre, isoft, epstn, dlam, d3D);
    
    for (i=0; i<6; i++) stress[i] = tau[i] + qn3D[i];
  }
/*----------------------------------------------------------------------*/
end:
/*Store values into Working-Array --> outside of this routine*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_plast_mises_3D */


/*----------------------------------------------------------------------*
 | program to establish local material law for 3D               sh 7/02 |
 | stress-strain law for isotropic material.                            |
 *----------------------------------------------------------------------*/
void mat_linel3D(DOUBLE ym,   /*young's modulus*/ 
                 DOUBLE pv,   /*poisson's ration*/
                 DOUBLE **d)  /*material matrix 3D*/
{
INT i,j;
DOUBLE d1,d2,d3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_linel3D");
#endif
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
d3=ym/((1.0 + pv)*2.0);
/*------------------------------------ set values in material-matrix ---*/
for (i=0; i<6; i++) for (j=0; j<6; j++) d[i][j] = 0.0;

d[0][0]=d1;
d[0][1]=d2;
d[0][2]=d2;

d[1][0]=d2;
d[1][1]=d1;
d[1][2]=d2;

d[2][0]=d2;
d[2][1]=d2;
d[2][2]=d1;

d[3][3]=d3;

d[4][4]=d3;

d[5][5]=d3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /*end of mat_linel3D*/




/*----------------------------------------------------------------------*
 |                                                        sh    7/02    |
 |   YIELD CRITERION  for  MISES PLASTICITY  -  3D                      |
 |           ----------------                                           |
 |        combined linear hardening law                                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
void yilcr_mises_lin(DOUBLE E, 
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
DOUBLE sx, sy, sz, sxy, sxz, syz;
DOUBLE sigym, hards, epstmax;
/*DOUBLE betah = 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("yilcr_mises_lin");
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
} /* end of yilcr_mises_lin */



/*----------------------------------------------------------------------*
 |                                                        sh    8/02    |
 |   radial return for elements with mises material model               |
 |   3D-Formulation with combined linear iso/kin. hardening             |
 |   --> see Simo & Hughes p.120f (plane_strain -> quasi 3D             |
 *----------------------------------------------------------------------*/
void radi_mises_lin(DOUBLE e, 
                    DOUBLE eh,
                    DOUBLE betah,
                    DOUBLE sigy,
                    DOUBLE vnu,
                    DOUBLE dia,
                    DOUBLE *sigma,
                    DOUBLE *qn,
                    INT    isoft,
                    DOUBLE *epstn,
                    DOUBLE *dlam)
{
/*----------------------------------------------------------------------*/
INT i;
INT isoft1 = 0;
INT nsoft  = 1;
DOUBLE half, ro23, q13, xsi1, xsi2, xsi3, xsi4, xsi5, xsi6, hard2, hards;
DOUBLE esig, epst, epstmax; 
DOUBLE J2, stps, g;
/*DOUBLE betah = 1.0;
/**/
DOUBLE fac;
DOUBLE x[6];
DOUBLE ft;
/**/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("radi_mises_lin");
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
} /* end of radi_mises_lin */

/*----------------------------------------------------------------------*
 |                                                        sh    8/02    |
 |   forms the elasto-plastic consistent tangent material tensor        |
 |   for combined, linear isotropic and kinematic hardening             |
 |                                                                      |
 |   see Simo & Huges: plane_strain: S.124 -> 3D                        |
 *----------------------------------------------------------------------*/
void mapl_mises_lin(DOUBLE e, 
                    DOUBLE eh,
                    DOUBLE betah,
                    DOUBLE sigy,
                    DOUBLE vnu,
                    DOUBLE dia,
                    DOUBLE *tau,   /*Praediktorspannungen*/
                    INT    isoft,
                    DOUBLE *epstn,
                    DOUBLE dlam,
                    DOUBLE **d)
{
/*----------------------------------------------------------------------*/
INT i, j;
INT nsoft  = 1;
DOUBLE x[6];
DOUBLE xsi1, xsi2, xsi3, xsi4, xsi5, xsi6, hards, g, epstmax, dum;
DOUBLE k,fac,fac1,fac2;
DOUBLE stps,J2;
DOUBLE A[6][6] = {1.,1.,1.,0.,0.,0.,
                  1.,1.,1.,0.,0.,0.,
                  1.,1.,1.,0.,0.,0.,
                  0.,0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,0.};

DOUBLE B[6][6] = { 2./3.,-1./3.,-1./3., 0.  , 0.  , 0. ,
                  -1./3., 2./3.,-1./3., 0.  , 0.  , 0. ,
                  -1./3.,-1./3., 2./3., 0.  , 0.  , 0. ,
                   0.   , 0.   , 0.   ,1./2., 0.  , 0. ,
                   0.   , 0.   , 0.   , 0.  ,1./2., 0. ,
                   0.   , 0.   , 0.   , 0.  , 0.  ,1./2.};
/*DOUBLE betah = 1.;
                  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mapl_mises_lin1");
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
} /* end of mapl_mises_lin */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/

