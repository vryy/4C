/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - mat_pl_dp_lin_main:      which calculates the constitutive matrix for 
                            a plasticity model based on the 'Drucker Prager' 
                            yield criterion, with a combined, isotropic and 
                            kinematic linear hardening law 
                            (see Dis. Menrath & Book: Simo & Hughes).
                            This routine is formulated in cartesian 
                            coordinate system, general 3D with the sorting 
                            [11,22,33,12,23,13]
 - mat_pl_dp_lin_radi:      which performs the radial return algorithm
 - mat_pl_dp_lin_mapl:      which calculates the elasto-plastic consistent
                            material matrix
 - mat_pl_dp_lin_radi_apex: which performs the radial return algorithm if
                            trial stresses are in the apex region
 - mat_pl_dp_lin_mapl_apex: which performs the radial return algorithm
                            for the apex region
 - mat_pl_dp_lin_preval:    calculates the gradient n (s/|s|) and the
                            norm of deviatoric stresses (|s|)

<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

/*! 
\addtogroup MAT 
*//*! @{ (documentation module open)*/ 


/*!----------------------------------------------------------------------
\brief consitutive matrix for 'Drucker Prager'-Plasticity Model with linear 
       hardening                                  

<pre>                                                            sh 09/03
This routine calculates the constitutive matrix and forces for a plasticity
model, based on the 'Drucker Prager' yield criterion with a combined isotropic
and kinematic, linear hardening law. (see Dis. Menrath & Simo)
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    ym        (i)  young's modulus 
\param  DOUBLE    pv        (i)  poisson's ration               
\param  DOUBLE    sigy      (i)  uniaxial yield stress         
\param  DOUBLE    eh        (i)  hardening modulus                 
\param  DOUBLE    gf        (i)  fracture energy
\param  DOUBLE    betah     (i)  controls the isotropic/kinematic hardening
\param  DOUBLE    alpha     (i)  coefficient of friction
\param  DOUBLE    stress[6] (o)  vector of stresses [11,22,33,12,23,13]
\param  DOUBLE    strain[6] (i)  actual strains from displacements  [11,22,33,12,23,13]
\param  DOUBLE  **d         (o)  constitutive matrix          
\param  INT      *iupd     (i/o) controls update of new stresses to wa         
\param  INT      *yip      (i/o) flag if global predictor step an if last step was plastic/elastic               
\param  DOUBLE   *epstn    (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE    sig[6]    (i)  stresses from WA  [11,22,33,12,23,13]
\param  DOUBLE    eps[6]    (i)  strains from WA  [11,22,33,12,23,13]
\param  DOUBLE    qn[6]     (i)  backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE    dia     (i)  internal length parameter from WA

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_mat_plast_dp()     [s9_mat_plast_dp.c]

*----------------------------------------------------------------------*/
void mat_pl_dp_lin_main(
             DOUBLE   ym,       
             DOUBLE   pv,       
             DOUBLE   sigy,     
             DOUBLE   eh,       
             DOUBLE   gf,       
             DOUBLE   betah,    
             DOUBLE   phi,      
             DOUBLE   stress[6],
             DOUBLE   strain[6],
             DOUBLE **d,        
             INT     *iupd,     
             INT     *yip,      
             DOUBLE  *epstn,    
             DOUBLE   sig[6],   
             DOUBLE   eps[6],   
             DOUBLE   qn[6],    
             DOUBLE   dia)      
{
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE q13,ro23;
DOUBLE yld;
DOUBLE dlam,ft1,ft2;
INT    isoft;
DOUBLE hards;
DOUBLE eta[6];          /*stress - backstress*/
DOUBLE sigma[6];        /*stress*/
DOUBLE delsig[6],deleps[6];
DOUBLE alpha1;
DOUBLE G,K;             /*shear/bulk modulus*/
DOUBLE trace;           /*trace of trial stresses*/
DOUBLE n[6];            /*Gradient in deviatoric direction of stresses (s/|s|)*/
DOUBLE norm;            /*Norm of n*/
DOUBLE k_eps;           /*uniaxial equivalent stresses -> isotropic hardening*/
DOUBLE epstmax;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_main");
#endif
/*----------------------------------------------------------------------*/
  *iupd=0;
/*---------------------------------------------- hardening/softening ---*/
 if(fabs(gf)>EPS10) /*softening*/
 {
   isoft = 1;
   eh = gf;
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
 }
 else /*hardening*/
 {
   isoft = 0;
   hards = (ym * eh) / (ym - eh);
 }

/*--------- calculate some needed prevalues -----------------------------*/
q13  = 1./3.;
ro23 = sqrt(2./3.);

phi    = phi * RAD;          
alpha1 = ro23 * ( q13 * (6.*sin(phi)) / (3.-sin(phi)) ); 

yld = sigy * (3.*cos(phi))/(3.-sin(phi));

K = ym / (3. - 6.* pv);
G = ym / (2. + 2.* pv);

ft2 = 1.0;   /*initialize for inverted cone*/

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
        }

        if(*yip==1)           /*elastic*/
        {
          *yip  = - *yip;
        }
        else                  /*plastic*/
        {
          dlam=0.;
          /*---- get norm (|s|) and gradient n for stress -----*/
          mat_pl_dp_lin_preval(stress,&norm,n);

          /*---- get Cep for global predictor step with dlam=0 ------*/
          mat_pl_dp_lin_mapl(hards,betah,alpha1,dlam,ym,pv,d,norm,n);

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
  
  for (i=0; i<6; i++) sigma[i] = sig[i] + delsig[i];
  for (i=0; i<6; i++) eta[i]   = sigma[i] - qn[i];


  /*---- get norm (|s|) and gradient n for eta (trial stresses) ---------*/
  mat_pl_dp_lin_preval(eta,&norm,n);

  /*-------- trace(sig - rho) --------------------------------------------*/
  trace = eta[0] + eta[1] + eta[2];
  /*-------- k(epstn)  ---------------------------------------------------*/
  /*    Beruecksichtigung des Verfestigungsverhaltens */
    if(isoft==0)
    {
      k_eps = ro23 * (yld + betah*hards* *epstn);
    }
    else
    {
      if(*epstn < epstmax)  k_eps = ro23 * (yld + betah*hards* *epstn);
      else k_eps = 0.01*sigy;    
    }

/*--- yield condition - Drucker Prager 3D - linear hardening/softening  ---*/
  ft1 = norm + alpha1 *trace - k_eps;                            /*main yield surface*/

  if (FABS(phi) > EPS5)   /*-- else == von Mieses yield criterion --*/
  ft2 = norm - 1./(3.*alpha1)*trace + k_eps/(3.*alpha1*alpha1);  /*inverted cone for DP*/


/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft1<EPS10) 
  {
    *yip = 1;
    for (i=0; i<6; i++) stress[i] = eta[i] + qn[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    *yip = 2;

    if (ft2 > EPS10) /*normal way of projection */
    {
        /* return -> new stresses, dlam, epstn */
        mat_pl_dp_lin_radi(hards,betah,epstn,G,K,alpha1,sigma,qn,&dlam,ft1,n);

        /* algorithmic elasto plastic tangent material matrix */
        for(i=0; i<6; i++) eta[i] = sigma[i] - qn[i];
        /*---- get norm (|s|) and gradient n for eta (updated stresses) ---------*/
        mat_pl_dp_lin_preval(eta,&norm,n);
        mat_pl_dp_lin_mapl(hards,betah,alpha1,dlam,ym,pv,d,norm,n);
    }
    else /*apex region, needs inverted cone*/
    {
        /* return -> new stresses, dlam, epstn */
        mat_pl_dp_lin_radi_apex(hards,betah,epstn,G,K,alpha1,sigma,qn,&dlam,ft1,ft2,n);

        /* algorithmic elasto plastic tangent material matrix for apex region -> inverted cone */
        for(i=0; i<6; i++) eta[i] = sigma[i] - qn[i];
        /*---- get norm (|s|) and gradient n for eta (updated stresses) ---------*/
        mat_pl_dp_lin_preval(eta,&norm,n);
        mat_pl_dp_lin_mapl_apex(hards,betah,alpha1,dlam,ym,pv,d,norm,n);
    }

    for (i=0; i<6; i++) stress[i] = sigma[i];
  }
/*----------------------------------------------------------------------*/
end:
/*Store values into Working-Array --> outside of this routine*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_main */


/*!----------------------------------------------------------------------
\brief return algorithm for 'Drucker Prager' material model

<pre>                                                            sh 09/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for the 'Drucker Prager' yield criterion with 
combined linear iso/kin. hardenig law  
</pre>
\param  DOUBLE    hards     (i)  (E*Eh)/(E-Eh)  
\param  DOUBLE    betah     (i)  controls the isotropic/kinematic hardening
\param  DOUBLE   *epstn    (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE    G         (i)  shear modulus  
\param  DOUBLE    K         (i)  bulk modulus  
\param  DOUBLE    alpha     (i)  coefficient of friction          
\param  DOUBLE    sigma[6] (i/o) trial stresses to be projected  
\param  DOUBLE    qn[6]    (i/o) backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE   *dlam      (o)  plastic multiplier 
\param  DOUBLE    ft_tr     (i)  yield criterion of trial state 
\param  DOUBLE    n[6]      (i)  gradient n of trial state 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_dp_lin_main()  [mat_pl_dp_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_dp_lin_radi(DOUBLE  hards, 
                        DOUBLE  betah,
                        DOUBLE *epstn,
                        DOUBLE  G,
                        DOUBLE  K,
                        DOUBLE  alpha,
                        DOUBLE  sigma[6],
                        DOUBLE  qn[6],
                        DOUBLE *dlam,
                        DOUBLE  ft_tr,
                        DOUBLE  n[6])
{
/*----------------------------------------------------------------------*/
INT    i;
DOUBLE ro23,q13;
DOUBLE H_k;
DOUBLE sig_m,s[6],qn_m,qn_dev[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_radi");
#endif
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;
H_k  = ro23*(1.-betah)*hards;

/*---- direct evaluation of new dlam -----------------------------------*/
*dlam = ft_tr / (2.*G + 9.*K*alpha*alpha + ro23*(1.+3.*alpha*alpha*(1.-betah))*hards);

/*-------- update of uniaxial plastic strain, stresses & backstress ----*/
*epstn = *epstn +  *dlam;

/*-- update of stresses and backstress   -------------------------------*/

/*--- hydrostatic part of sigma ------------------*/
sig_m = q13 * (sigma[0] + sigma[1] + sigma[2]);
sig_m = sig_m - 3.* K * *dlam * alpha;

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;

for (i=0; i<6; i++) s[i] = s[i] - 2.* G * *dlam * n[i];

sigma[0] = sig_m + s[0];
sigma[1] = sig_m + s[1];
sigma[2] = sig_m + s[2];
sigma[3] =         s[3];
sigma[4] =         s[4];
sigma[5] =         s[5];

/*--- hydrostatic part of backstress vector ------------------*/
qn_m = q13 * (qn[0] + qn[1] + qn[2]);
qn_m = qn_m + *dlam * H_k * alpha;

/*--- deviatoric part of backstress vector -------------------*/
qn_dev[0]= q13 * (2. * qn[0] - qn[1] - qn[2] );
qn_dev[1]= q13 * (2. * qn[1] - qn[0] - qn[2] );
qn_dev[2]= q13 * (2. * qn[2] - qn[0] - qn[1] );
qn_dev[3]= qn[3] ;
qn_dev[4]= qn[4] ;
qn_dev[5]= qn[5] ;

for (i=0; i<6; i++) qn_dev[i] = qn_dev[i] + *dlam * H_k * n[i];

qn[0] = qn_m + qn_dev[0];
qn[1] = qn_m + qn_dev[1];
qn[2] = qn_m + qn_dev[2];
qn[3] =        qn_dev[3];
qn[4] =        qn_dev[4];
qn[5] =        qn_dev[5];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_radi*/


/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 von 'Drucker Prager' plasticity                                   

<pre>                                                            sh 09/03
This routine calculates the elasto-plastic consistent material tangent
for the 'Drucker Prager plasticity with combined, linear iso/kin. hardening
law. 
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    hards     (i)  (E*Eh)/(E-Eh)  
\param  DOUBLE    betah     (i)  controls the isotropic/kinematic hardening
\param  DOUBLE    alpha     (i)  coefficient of friction          
\param  DOUBLE    dlam      (i)  plastic multiplier 
\param  DOUBLE    e         (i)  youngs modulus  
\param  DOUBLE    vnu       (i)  poisson's ration  
\param  DOUBLE  **d         (o)  material matrix to be calculated   
                                 i: elastic material matric (not needed)
                                 o: konsistent elasto-plastic material matrix          
\param  DOUBLE    norm      (i)  norm of deviatoric stresses |s| (projected) 
\param  DOUBLE    n[6]      (i)  gradient n (s/|s|)  (projected)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_dp_lin_main()     [mat_pl_dp_lin.c]

*-----------------------------------------------------------------------*/
void mat_pl_dp_lin_mapl(DOUBLE   hards,
                        DOUBLE   betah,
                        DOUBLE   alpha,
                        DOUBLE   dlam,
                        DOUBLE   e,
                        DOUBLE   vnu,
                        DOUBLE **d,
                        DOUBLE   norm,
                        DOUBLE   n[6])
{
/*----------------------------------------------------------------------*/
INT    i,j, cc, irc;
INT    isix=6, ione=1;
DOUBLE dzero = 0.;
DOUBLE ro23,q13,q23;
DOUBLE fac;
DOUBLE df2dss[36];
DOUBLE dfds[6];
DOUBLE dfdk;
DOUBLE hm[6][6],CC[36], CI_vec[36];
DOUBLE nenner, zaehler1[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_mapl");
#endif
/*----------------------------------------------------------------------*/

/*initialize needed matrices*/
nenner = 0.0;
for (i=0; i<6; i++)  zaehler1[i]     = 0.0;
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q13  = 1./3.;
q23  = 2./3.;
/*----------------------------------------------------------------------*/

/*- get INVERSE of linear elastic isotropic material tensor --- Cel_-1 -*/
mat_el_iso_inv(e, vnu, CI_vec); /* This routine returns a vector [36] --*/

/*------  df2dss = 1/|n| * (P_dev - n dyad n) --------------------------*/
df2dss[ 0] = ( q23 - n[0] * n[0]); 
df2dss[ 1] = (-q13 - n[0] * n[1]); 
df2dss[ 2] = (-q13 - n[0] * n[2]); 
df2dss[ 3] = (     - n[0] * n[3]); 
df2dss[ 4] = (     - n[0] * n[4]); 
df2dss[ 5] = (     - n[0] * n[5]); 

df2dss[ 6] = (-q13 - n[1] * n[0]); 
df2dss[ 7] = ( q23 - n[1] * n[1]); 
df2dss[ 8] = (-q13 - n[1] * n[2]); 
df2dss[ 9] = (     - n[1] * n[3]); 
df2dss[10] = (     - n[1] * n[4]); 
df2dss[11] = (     - n[1] * n[5]); 

df2dss[12] = (-q13 - n[2] * n[0]); 
df2dss[13] = (-q13 - n[2] * n[1]); 
df2dss[14] = ( q23 - n[2] * n[2]); 
df2dss[15] = (     - n[2] * n[3]); 
df2dss[16] = (     - n[2] * n[4]); 
df2dss[17] = (     - n[2] * n[5]); 

df2dss[18] = (     - n[3] * n[0]); 
df2dss[19] = (     - n[3] * n[1]); 
df2dss[20] = (     - n[3] * n[2]); 
df2dss[21] = ( 2.  - n[3] * n[3]); 
df2dss[22] = (     - n[3] * n[4]); 
df2dss[23] = (     - n[3] * n[5]); 

df2dss[24] = (     - n[4] * n[0]); 
df2dss[25] = (     - n[4] * n[1]); 
df2dss[26] = (     - n[4] * n[2]); 
df2dss[27] = (     - n[4] * n[3]); 
df2dss[28] = ( 2.  - n[4] * n[4]); 
df2dss[29] = (     - n[4] * n[5]); 

df2dss[30] = (     - n[5] * n[0]); 
df2dss[31] = (     - n[5] * n[1]); 
df2dss[32] = (     - n[5] * n[2]); 
df2dss[33] = (     - n[5] * n[3]); 
df2dss[34] = (     - n[5] * n[4]); 
df2dss[35] = ( 2.  - n[5] * n[5]); 

/*------  Cel_-1 + dlam*df2dss -------------------------------------------*/
fac = dlam / norm;
for (i=0; i<36; i++) CI_vec[i] += fac * df2dss[i];

/*------ Theta(CC/hm) = ( Cel_-1 + dlam*df2dss )-1 -----------------------*/
c1inv6(CI_vec,CC,&irc); /*fortran routine, uses a vector*/

/*----------- write vector back to matrix -------------------*/
cc=0;
for (i=0; i<6; i++) for (j=0; j<6; j++) hm[i][j] = CC[cc++];

/*------ dfdk = -sqrt(2/3)*(1+3*alpha^2*(1-betah))*H -------------------------*/
dfdk = -ro23*(1. + 3.*alpha*alpha*(1.-betah))*hards;

/*------ dfds[6] = n[6] + alpha*1 -------------------------------------------*/
dfds[0] = n[0] + alpha;
dfds[1] = n[1] + alpha;
dfds[2] = n[2] + alpha;
dfds[3] = n[3] * 2.;
dfds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
dfds[5] = n[5] * 2.;

/*----- nenner = dfds:hm:dfds - dfdk --------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) nenner += dfds[j]*hm[j][i]*dfds[i];
nenner = nenner -dfdk;

/*------- zaehler1 = hm * dfds --------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) zaehler1[i] += hm[i][j]*dfds[j];

/*---- Cep = hm - [(hm:dfds)dyad(hm:dfds)]/(dfds:hm:dfds - dfdk) ---------*/
fac = 1./nenner;

for (i=0; i<6; i++) for (j=0; j<6; j++) 
{
      d[i][j] = hm[i][j] - fac * zaehler1[i] * zaehler1[j]; 
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_mapl */


/*!----------------------------------------------------------------------
\brief return algorithm for 'Drucker Prager' material model for the
inverted cone

<pre>                                                            sh 09/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorithm for the 'Drucker Prager' yield criterion with 
combined linear iso/kin. hardenig law  
</pre>
\param  DOUBLE    hards     (i)  (E*Eh)/(E-Eh)  
\param  DOUBLE    betah     (i)  controls the isotropic/kinematic hardening
\param  DOUBLE   *epstn    (i/o) uniaxial equivalent strain -> WA
\param  DOUBLE    G         (i)  shear modulus  
\param  DOUBLE    K         (i)  bulk modulus 
\param  DOUBLE    alpha1    (i)  coefficient of friction of main yield surface         
\param  DOUBLE    sigma[6] (i/o) trial stresses to be projected  
\param  DOUBLE    qn[6]    (i/o) backstress vector from WA  [11,22,33,12,23,13]
\param  DOUBLE   *dlam      (o)  plastic multiplier (dlam1 + dlam2 ; of main yield surface + inverted cone)
\param  DOUBLE    ft1       (i)  yilcr at trial state of main yield surface
\param  DOUBLE    ft2       (i)  yilcr at trial state of inverted cone

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_dp_lin_main()  [mat_pl_dp_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_dp_lin_radi_apex(DOUBLE  hards, 
                             DOUBLE  betah,
                             DOUBLE *epstn,
                             DOUBLE  G,
                             DOUBLE  K,
                             DOUBLE  alpha1,
                             DOUBLE  sigma[6],
                             DOUBLE  qn[6],
                             DOUBLE *dlam,
                             DOUBLE  ft1,
                             DOUBLE  ft2,
                             DOUBLE  n[6])
{
/*----------------------------------------------------------------------*/
INT    i;
DOUBLE ro23,q13;
DOUBLE alpha2;
DOUBLE dlam1,dlam2;
DOUBLE d11,d12,d21,d22,det;
DOUBLE H_k;
DOUBLE sig_m,s[6],qn_m,qn_dev[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_radi_apex");
#endif
/*----------------------------------------------------------------------*/
ro23   = sqrt(2./3.);
q13    = 1./3.;
alpha2 = - 1./(3.*alpha1);
H_k    = ro23*(1.-betah)*hards;

/*------- direct evaluation of dlam1 and dlam2 ---- (lin. hardening) ---*/
d11 = 2.*G + 9.*K*alpha1*alpha1 + ro23*(1.+3.*alpha1*alpha1)*(1.-betah)*hards + ro23*betah*hards;
d12 = 2.*G - 3.*K;
d21 = 2.*G - 3.*K - 1./(3.*alpha1*alpha1)*ro23*betah*hards;
d22 = 2.*G + K/(alpha1*alpha1);

det = d11*d22 - d12*d21;

/*---- direct evaluation of new dlam1 & dlam2 --------------------------*/
dlam1 = 1./det * (+d22*ft1 - d12*ft2);
dlam2 = 1./det * (-d21*ft1 + d11*ft2);

/*---- check for inactive yield surfaces -------------------------------*/
/*---- NOTE: projection in the apex region implies that both yield -----*/
/*----       surfaces are active during return process             -----*/
if   (dlam1 < 0.0) dserror("dlam1 < 0 in mat_pl_dp_lin -> radi_apex");
if   (dlam2 > 0.0) dserror("dlam2 > 0 in mat_pl_dp_lin -> radi_apex");

/*---- new dlam for mapl_apex  -----------------------------------------*/
*dlam = dlam1 + dlam2; 

/*-------- update of uniaxial plastic strain, stresses & backstress ----*/
*epstn = *epstn +  dlam1;


/*-- update of stresses and backstress   -------------------------------*/

/*--- hydrostatic part of sigma ------------------*/
sig_m = q13 * (sigma[0] + sigma[1] + sigma[2]);
sig_m = sig_m - 3.* K * dlam1 * alpha1
              - 3.* K * dlam2 * alpha2;

/*--- deviatoric part of sigma -------------------*/
s[0]= q13 * (2. * sigma[0] - sigma[1] - sigma[2] );
s[1]= q13 * (2. * sigma[1] - sigma[0] - sigma[2] );
s[2]= q13 * (2. * sigma[2] - sigma[0] - sigma[1] );
s[3]= sigma[3] ;
s[4]= sigma[4] ;
s[5]= sigma[5] ;

for (i=0; i<6; i++) s[i] = s[i] - 2.* G * dlam1 * n[i]
                                - 2.* G * dlam2 * n[i];

sigma[0] = sig_m + s[0];
sigma[1] = sig_m + s[1];
sigma[2] = sig_m + s[2];
sigma[3] =         s[3];
sigma[4] =         s[4];
sigma[5] =         s[5];


/*--- NOTE: H_k is only defined for main yield surface -------*/
/*--- hydrostatic part of backstress vector ------------------*/
qn_m = q13 * (qn[0] + qn[1] + qn[2]);
qn_m = qn_m + dlam1 * H_k * alpha1;

/*--- deviatoric part of backstress vector -------------------*/
qn_dev[0]= q13 * (2. * qn[0] - qn[1] - qn[2] );
qn_dev[1]= q13 * (2. * qn[1] - qn[0] - qn[2] );
qn_dev[2]= q13 * (2. * qn[2] - qn[0] - qn[1] );
qn_dev[3]= qn[3] ;
qn_dev[4]= qn[4] ;
qn_dev[5]= qn[5] ;

for (i=0; i<6; i++) qn_dev[i] = qn_dev[i] + dlam1 * H_k * n[i];

qn[0] = qn_m + qn_dev[0];
qn[1] = qn_m + qn_dev[1];
qn[2] = qn_m + qn_dev[2];
qn[3] =        qn_dev[3];
qn[4] =        qn_dev[4];
qn[5] =        qn_dev[5];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_radi_apex */


/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 von 'Drucker Prager' plasticity for trial stresses in the apex region
 --> inverted cone                                  

<pre>                                                            sh 09/03
This routine calculates the elasto-plastic consistent material tangent
for the 'Drucker Prager plasticity with combined, linear iso/kin. hardening
law. --> inverted cone
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    hards     (i)  (E*Eh)/(E-Eh)  
\param  DOUBLE    betah     (i)  controls the isotropic/kinematic hardening
\param  DOUBLE    alpha1    (i)  coefficient of friction for main yield surface          
\param  DOUBLE    dlam      (i)  plastic multiplier (dlam1 + dlam2 ; of main yield surface + inverted cone)
\param  DOUBLE    e         (i)  youngs modulus  
\param  DOUBLE    vnu       (i)  poisson's ratio  
\param  DOUBLE  **d         (o)  material matrix to be calculated   
                                 i: elastic material matric (not needed)
                                 o: konsistent elasto-plastic material matrix          
\param  DOUBLE    norm      (i)  norm of deviatoric stresses (|s|) (projected)  
\param  DOUBLE    n[6]      (i)  gradient n (s/|s|) (projected) 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_dp_lin_main()     [mat_pl_dp_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_dp_lin_mapl_apex(DOUBLE   hards,
                             DOUBLE   betah,
                             DOUBLE   alpha1,
                             DOUBLE   dlam,
                             DOUBLE   e,
                             DOUBLE   vnu,
                             DOUBLE **d,
                             DOUBLE   norm, 
                             DOUBLE   n[6]) 
{
/*----------------------------------------------------------------------*/
INT    i,j,k,cc,irc;
INT    isix=6, ione=1;
DOUBLE dzero = 0.;
DOUBLE alpha2;
DOUBLE det;
DOUBLE ro23,q12,q13,q23;
DOUBLE df1ds[6], df2ds[6];
DOUBLE E[2][2],E_help[2][2],U[6][2];
DOUBLE hm[6][6];
DOUBLE H1[6][2], H2[6][2];
DOUBLE HM2[6][6];
DOUBLE CC[36], CI_vec[36];
DOUBLE fac;
DOUBLE df2dss[36];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_mapl_apex");
#endif
/*----------------------------------------------------------------------*/
ro23 = sqrt(2./3.);
q12  = 1./2.;
q13  = 1./3.;
q23  = 2./3.;
alpha2 = - 1./(3.*alpha1);
/*----------------------------------------------------------------------*/

/*initialize needed matrices*/
for (i=0; i<2; i++) for (j=0; j<2; j++) E_help[i][j] = 0.0;
for (i=0; i<6; i++) for (j=0; j<2; j++) H1[i][j]     = 0.0;
for (i=0; i<6; i++) for (j=0; j<2; j++) H2[i][j]     = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++) HM2[i][j]    = 0.0;
/*----------------------------------------------------------------------*/

/*- get INVERSE of linear elastic isotropic material tensor --- Cel_-1 -*/
mat_el_iso_inv(e, vnu, CI_vec); /* This routine returns a vector [36] --*/

/*------  df2dss = 1/|n| * (P_dev - n dyad n) --------------------------*/
df2dss[ 0] = ( q23 - n[0] * n[0]); 
df2dss[ 1] = (-q13 - n[0] * n[1]); 
df2dss[ 2] = (-q13 - n[0] * n[2]); 
df2dss[ 3] = (     - n[0] * n[3]); 
df2dss[ 4] = (     - n[0] * n[4]); 
df2dss[ 5] = (     - n[0] * n[5]); 

df2dss[ 6] = (-q13 - n[1] * n[0]); 
df2dss[ 7] = ( q23 - n[1] * n[1]); 
df2dss[ 8] = (-q13 - n[1] * n[2]); 
df2dss[ 9] = (     - n[1] * n[3]); 
df2dss[10] = (     - n[1] * n[4]); 
df2dss[11] = (     - n[1] * n[5]); 

df2dss[12] = (-q13 - n[2] * n[0]); 
df2dss[13] = (-q13 - n[2] * n[1]); 
df2dss[14] = ( q23 - n[2] * n[2]); 
df2dss[15] = (     - n[2] * n[3]); 
df2dss[16] = (     - n[2] * n[4]); 
df2dss[17] = (     - n[2] * n[5]); 

df2dss[18] = (     - n[3] * n[0]); 
df2dss[19] = (     - n[3] * n[1]); 
df2dss[20] = (     - n[3] * n[2]); 
df2dss[21] = ( 2.  - n[3] * n[3]); 
df2dss[22] = (     - n[3] * n[4]); 
df2dss[23] = (     - n[3] * n[5]); 

df2dss[24] = (     - n[4] * n[0]); 
df2dss[25] = (     - n[4] * n[1]); 
df2dss[26] = (     - n[4] * n[2]); 
df2dss[27] = (     - n[4] * n[3]); 
df2dss[28] = ( 2.  - n[4] * n[4]); 
df2dss[29] = (     - n[4] * n[5]); 

df2dss[30] = (     - n[5] * n[0]); 
df2dss[31] = (     - n[5] * n[1]); 
df2dss[32] = (     - n[5] * n[2]); 
df2dss[33] = (     - n[5] * n[3]); 
df2dss[34] = (     - n[5] * n[4]); 
df2dss[35] = ( 2.  - n[5] * n[5]); 

/*------  Cel_-1 + dlam*df2dss -------------------------------------------*/
fac = dlam / norm;
for (i=0; i<36; i++) CI_vec[i] += fac * df2dss[i];

/*------ Theta(CC/hm) = ( Cel_-1 + dlam*df2dss )-1 -----------------------*/
c1inv6(CI_vec,CC,&irc); /*fortran routine, uses a vector*/
/*----------- write vector back to matrix -------------------*/
cc=0;
for (i=0; i<6; i++) for (j=0; j<6; j++) hm[i][j] = CC[cc++];


/*----- Matrix of plastic Modules --> see Menrath (A.8)  ---------------------*/
E[0][0] = -ro23*(1.+3.*alpha1*alpha1*(1.-betah))*hards;
E[0][1] = 0.0;
E[1][0] = (ro23*betah*hards)/(3.*alpha1*alpha1);
E[1][1] = 0.0;

/*------ df1ds[6] = n[6] + alpha1*1 -------------------------------------------*/
df1ds[0] = n[0] + alpha1;
df1ds[1] = n[1] + alpha1;
df1ds[2] = n[2] + alpha1;
df1ds[3] = n[3] * 2.;
df1ds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
df1ds[5] = n[5] * 2.;

/*------ df2ds[6] = n[6] + alpha2*1 -------------------------------------------*/
df2ds[0] = n[0] + alpha2;
df2ds[1] = n[1] + alpha2;
df2ds[2] = n[2] + alpha2;
df2ds[3] = n[3] * 2.;
df2ds[4] = n[4] * 2.;  /*multiply with 2. for vector/matrix notation !!*/
df2ds[5] = n[5] * 2.;

/*----- E_help = UT*hm*U --------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[0][0] += df1ds[j] * hm[j][i] * df1ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[0][1] += df1ds[j] * hm[j][i] * df2ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[1][0] += df2ds[j] * hm[j][i] * df1ds[i];
for (i=0; i<6; i++) for (j=0; j<6; j++) E_help[1][1] += df2ds[j] * hm[j][i] * df2ds[i];
              
/*------ E_help = E + UT*hm*U ------------------------------*/
E_help[0][0] += E[0][0];
E_help[0][1] += E[0][1];
E_help[1][0] += E[1][0];
E_help[1][1] += E[1][1];

/*------------ E = ( E + UT*hm*U )-1 -----------------------*/
det = E_help[0][0]*E_help[1][1] - E_help[0][1]*E_help[1][0];
det = 1./det;

E[0][0] =  det * E_help[1][1];
E[0][1] = -det * E_help[0][1];
E[1][0] = -det * E_help[1][0];
E[1][1] =  det * E_help[0][0];

/*-----------------------------------------------------------*/
U[0][0] = df1ds[0];
U[1][0] = df1ds[1];
U[2][0] = df1ds[2];
U[3][0] = df1ds[3];
U[4][0] = df1ds[4];
U[5][0] = df1ds[5];

U[0][1] = df2ds[0];
U[1][1] = df2ds[1];
U[2][1] = df2ds[2];
U[3][1] = df2ds[3];
U[4][1] = df2ds[4];
U[5][1] = df2ds[5];

/*--------- H1[6][2] = hm[6][6]*U[6][2] ---------------------*/
for (j=0; j<6; j++) for (k=0; k<6; k++) for (i=0; i<2; i++) H1[j][i] += hm[j][k] * U[k][i]; 

/*--------- H2[6][2] = hm[6][6]*U[6][2]*E[2][2] = H1[6][2]*E[2][2] -----*/
for (j=0; j<6; j++) for (k=0; k<2; k++) for (i=0; i<2; i++) H2[j][i] += H1[j][k] * E[k][i]; 

/*--------- HM2[6][6] = hm[6][6]*U[6][2]*E[2][2]*UT[2][6]*hm[6][6] = H2[6][2] * H1T[2][6] ----*/
for (j=0; j<6; j++) for (k=0; k<2; k++) for (i=0; i<6; i++) HM2[j][i] += H2[j][k] * H1[i][k]; 

/*--- Cep = hm - HM2 ---------------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) d[i][j] = hm[i][j] - HM2[i][j]; 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_mapl_apex */


/*!----------------------------------------------------------------------
\brief calculates some neede values for Drucker-Prager Model with
combined linear iso/kin hardening                                  

<pre>                                                            sh 09/03
This routine calculates the 
- norm of deviatoric stresses
- the gradient n (dfds)
for the 'Drucker Prager plasticity with combined, linear iso/kin. hardening
law. 
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    sig[6]  (i/o)  stresses (input, could be modified due 
                                 to numerical aspects if |s| is close to ZERO  
\param  DOUBLE   *norm     (o)   norm of deviatoric stresses |s|
\param  DOUBLE    n[6]     (o)   gradient n of deviatoric stresses (s/|s|)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_dp_lin_main()     [mat_pl_dp_lin.c]

*----------------------------------------------------------------------*/
void mat_pl_dp_lin_preval(DOUBLE   sig[6], 
                          DOUBLE  *norm,
                          DOUBLE   n[6])
{
/*----------------------------------------------------------------------*/
DOUBLE q13;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_dp_lin_preval");
#endif
/*----------------------------------------------------------------------*/
q13  = 1./3.;

/*-------- check for numerical problems  --> |s| = ZERO ----------------*/
if (FABS(sig[0]-sig[1]) < EPS8  &&    
    FABS(sig[0]-sig[2]) < EPS8  &&     /*sig[0]=sig[1]=sig[2]*/
    FABS(sig[1]-sig[2]) < EPS8)
{
   /*-- modify the stress vector slightly to avoid numerical problems --*/
   if(FABS(sig[2]) < EPS8)  sig[2] = EPS5;
   else                     sig[2] = 1.001*sig[2];
}    

/*-------- |s| ---------------------------------------------------------*/
*norm = sqrt(2./3.*(sig[0]*sig[0] + sig[1]*sig[1] + sig[2]*sig[2] -
                    sig[0]*sig[1] - sig[0]*sig[2] - sig[1]*sig[2])
              + 2.*(sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]));

/*-------- components of gradient n ------------------------------------*/
n[0]= q13 * (2. * sig[0] - sig[1] - sig[2] ) / *norm;
n[1]= q13 * (2. * sig[1] - sig[0] - sig[2] ) / *norm;
n[2]= q13 * (2. * sig[2] - sig[0] - sig[1] ) / *norm;
n[3]= sig[3] / *norm;
n[4]= sig[4] / *norm;
n[5]= sig[5] / *norm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_dp_lin_preval */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/

