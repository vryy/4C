/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - mat_pl_hoff_main: which calculates the constitutive matrix for a
                     anisotropic-plasticity model based on the
                     Hoffman criterion (see Dis. Hoermann p.57ff)
                     this routine is formulated in cartesian coordinate
                     system, general 3D with the sorting [11,22,33,12,23,13]
- mat_pl_hoff_homa: which calculates the coupling matrix P and the
                    coupling vector Q
- mat_pl_hoff_mapl: which calculates the elasto-plastic consistent
                    material matrix
- mat_pl_hoff_yilcr: which calculates the hoffman yield criterion
- mat_pl_hoff_radi: which performs the radial return algorithm
- mat_pl_hoff_hard: which calculates the necessary hardening variables
- mat_pl_hoff_serv: which calculates some needed prevalues

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
\brief consitutive matrix for von Anisotropic-Plasticity,based on the 
Hoffman-criterion (see Dis. Hoermann p.57f)                                    

<pre>                                                            sh 03/03
This routine calculates the constitutive matrix and forces for an 
Anisotropic-Plasticity model with a linear hardening law, based on the
Hoffman-criterion. 
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  PL_HOFF   *mat       (i)  material properties for hoffman material 
\param  DOUBLE     stress[6] (i)  vector of stresses [11,22,33,12,23,13]               
\param  DOUBLE   **C         (o)  constitutive matrix          
\param  INT       *iupd     (i/o) controls update of new stresses to wa         
\param  INT       *yip      (i/o) flag if global predictor step an if last step was plastic/elastic                 
\param  DOUBLE    *dhard    (i/o) hardening value from WA
\param  DOUBLE     strain[6] (i)  actual strains from displacements  [11,22,33,12,23,13]
\param  DOUBLE     sig[6]    (i)  stresses from WA  [11,22,33,12,23,13]
\param  DOUBLE     eps[6]    (i)  strains from WA  [11,22,33,12,23,13]
\param  DOUBLE     dkappa[6] (i)  
\param  DOUBLE     gamma[6]  (i)  
\param  DOUBLE     rkappa[9] (i)  

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_mat_plast_hoff()     [s9_mat_plast_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_main(
                 PL_HOFF   *mat,      /*!< material properties          */
                 DOUBLE     stress[6],/*!< vector of stresses [11,22,33,12,23,13]  */
                 DOUBLE   **C,        /*!< constitutive matrix          */
                 INT       *iupd,     /*!< controls update of new stresses to wa */
                 INT       *yip,      /*!< from WA*/
                 DOUBLE    *dhard,    /*!< from WA*/
                 DOUBLE     strain[6],/*!< actual strains from displacements*/
                 DOUBLE     sig[6],   /*!< stresses from WA*/
                 DOUBLE     eps[6],   /*!< strains from WA*/
                 DOUBLE     dkappa[6],    
                 DOUBLE     gamma[6],
                 DOUBLE     rkappa[9])    
{
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE deleps[6], delsig[6];
DOUBLE tau[6];                /*predictor stress*/
DOUBLE signew[6];
DOUBLE P[6][6],P_hat[6][6];
DOUBLE Q[6],   Q_hat[6];
DOUBLE dlam;
DOUBLE ft;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_main");
#endif
/*----------------------------------------------------------------------*/
  *iupd=0;
/*-------- calculate the P-Matrix and the Q-Vector ---------------------*/
mat_pl_hoff_homa(mat, P, Q, rkappa);

/*----- copy P and Q to P_hat and Q_hat for projection -----------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) P_hat[i][j] = P[i][j]; 
for (i=0; i<6; i++)                     Q_hat[i]    = Q[i]; 

/*----- get elastic orthotropic material tensor ------------------------*/
mat_el_orth(mat->emod1 ,
            mat->emod2 ,
            mat->emod3 ,
            mat->xnue12,
            mat->xnue23,
            mat->xnue13,
            mat->gmod12,
            mat->gmod23,
            mat->gmod13,
            C);

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
          mat_pl_hoff_mapl( stress, C, P_hat , Q_hat , dlam, dhard, gamma); 
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

  /*deleps[3] = 2. * deleps[3];*/   /*write as vector in xx_call_mat*/
  /*deleps[4] = 2. * deleps[4];*/
  /*deleps[5] = 2. * deleps[5];*/

  for (i=0; i<6; i++) delsig[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) delsig[i] += C[i][j] * deleps[j];
  
  for (i=0; i<6; i++) tau[i] = sig[i] + delsig[i];

/*---------yield condition - Hoffman criterion  ------------------------*/
  ft = mat_pl_hoff_yilcr(tau, P, Q, mat->uniax);
  
/*------------- state of stress within yield surface - E L A S T I C ---*/
  if (ft<EPS10) 
  {
    *yip = 1;
    for (i=0; i<6; i++) stress[i] = tau[i];
  }
/*------------ state of stress outside yield surface - P L A S T I C ---*/
  else 
  {
    *yip = 2;
    /* return -> new stresses, dlam, epstn */
    mat_pl_hoff_radi(mat,tau,C,P,Q,rkappa,dhard,&dlam,signew,P_hat,Q_hat,gamma,dkappa);
    /* algorithmic elasto plastic tangent material matrix */
    mat_pl_hoff_mapl( signew, C, P_hat , Q_hat , dlam, dhard, gamma); 

    for (i=0; i<6; i++) stress[i] = signew[i];
  }
/*----------------------------------------------------------------------*/
end:
/*Store values into Working-Array --> outside of this routine*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_hoff_main */



/*!----------------------------------------------------------------------
\brief calculates needed Matrices P and Q for anisotropic plasticity
model                                   

<pre>                                                            sh 03/03
This routine calculates the neede Matrices P and Q for anisotropic 
plasticity model, based on the hoffman yield criterion
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  PL_HOFF   *mat       (i)  material properties for hoffman material  
\param  DOUBLE     P[6][6]   (o)  coupling matrix P             
\param  DOUBLE     Q[6]      (o)  coupling vector Q         
\param  DOUBLE     rkappa[9] (i)  vector of equivalent strains 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_homa(PL_HOFF   *mat,       /*!< material properties          */
                      DOUBLE     P[6][6],   /*!< */
                      DOUBLE     Q[6],      /*!< */
                      DOUBLE     rkappa[9]) /*!< */
{
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE uniax;
DOUBLE sig[9];
DOUBLE sig_I[9];
DOUBLE sig_H[9];
DOUBLE H[9];
DOUBLE alfa[9];
DOUBLE c1,c2,c3,c1c2c3;

#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_homa");
#endif
/*--------- initialize P,Q ---------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++)  P[i][j] = 0.0;
for (i=0; i<6; i++)                      Q[i]    = 0.0;
/*--------- get material properties ------------------------------------*/
uniax = mat->uniax;

sig_I[0] = mat->s11T; 
sig_I[1] = mat->s11C; 
sig_I[2] = mat->s22T; 
sig_I[3] = mat->s22C; 
sig_I[4] = mat->s33T; 
sig_I[5] = mat->s33C; 
sig_I[6] = mat->s12;
sig_I[7] = mat->s23;
sig_I[8] = mat->s13;

sig_H[0] = mat->sh11T;
sig_H[1] = mat->sh11C;
sig_H[2] = mat->sh22T;
sig_H[3] = mat->sh22C;
sig_H[4] = mat->sh33T;
sig_H[5] = mat->sh33C;
sig_H[6] = mat->sh12; 
sig_H[7] = mat->sh23; 
sig_H[8] = mat->sh13; 

H[0] = mat->ha11T;
H[1] = mat->ha11C;
H[2] = mat->ha22T;
H[3] = mat->ha22C;
H[4] = mat->ha33T;
H[5] = mat->ha33C;
H[6] = mat->ha12; 
H[7] = mat->ha23; 
H[8] = mat->ha13; 

/*--------- updated yield stress ---------------------------------------*/
sig[0] = sig_I[0] + (sig_H[0]-sig_I[0])/(H[0]) * rkappa[0];
sig[1] = sig_I[1] + (sig_H[1]-sig_I[1])/(H[1]) * rkappa[1];
sig[2] = sig_I[2] + (sig_H[2]-sig_I[2])/(H[2]) * rkappa[2];
sig[3] = sig_I[3] + (sig_H[3]-sig_I[3])/(H[3]) * rkappa[3];
sig[4] = sig_I[4] + (sig_H[4]-sig_I[4])/(H[4]) * rkappa[4];
sig[5] = sig_I[5] + (sig_H[5]-sig_I[5])/(H[5]) * rkappa[5];
sig[6] = sig_I[6] + (sig_H[6]-sig_I[6])/(H[6]) * rkappa[6];
sig[7] = sig_I[7] + (sig_H[7]-sig_I[7])/(H[7]) * rkappa[7];
sig[8] = sig_I[8] + (sig_H[8]-sig_I[8])/(H[8]) * rkappa[8];

alfa[0] =  uniax * uniax *(sig[1]-sig[0])/(sig[0] * sig[1]);
alfa[1] =  uniax * uniax *(sig[3]-sig[2])/(sig[2] * sig[3]);
alfa[2] =  uniax * uniax *(sig[5]-sig[4])/(sig[4] * sig[5]);
alfa[3] =  uniax * uniax *(1./(sig[1] * sig[0]) + 1./(sig[3] * sig[2]) - 1./(sig[5] * sig[4]));
alfa[4] =  uniax * uniax *(1./(sig[3] * sig[2]) + 1./(sig[5] * sig[4]) - 1./(sig[1] * sig[0]));
alfa[5] =  uniax * uniax *(1./(sig[5] * sig[4]) + 1./(sig[1] * sig[0]) - 1./(sig[3] * sig[2]));
alfa[6] = (uniax * uniax)/(3. * sig[6] * sig[6]);
alfa[7] = (uniax * uniax)/(3. * sig[7] * sig[7]);
alfa[8] = (uniax * uniax)/(3. * sig[8] * sig[8]);

/*--------- check convecity of yield surface ---------------------------------------*/
/* see Pankaj et al. (1999), 'Convexity studies of two anisotropic                  */
/*                            yield criteria in principal stress                    */
/*                            space', Engineering Computations Vol 16               */
/*                            No. 2, pp 215-229                                     */
c1 = alfa[4]/(uniax * uniax);
c2 = alfa[5]/(uniax * uniax);
c3 = alfa[3]/(uniax * uniax);
c1c2c3 = c1 * c2 + c2 * c3 + c3 * c1;
if (c1c2c3 < 0.0)
{
   printf("--- Convecity of Yield Function is not given with -----\n");
   printf("--- the current Set of Material Parameters ------------\n");
   printf("-------------------------------------------------------\n");
   printf("--- see Pankaj et al. (1999) --------------------------\n");
   printf("--- 'Convecity studies of two anisotropic yield -------\n");
   printf("---  criteria in principal stress space ---------------\n");
   printf("--- Engineering Computatins Vol 16, No. 2, pp 215-229 -\n");
   printf("-------------------------------------------------------\n");
}

/*--------- fill P-Matrix and Q-Vector ---------------------------------*/
P[0][0] = alfa[3]+alfa[5];
P[1][1] = alfa[3]+alfa[4];
P[2][2] = alfa[4]+alfa[5];
P[3][3] = 6. * alfa[6];
P[4][4] = 6. * alfa[7];
P[5][5] = 6. * alfa[8];
P[0][1] = -alfa[3];
P[0][2] = -alfa[5];
P[2][1] = -alfa[4];
P[1][0] = P[0][1];
P[2][0] = P[0][2];
P[1][2] = P[2][1];

Q[0] = alfa[0];
Q[1] = alfa[1];
Q[2] = alfa[2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_hoff_homa */

/*!----------------------------------------------------------------------
\brief calculates the elasto-plastic consistent material tangent for
 anisotropic plasticity model                                   

<pre>                                                            sh 03/03
This routine calculates the elasto-plastic consistent material tangent
for the anisotropic plasticity model which is based on the hoffman
yield criterion, combined with a anisotropic, linear hardening law
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE     stress[6] (i)  current stresses  
\param  DOUBLE   **C        (i/o) material matrix to be calculated             
\param  DOUBLE     P[6][6]  (o)  coupling matrix P         
\param  DOUBLE     Q[6]     (i)  coupling vector Q 
\param  DOUBLE     dlam     (i)  plastic multiplier 
\param  DOUBLE    *dhard    (i)  hardening value from WA 
\param  DOUBLE     gamma[6] (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_mapl(DOUBLE     stress[6],   /*!< current stresses          */
                      DOUBLE   **C,           /*!< material matrix to be calculated */
                      DOUBLE     P[6][6],     /*!< coupling matrix P */
                      DOUBLE     Q[6],        /*!< coupling vector Q */
                      DOUBLE     dlam,        /*!< plastic multiplier */
                      DOUBLE    *dhard,       /*!< hardening value from WA*/
                      DOUBLE     gamma[6])    /*!< from WA*/ 
{
INT    i,j, cc, irc;
INT    isix=6, ione=1;
DOUBLE dzero = 0.;
DOUBLE xhix;
DOUBLE HIa[6][6];
DOUBLE ma66h[36],CI[36], H[36], HI[36], P_vec[36];
DOUBLE PS[6], PSQ[6], X[6], X1[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_mapl");
#endif
/*--------- H = CI + dlam * P ------------------------------------------*/
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = C[i][j];  
  c1inv6(ma66h,CI,&irc); /*fortran routine, uses a vector*/

  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) P_vec[cc++] = P[i][j];  

  for (i=0; i<36; i++)      H[i] = CI[i] + dlam * P_vec[i];

/*--------- inverse of H ----------------------------------------------*/
  c1inv6(H,HI,&irc); /*fortran routine, uses a vector*/

/*--------- dGdS ------------------------------------------------------*/
  c1ab (P_vec,stress,PS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++)      PSQ[i] = PS[i] + Q[i];

  c1ab (HI,PSQ,X,&isix,&isix,&ione,&dzero);
  c1ab (HI,gamma,X1,&isix,&isix,&ione,&dzero);

  xhix = 0.;
  for (i=0; i<6; i++) xhix += PSQ[i] * X1[i];

  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) HIa[i][j] = HI[cc++];
  
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
        C[i][j] = HIa[i][j] - X1[i] * X[j] / (xhix - (*dhard));
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mat_pl_hoff_mapl */


/*!----------------------------------------------------------------------
\brief hoffman yield criterion

<pre>                                                            sh 03/03
This function calculates the hoffman yield criterion
f = 1/2 * (ST*P*S)  +  (ST*Q) - (uniax*uniax) 
</pre>
\param  DOUBLE  stress[6] (i)  stresses  
\param  DOUBLE  P[6][6]   (i)  coupling matrix P             
\param  DOUBLE  Q[6]      (i)  coupling vector Q 
\param  DOUBLE  uniax     (i)  uniaxial yield stress 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
DOUBLE mat_pl_hoff_yilcr(DOUBLE   stress[6], 
                         DOUBLE   P[6][6],
                         DOUBLE   Q[6],
                         DOUBLE   uniax)
{
INT    i,j;
DOUBLE ft = 0.;
DOUBLE PS[6];
DOUBLE SPS = 0.;
DOUBLE QS = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_yilcr");
#endif
/*---------------- initialize  -----------------------------------------*/
  for (i=0; i<6; i++) PS[i] = 0.0;
/*----------------------------------------------------------------------*/

  for (i=0; i<6; i++) for (j=0; j<6; j++) PS[i] += P[i][j] * stress[j];
  
  for (i=0; i<6; i++)  SPS += stress[i] * PS[i];
  for (i=0; i<6; i++)  QS  += Q[i] * stress[i];
  
  ft = 0.5 * SPS + QS - (uniax * uniax);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ft;
} /* end of mat_pl_hoff_yilcr */


/*!----------------------------------------------------------------------
\brief return algorithm for anisotropic plastic material model

<pre>                                                            sh 03/03
This routine projects the trial stresses back to the yield surface. In 
here is the return algorith for the anisotropic plastic material model, 
using a hoffman yield criterion and a anisotropic linear hardening law.
</pre>
\param  PL_HOFF   *mat         (i)  material properties  
\param  DOUBLE     sigtr[6]    (i)  trial stresses           
\param  DOUBLE   **C           (i)  material matrix 
\param  DOUBLE     P[6][6]     (i)  coupling matrix P for yield function 
\param  DOUBLE     Q[6]        (i)  coupling vector Q for yield function
\param  DOUBLE     rkappal[9]  (i)  hardening values             
\param  DOUBLE    *dhard       (i)  hardening variable 
\param  DOUBLE    *dlam        (o)  plastic multiplier 
\param  DOUBLE     signew[6]   (o)  projected stresses  
\param  DOUBLE     P_hat[6][6] (i)  coupling matrix P_hat for plastic potential           
\param  DOUBLE     Q_hat[6]    (i)  coupling vector Q_hat for plastic potential
\param  DOUBLE     gamma[6]    (i)   
\param  DOUBLE     dkappa[6]   (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_radi(PL_HOFF   *mat,
                      DOUBLE     sigtr[6], 
                      DOUBLE   **C,
                      DOUBLE     P[6][6],
                      DOUBLE     Q[6],
                      DOUBLE     rkappal[9],
                      DOUBLE    *dhard,
                      DOUBLE    *dlam,
                      DOUBLE     signew[6],
                      DOUBLE     P_hat[6][6],
                      DOUBLE     Q_hat[6],
                      DOUBLE     gamma[6],
                      DOUBLE     dkappa[6])
{
INT      i,j, cc, irc;
INT      iter=0;
INT      inine=9, isix=6, ione=1;
DOUBLE   dzero = 0.;
DOUBLE   ft;
DOUBLE   sig[6];
DOUBLE   rkappa[9];
DOUBLE   dfdl;
DOUBLE   dfdl1;
DOUBLE   dfdl2;
DOUBLE   H[6][6];
DOUBLE   HICPS[6], TKTL[6], HELP[6], ETA[9], RMK[9], DKDS[6][6], DKDL[6];
DOUBLE   CQ[6], PS[6], CPS[6], DGDS[6], DFDS[6], DSDL[6], gamma2[6]; 
DOUBLE   C_vec[36], P_hat_vec[36], CP_vec[36], ma66h[36], HI_vec[36];
DOUBLE   P_vec[36], RM[54], RM_DKDS[54], DKDS_vec[36];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_radi");
#endif
/*----------------------------------------------------------------------*/
  *dlam = 0.;
  for (i=0; i<6; i++)   signew[i] = 0.0;
  for (i=0; i<9; i++)   rkappa[i] = 0.0;
  for (i=0; i<54; i++)  RM[i]     = 0.0;

/*------- matrix -> vector ---------------------------------------------*/
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) C_vec[cc++]     = C[i][j];  
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) P_hat_vec[cc++] = P_hat[i][j];  
/*------- iteration scheme to obtain increment of plastic multiplier ---*/
L500:
  ++iter;
/*-------calculation of new stresses: S = S* - dlam*C*dG/dS ------------*/
/*                           [I+DLAMBDA*P*C]^(-1)*(STRIAL-DLAMBDA*C*Q)  */
/*                      HIER:    H^1*(STRIAL-DLAMBDA*C*Q)               */

  c1ab (C_vec,Q_hat,CQ,&isix,&isix,&ione,&dzero);
  c1ab (C_vec,P_hat_vec,CP_vec,&isix,&isix,&isix,&dzero);

  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) H[i][j] = *dlam * CP_vec[cc++];
  H[0][0] = 1. + H[0][0];
  H[1][1] = 1. + H[1][1];
  H[2][2] = 1. + H[2][2];
  H[3][3] = 1. + H[3][3];
  H[4][4] = 1. + H[4][4];
  H[5][5] = 1. + H[5][5];

  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = H[i][j];  
  c1inv6(ma66h,HI_vec,&irc); /*fortran routine, uses a vector*/

  for (i=0; i<6; i++)  sig[i] = sigtr[i] - *dlam * CQ[i]; 

  c1ab (HI_vec,sig,signew,&isix,&isix,&ione,&dzero);

/*------------------- dGdS = dG/dS ------------------------------------*/
  c1ab (P_hat_vec,signew,PS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) DGDS[i] = PS[i] + Q_hat[i];  

/*------------------- hardening law -----------------------------------*/
  mat_pl_hoff_hard(DGDS,dlam,rkappal,rkappa,dkappa,iter,RM);

/*------------ caluclate material matrix due to updated kappa ---------*/
  mat_pl_hoff_homa(mat, P, Q, rkappa);
  /*------- matrix -> vector ------------------------------------------*/
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) P_vec[cc++]     = P[i][j];  

/*---------yield condition - Hoffman criterion  -----------------------*/
  ft = mat_pl_hoff_yilcr(signew, P, Q, mat->uniax);

/*--------- calculate some needed values  -----------------------------*/
  mat_pl_hoff_serv(DKDL,mat,rkappa,P_hat,Q_hat,signew,ETA,dlam,DKDS);

/*--------- derivative of yield criteria with respect to dlam  --------*/
  /*------------------------------------------------  DFDS=DPHI/DSIGMA */
  c1ab (P_vec,signew,PS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) DFDS[i] = PS[i] + Q[i];  
  /*---------------------------------------------- DSDL=DSIGMA/DLAMBDA */
  c1ab (CP_vec,sig,CPS,&isix,&isix,&ione,&dzero);

  c1ab (HI_vec,CPS,HICPS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) CPS[i] = HICPS[i] + CQ[i];  
  c1ab (HI_vec,CPS,DSDL,&isix,&isix,&ione,&dzero);

  for (i=0; i<6; i++) DSDL[i] = -DSDL[i];  
  /*------------------------------------------------ DFDL=DPHI/DLAMBDA */
     /*----- CALCULATE FIRST PART OF DPHI/DLAMBDA */
  dfdl1 = 0.0;
  for (i=0; i<6; i++) dfdl1 += DFDS[i] * DSDL[i];  
     /*---- CALCULATE SECOND PART OF DPHI/DLAMBDA TKTL = DKDL + DKDS * DSDL */
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) DKDS_vec[cc++]     = DKDS[i][j];  

  c1ab (DKDS_vec,DSDL,HELP,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) TKTL[i] = DKDL[i] + HELP[i];  
            /*------------ ETA*M*TKTL ---*/
  c1ab (RM,TKTL,RMK,&inine,&isix,&ione,&dzero);
  dfdl2 = 0.0;
  for (i=0; i<9; i++) dfdl2 += ETA[i] * RMK[i];  

  dfdl = dfdl1 + dfdl2;

/*------------------------------------------------ new dlam  --------*/
  *dlam = *dlam - (ft/dfdl);

/* check convergence */
  if (fabs(ft/mat->uniax) > EPS8) 
  {
     if (iter > 50) 
     {
              dserror("local iteration in 'radi_hoff' exceeds limit");
     }
     goto L500; 
  }
/*------------------------------------------------- update of sigma ----*/
  c1ab (C_vec,Q_hat,CQ,&isix,&isix,&ione,&dzero);
  c1ab (C_vec,P_hat_vec,CP_vec,&isix,&isix,&isix,&dzero);

  cc=0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) H[i][j] = *dlam * CP_vec[cc++];
  H[0][0] = 1. + H[0][0];
  H[1][1] = 1. + H[1][1];
  H[2][2] = 1. + H[2][2];
  H[3][3] = 1. + H[3][3];
  H[4][4] = 1. + H[4][4];
  H[5][5] = 1. + H[5][5];

  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) ma66h[cc++] = H[i][j];  
  c1inv6(ma66h,HI_vec,&irc); /*fortran routine, uses a vector*/

  for (i=0; i<6; i++)  sig[i] = sigtr[i] - *dlam * CQ[i]; 

  c1ab (HI_vec,sig,signew,&isix,&isix,&ione,&dzero);

/*--------------- update of dGdS = dG/dS ------------------------------*/
  c1ab (P_hat_vec,signew,PS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) DGDS[i] = PS[i] + Q_hat[i];  

/*----------------update of rkappe -> hardening law -------------------*/
  mat_pl_hoff_hard(DGDS,dlam,rkappal,rkappa,dkappa,iter,RM);
  for (i=0; i<9; i++)  rkappal[i] = rkappa[i]; 

/*---- caluclate updated material matrix due to updated kappa ---------*/
  mat_pl_hoff_homa(mat, P, Q, rkappa);
  /*------- matrix -> vector ------------------------------------------*/
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) P_vec[cc++]     = P[i][j];  

/*------------------ update of DKDS, DKDL  ----------------------------*/
  mat_pl_hoff_serv(DKDL,mat,rkappa,P_hat,Q_hat,signew,ETA,dlam,DKDS);

/*----------------------------------------------- update of DFDS ------*/
  c1ab (P_vec,signew,PS,&isix,&isix,&ione,&dzero);
  for (i=0; i<6; i++) DFDS[i] = PS[i] + Q[i];  

/*--------------------------------------- calculate dhard for Cep -----*/
  c1ab (RM,DKDL,RMK,&inine,&isix,&ione,&dzero);
  *dhard = 0.;
  for (i=0; i<9; i++) *dhard += ETA[i] * RMK[i];  

/*--------------------------------------- calculate gamma for Cep -----*/
  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<6; i++) DKDS_vec[cc++]     = DKDS[i][j];  

  c1ab (RM,DKDS_vec,RM_DKDS,&inine,&isix,&isix,&dzero);
  c1ab (ETA,RM_DKDS,gamma2,&ione,&inine,&isix,&dzero);
  for (i=0; i<6; i++) gamma[i] = DFDS[i] + gamma2[i];  



/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ;
} /* end of mat_pl_hoff_radi */


/*!----------------------------------------------------------------------
\brief hardening law for anisotropic plastic material

<pre>                                                            sh 03/03
anisotropic linear hardening law for anisotropic plastic material model
based on the hoffman yield criterion
</pre>
\param  DOUBLE     DGDS[6]        (i)    
\param  DOUBLE    *dlam           (i)  plastic multiplier         
\param  DOUBLE     rkappal[9]     (i)  
\param  DOUBLE     rkappa[9]     (i/o) to be modified -> update 
\param  DOUBLE     dkappa_hat[6]  (i)   
\param  INT        iter           (i)  number of local iteration in "rad" 
\param  DOUBLE     RM[54]         (o)  mapping matrix 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_hard(DOUBLE     DGDS[6],
                      DOUBLE    *dlam,
                      DOUBLE     rkappal[9],
                      DOUBLE     rkappa[9],
                      DOUBLE     dkappa_hat[6],
                      INT        iter,
                      DOUBLE     RM[54])
{
INT    i,j, cc;
INT    f1,f2,f3,f4,f5,f6;
INT    inine=9, isix=6, ione=1;
DOUBLE dzero = 0.;
DOUBLE dkappa[9];
DOUBLE RMa[9][6];
DOUBLE epspl1[6], epspl2[6];
DOUBLE norm = 0.;
DOUBLE deps = 0.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_hard");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<54; i++)                    RM[i]     = 0.0;
for (i=0; i<9; i++) for (j=0; j<6; j++) RMa[i][j] = 0.0;

/*------------------------------------------ skalar hardening parameter */
for (i=0; i<6; i++)  epspl1[i] = *dlam * DGDS[i];
for (i=0; i<6; i++)  epspl2[i] = *dlam * DGDS[i];

epspl2[3] = epspl2[3] / 2.;
epspl2[4] = epspl2[4] / 2.;
epspl2[5] = epspl2[5] / 2.;

for (i=0; i<6; i++)  deps += epspl1[i] * epspl2[i];
deps = (2./3.) * deps;
deps = sqrt(deps);

/*------------------------------------------- NORM OF DPOTENTIAL/DSIGMA */       
for (i=0; i<6; i++)  norm += DGDS[i] * DGDS[i];
norm = sqrt(norm);
/*--------------------------------------- VECTOR OF INTERNAL PARAMETERS */      
for (i=0; i<6; i++)  dkappa_hat[i]= DGDS[i] * deps/norm;

if (iter == 1) for (i=0; i<6; i++) dkappa_hat[i] = DGDS[i];

if (dkappa_hat[0] > 0)
{
  f1 = 1;
  f2 = 0;
}
else
{
  f1 = 0;
  f2 = -1;
}
if (dkappa_hat[1] > 0)
{
  f3 = 1;
  f4 = 0;
}
else
{
  f3 = 0;
  f4 = -1;
}
if (dkappa_hat[2] > 0)
{
  f5 = 1;
  f6 = 0;
}
else
{
  f5 = 0;
  f6 = -1;
}

RMa[0][0] = f1;
RMa[1][0] = f2;
RMa[2][1] = f3;
RMa[3][1] = f4;
RMa[4][2] = f5;
RMa[5][2] = f6;

if (dkappa_hat[3] > 0)
{
  RMa[6][3] = 1;
}
else
{
  RMa[6][3] = -1;
}
if (dkappa_hat[4] > 0)
{
  RMa[7][4] = 1;
}
else
{
  RMa[7][4] = -1;
}
if (dkappa_hat[5] > 0)
{
  RMa[8][5] = 1;
}
else
{
  RMa[8][5] = -1;
}

if (iter == 1) for (i=0; i<6; i++) dkappa_hat[i] = 0.;

  cc = 0;
  for (j=0; j<6; j++) for (i=0; i<9; i++) RM[cc++] = RMa[i][j];  
  c1ab (RM,dkappa_hat,dkappa,&inine,&isix,&ione,&dzero);

  for (i=0; i<9; i++)  rkappa[i]= rkappal[i] + dkappa[i];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ;
} /* end of mat_pl_hoff_hard */


/*!----------------------------------------------------------------------
\brief some needed prevalues for the anisotropic plastic model

<pre>                                                            sh 03/03
this routine calculates some needed prevalues for anisotropic plastic 
material model based on the hoffman yield criterion
</pre>
\param  DOUBLE     DKDL[6]       (o)  partial differentiation of kappa do dlam 
\param  PL_HOFF   *mat           (i)  material properties
\param  DOUBLE     rkappa[9]     (i) 
\param  DOUBLE     P[6][6]       (i)  coupling matrix P for yield function 
\param  DOUBLE     Q[6]          (i)  coupling vector Q for yield function
\param  DOUBLE     sig[6]        (i)  stresses
\param  DOUBLE     eta[9]        (o)  partial differentiation of yield function to kappa
\param  DOUBLE     DKDS[6][6]    (o)  partial differentiation of kappa to stress         

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: mat_pl_hoff_main()     [mat_pl_hoff.c]

*----------------------------------------------------------------------*/
void mat_pl_hoff_serv(DOUBLE     DKDL[6],
                      PL_HOFF   *mat,
                      DOUBLE     rkappa[9],
                      DOUBLE     P[6][6],
                      DOUBLE     Q[6],
                      DOUBLE     sig[6],
                      DOUBLE     eta[9],
                      DOUBLE    *dlam,
                      DOUBLE     DKDS[6][6])
{
INT     i,j,k,m,n, cc;
INT     isix=6, ione=1;
DOUBLE  dzero = 0.;
DOUBLE  nenner, zaehler, faktor;
DOUBLE  uniax;
DOUBLE  SIG1, SIG2, SIG3, SIG4, SIG5, SIG6, SIG7, SIG8, SIG9;
DOUBLE  sig_I[9];
DOUBLE  sig_H[9];
DOUBLE  H[9];
DOUBLE  DPDK[6][6][9], DPDA[6][6][9];
DOUBLE  DQDK[6][9], DQDA[6][9], DADS[9][9], DSDK[9][9];
DOUBLE  PS[6], PSQ[6], PSQ1[6];
DOUBLE  P_vec[36];
DOUBLE  XX[6][6],XX1[6][6], T[6][6]; 
DOUBLE  XX_vec[36], XX1_vec[36], T_vec[36], XP_vec[36], XX1T_vec[36];
DOUBLE  DKDS_vec[36];
DOUBLE  DPDK_RED[6][6], DQDK_RED[6];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_pl_hoff_serv");
#endif
/*--------- initialize  ------------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) for (k=0; k<9; k++) DPDK[i][j][k] = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++) for (k=0; k<9; k++) DPDA[i][j][k] = 0.0;
for (i=0; i<6; i++) for (j=0; j<9; j++)                     DQDK[i][j]    = 0.0;
for (i=0; i<6; i++) for (j=0; j<9; j++)                     DQDA[i][j]    = 0.0;
for (i=0; i<9; i++) for (j=0; j<9; j++)                     DADS[i][j]    = 0.0;
for (i=0; i<9; i++)                                         sig_I[i]      = 0.0;
for (i=0; i<9; i++)                                         sig_H[i]      = 0.0;
for (i=0; i<9; i++)                                         H[i]          = 0.0;
for (i=0; i<9; i++)                                         eta[i]        = 0.0;
for (i=0; i<9; i++) for (j=0; j<9; j++)                     DSDK[i][j]    = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++)                     XX[i][j]      = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++)                     XX1[i][j]     = 0.0;
for (i=0; i<6; i++) for (j=0; j<6; j++)                     T[i][j]       = 0.0;
/*--------- get material properties ------------------------------------*/
uniax = mat->uniax;

sig_I[0] = mat->s11T; 
sig_I[1] = mat->s11C; 
sig_I[2] = mat->s22T; 
sig_I[3] = mat->s22C; 
sig_I[4] = mat->s33T; 
sig_I[5] = mat->s33C; 
sig_I[6] = mat->s12;
sig_I[7] = mat->s23;
sig_I[8] = mat->s13;

sig_H[0] = mat->sh11T;
sig_H[1] = mat->sh11C;
sig_H[2] = mat->sh22T;
sig_H[3] = mat->sh22C;
sig_H[4] = mat->sh33T;
sig_H[5] = mat->sh33C;
sig_H[6] = mat->sh12; 
sig_H[7] = mat->sh23; 
sig_H[8] = mat->sh13; 

H[0] = mat->ha11T;
H[1] = mat->ha11C;
H[2] = mat->ha22T;
H[3] = mat->ha22C;
H[4] = mat->ha33T;
H[5] = mat->ha33C;
H[6] = mat->ha12; 
H[7] = mat->ha23; 
H[8] = mat->ha13; 

/*--------- updated yield stress ---------------------------------------*/
SIG1 = sig_I[0] + (sig_H[0]-sig_I[0])/(H[0]) * rkappa[0];
SIG2 = sig_I[1] + (sig_H[1]-sig_I[1])/(H[1]) * rkappa[1];
SIG3 = sig_I[2] + (sig_H[2]-sig_I[2])/(H[2]) * rkappa[2];
SIG4 = sig_I[3] + (sig_H[3]-sig_I[3])/(H[3]) * rkappa[3];
SIG5 = sig_I[4] + (sig_H[4]-sig_I[4])/(H[4]) * rkappa[4];
SIG6 = sig_I[5] + (sig_H[5]-sig_I[5])/(H[5]) * rkappa[5];
SIG7 = sig_I[6] + (sig_H[6]-sig_I[6])/(H[6]) * rkappa[6];
SIG8 = sig_I[7] + (sig_H[7]-sig_I[7])/(H[7]) * rkappa[7];
SIG9 = sig_I[8] + (sig_H[8]-sig_I[8])/(H[8]) * rkappa[8];

/*-------------------------------------------- DP/DALPHAJ DPDA(6,6,J) --*/ 
DPDA[0][0][3] =  1.;
DPDA[0][1][3] = -1.;
DPDA[1][0][3] = -1.;
DPDA[1][1][3] =  1.;
       
DPDA[1][1][4] =  1.;
DPDA[1][2][4] = -1.;
DPDA[2][1][4] = -1.;
DPDA[2][2][4] =  1.;
       
DPDA[0][0][5] =  1.;
DPDA[0][2][5] = -1.;
DPDA[2][0][5] = -1.;
DPDA[2][2][5] =  1.;
       
DPDA[3][3][6] =  6.;
       
DPDA[4][4][7] =  6.;
       
DPDA[5][5][8] =  6.;

/*----------------------------------------------- DALPHAJ/DSIGMAK ------*/       
DADS[0][0] = -uniax*uniax/(SIG2*SIG1)*(1.+(SIG2-SIG1)/SIG1);
DADS[0][1] = -uniax*uniax/(SIG2*SIG1)*(-1.+(SIG2-SIG1)/(SIG2));

DADS[1][2] = -uniax*uniax/(SIG4*SIG3)*(1.+(SIG4-SIG3)/(SIG3));
DADS[1][3] = -uniax*uniax/(SIG4*SIG3)*(-1.+(SIG4-SIG3)/(SIG4));

DADS[2][4] = -uniax*uniax/(SIG5*SIG6)*(1.+(SIG6-SIG5)/(SIG5));
DADS[2][5] = -uniax*uniax/(SIG5*SIG6)*(-1.+(SIG6-SIG5)/(SIG6));

DADS[3][0] = -uniax*uniax/(SIG2*SIG1*SIG1);
DADS[3][1] = -uniax*uniax/(SIG2*SIG2*SIG1);
DADS[3][2] = -uniax*uniax/(SIG4*SIG3*SIG3);
DADS[3][3] = -uniax*uniax/(SIG4*SIG4*SIG3);
DADS[3][4] =  uniax*uniax/(SIG6*SIG5*SIG5);
DADS[3][5] =  uniax*uniax/(SIG6*SIG6*SIG5);

DADS[4][0] =  uniax*uniax/(SIG2*SIG1*SIG1);
DADS[4][1] =  uniax*uniax/(SIG2*SIG2*SIG1);
DADS[4][2] = -uniax*uniax/(SIG4*SIG3*SIG3);
DADS[4][3] = -uniax*uniax/(SIG4*SIG4*SIG3);
DADS[4][4] = -uniax*uniax/(SIG6*SIG5*SIG5);
DADS[4][5] = -uniax*uniax/(SIG6*SIG6*SIG5);

DADS[5][0] = -uniax*uniax/(SIG2*SIG1*SIG1);
DADS[5][1] = -uniax*uniax/(SIG2*SIG2*SIG1);
DADS[5][2] =  uniax*uniax/(SIG4*SIG3*SIG3);
DADS[5][3] =  uniax*uniax/(SIG4*SIG4*SIG3);
DADS[5][4] = -uniax*uniax/(SIG6*SIG5*SIG5);
DADS[5][5] = -uniax*uniax/(SIG6*SIG6*SIG5);

DADS[6][6] = -2*uniax*uniax/(3*SIG7*SIG7*SIG7);

DADS[7][7] = -2*uniax*uniax/(3*SIG8*SIG8*SIG8);

DADS[8][8] = -2*uniax*uniax/(3*SIG9*SIG9*SIG9);
      
/*-------------------------------------------------- DSIGMAK/DKAPPAI ---*/       
for (i=0; i<9; i++)  DSDK[i][i] = (sig_H[i]-sig_I[i])/(H[i]);

/*----------------------------------------- DQ/DALPHAJ   DQDA(6,J) -----*/
DQDA[0][0] = 1.;
DQDA[1][1] = 1.;
DQDA[2][2] = 1.;

/*===================================================================== */
/*---------------------------------------------DP/DKAPPA=DPDA*DADS*DSDK */
/*===================================================================== */
for (i=0; i<9; i++)
   for (j=0; j<6; j++)
      for (n=0; n<6; n++)
         for (k=3; k<9; k++)
            for (m=0; m<9; m++) 
                 DPDK[j][n][i] = DPDK[j][n][i] + DPDA[j][n][k] * DADS[k][m] * DSDK[m][i];

/*===================================================================== */
/*---------------------------------------------DQ/DKAPPA=DQDA*DADS*DSDK */
/*===================================================================== */
for (i=0; i<9; i++)
   for (j=0; j<6; j++)
      for (n=0; n<3; n++)
         for (m=0; m<9; m++) 
                 DQDK[j][i] = DQDK[j][i] + DQDA[j][n] * DADS[n][m] * DSDK[m][i];

/*===================================================================== */
/*------------------------------------------------------ DKAPPA/DLAMBDA */     
/*===================================================================== */
cc = 0;
for (j=0; j<6; j++) for (i=0; i<6; i++) P_vec[cc++] = P[i][j];  

c1ab (P_vec,sig,PS,&isix,&isix,&ione,&dzero);
for (i=0; i<6; i++)      PSQ[i] = PS[i] + Q[i];

nenner = 0.0;
for (i=0; i<6; i++) nenner += PSQ[i] * PS[i];
nenner = sqrt(nenner);

for (i=0; i<6; i++) PSQ1[i] = PSQ[i];
PSQ1[3] = 0.5 * PSQ1[3];
PSQ1[4] = 0.5 * PSQ1[4];
PSQ1[5] = 0.5 * PSQ1[5];

zaehler = 0.0;
for (i=0; i<6; i++) zaehler += PSQ[i] * PSQ1[i];
zaehler = 2./3. * zaehler;
zaehler = sqrt(zaehler);

for (i=0; i<6; i++) DKDL[i] = PSQ[i] * zaehler/nenner;

/*===================================================================== */
/*------------------------------------------------------- DKAPPA/DSIGMA */     
/*===================================================================== */
for (i=0; i<6; i++) for(j=0; j<6; j++) XX[i][j] = PSQ[i] * PSQ[j];

cc = 0;
for (j=0; j<6; j++) for (i=0; i<6; i++) XX_vec[cc++] = XX[i][j];  

c1ab (XX_vec,P_vec,XP_vec,&isix,&isix,&isix,&dzero);

nenner = 1./(nenner * nenner);

for (i=0; i<36; i++) DKDS_vec[i] = P_vec[i] - nenner * XP_vec[i];

for (i=0; i<6; i++) for(j=0; j<6; j++) XX1[i][j] = XX[i][j];
cc = 0;
for (j=0; j<6; j++) for (i=0; i<6; i++) XX1_vec[cc++] = XX1[i][j];  

T[0][0] = 2./3.;
T[1][1] = 2./3.;
T[2][2] = 2./3.;
T[3][3] = 1./3.;
T[4][4] = 1./3.;
T[5][5] = 1./3.;

cc = 0;
for (j=0; j<6; j++) for (i=0; i<6; i++) T_vec[cc++] = T[i][j];  

c1ab (XX1_vec,T_vec,XX1T_vec,&isix,&isix,&isix,&dzero);
c1ab (XX1T_vec,P_vec,XP_vec,&isix,&isix,&isix,&dzero);

faktor = 1./(zaehler * zaehler);

for (i=0; i<36; i++) XP_vec[i] = faktor * XP_vec[i];

for (i=0; i<36; i++) DKDS_vec[i] = DKDS_vec[i] + XP_vec[i];
for (i=0; i<36; i++) DKDS_vec[i] = *dlam * zaehler * DKDS_vec[i];

nenner = 0.0;
for (i=0; i<6; i++) nenner += PSQ[i] * PSQ[i];
nenner = sqrt(nenner);
nenner = 1./nenner;
for (i=0; i<36; i++) DKDS_vec[i] = nenner * DKDS_vec[i];

cc = 0;
for (j=0; j<6; j++) for (i=0; i<6; i++) DKDS[i][j] = DKDS_vec[cc++];  

/*===================================================================== */
/*----------------------------------------------------------------- ETA */
/*   ETAI = 1/2*SIG*DPDKI*SIG + SIG*DQDKI                               */
/*===================================================================== */
for (i=0; i<9; i++) 
{
  for (j=0; j<6; j++)
  {
    for (k=0; k<6; k++) DPDK_RED[j][k] = DPDK[j][k][i]; 
    DQDK_RED[j] = DQDK[j][i];
  } 
  eta[i] = mat_pl_hoff_yilcr(sig, DPDK_RED, DQDK_RED, 0.0);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return ;
} /* end of mat_pl_hoff_serv */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/
