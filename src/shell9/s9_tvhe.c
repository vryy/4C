/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_tvhe: re-calculates the metrics (geom. linear/nonlinear)
 - s9_tvhe_lin: re-calculates the metrics (geom. linear) -> not used

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief re-calculates the metrics (geom. linear/nonlinear)                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine re-calculates the metrics (geom. linear/nonlinear) to take 
the neglection of the in theta_3 quadratic parts into account. The 
metics 'gmkovc' and 'gmkonc' are re-calculated.
</pre>
\param  DOUBLE   **gmkovr   (i)  kovariant metric at specified location in shell body (ref. config.)
\param  DOUBLE   **gmkovc  (i/o) kovariant metric at specified location in shell body (cur. config.)
\param  DOUBLE   **gmkonr   (o)  kontravariant metric at specified location in shell body (ref. config.)
\param  DOUBLE   **gmkonc  (i/o) kontravariant metric at specified location in shell body (cur. config.)
\param  DOUBLE   **gkovr    (i)  kovariant basis vectors at specified location in shell body (ref. config.)
\param  DOUBLE   **gkovc    (i)  kovariant basis vectors at specified location in shell body (cur. config.)
\param  DOUBLE    *detr     (o)  determinant of gkovr (ref. config.)
\param  DOUBLE    *detc     (o)  determinant of gkovc (cur. config.)
\param  DOUBLE  ***amkovc   (i)  kovariant metric in reference layer of each kin. layer (cur. config.)
\param  DOUBLE  ***amkovr   (i)  kovariant metric in reference layer of each kin. layer (ref. config.)
\param  DOUBLE  ***akovc    (i)  kovariant basis vectors in reference layer of each kin. layer (cur. config.)
\param  DOUBLE  ***akovr    (i)  kovariant basis vectors in reference layer of each kin. layer (ref. config.)
\param  DOUBLE  ***a3kvpc   (i)  partial derivatives of a3_L of each kin. Layer (cur. config.)
\param  DOUBLE  ***a3kvpr   (i)  partial derivatives of a3_L of each kin. Layer (ref. config.)
\param  DOUBLE     e3       (i)  local theta3 of actual mat. layer
\param  INT        kintyp   (i)  kintyp=0: geo_lin; =1: upd_lagr; =2: tot_lagr 
\param  DOUBLE     h        (i)  total thickness of this element
\param  DOUBLE    *klayhgt  (i)  hight of kin layer in % of total thickness of shell 
\param  DOUBLE    *mlayhgt  (i)  hight of mat layer in % of adjacent kin layer 
\param  INT        num_klay (i)  number of kin layers to this element  
\param  INT        num_mlay (i)  number of mat layers to this kin layer 
\param  INT        klay     (i)  actual kin layer 
\param  INT        mlay     (i)  actual mat layer of this kin layer 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]
                             s9_stress()     [s9_stress.c]

*----------------------------------------------------------------------*/
void s9_tvhe(DOUBLE  **gmkovr,
             DOUBLE  **gmkovc,
             DOUBLE  **gmkonr,
             DOUBLE  **gmkonc,
             DOUBLE  **gkovr,
             DOUBLE  **gkovc,
             DOUBLE   *detr,
             DOUBLE   *detc,
             DOUBLE ***amkovc,
             DOUBLE ***amkovr,
             DOUBLE ***akovc,
             DOUBLE ***akovr,
             DOUBLE ***a3kvpc,
             DOUBLE ***a3kvpr,
             DOUBLE    e3,
             INT       kintyp,    /* typ of kinematic formulation */
             DOUBLE    h,         /* total thickness of this element */
             DOUBLE   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT       num_klay,  /* number of kin layers to this element */
             INT       num_mlay,  /* number of mat layers to this kin layer */
             INT       klay,      /* actual kin layer */
             INT       mlay,      /* actual mat layer of this kin layer */
             DOUBLE    condfac)
{
INT i,j,k,jlay;
DOUBLE b11c=0.0;
DOUBLE b12c=0.0;
DOUBLE b21c=0.0;
DOUBLE b22c=0.0;
DOUBLE b31c=0.0;
DOUBLE b32c=0.0;

DOUBLE b11r=0.0;
DOUBLE b12r=0.0;
DOUBLE b21r=0.0;
DOUBLE b22r=0.0;
DOUBLE b31r=0.0;
DOUBLE b32r=0.0;

DOUBLE deltah, h_mlay, h_kl;
DOUBLE zeta_kl,zeta,det_dummy;
DOUBLE heps[3][3];

#ifdef DEBUG 
dstrc_enter("s9_tvhe");
#endif
/*----------------------------------------------------------------------*/
if (kintyp > 0) /*geom. nonlinear*/
{
  gmkovc[0][0] = gmkovr[0][0] + (amkovc[0][0][klay]-amkovr[0][0][klay]);
  gmkovc[1][1] = gmkovr[1][1] + (amkovc[1][1][klay]-amkovr[1][1][klay]);
  gmkovc[2][2] = gmkovr[2][2] + (amkovc[2][2][klay]-amkovr[2][2][klay]);
  gmkovc[0][1] = gmkovr[0][1] + (amkovc[0][1][klay]-amkovr[0][1][klay]);
  gmkovc[0][2] = gmkovr[0][2] + (amkovc[0][2][klay]-amkovr[0][2][klay]);
  gmkovc[1][2] = gmkovr[1][2] + (amkovc[1][2][klay]-amkovr[1][2][klay]);
  /*----------------------------------------------------------------------*/
  /*- calculate zeta_kl of kinematic layer due to local zeta_ml of material layer -> old s9tmtr --*/
  h_kl   = (klayhgt[klay]/100.)*h;              /* absolute hight of the actual kinematic layer */
  deltah = (mlayhgt[mlay]/100.)*h_kl;           /* absolute hight of the actual material layer */
  h_mlay   = 0.0;
  for (i=0; i<=mlay; i++)                        /* sum of the absolute hights of the material layers  */
  {                                             /* within the actual kinematic layer up to the actual */
     h_mlay += (mlayhgt[i]/100.)*h_kl;          /* material layer */
  }
  zeta_kl = -1. + (-deltah*(1.-e3)+2.*h_mlay)/h_kl;  /* equals the theta3 value of the act. kinematic layer, see equation (5.45) in Dis. Braun */

  /*--- loop over all kinematic layers to get the continuity coeficients (Pic. 7.4 on p. 116 in Dis. Braun)-> old s9mmtr ------*/
  for (jlay=0; jlay<num_klay; jlay++)
  {
     zeta = zeta_kl;
     zeta = s9con(zeta,num_klay,klay,jlay,condfac);
     if (zeta != 0.0)  /*---- interpolation -------- g1,g2 (kov.) -------*/
     {
        b11c=0.0;
        b12c=0.0;
        b21c=0.0;
        b22c=0.0;
        b31c=0.0;
        b32c=0.0;

        b11r=0.0;
        b12r=0.0;
        b21r=0.0;
        b22r=0.0;
        b31r=0.0;
        b32r=0.0;
        
        for (i=0; i<3; i++) b11c += akovc[i][0][0]   *a3kvpc[i][0][jlay];
        for (i=0; i<3; i++) b12c += akovc[i][0][0]   *a3kvpc[i][1][jlay];
        for (i=0; i<3; i++) b21c += akovc[i][1][0]   *a3kvpc[i][0][jlay];
        for (i=0; i<3; i++) b22c += akovc[i][1][0]   *a3kvpc[i][1][jlay];
        for (i=0; i<3; i++) b31c += akovc[i][2][klay]*a3kvpc[i][0][jlay];
        for (i=0; i<3; i++) b32c += akovc[i][2][klay]*a3kvpc[i][1][jlay];

        for (i=0; i<3; i++) b11r += akovr[i][0][0]   *a3kvpr[i][0][jlay];
        for (i=0; i<3; i++) b12r += akovr[i][0][0]   *a3kvpr[i][1][jlay];
        for (i=0; i<3; i++) b21r += akovr[i][1][0]   *a3kvpr[i][0][jlay];
        for (i=0; i<3; i++) b22r += akovr[i][1][0]   *a3kvpr[i][1][jlay];
        for (i=0; i<3; i++) b31r += akovr[i][2][klay]*a3kvpr[i][0][jlay];
        for (i=0; i<3; i++) b32r += akovr[i][2][klay]*a3kvpr[i][1][jlay];

        gmkovc[0][0] = gmkovc[0][0] + zeta*2.0*(b11c-b11r);
        gmkovc[1][1] = gmkovc[1][1] + zeta*2.0*(b22c-b22r);
        gmkovc[2][2] = gmkovc[2][2] ;
        gmkovc[0][1] = gmkovc[0][1] + zeta*(b21c+b12c-b21r-b12r);
        gmkovc[0][2] = gmkovc[0][2] + zeta*(b31c-b31r);
        gmkovc[1][2] = gmkovc[1][2] + zeta*(b32c-b32r);
     }
  } /*================================ end loop over all layers ========*/
}/*======== end geom. nonlinear ========================================*/
else /*geom. linear*/
{
  for (i=0; i<3; i++)
  for (j=0; j<3; j++)
  {
     heps[i][j]=0.0;
     for (k=0; k<3; k++)  heps[i][j] += gkovc[k][i] * gkovr[k][j];
  }
  /*--------------------------------------------------------------------*/
  gmkovc[0][0] = 2.0 * heps[0][0] - gmkovr[0][0];
  gmkovc[1][1] = 2.0 * heps[1][1] - gmkovr[1][1];
  gmkovc[2][2] = 2.0 * heps[2][2] - gmkovr[2][2];

  gmkovc[0][1] = heps[0][1]+heps[1][0] - gmkovr[0][1];
  gmkovc[0][2] = heps[0][2]+heps[2][0] - gmkovr[0][2];
  gmkovc[1][2] = heps[1][2]+heps[2][1] - gmkovr[1][2];
}/*=========== end geom. linear ========================================*/ 

gmkovc[2][0] = gmkovc[0][2]; 
gmkovc[2][1] = gmkovc[1][2]; 
gmkovc[1][0] = gmkovc[0][1]; 

/*----------------------------------------------------------------------*/
math_array_copy(gmkovr,3,3,gmkonr);
math_inv3(gmkonr,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detr = sqrt(det_dummy);

math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tvhe */



/*!----------------------------------------------------------------------
\brief re-calculates the metrics (geom. linear)                                      

<pre>                                                          m.gee 6/01              
This routine re-calculates the metrics (geom. linear) to take 
the neglection of the in theta_3 quadratic parts into account. The 
metics 'gmkovc' and 'gmkonc' are re-calculated. This Routine is not used
as it is part of 's9_tvhe'!
</pre>
\param  DOUBLE   **gmkovr   (i)  kovariant metric at specified location in shell body (ref. config.)
\param  DOUBLE   **gmkovc  (i/o) kovariant metric at specified location in shell body (cur. config.)
\param  DOUBLE   **gmkonr   (o)  kontravariant metric at specified location in shell body (ref. config.)
\param  DOUBLE   **gmkonc  (i/o) kontravariant metric at specified location in shell body (cur. config.)
\param  DOUBLE   **gkovr    (i)  kovariant basis vectors at specified location in shell body (ref. config.)
\param  DOUBLE   **gkovc    (i)  kovariant basis vectors at specified location in shell body (cur. config.)
\param  DOUBLE    *detr     (o)  determinant of gkovr (ref. config.)
\param  DOUBLE    *detc     (o)  determinant of gkovc (cur. config.)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ---------

*----------------------------------------------------------------------*/
void s9_tvhe_lin(DOUBLE  **gmkovr,
                 DOUBLE  **gmkovc,
                 DOUBLE  **gmkonr,
                 DOUBLE  **gmkonc,
                 DOUBLE  **gkovr,
                 DOUBLE  **gkovc,
                 DOUBLE   *detr,
                 DOUBLE   *detc)
{
INT i,j,k;

DOUBLE heps[3][3];
DOUBLE det_dummy;


#ifdef DEBUG 
dstrc_enter("s9_tvhe_lin");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   heps[i][j]=0.0;
   for (k=0; k<3; k++)  heps[i][j] += gkovc[k][i] * gkovr[k][j];
}
/*----------------------------------------------------------------------*/
gmkovc[0][0] = 2.0 * heps[0][0] - gmkovr[0][0];
gmkovc[1][1] = 2.0 * heps[1][1] - gmkovr[1][1];
gmkovc[2][2] = 2.0 * heps[2][2] - gmkovr[2][2];

gmkovc[0][1] = heps[0][1]+heps[1][0] - gmkovr[0][1];
gmkovc[0][2] = heps[0][2]+heps[2][0] - gmkovr[0][2];
gmkovc[1][2] = heps[1][2]+heps[2][1] - gmkovr[1][2];

gmkovc[1][0] = gmkovc[0][1];
gmkovc[2][0] = gmkovc[0][2];
gmkovc[2][1] = gmkovc[1][2];
/*----------------------------------------------------------------------*/
math_array_copy(gmkovr,3,3,gmkonr);
math_inv3(gmkonr,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detr = sqrt(det_dummy);

math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tvhe_lin */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
