/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_tmtr: which calculates the metrics within the shell body (theta_3
            is variable)

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the metrics within the shell body for variable theta_3                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the metrics (gkov, gkon, gmkov, gmkon) within the
shell body for a variable theta_3 (thickness direction). For a multilayer
formulation, the e3-value has to be transformed into a 'unity'-thickness
format so that the gaussian integration is valid again. For more detail
see Dis. Braun (p.77 -> material layer & p.115f for kinematic layers)
</pre>
\param  DOUBLE     e3      (i)  zeta_3_L (within one material layer)
\param  DOUBLE   **gkov    (o)  kovariant basis vectors at specified location in shell body
\param  DOUBLE   **gkon    (o)  kontravariant basis vectors at specified location in shell body
\param  DOUBLE   **gmkov   (o)  kovariant metric at specified location in shell body
\param  DOUBLE   **gmkon   (o)  kontravariant metric at specified location in shell body
\param  DOUBLE    *det     (o)  determinant of gkov (-> Jacobi determinant)
\param  DOUBLE  ***akov    (i)  kovariant basis vectors in reference layer of each kinematic layer
\param  DOUBLE  ***a3kvp   (i)  partial derivatives of a3_L for each kinematic layer
\param  DOUBLE     h       (i)  total thickness of this element
\param  DOUBLE    *klayhgt (i)  hight of kin layer in % of total thickness of shell 
\param  DOUBLE    *mlayhgt (i)  hight of mat layer in % of adjacent kin layer 
\param  INT        num_klay(i)  number of kin layers to this element  
\param  INT        klay    (i)  actual kin layer 
\param  INT        mlay    (i)  actual mat layer of this kin layer 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]
                             s9_stress()     [s9_stress.c]
                             s9jaco()        [s9_jaco.c]

*----------------------------------------------------------------------*/
void s9_tmtr(DOUBLE     e3,
             DOUBLE   **gkov,
             DOUBLE   **gkon,
             DOUBLE   **gmkov,
             DOUBLE   **gmkon,
             DOUBLE    *det,
             DOUBLE  ***akov,
             DOUBLE  ***a3kvp,
             DOUBLE     h,           /* total thickness of this element */
             DOUBLE    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             DOUBLE    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             INT        num_klay,    /* number of kin layers to this element */  
             INT        klay,        /* actual kin layer */
             INT        mlay,        /* actual mat layer of this kin layer */
             DOUBLE     condfac)
{
INT    i,j,k,idim,jlay;
DOUBLE deltah, h_mlay, h_kl;
DOUBLE zeta_kl,zeta,det_dummy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_tmtr");
#endif
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

/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
{
   gkov[i][0] = 0.0;
   gkov[i][1] = 0.0;
}
/*--- loop over all kinematic layers to get the continuity coeficients (Pic. 7.4 on p. 116 in Dis. Braun)-> old s9mmtr ------*/
for (jlay=0; jlay<num_klay; jlay++)
{
   zeta = zeta_kl;
   zeta = s9con(zeta,num_klay,klay,jlay,condfac);
   if (zeta != 0.0)  /*---- interpolation -------- g1,g2 (kov.) -------*/
   {
      for(idim=0; idim<3; idim++)
      {
        gkov[idim][0] += zeta * a3kvp[idim][0][jlay];
        gkov[idim][1] += zeta * a3kvp[idim][1][jlay];
      }
   }
} /*================================ end loop over all layers ========*/

for (idim=0; idim<3; idim++)
{
   gkov[idim][0] += akov[idim][0][0];
   gkov[idim][1] += akov[idim][1][0];
}
/*------------------------- interpolation ---------- g3 --------------*/
for (idim=0; idim<3; idim++) gkov[idim][2] = akov[idim][2][klay];
/*--------------- kontravariant basis vectors g1,g2,g3 (inverse of kov) */
   math_array_copy(gkov,3,3,gkon);
   math_inv3(gkon,det);
   math_tran(gkon,3);
/*--------------------------------------------- kovariant metrik tensor */
for (i=0; i<3; i++)
{
   for (j=i; j<3; j++)
   {
      gmkov[i][j]=0.0;
      for (k=0; k<3; k++) 
      gmkov[i][j] += gkov[k][i]*gkov[k][j];
   }
}   
      gmkov[1][0] = gmkov[0][1];
      gmkov[2][0] = gmkov[0][2];
      gmkov[2][1] = gmkov[1][2];
/*----------------------------------------- kontravariant metrik tensor */
    math_array_copy(gmkov,3,3,gmkon);
    math_inv3(gmkon,&det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tmtr */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
