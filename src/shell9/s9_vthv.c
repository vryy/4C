/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_vthv: which calculates the modifications to the metrics due to EAS 
           (additional strains Etilde)


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief modified metrics due to EAS                                      

<pre>                     m.gee 11/01              modified by   sh 02/03
This routine calculates the modification to the metrics in the shell body
due to EAS.
</pre>
\param  DOUBLE   **gmkovc   (o)  kovariant metric at specified location in shell body (cur. config.)
\param  DOUBLE   **gmkonc   (o)  kontravariant metric at specified location in shell body (cur. config.)
\param  DOUBLE    *epsh     (i)  strains from EAS (E_tilde) of this kin. layer
\param  DOUBLE    *detc     (o)  determinant of gkovc (cur. config.)
\param  DOUBLE     e3       (i)  zeta_3_L (within one material layer)
\param  DOUBLE     h        (i)  total thickness of this element
\param  DOUBLE    *klayhgt  (i)  hight of kin layer in % of total thickness of shell 
\param  DOUBLE    *mlayhgt  (i)  hight of mat layer in % of adjacent kin layer 
\param  INT        num_klay (i)  number of kin layers to this element  
\param  INT        klay     (i)  actual kin layer 
\param  INT        mlay     (i)  actual mat layer of this kin layer 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_vthv(DOUBLE **gmkovc,
             DOUBLE **gmkonc,
             DOUBLE  *epsh,
             DOUBLE  *detc,
             DOUBLE   e3,
             DOUBLE   h,         /* total thickness of this element */
             DOUBLE  *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE  *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT      num_klay,  /* number of kin layers to this element */
             INT      klay,      /* actual kin layer */
             INT      mlay,      /* actual mat layer of this kin layer */
             DOUBLE   condfac)
{
INT i;
DOUBLE deltah, h_mlay, h_kl;
DOUBLE zeta_kl,zeta,det_dummy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_vthv");
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

zeta = zeta_kl;
zeta = s9con(zeta,num_klay,klay,klay,condfac);

  gmkovc[0][0] = gmkovc[0][0] + 2.0 * (epsh[0]+zeta*epsh[6])    ;
  gmkovc[1][0] = gmkovc[1][0] +       (epsh[1]+zeta*epsh[7])    ;
  gmkovc[2][0] = gmkovc[2][0] +       (epsh[3]+zeta*epsh[9])    ;
  gmkovc[1][1] = gmkovc[1][1] + 2.0 * (epsh[2]+zeta*epsh[8])    ;
  gmkovc[2][1] = gmkovc[2][1] +       (epsh[4]+zeta*epsh[10])   ;
  gmkovc[2][2] = gmkovc[2][2] + 2.0 * (epsh[5]+zeta*epsh[11])   ;

gmkovc[0][2] = gmkovc[2][0];
gmkovc[1][2] = gmkovc[2][1];
gmkovc[0][1] = gmkovc[1][0];
/*----------------------------------------------------------------------*/
math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = -det_dummy;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_vthv */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
