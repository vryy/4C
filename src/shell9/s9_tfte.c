/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9_tfte: calculates forces (stress resultants -> "Schnittgroessen"
            from stresses for each kinematic layer
            NOTE: not in use!! 

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates forces from stresses                                      

<pre>                     m.gee 12/01             modified by    sh 01/03
This routine calculates forces (stress resultants -> "Schnittgroessen")
from stresses.
NOTE: this routine is not in use right now, as the interpretation of 
      these results is not very easy -> rather stresses than forces
      for multilayer shell element!
......................................................................
.    BERECHNEN DER KRAEFTE UND MOMENTE AUS SPANNUNGEN                .
.    STRESS    -  SPANNUNGEN (KONTRAVARIANT)                         .
.    GMKOV     -  KOVARIANTER METRIKTENSOR  (SCHALENRAUM)            .
.    AMKOV     -  KOVARIANTER METRIKTENSOR  (SCHALENMITTELFLAECHE)   .
.    GMKON     -  KONTRAVARIANTER METRIKTENSOR (SCHALENRAUM)         .
.    H         -  DICKE AM GAUSSPUNKT                                .
.    E3        -  NORMIERTE DICKENKOORDINATE                         .
.    DETSM     -  'FLAECHENELEMENT' DER SMF                          .
.    DETSR     -      "             DES SCHALENRAUMS                 .
......................................................................
</pre>
\param  DOUBLE   **force   (o) forces to be calculated
\param  INT        ngauss  (i) ID of actual GP
\param  DOUBLE    *stress  (i) stresses at GP

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: not in use!!

*----------------------------------------------------------------------*/
void s9_tfte(DOUBLE  **force,
             INT       ngauss,
             DOUBLE   *stress,
             DOUBLE  **gkov,
             DOUBLE ***akon,
             DOUBLE ***amkov,
             DOUBLE    h,
             DOUBLE    e3,
             DOUBLE   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT       klay,      /* actual kin layer */
             INT       mlay,      /* actual mat layer of this kin layer */
             DOUBLE    fact,
             DOUBLE    detsm,
             DOUBLE    detsr)
{
INT    i,j;
DOUBLE s11,s12,s21,s13,s31,s22,s23,s32,s33;
DOUBLE xu;
DOUBLE hh;
DOUBLE wgt;
DOUBLE wgthe3;
DOUBLE xu11=0.0;
DOUBLE xu21=0.0;
DOUBLE xu31=0.0;
DOUBLE xu12=0.0;
DOUBLE xu22=0.0;
DOUBLE xu32=0.0;
DOUBLE xu13=0.0;
DOUBLE xu23=0.0;
DOUBLE xu33=0.0;

DOUBLE deltah, h_mlay, h_kl,zeta_kl;

#ifdef DEBUG 
dstrc_enter("s9_tfte");
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

s11 = stress[0];
s12 = stress[1];
s21 = s12;
s22 = stress[2];
s13 = stress[3]*2;
s31 = s13;
s23 = stress[4]*2;
s32 = s23;
s33 = stress[5]*4;
/*----------------------------------------------------------------------*/
xu     = (detsr/detsm);                /*shell shifter*/
hh     = sqrt(amkov[2][2][klay])*0.5;  /*-> norm of a3! */
wgt    = hh * fact * xu;
wgthe3 = wgt * zeta_kl * hh;
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++) xu11 += gkov[i][0]*akon[i][0][klay];
for (i=0; i<3; i++) xu12 += gkov[i][0]*akon[i][1][klay];
for (i=0; i<3; i++) xu13 += gkov[i][0]*akon[i][2][klay];
for (i=0; i<3; i++) xu21 += gkov[i][1]*akon[i][0][klay];
for (i=0; i<3; i++) xu22 += gkov[i][1]*akon[i][1][klay];
for (i=0; i<3; i++) xu23 += gkov[i][1]*akon[i][2][klay];
for (i=0; i<3; i++) xu31 += gkov[i][2]*akon[i][0][klay];
for (i=0; i<3; i++) xu32 += gkov[i][2]*akon[i][1][klay];
for (i=0; i<3; i++) xu33 += gkov[i][2]*akon[i][2][klay];
/*----------------------------------------------------------------------*/
/*..................................................................N11 */
force[0][ngauss]  += (s11*xu11+s12*xu21+s13*xu31) * wgt;
/*..................................................................N12 */
force[2][ngauss]  += (s11*xu12+s12*xu22+s13*xu32) * wgt;
/*..................................................................N21 */
force[8][ngauss]  += (s21*xu11+s22*xu21+s23*xu31) * wgt;
/*..................................................................N22 */
force[1][ngauss]  += (s21*xu12+s22*xu22+s23*xu32) * wgt;
/*...............................................................N13/Q1 */
force[3][ngauss]  += (s11*xu13+s12*xu23+s13*xu33) * wgt;
/*...............................................................N23/Q2 */
force[4][ngauss]  += (s21*xu13+s22*xu23+s23*xu33) * wgt;
/*..................................................................N31 */
force[16][ngauss] += (s31*xu11+s32*xu21+s33*xu31) * wgt;
/*..................................................................N32 */
force[17][ngauss] += (s31*xu12+s32*xu22+s33*xu32) * wgt;
/*..................................................................N33 */
force[9][ngauss]  += (s31*xu13+s32*xu23+s33*xu33) * wgt;
/*..................................................................M11 */
force[5][ngauss]  += (s11*xu11+s12*xu21+s13*xu31) * wgthe3;
/*..................................................................M12 */
force[7][ngauss]  += (s11*xu12+s12*xu22+s13*xu32) * wgthe3;
/*..................................................................M22 */
force[6][ngauss]  += (s21*xu12+s22*xu22+s23*xu32) * wgthe3;
/*..................................................................M13 */
force[10][ngauss] += (s11*xu13+s12*xu23+s13*xu33) * wgthe3;
/*..................................................................M23 */
force[11][ngauss] += (s21*xu13+s22*xu23+s23*xu33) * wgthe3;
/*..................................................................M31 */
force[12][ngauss] += (s31*xu11+s32*xu21+s33*xu31) * wgthe3;
/*..................................................................M32 */
force[13][ngauss] += (s31*xu12+s32*xu22+s33*xu32) * wgthe3;
/*..................................................................M21 */
force[14][ngauss] += (s21*xu11+s22*xu21+s23*xu31) * wgthe3;
/*..................................................................M33 */
force[15][ngauss] += (s31*xu13+s32*xu23+s33*xu33) * wgthe3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tfte */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
