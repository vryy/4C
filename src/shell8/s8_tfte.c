/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates forces from stresses                       m.gee 12/01   |
C---------------------------------------------------------------------
C!    BERECHNEN DER KRAEFTE UND MOMENTE AUS SPANNUNGEN
C!    STRESS    -  SPANNUNGEN (KONTRAVARIANT)
C!    GMKOV     -  KOVARIANTER METRIKTENSOR  (SCHALENRAUM)
C!    AMKOV     -  KOVARIANTER METRIKTENSOR  (SCHALENMITTELFLAECHE)
C!    GMKON     -  KONTRAVARIANTER METRIKTENSOR (SCHALENRAUM)
C!    H         -  DICKE AM GAUSSPUNKT
C!    E3        -  NORMIERTE DICKENKOORDINATE
C!    DETSM     -  'FLAECHENELEMENT' DER SMF
C!    DETSR     -      "             DES SCHALENRAUMS
C---------------------------------------------------------------------
 *----------------------------------------------------------------------*/
void s8_tfte(DOUBLE **force,
             INT      ngauss,
             DOUBLE  *stress,
             DOUBLE **gkov,
             DOUBLE **akon,
             DOUBLE **gmkov,
             DOUBLE **gmkon,
             DOUBLE **amkov,
             DOUBLE **amkon,
             DOUBLE   h,
             DOUBLE   e3,
             DOUBLE   fact,
             DOUBLE   detsm,
             DOUBLE   detsr)
{
INT    i;
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
#ifdef DEBUG
dstrc_enter("s8_tfte");
#endif
/*----------------------------------------------------------------------*/
s11 = stress[0];
s12 = stress[1];
s21 = s12;
s13 = stress[2];
s31 = s13;
s22 = stress[3];
s23 = stress[4];
s32 = s23;
s33 = stress[5];
/*----------------------------------------------------------------------*/
xu     = detsr/detsm;
hh     = sqrt(amkov[2][2]);
wgt    = hh * fact * xu;
wgthe3 = wgt * e3 * hh;
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++) xu11 += gkov[i][0]*akon[i][0];
for (i=0; i<3; i++) xu21 += gkov[i][0]*akon[i][1];
for (i=0; i<3; i++) xu31 += gkov[i][0]*akon[i][2];
for (i=0; i<3; i++) xu12 += gkov[i][1]*akon[i][0];
for (i=0; i<3; i++) xu22 += gkov[i][1]*akon[i][1];
for (i=0; i<3; i++) xu32 += gkov[i][1]*akon[i][2];
for (i=0; i<3; i++) xu13 += gkov[i][2]*akon[i][0];
for (i=0; i<3; i++) xu23 += gkov[i][2]*akon[i][1];
for (i=0; i<3; i++) xu33 += gkov[i][2]*akon[i][2];
/*----------------------------------------------------------------------*/
/*..................................................................N11 */
force[0][ngauss]  += (s11*xu11+s12*xu21+s13*xu31) * wgt;
/*..................................................................N12 */
force[2][ngauss]  += (s11*xu12+s12*xu22+s13*xu32) * wgt;
/*..................................................................N21 */
force[8][ngauss]  += (s21*xu11+s22*xu21+s23*xu31) * wgt;
/*..................................................................N22 */
force[1][ngauss]  += (s21*xu12+s22*xu22+s23*xu32) * wgt;
/*..................................................................Q1 */
force[3][ngauss]  += (s11*xu13+s12*xu23+s13*xu33) * wgt;
/*..................................................................Q2 */
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
} /* end of s8_tfte */
#endif
