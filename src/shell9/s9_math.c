/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_Msort_bs9:   which rearranges the constitutive matrix from "brick" to "shell9"
 - s9_Msort_s9b:   which rearranges the constitutive matrix from "shell9" to "brick"

 - s9_Vsort_bs9:   sorting of 'stress'/'strain'-vector from "brick" to "shell9"
 - s9_Vsort_s9b:   sorting of 'stress'/'strain'-vector from "shell9" to "brick"

 - s9shear_cor:    which makes some modifications to the constitutive matrix
                   due to some shear correction aspects

 - s9_rot:         which calculates the basis vectors of the material/laminat
                   coordinate system according to the cartesian coord. sys.
 - s9_Teps:        which calculates the transformation matrix T_sup_epsilon
 - s9_Ecama:       which transforms strains from cartesian to material/laminat coord. sys
 - s9_Smaca:       which transforms stresses from material/laminat to cartesian coord. sys
 - s9_Cmaca:       which transforms the constitutive matrix from material/laminat to
                   cartesian coord. sys


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief sorts the constitutive matrix according to shell9,... conventions
       sorts from "brick" to "shell9"

<pre>                                                            sh 05/03
This routine rearranges the constitutive matrix due to the order
convention of different elements and material laws, e.g.
brick:  [11,22,33,12,23,13]
shell8: [11,12,13,22,23,33]
shell9: [11,12,22,13,23,33] ...
</pre>
\param  DOUBLE **M      (i/o) matrix do be rearranged ([6][6])

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Msort_bs9(DOUBLE **M)
{
INT    i,j;
DOUBLE C_HELP[6][6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Msort_bs9");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) C_HELP[i][j] = M[i][j];

   M[0][0] =  C_HELP[0][0];
   M[0][1] =  C_HELP[0][3];
   M[0][2] =  C_HELP[0][1];
   M[0][3] =  C_HELP[0][5];
   M[0][4] =  C_HELP[0][4];
   M[0][5] =  C_HELP[0][2];

   M[1][0] =  C_HELP[3][0];
   M[1][1] =  C_HELP[3][3];
   M[1][2] =  C_HELP[3][1];
   M[1][3] =  C_HELP[3][5];
   M[1][4] =  C_HELP[3][4];
   M[1][5] =  C_HELP[3][2];

   M[2][0] =  C_HELP[1][0];
   M[2][1] =  C_HELP[1][3];
   M[2][2] =  C_HELP[1][1];
   M[2][3] =  C_HELP[1][5];
   M[2][4] =  C_HELP[1][4];
   M[2][5] =  C_HELP[1][2];

   M[3][0] =  C_HELP[5][0];
   M[3][1] =  C_HELP[5][3];
   M[3][2] =  C_HELP[5][1];
   M[3][3] =  C_HELP[5][5];
   M[3][4] =  C_HELP[5][4];
   M[3][5] =  C_HELP[5][2];

   M[4][0] =  C_HELP[4][0];
   M[4][1] =  C_HELP[4][3];
   M[4][2] =  C_HELP[4][1];
   M[4][3] =  C_HELP[4][5];
   M[4][4] =  C_HELP[4][4];
   M[4][5] =  C_HELP[4][2];

   M[5][0] =  C_HELP[2][0];
   M[5][1] =  C_HELP[2][3];
   M[5][2] =  C_HELP[2][1];
   M[5][3] =  C_HELP[2][5];
   M[5][4] =  C_HELP[2][4];
   M[5][5] =  C_HELP[2][2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Msort_bs9 */


/*!----------------------------------------------------------------------
\brief sorts the constitutive matrix according to shell9,... conventions
       sorts from "shell9" to "brick"

<pre>                                                            sh 05/03
This routine rearranges the constitutive matrix due to the order
convention of different elements and material laws, e.g.
brick:  [11,22,33,12,23,13]
shell8: [11,12,13,22,23,33]
shell9: [11,12,22,13,23,33] ...
</pre>
\param  DOUBLE **M      (i/o) matrix do be rearranged ([6][6])

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Msort_s9b(DOUBLE **M)
{
INT    i,j;
DOUBLE C_HELP[6][6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Msort_s9b");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<6; i++) for (j=0; j<6; j++) C_HELP[i][j] = M[i][j];

   M[0][0] =  C_HELP[0][0];
   M[0][3] =  C_HELP[0][1];
   M[0][1] =  C_HELP[0][2];
   M[0][5] =  C_HELP[0][3];
   M[0][4] =  C_HELP[0][4];
   M[0][2] =  C_HELP[0][5];

   M[3][0] =  C_HELP[1][0];
   M[3][3] =  C_HELP[1][1];
   M[3][1] =  C_HELP[1][2];
   M[3][5] =  C_HELP[1][3];
   M[3][4] =  C_HELP[1][4];
   M[3][2] =  C_HELP[1][5];

   M[1][0] =  C_HELP[2][0];
   M[1][3] =  C_HELP[2][1];
   M[1][1] =  C_HELP[2][2];
   M[1][5] =  C_HELP[2][3];
   M[1][4] =  C_HELP[2][4];
   M[1][2] =  C_HELP[2][5];

   M[5][0] =  C_HELP[3][0];
   M[5][3] =  C_HELP[3][1];
   M[5][1] =  C_HELP[3][2];
   M[5][5] =  C_HELP[3][3];
   M[5][4] =  C_HELP[3][4];
   M[5][2] =  C_HELP[3][5];

   M[4][0] =  C_HELP[4][0];
   M[4][3] =  C_HELP[4][1];
   M[4][1] =  C_HELP[4][2];
   M[4][5] =  C_HELP[4][3];
   M[4][4] =  C_HELP[4][4];
   M[4][2] =  C_HELP[4][5];

   M[2][0] =  C_HELP[5][0];
   M[2][3] =  C_HELP[5][1];
   M[2][1] =  C_HELP[5][2];
   M[2][5] =  C_HELP[5][3];
   M[2][4] =  C_HELP[5][4];
   M[2][2] =  C_HELP[5][5];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Msort_s9b */


/*!----------------------------------------------------------------------
\brief sorts the vectors 'strain'/'stress' according to shell9,... conventions
       sorts form "brick" to "shell9"

<pre>                                                            sh 05/03
This routine rearranges the vectors 'strain'/'stress' due to the order
convention of different elements and material laws, e.g.
brick:  [11,22,33,12,23,13]
shell8: [11,12,13,22,23,33]
shell9: [11,12,22,13,23,33] ...
</pre>
\param  DOUBLE *vec     (i/o) matrix do be rearranged ([6])

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Vsort_bs9(DOUBLE vec[6])
{
INT    i;
DOUBLE vec_help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Vsort_bs9");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<6; i++) vec_help[i] = vec[i];

   vec[0] =  vec_help[0];
   vec[1] =  vec_help[3];
   vec[2] =  vec_help[1];
   vec[3] =  vec_help[5];
   vec[4] =  vec_help[4];
   vec[5] =  vec_help[2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Vsort_bs9 */


/*!----------------------------------------------------------------------
\brief sorts the vectors 'strain'/'stress' according to shell9,... conventions
       sorts from "shell9" to "brick"

<pre>                                                            sh 05/03
This routine rearranges the vectors 'strain'/'stress' due to the order
convention of different elements and material laws, e.g.
brick:  [11,22,33,12,23,13]
shell8: [11,12,13,22,23,33]
shell9: [11,12,22,13,23,33] ...
</pre>
\param  DOUBLE *vec     (i/o) matrix do be rearranged ([6])

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Vsort_s9b(DOUBLE vec[6])
{
INT    i;
DOUBLE vec_help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Vsort_s9b");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<6; i++) vec_help[i] = vec[i];

   vec[0] =  vec_help[0];
   vec[3] =  vec_help[1];
   vec[1] =  vec_help[2];
   vec[5] =  vec_help[3];
   vec[4] =  vec_help[4];
   vec[2] =  vec_help[5];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Vsort_s9b */

/*!----------------------------------------------------------------------
\brief make some modification to material matrix due to shear correction
       factors

<pre>                                                            sh 02/03
This routine modifies the constitutive matrix with the shear correction
factor 'xsi' for [13] and [23]
</pre>
\param  DOUBLE  **C     (i/o) consitutive matrix to be modified
\param  DOUBLE    xsi    (i)  shear correction factor

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9shear_cor(DOUBLE **C, DOUBLE xsi)
{
#ifdef DEBUG
dstrc_enter("s9shear_cor");
#endif
/*----------------------------------------------------------------------*/
C[3][3] = C[3][3]/xsi;
C[3][4] = C[3][4]/xsi;

C[4][3] = C[4][3]/xsi;
C[4][4] = C[4][4]/xsi;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9shear_cor */


/*!----------------------------------------------------------------------
\brief calculates the basis vectors of a orthonormal material coordinate
system laying tangential to the shell surface

<pre>                                                            sh 09/03
theorie: Dis. Braun Kap. 5.4 and Anhang B
</pre>
\param  DOUBLE    phi       (i)  rotation angle
\param  INT       rot_axis  (i)  1=global x, 2=global y, 3=global z
\param  DOUBLE    gkov      (i)  kovariant basis vectors
\param  DOUBLE    g1,2,3[3] (o)  basis vectors of orthonormal material
                                 coordinate system tangential to shell surface

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_rot(DOUBLE phi,   INT    rot_axis, DOUBLE **gkov,
            DOUBLE g1[3], DOUBLE g2[3],    DOUBLE g3[3])
{
DOUBLE  rad;
DOUBLE  i1[3], j1[3];          /*i1 = i1', j1 = i1''*/
DOUBLE  i2[3], j2[3];          /*i2 = i2', j2 = i2''*/
DOUBLE  i3[3], j3[3];          /*i3 = i3', j3 = i3''*/
DOUBLE  norm;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_rot");
#endif
/*----------------------------------------------------------------------*/
rad = RAD * phi;

g3[0] = gkov[0][2];
g3[1] = gkov[1][2];
g3[2] = gkov[2][2];

if (rot_axis==1)
{
   /*calculate i2' */
   i2[0] = 0.0;
   i2[1] = cos(rad);
   i2[2] = sin(rad);

   /*calculate j2=i2''=i2' x g3*/
   j2[0] = i2[1]*g3[2] - i2[2]*g3[1];
   j2[1] = i2[2]*g3[0] - i2[0]*g3[2];
   j2[2] = i2[0]*g3[1] - i2[1]*g3[0];

   if((j2[0]*j2[0]+j2[1]*j2[1]+j2[2]*j2[2])<EPS12)  /*Singularity -> B.7*/
   {
      g2[0] =  0.0;
      g2[1] = -g3[2];
      g2[2] =  g3[1];

      g1[0] = g2[1]*g3[2] - g2[2]*g3[1];
      g1[1] = g2[2]*g3[0] - g2[0]*g3[2];
      g1[2] = g2[0]*g3[1] - g2[1]*g3[0];
   }
   else  /*Normal -> B.6*/
   {
      g1[0] = g3[1]*j2[2] - g3[2]*j2[1];
      g1[1] = g3[2]*j2[0] - g3[0]*j2[2];
      g1[2] = g3[0]*j2[1] - g3[1]*j2[0];

      g2[0] = g3[1]*g1[2] - g3[2]*g1[1];
      g2[1] = g3[2]*g1[0] - g3[0]*g1[2];
      g2[2] = g3[0]*g1[1] - g3[1]*g1[0];
   }
}
else if (rot_axis==2)
{
   /*calculate i3' */
   i3[0] = sin(rad);
   i3[1] = 0.0;
   i3[2] = cos(rad);

   /*calculate j3=i3''=i3' x g3*/
   j3[0] = i3[1]*g3[2] - i3[2]*g3[1];
   j3[1] = i3[2]*g3[0] - i3[0]*g3[2];
   j3[2] = i3[0]*g3[1] - i3[1]*g3[0];

   if((j3[0]*j3[0]+j3[1]*j3[1]+j3[2]*j3[2])<EPS12)  /*Singularity -> B.5*/
   {
      g2[0] =  g3[2];
      g2[1] =  0.0;
      g2[2] = -g3[0];

      g1[0] = g2[1]*g3[2] - g2[2]*g3[1];
      g1[1] = g2[2]*g3[0] - g2[0]*g3[2];
      g1[2] = g2[0]*g3[1] - g2[1]*g3[0];
   }
   else  /*Normal -> B.4*/
   {
      g1[0] = g3[1]*j3[2] - g3[2]*j3[1];
      g1[1] = g3[2]*j3[0] - g3[0]*j3[2];
      g1[2] = g3[0]*j3[1] - g3[1]*j3[0];

      g2[0] = g3[1]*g1[2] - g3[2]*g1[1];
      g2[1] = g3[2]*g1[0] - g3[0]*g1[2];
      g2[2] = g3[0]*g1[1] - g3[1]*g1[0];
   }
}
else if (rot_axis==3)
{
   /*calculate i1' */
   i1[0] = cos(rad);
   i1[1] = sin(rad);
   i1[2] = 0.0;

   /*calculate j1=i1''=i1' x g3*/
   j1[0] = i1[1]*g3[2] - i1[2]*g3[1];
   j1[1] = i1[2]*g3[0] - i1[0]*g3[2];
   j1[2] = i1[0]*g3[1] - i1[1]*g3[0];

   if((j1[0]*j1[0]+j1[1]*j1[1]+j1[2]*j1[2])<EPS12)  /*Singularity -> B.2*/
   {
      g2[0] = -g3[1];
      g2[1] =  g3[0];
      g2[2] =  0.0;

      g1[0] = g2[1]*g3[2] - g2[2]*g3[1];
      g1[1] = g2[2]*g3[0] - g2[0]*g3[2];
      g1[2] = g2[0]*g3[1] - g2[1]*g3[0];
   }
   else  /*Normal -> B.1*/
   {
      g1[0] = g3[1]*j1[2] - g3[2]*j1[1];
      g1[1] = g3[2]*j1[0] - g3[0]*j1[2];
      g1[2] = g3[0]*j1[1] - g3[1]*j1[0];

      g2[0] = g3[1]*g1[2] - g3[2]*g1[1];
      g2[1] = g3[2]*g1[0] - g3[0]*g1[2];
      g2[2] = g3[0]*g1[1] - g3[1]*g1[0];
   }
}
/*norm the basis vectors*/
norm = sqrt(g1[0]*g1[0]+g1[1]*g1[1]+g1[2]*g1[2]);
norm = 1./norm;
g1[0] *= norm;
g1[1] *= norm;
g1[2] *= norm;

norm = sqrt(g2[0]*g2[0]+g2[1]*g2[1]+g2[2]*g2[2]);
norm = 1./norm;
g2[0] *= norm;
g2[1] *= norm;
g2[2] *= norm;

norm = sqrt(g3[0]*g3[0]+g3[1]*g3[1]+g3[2]*g3[2]);
norm = 1./norm;
g3[0] *= norm;
g3[1] *= norm;
g3[2] *= norm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_rot */


/*!----------------------------------------------------------------------
\brief calculates transformation Matrix T_sup_epsilon
see Altenbach, Altenbach & Rikards p.29
"Einfuehrung in die Mechanik der Laminat- und Sandwichtragwerke"

<pre>                                                            sh 09/03
T_sup_epsilon is needed to transform strains, stresses and C from
cartesian coordinate system to material/laminat coordinate system

eps' =  T_eps * eps
C    = (T_eps)T * C' * T_eps
sig  = (T_eps)T * sig'
</pre>
\param  DOUBLE    g1,2,3[3] (i)  basis vectors of orthonormal material
                                 coordinate system tangential to shell surface
\param  DOUBLE    T[6][6]   (o)  transformation Matrix T_sup_epsilon

\warning      be careful about sorting [11,22,33,12,23,13]
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Teps(DOUBLE g1[3], DOUBLE g2[3], DOUBLE g3[3], DOUBLE T[6][6])
{
DOUBLE  R11,R12,R13,R21,R22,R23,R31,R32,R33;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Teps");
#endif
/*----------------------------------------------------------------------*/
/*-- basis vectors of material coord. sys. according to cartesian coord. sys. --*/
 R11 = g1[0];
 R12 = g1[1];
 R13 = g1[2];
 R21 = g2[0];
 R22 = g2[1];
 R23 = g2[2];
 R31 = g3[0];
 R32 = g3[1];
 R33 = g3[2];

/*get transformation matrix according to: Altenbach,Altenbach,Rikards:
  Einfuehrung in die Mechanik der Laminat- und Sandwichtragwerke" p.29*/
/*Sorting: [11,22,33,12,23,13]*/
 /*1.Zeile*/
 T[0][0] = R11*R11;
 T[0][1] = R12*R12;
 T[0][2] = R13*R13;
 T[0][3] = R11*R12;
 T[0][4] = R12*R13;
 T[0][5] = R11*R13;
 /*2.Zeile*/
 T[1][0] = R21*R21;
 T[1][1] = R22*R22;
 T[1][2] = R23*R23;
 T[1][3] = R21*R22;
 T[1][4] = R22*R23;
 T[1][5] = R21*R23;
 /*3.Zeile*/
 T[2][0] = R31*R31;
 T[2][1] = R32*R32;
 T[2][2] = R33*R33;
 T[2][3] = R31*R32;
 T[2][4] = R32*R33;
 T[2][5] = R31*R33;
 /*4.Zeile*/
 T[3][0] = 2.*R11*R21;
 T[3][1] = 2.*R12*R22;
 T[3][2] = 2.*R13*R23;
 T[3][3] = R11*R22 + R12*R21;
 T[3][4] = R12*R23 + R13*R22;
 T[3][5] = R11*R23 + R13*R21;
 /*5.Zeile*/
 T[4][0] = 2.*R21*R31;
 T[4][1] = 2.*R22*R32;
 T[4][2] = 2.*R23*R33;
 T[4][3] = R21*R32 + R22*R31;
 T[4][4] = R22*R33 + R23*R32;
 T[4][5] = R21*R33 + R23*R31;
 /*6.Zeile*/
 T[5][0] = 2.*R11*R31;
 T[5][1] = 2.*R12*R32;
 T[5][2] = 2.*R13*R33;
 T[5][3] = R11*R32 + R12*R31;
 T[5][4] = R12*R33 + R13*R32;
 T[5][5] = R11*R33 + R13*R31;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Teps */


/*!----------------------------------------------------------------------
\brief transforms strains from cartesian coordinate system to material/
laminat coordinate system
see Altenbach, Altenbach & Rikards p.29
"Einfuehrung in die Mechanik der Laminat- und Sandwichtragwerke"

<pre>                                                            sh 09/03
eps' =  T_eps * eps
</pre>
\param  DOUBLE    E[6]     (i/0)  strains to be transformed
\param  DOUBLE    T[6][6]   (i)  transformation Matrix T_sup_epsilon

\warning      be careful about sorting [11,22,33,12,23,13]
              E={ e11, e22, e33, 2*e12, 2*e23, 2*e13 }
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Ecama(DOUBLE E[6], DOUBLE T[6][6])
{
INT     i,j;
DOUBLE  E_help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Ecama");
#endif
/*----------------------------------------------------------------------*/
E_help[0] = E[0];
E_help[1] = E[1];
E_help[2] = E[2];
E_help[3] = E[3];
E_help[4] = E[4];
E_help[5] = E[5];

E[0] = 0.0;
E[1] = 0.0;
E[2] = 0.0;
E[3] = 0.0;
E[4] = 0.0;
E[5] = 0.0;

/*  eps' = T_eps * eps   */
for (i=0; i<6; i++) for (j=0; j<6; j++) E[i] += T[i][j] * E_help[j];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Ecama */



/*!----------------------------------------------------------------------
\brief transforms stresses from material/laminat coordinate system to
cartesian coordinate system
see Altenbach, Altenbach & Rikards p.29
"Einfuehrung in die Mechanik der Laminat- und Sandwichtragwerke"

<pre>                                                            sh 09/03
sig = (T_eps)T * sig'
</pre>
\param  DOUBLE    S[6]     (i/0)  stresses to be transformed
\param  DOUBLE    T[6][6]   (i)  transformation Matrix T_sup_epsilon

\warning      be careful about sorting [11,22,33,12,23,13]
              E={ e11, e22, e33, 2*e12, 2*e23, 2*e13 }
              S={ s11, s22, s33,   s12,   s23,   s13 }
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Smaca(DOUBLE S[6], DOUBLE T[6][6])
{
INT     i,j;
DOUBLE  S_help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Smaca");
#endif
/*----------------------------------------------------------------------*/
S_help[0] = S[0];
S_help[1] = S[1];
S_help[2] = S[2];
S_help[3] = S[3];
S_help[4] = S[4];
S_help[5] = S[5];

S[0] = 0.0;
S[1] = 0.0;
S[2] = 0.0;
S[3] = 0.0;
S[4] = 0.0;
S[5] = 0.0;

/*  sig = (T_eps)T * sig'   */
for (i=0; i<6; i++) for (j=0; j<6; j++) S[i] += T[j][i] * S_help[j];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Smaca */



/*!----------------------------------------------------------------------
\brief transforms the constitutive matrix from orthonormal material coord.
system to global cartesian coord. system
see Altenbach, Altenbach & Rikards p.29
"Einfuehrung in die Mechanik der Laminat- und Sandwichtragwerke"

<pre>                                                            sh 09/03
C = (T_eps)T * C' * T_eps
</pre>
\param  DOUBLE  **C        (i/0) constitutive Matrix to be transformed
\param  DOUBLE    T[6][6]   (i)  transformation Matrix T_sup_epsilon

\warning      be careful about sorting [11,22,33,12,23,13]
              E={ e11, e22, e33, 2*e12, 2*e23, 2*e13 }
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Cmaca(DOUBLE **C, DOUBLE T[6][6])
{
INT     i,j,k;
DOUBLE  C_help[6][6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Cmaca");
#endif
/*----------------------------------------------------------------------*/
/* initialize C_help*/
for (i=0; i<6; i++) for (j=0; j<6; j++) C_help[i][j] = 0.0;

/*  C_help = (T_eps)T * C'   */
for (i=0; i<6; i++)
    for (j=0; j<6; j++)
        for (k=0; k<6; k++) C_help[i][j] += T[k][i] * C[k][j];


/* initialize C*/
for (i=0; i<6; i++) for (j=0; j<6; j++) C[i][j] = 0.0;

/* C = (T_eps)T * C' * T_eps -> C = C_help * T_eps   */
for (i=0; i<6; i++)
    for (j=0; j<6; j++)
        for (k=0; k<6; k++) C[i][j] += C_help[i][k] * T[k][j];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Cmaca */




/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

#endif
