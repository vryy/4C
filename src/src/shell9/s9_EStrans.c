/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_Ecacu: transforms strains (E) from cartesian to curvilinear
 - s9_Ecuca: transforms strains (E) from curvilinear to cartesian
 - s9_Scacu: transforms stresses (S) from cartesian to curvilinear
 - s9_Scuca: transforms stresses (S) from curvilinear to cartesian

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
\brief transforms strains (E) from cartesian to curvilinear coordinate system

<pre>                                                            sh 05/03
This routine transforms the strains (E) from cartesian to curvilinear
coordinate system. This is very useful in order to use general material laws,
which are written in orthonormal coordinate systems.
</pre>
\param  DOUBLE  *str    (i/o) stresses or strains to be modified
\param  DOUBLE **gkov    (i)  kovariant basis vectors

\warning Always call this routine with the kovariant basis vectors (gkov).
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Ecacu (DOUBLE E[6], DOUBLE **gkov)
{
DOUBLE   help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Ecacu");
#endif
/*----------------------------------------------------------------------*/
help[0] = 0.0;
help[1] = 0.0;
help[2] = 0.0;
help[3] = 0.0;
help[4] = 0.0;
help[5] = 0.0;

help[0] =      gkov[0][0] * gkov[0][0] * E[0] +
          2. * gkov[1][0] * gkov[0][0] * E[1] +
          2. * gkov[2][0] * gkov[0][0] * E[3] +
               gkov[1][0] * gkov[1][0] * E[2] +
          2. * gkov[2][0] * gkov[1][0] * E[4] +
               gkov[2][0] * gkov[2][0] * E[5] ;

help[2] =      gkov[0][1] * gkov[0][1] * E[0] +
          2. * gkov[1][1] * gkov[0][1] * E[1] +
          2. * gkov[2][1] * gkov[0][1] * E[3] +
               gkov[1][1] * gkov[1][1] * E[2] +
          2. * gkov[2][1] * gkov[1][1] * E[4] +
               gkov[2][1] * gkov[2][1] * E[5] ;

help[5] =      gkov[0][2] * gkov[0][2] * E[0] +
          2. * gkov[1][2] * gkov[0][2] * E[1] +
          2. * gkov[2][2] * gkov[0][2] * E[3] +
               gkov[1][2] * gkov[1][2] * E[2] +
          2. * gkov[2][2] * gkov[1][2] * E[4] +
               gkov[2][2] * gkov[2][2] * E[5] ;

help[1] = gkov[0][0] * gkov[0][1] * E[0] +
          gkov[1][0] * gkov[0][1] * E[1] +
          gkov[2][0] * gkov[0][1] * E[3] +
          gkov[0][0] * gkov[1][1] * E[1] +
          gkov[1][0] * gkov[1][1] * E[2] +
          gkov[2][0] * gkov[1][1] * E[4] +
          gkov[0][0] * gkov[2][1] * E[3] +
          gkov[1][0] * gkov[2][1] * E[4] +
          gkov[2][0] * gkov[2][1] * E[5] ;

help[3] = gkov[0][0] * gkov[0][2] * E[0] +
          gkov[1][0] * gkov[0][2] * E[1] +
          gkov[2][0] * gkov[0][2] * E[3] +
          gkov[0][0] * gkov[1][2] * E[1] +
          gkov[1][0] * gkov[1][2] * E[2] +
          gkov[2][0] * gkov[1][2] * E[4] +
          gkov[0][0] * gkov[2][2] * E[3] +
          gkov[1][0] * gkov[2][2] * E[4] +
          gkov[2][0] * gkov[2][2] * E[5] ;

help[4] = gkov[0][1] * gkov[0][2] * E[0] +
          gkov[1][1] * gkov[0][2] * E[1] +
          gkov[2][1] * gkov[0][2] * E[3] +
          gkov[0][1] * gkov[1][2] * E[1] +
          gkov[1][1] * gkov[1][2] * E[2] +
          gkov[2][1] * gkov[1][2] * E[4] +
          gkov[0][1] * gkov[2][2] * E[3] +
          gkov[1][1] * gkov[2][2] * E[4] +
          gkov[2][1] * gkov[2][2] * E[5] ;

E[0] = help[0];
E[1] = help[1];
E[2] = help[2];
E[3] = help[3];
E[4] = help[4];
E[5] = help[5];

/*----------------------------------------------------------------
This is better to understand, but not as fast

A[0][0] = E[0];
A[1][1] = E[2];
A[2][2] = E[5];
A[0][1] = E[1];
A[1][0] = E[1];
A[0][2] = E[3];
A[2][0] = E[3];
A[1][2] = E[4];
A[2][1] = E[4];

B[0][0] =  0.0;
B[1][1] =  0.0;
B[2][2] =  0.0;
B[0][1] =  0.0;
B[0][2] =  0.0;
B[1][2] =  0.0;

for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   B[0][0] += gkov[k][0] * gkov[l][0] * A[k][l];
   B[1][1] += gkov[k][1] * gkov[l][1] * A[k][l];
   B[2][2] += gkov[k][2] * gkov[l][2] * A[k][l];
   B[0][1] += gkov[k][0] * gkov[l][1] * A[k][l];
   B[0][2] += gkov[k][0] * gkov[l][2] * A[k][l];
   B[1][2] += gkov[k][1] * gkov[l][2] * A[k][l];
}

E[0] = B[0][0];
E[2] = B[1][1];
E[5] = B[2][2];
E[1] = B[0][1];
E[3] = B[0][2];
E[4] = B[1][2];
------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Ecacu */



/*!----------------------------------------------------------------------
\brief transforms strains (E) from curvilinear to cartesian coordinate system

<pre>                                                            sh 05/03
This routine transforms the strains (E) from curvilinear to cartesian
coordinate system. This is very useful in order to use general material laws,
which are written in orthonormal coordinate systems.
</pre>
\param  DOUBLE  *str    (i/o) stresses or strains to be modified
\param  DOUBLE **gkon    (i)  kontravariant basis vectors

\warning Always call this routine with the kontravariant basis vectors (gkon).
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Ecuca (DOUBLE E[6], DOUBLE **gkon)
{
DOUBLE   help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Ecuca");
#endif
/*----------------------------------------------------------------------*/
help[0] = 0.0;
help[1] = 0.0;
help[2] = 0.0;
help[3] = 0.0;
help[4] = 0.0;
help[5] = 0.0;

help[0] =      gkon[0][0] * gkon[0][0] * E[0] +
          2. * gkon[0][1] * gkon[0][0] * E[1] +
          2. * gkon[0][2] * gkon[0][0] * E[3] +
               gkon[0][1] * gkon[0][1] * E[2] +
          2. * gkon[0][2] * gkon[0][1] * E[4] +
               gkon[0][2] * gkon[0][2] * E[5] ;

help[2] =      gkon[1][0] * gkon[1][0] * E[0] +
          2. * gkon[1][1] * gkon[1][0] * E[1] +
          2. * gkon[1][2] * gkon[1][0] * E[3] +
               gkon[1][1] * gkon[1][1] * E[2] +
          2. * gkon[1][2] * gkon[1][1] * E[4] +
               gkon[1][2] * gkon[1][2] * E[5] ;

help[5] =      gkon[2][0] * gkon[2][0] * E[0] +
          2. * gkon[2][1] * gkon[2][0] * E[1] +
          2. * gkon[2][2] * gkon[2][0] * E[3] +
               gkon[2][1] * gkon[2][1] * E[2] +
          2. * gkon[2][2] * gkon[2][1] * E[4] +
               gkon[2][2] * gkon[2][2] * E[5] ;

help[1] = gkon[0][0] * gkon[1][0] * E[0] +
          gkon[0][1] * gkon[1][0] * E[1] +
          gkon[0][2] * gkon[1][0] * E[3] +
          gkon[0][0] * gkon[1][1] * E[1] +
          gkon[0][1] * gkon[1][1] * E[2] +
          gkon[0][2] * gkon[1][1] * E[4] +
          gkon[0][0] * gkon[1][2] * E[3] +
          gkon[0][1] * gkon[1][2] * E[4] +
          gkon[0][2] * gkon[1][2] * E[5] ;

help[3] = gkon[0][0] * gkon[2][0] * E[0] +
          gkon[0][1] * gkon[2][0] * E[1] +
          gkon[0][2] * gkon[2][0] * E[3] +
          gkon[0][0] * gkon[2][1] * E[1] +
          gkon[0][1] * gkon[2][1] * E[2] +
          gkon[0][2] * gkon[2][1] * E[4] +
          gkon[0][0] * gkon[2][2] * E[3] +
          gkon[0][1] * gkon[2][2] * E[4] +
          gkon[0][2] * gkon[2][2] * E[5] ;

help[4] = gkon[1][0] * gkon[2][0] * E[0] +
          gkon[1][1] * gkon[2][0] * E[1] +
          gkon[1][2] * gkon[2][0] * E[3] +
          gkon[1][0] * gkon[2][1] * E[1] +
          gkon[1][1] * gkon[2][1] * E[2] +
          gkon[1][2] * gkon[2][1] * E[4] +
          gkon[1][0] * gkon[2][2] * E[3] +
          gkon[1][1] * gkon[2][2] * E[4] +
          gkon[1][2] * gkon[2][2] * E[5] ;

E[0] = help[0];
E[1] = help[1];
E[2] = help[2];
E[3] = help[3];
E[4] = help[4];
E[5] = help[5];

/*----------------------------------------------------------------
This is better to understand, but not as fast

A[0][0] = E[0];
A[1][1] = E[2];
A[2][2] = E[5];
A[0][1] = E[1];
A[1][0] = E[1];
A[0][2] = E[3];
A[2][0] = E[3];
A[1][2] = E[4];
A[2][1] = E[4];

B[0][0] =  0.0;
B[1][1] =  0.0;
B[2][2] =  0.0;
B[0][1] =  0.0;
B[0][2] =  0.0;
B[1][2] =  0.0;

for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   B[0][0] += gkon[0][k] * gkon[0][l] * A[k][l];
   B[1][1] += gkon[1][k] * gkon[1][l] * A[k][l];
   B[2][2] += gkon[2][k] * gkon[2][l] * A[k][l];
   B[0][1] += gkon[0][k] * gkon[1][l] * A[k][l];
   B[0][2] += gkon[0][k] * gkon[2][l] * A[k][l];
   B[1][2] += gkon[1][k] * gkon[2][l] * A[k][l];
}

E[0] = B[0][0];
E[2] = B[1][1];
E[5] = B[2][2];
E[1] = B[0][1];
E[3] = B[0][2];
E[4] = B[1][2];

------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Ecuca */

/*!----------------------------------------------------------------------
\brief transforms stresses (S) from cartesian to curvilinear coordinate system

<pre>                                                            sh 05/03
This routine transforms the stresses (S) from cartesian to curvilinear
coordinate system. This is very useful in order to use general material laws,
which are written in orthonormal coordinate systems.
</pre>
\param  DOUBLE  *str    (i/o) stresses or strains to be modified
\param  DOUBLE **gkon    (i)  kontravariant basis vectors

\warning Always call this routine with the kontravariant basis vectors (gkon).
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Scacu (DOUBLE S[6], DOUBLE **gkon)
{
DOUBLE   help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Scacu");
#endif
/*----------------------------------------------------------------------*/
help[0] = 0.0;
help[1] = 0.0;
help[2] = 0.0;
help[3] = 0.0;
help[4] = 0.0;
help[5] = 0.0;

help[0] =      gkon[0][0] * gkon[0][0] * S[0] +
          2. * gkon[1][0] * gkon[0][0] * S[1] +
          2. * gkon[2][0] * gkon[0][0] * S[3] +
               gkon[1][0] * gkon[1][0] * S[2] +
          2. * gkon[2][0] * gkon[1][0] * S[4] +
               gkon[2][0] * gkon[2][0] * S[5] ;

help[2] =      gkon[0][1] * gkon[0][1] * S[0] +
          2. * gkon[1][1] * gkon[0][1] * S[1] +
          2. * gkon[2][1] * gkon[0][1] * S[3] +
               gkon[1][1] * gkon[1][1] * S[2] +
          2. * gkon[2][1] * gkon[1][1] * S[4] +
               gkon[2][1] * gkon[2][1] * S[5] ;

help[5] =      gkon[0][2] * gkon[0][2] * S[0] +
          2. * gkon[1][2] * gkon[0][2] * S[1] +
          2. * gkon[2][2] * gkon[0][2] * S[3] +
               gkon[1][2] * gkon[1][2] * S[2] +
          2. * gkon[2][2] * gkon[1][2] * S[4] +
               gkon[2][2] * gkon[2][2] * S[5] ;

help[1] = gkon[0][0] * gkon[0][1] * S[0] +
          gkon[1][0] * gkon[0][1] * S[1] +
          gkon[2][0] * gkon[0][1] * S[3] +
          gkon[0][0] * gkon[1][1] * S[1] +
          gkon[1][0] * gkon[1][1] * S[2] +
          gkon[2][0] * gkon[1][1] * S[4] +
          gkon[0][0] * gkon[2][1] * S[3] +
          gkon[1][0] * gkon[2][1] * S[4] +
          gkon[2][0] * gkon[2][1] * S[5] ;

help[3] = gkon[0][0] * gkon[0][2] * S[0] +
          gkon[1][0] * gkon[0][2] * S[1] +
          gkon[2][0] * gkon[0][2] * S[3] +
          gkon[0][0] * gkon[1][2] * S[1] +
          gkon[1][0] * gkon[1][2] * S[2] +
          gkon[2][0] * gkon[1][2] * S[4] +
          gkon[0][0] * gkon[2][2] * S[3] +
          gkon[1][0] * gkon[2][2] * S[4] +
          gkon[2][0] * gkon[2][2] * S[5] ;

help[4] = gkon[0][1] * gkon[0][2] * S[0] +
          gkon[1][1] * gkon[0][2] * S[1] +
          gkon[2][1] * gkon[0][2] * S[3] +
          gkon[0][1] * gkon[1][2] * S[1] +
          gkon[1][1] * gkon[1][2] * S[2] +
          gkon[2][1] * gkon[1][2] * S[4] +
          gkon[0][1] * gkon[2][2] * S[3] +
          gkon[1][1] * gkon[2][2] * S[4] +
          gkon[2][1] * gkon[2][2] * S[5] ;

S[0] = help[0];
S[1] = help[1];
S[2] = help[2];
S[3] = help[3];
S[4] = help[4];
S[5] = help[5];


/*----------------------------------------------------------------
This is better to understand, but not as fast

A[0][0] = S[0];
A[1][1] = S[2];
A[2][2] = S[5];
A[0][1] = S[1];
A[1][0] = S[1];
A[0][2] = S[3];
A[2][0] = S[3];
A[1][2] = S[4];
A[2][1] = S[4];

B[0][0] =  0.0;
B[1][1] =  0.0;
B[2][2] =  0.0;
B[0][1] =  0.0;
B[0][2] =  0.0;
B[1][2] =  0.0;

for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   B[0][0] += gkon[k][0] * gkon[l][0] * A[k][l];
   B[1][1] += gkon[k][1] * gkon[l][1] * A[k][l];
   B[2][2] += gkon[k][2] * gkon[l][2] * A[k][l];
   B[0][1] += gkon[k][0] * gkon[l][1] * A[k][l];
   B[0][2] += gkon[k][0] * gkon[l][2] * A[k][l];
   B[1][2] += gkon[k][1] * gkon[l][2] * A[k][l];
}

S[0] = B[0][0];
S[2] = B[1][1];
S[5] = B[2][2];
S[1] = B[0][1];
S[3] = B[0][2];
S[4] = B[1][2];

------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Scacu */


/*!----------------------------------------------------------------------
\brief transforms stresses (S) from curvilinear to cartesian coordinate system

<pre>                                                            sh 05/03
This routine transforms the stresses (S) from curvilinear to cartesian
coordinate system. This is very useful in order to use general material laws,
which are written in orthonormal coordinate systems.
</pre>
\param  DOUBLE  *str    (i/o) stresses or strains to be modified
\param  DOUBLE **gkov    (i)  kovariant basis vectors

\warning Always call this routine with the kovariant basis vectors (gkov).
\return void
\sa calling: ---; called by:

*-----------------------------------------------------------------------*/
void s9_Scuca (DOUBLE S[6], DOUBLE **gkov)
{
DOUBLE   help[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_Scuca");
#endif
/*----------------------------------------------------------------------*/
help[0] = 0.0;
help[1] = 0.0;
help[2] = 0.0;
help[3] = 0.0;
help[4] = 0.0;
help[5] = 0.0;

help[0] =      gkov[0][0] * gkov[0][0] * S[0] +
          2. * gkov[0][1] * gkov[0][0] * S[1] +
          2. * gkov[0][2] * gkov[0][0] * S[3] +
               gkov[0][1] * gkov[0][1] * S[2] +
          2. * gkov[0][2] * gkov[0][1] * S[4] +
               gkov[0][2] * gkov[0][2] * S[5] ;

help[2] =      gkov[1][0] * gkov[1][0] * S[0] +
          2. * gkov[1][1] * gkov[1][0] * S[1] +
          2. * gkov[1][2] * gkov[1][0] * S[3] +
               gkov[1][1] * gkov[1][1] * S[2] +
          2. * gkov[1][2] * gkov[1][1] * S[4] +
               gkov[1][2] * gkov[1][2] * S[5] ;

help[5] =      gkov[2][0] * gkov[2][0] * S[0] +
          2. * gkov[2][1] * gkov[2][0] * S[1] +
          2. * gkov[2][2] * gkov[2][0] * S[3] +
               gkov[2][1] * gkov[2][1] * S[2] +
          2. * gkov[2][2] * gkov[2][1] * S[4] +
               gkov[2][2] * gkov[2][2] * S[5] ;

help[1] = gkov[0][0] * gkov[1][0] * S[0] +
          gkov[0][1] * gkov[1][0] * S[1] +
          gkov[0][2] * gkov[1][0] * S[3] +
          gkov[0][0] * gkov[1][1] * S[1] +
          gkov[0][1] * gkov[1][1] * S[2] +
          gkov[0][2] * gkov[1][1] * S[4] +
          gkov[0][0] * gkov[1][2] * S[3] +
          gkov[0][1] * gkov[1][2] * S[4] +
          gkov[0][2] * gkov[1][2] * S[5] ;

help[3] = gkov[0][0] * gkov[2][0] * S[0] +
          gkov[0][1] * gkov[2][0] * S[1] +
          gkov[0][2] * gkov[2][0] * S[3] +
          gkov[0][0] * gkov[2][1] * S[1] +
          gkov[0][1] * gkov[2][1] * S[2] +
          gkov[0][2] * gkov[2][1] * S[4] +
          gkov[0][0] * gkov[2][2] * S[3] +
          gkov[0][1] * gkov[2][2] * S[4] +
          gkov[0][2] * gkov[2][2] * S[5] ;

help[4] = gkov[1][0] * gkov[2][0] * S[0] +
          gkov[1][1] * gkov[2][0] * S[1] +
          gkov[1][2] * gkov[2][0] * S[3] +
          gkov[1][0] * gkov[2][1] * S[1] +
          gkov[1][1] * gkov[2][1] * S[2] +
          gkov[1][2] * gkov[2][1] * S[4] +
          gkov[1][0] * gkov[2][2] * S[3] +
          gkov[1][1] * gkov[2][2] * S[4] +
          gkov[1][2] * gkov[2][2] * S[5] ;

S[0] = help[0];
S[1] = help[1];
S[2] = help[2];
S[3] = help[3];
S[4] = help[4];
S[5] = help[5];

/*----------------------------------------------------------------
This is better to understand, but not as fast

A[0][0] = S[0];
A[1][1] = S[2];
A[2][2] = S[5];
A[0][1] = S[1];
A[1][0] = S[1];
A[0][2] = S[3];
A[2][0] = S[3];
A[1][2] = S[4];
A[2][1] = S[4];

B[0][0] =  0.0;
B[1][1] =  0.0;
B[2][2] =  0.0;
B[0][1] =  0.0;
B[0][2] =  0.0;
B[1][2] =  0.0;

for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   B[0][0] += gkov[0][k] * gkov[0][l] * A[k][l];
   B[1][1] += gkov[1][k] * gkov[1][l] * A[k][l];
   B[2][2] += gkov[2][k] * gkov[2][l] * A[k][l];
   B[0][1] += gkov[0][k] * gkov[1][l] * A[k][l];
   B[0][2] += gkov[0][k] * gkov[2][l] * A[k][l];
   B[1][2] += gkov[1][k] * gkov[2][l] * A[k][l];
}

S[0] = B[0][0];
S[2] = B[1][1];
S[5] = B[2][2];
S[1] = B[0][1];
S[3] = B[0][2];
S[4] = B[1][2];

------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_Scuca */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
#endif
