/*!-----------------------------------------------------------------------
\file
\brief contains the routines 'if_ke' which calculates the usual
stiffness matrix for an interface element and 'if_fint' which calculates
the internal forces for an interface element
<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122

*-----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*!
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates usual stiffness matrix for small strains

<pre>                                                              ah 05/03
This routine calculates usual stiffness matrix for small strains

</pre>
\param   iel       INT        (I) number of element nodes
\param   flag      INT        (I) flag
\param **stiff     DOUBLE     (O) element stiffness matrix
\param **bop       DOUBLE     (I) B-operator
\param **Q         DOUBLE     (I) material tangent
\param   fac       DOUBLE     (I) integration factor

\warning There is nothing special to this routine
\return void
\sa calling:  ---;
    caled by: ifstatic_ke();

*----------------------------------------------------------------------*/
void if_ke(INT       iel,
           INT       flag,
           DOUBLE  **stiff,
           DOUBLE  **bop,
           DOUBLE  **Q,
           DOUBLE    fac)
{
INT            i,j,k,l,m;
DOUBLE         dum;
DOUBLE         CB[2];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("if_ke");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<2*iel; j++)
{
  for (k=0; k<2; k++)
  {
   CB[k] = 0.0 ;
    for (l=0; l<2; l++)
    {
    CB[k] = CB[k] + Q[k][l]*bop[l][j]*fac ;
    }
  }
  for (i=0; i<2*iel; i++)
  {
/*---------------------------------- dum = BT C B * fac of this GP ---*/
    dum = 0.0 ;
    for (m=0; m<2; m++)
    {
     dum = dum + bop[m][i]*CB[m] ;
    }
/*---------------------------------------- sum over gaussian points ---*/
     stiff[i][j] = stiff[i][j] + dum ;
  }
}
/*---- psoudo-stiffness for main-diagonal terms of eliminated nodes ---*/
# if 0
if (iel == 8)
{
 if(flag == 1)
 {
 stiff[10][10]=1000;
 stiff[11][11]=1000;
 stiff[14][14]=1000;
 stiff[15][15]=1000;
 }
 else if(flag ==2)
 {
 stiff[8][8]=1000;
 stiff[9][9]=1000;
 stiff[12][12]=1000;
 stiff[13][13]=1000;
 }
}
# endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of if_ke */

/*!----------------------------------------------------------------------
\brief calculates internal forces

<pre>                                                             ah 05/03
This routine calculates the interfal forces
</pre>

\param   iel       INT        (I) number of element nodes
\param  *T         DOUBLE     (I) stresses
\param   fac       DOUBLE     (I) integration factor
\param **bop       DOUBLE     (I) B-operator
\param  *fint      DOUBLE     (O) element internal force

\warning There is nothing special to this routine
\return void
\sa calling:  ---;
    caled by: ifstatic_ke();

*----------------------------------------------------------------------*/
void if_fint(INT      iel,
             DOUBLE  *T,
             DOUBLE   fac,
             DOUBLE **bop,
             DOUBLE  *fint)
{
/*----------------------------------------------------------------------*/
INT i;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("if_fint");
#endif
/*----------------------------- updated lagrange or geometric linear ---*/
  for (i=0; i<2*iel; i++)
  {
    fint[i] += bop[0][i]*T[0]*fac + bop[1][i]*T[1]*fac;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of if_fint */




/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
#endif
