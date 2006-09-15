/*======================================================================*/
/*!
\file
\brief Coordinate transformations

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Transform local/global heat flux and temperature gradient

Local 3-direction is zero

\author bborn
\date 03/06
*/
void th2_trafo_lss(DOUBLE     *a,
                   DOUBLE    **G,
                   DOUBLE    **GI,
                   INT         it)
{
/*----------------------------------------------------------------------*/
INT numeps=4;
INT i;
DOUBLE ga[4];
#ifdef DEBUG
dstrc_enter("th2_trafo_lss");
#endif
 /*----------------------------------------------------------------*
 /  a  ... vector to be transformed                       (i/o)    |
 |  G  ... elements of transformation matrix                (i)    |
 |  GI ... inverse of G                                     (i)    |
 |  it	=	0	Sig(global)  --->  Sig(local)       (i)    |
 |  it	=	1	Sig(local)   --->  Sig(global)      (i)    |
 |  it	=	2	Eps(global)  --->  Eps(local)       (i)    |
 |  it	=	3	Eps(local)   --->  Eps(global)      (i)    |
 |----------------------------------------------------------------*/
switch (it) {
	case 0: goto L1;
	case 1: goto L2;
	case 2: goto L3;
	case 3: goto L4;
}


 /*----------------------------------------------------------------*
 |  Sig(global)  --->  Sig(local)  S(loc) = G(Inv) * S(glo)       |
 |----------------------------------------------------------------*/
L1:
    math_matvecdense(ga,GI,a,numeps,numeps,0,1.);
    goto L100;

 /*----------------------------------------------------------------*
 |  Sig(local)   --->  Sig(global) S(glo) = G * S(loc)             |
 |----------------------------------------------------------------*/
L2:
    math_matvecdense(ga,G,a,numeps,numeps,0,1.);
    goto L100;

 /*----------------------------------------------------------------*
 |  Eps(global)  --->  Eps(local)  E(loc) = H(Inv) * E(glo)        |
 |                                          H(Inv) = G(Trans)      |
 |----------------------------------------------------------------*/
L3:
    math_mattrnvecdense(ga,G,a,numeps,numeps,0,1.);
    goto L100;

 /*----------------------------------------------------------------*
 |  Eps(local)   --->  Eps(global) E(loc) = H * E(glo)             |
 |                                          H = G(Inv)(Trans)      |
 |----------------------------------------------------------------*/
L4:
    math_mattrnvecdense(ga,GI,a,numeps,numeps,0,1.);

 /*----------------------------------------------------------------*
 |  Copy results of operation to vector a                          |
 |----------------------------------------------------------------*/
L100:

    for (i=0; i<numeps; i++)
    {
    	a[i] = ga[i];
    }
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of th2_trafo_lss */
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*!
\brief Set Transformation Matrices G and GI

\author bborn
\date 03/06
*/
void th2_trafo_sett(DOUBLE    **A,
                    DOUBLE    **B,
                    DOUBLE    **C)
{
  /*--------------------------------------------------------------------*/
  INT numeps=4;
  INT l;
  DOUBLE a11=0.0;
  DOUBLE a12=0.0;
  DOUBLE a21=0.0;
  DOUBLE a22=0.0;
  DOUBLE zero=0.;
  DOUBLE one=1.;
  INT i,j;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_trafo_sett");
#endif

  /*--------------------------------------------------------------------*/
  /*----------------------------------------------------------------*
  |  A  ... matrix of direction cosines                     (i)     |
  |  B  ... G-Matrix      (L=1)                             (o)     |
  |  C  ... inverse of G  (L=2)                             (o)     |
  *----------------------------------------------------------------*/

  for (l=1; l<3; l++)
  {
	if (l==1)
	{
	a11=A[0][0];
	a21=A[0][1];
	a12=A[1][0];
	a22=A[1][1];
	}
	else if (l==2)
	{
	a11=A[0][0];
	a21=A[1][0];
	a12=A[0][1];
	a22=A[1][1];
	}

	B[0][0]=a11*a11;
	B[1][0]=a12*a12;
	B[2][0]=a11*a12;
	B[3][0]=zero;

	B[0][1]=a21*a21;
	B[1][1]=a22*a22;
	B[2][1]=a21*a22;
	B[3][1]=zero;

	B[0][2]=a11*a21*2.;
	B[1][2]=a12*a22*2.;
	B[2][2]=a11*a22+a21*a12;
	B[3][2]=zero;

	B[0][3]=zero;
	B[1][3]=zero;
	B[2][3]=zero;
	B[3][3]=one;

	if (l==1)
	{
	for (i=0; i<numeps; i++)
	{
	   for (j=0; j<numeps; j++)
	   {
	      C[j][i]=B[j][i];
	   }
        }
	}
}

  /*--------------------------------------------------------------------*/ 
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_trafo_sett */





/*----------------------------------------------------------------------*/
/*!
\brief Calculate Transformation Matrices G and G(Inv)

\param  **xjm       DOUBLE   (i)    Jacobi matrix
\param  **G         DOUBLE   (o)    transformation matrix 
\param  **GI        DOUBLE   (o)    transformation matrix 
\param  **dum       DOUBLE   (io)   working array
                                      (reasing behind this is totally
                                      unclear)

\author bborn
\date 03/06
*/
void th2_trafo_tram(DOUBLE   **xjm,
                    DOUBLE    **G,
                    DOUBLE    **GI,
                    DOUBLE    **dum)
{
  /*--------------------------------------------------------------------*/
  const INT ndim = 2;  /* genprob.ndim; */ /* must be 2 */
  DOUBLE xr[ndim];
  DOUBLE xs[ndim];
  DOUBLE xc[2];
  DOUBLE x3n;
  DOUBLE one=1.;
  DOUBLE length;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_trafo_tram");
#endif

  /*--------------------------------------------------------------------*/
  /*----------------------------------------------------------------*
  |  Local Axes : x = Tangent at r-axes                             |
  |               y = cross product z*x                             |
  |               z = normal at surface                             |
  |-----------------------------------------------------------------|
  |  xjm --- Jacobian matrix r,s-direction                          |
  |  G   --- Transformation Matrix S(Global) = G * S(Local)         |
  |  GI  --- Inverse of G          S(Local)  = G * S(Global)        |
  *----------------------------------------------------------------*/

 
  /*----------------------------------------------------------------*
  |  Calculate local coordinate axes                                |
  *----------------------------------------------------------------*/
  xr[0]=xjm[0][0];
  xr[1]=xjm[0][1];
  xs[0]=xjm[1][0];
  xs[1]=xjm[1][1];

  x3n=one;

  xc[0]=-xr[1];
  xc[1]=xr[0];

  /*----------------------------------------------------------------*
  |  Reduce to unit vectors uxr, uxc                                |
  *----------------------------------------------------------------*/
  math_unvc(&length,xr,ndim);
  math_unvc(&length,xc,ndim);

  /*----------------------------------------------------------------*
  |  Set matrix of direction cosines                                |
  *----------------------------------------------------------------*/
  dum[0][0]=xr[0];
  dum[0][1]=xr[1];
  dum[1][0]=xc[0];
  dum[1][1]=xc[1];

  /*----------------------------------------------------------------*
  |  Calculate final transformation matrix G and G(Inv)             |
  *----------------------------------------------------------------*/
  th2_trafo_sett(dum,G,GI);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
return;
}  /* end of th2_trafo_tram */


/*======================================================================*/
#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
