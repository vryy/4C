/*!----------------------------------------------------------------------
\file
\brief evaluate bubble parts of large-scale matrices for fluid2
<pre>                                                       gravem 07/03

Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Kvv for fluid2

<pre>                                                       gravem 07/03


In this routine, the bubble part of the matrix Kvv is calculated.

NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  element stiffne matrix
\param  *velint    DOUBLE	   (i)    velocity at int point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param  *vbubint   DOUBLE	   (i)    velocity bubble functions
\param **vbubderxy DOUBLE	   (i)    global deriv. of vel. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbkvv(DOUBLE         **estif,
		DOUBLE          *velint,
		DOUBLE         **vderxy,
		DOUBLE          *funct,
		DOUBLE         **derxy,
		DOUBLE          *vbubint,
		DOUBLE         **vbubderxy,
		DOUBLE           fac,
		DOUBLE           visc,
		INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow,icol,irn,icn;
DOUBLE  con,aux,beta,divv;

#ifdef DEBUG
dstrc_enter("f2_calbkvv");
#endif

/*--------------------------------------------------------------------- */
fdyn = alldyn[genprob.numff].fdyn;
con=fac*visc;

if (fdyn->vite==0)
{
/*----------------------------------------------------------------------*
   Calculate velocity bubble part of matrix K:
    /
   |  nue * grad(v) : grad(u_bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][icol]     += con*derxy[0][irn]*vbubderxy[0][icn];
      estif[irow][icol+1]   += con*derxy[1][irn]*vbubderxy[0][icn];
      estif[irow+1][icol]   += con*derxy[0][irn]*vbubderxy[1][icn];
      estif[irow+1][icol+1] += con*derxy[1][irn]*vbubderxy[1][icn];
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */
}
else
{
/*----------------------------------------------------------------------*
   Calculate velocity bubble part of matrix K:
    /
   |  2 * nue * eps(v) : eps(u_bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][icol]     += con*(TWO*derxy[0][irn]*vbubderxy[0][icn] \
                                    + derxy[1][irn]*vbubderxy[1][icn]);
      estif[irow][icol+1]   += con*(    derxy[1][irn]*vbubderxy[0][icn]);
      estif[irow+1][icol]   += con*(    derxy[0][irn]*vbubderxy[1][icn]);
      estif[irow+1][icol+1] += con*(TWO*derxy[1][irn]*vbubderxy[1][icn] \
                                    + derxy[0][irn]*vbubderxy[0][icn]);
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*
   Calculate velocity bubble part of matrix Nc(u)
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  v * u_old * grad(u_bub)     d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      aux = (velint[0]*vbubderxy[0][icn] \
            +velint[1]*vbubderxy[1][icn])*funct[irn]*fac;
      estif[irow][icol]     += aux;
      estif[irow+1][icol+1] += aux;
      irow += 2;
   } /* end loop over irn */
   icol += 2;
} /* end loop over icn */

if (fdyn->conte!=0)
{
/*----------------------------------------------------------------------*
    /
   | beta * v * u_bub * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  if (fdyn->conte==1) beta = ONE;
  else beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1];
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    aux = beta*vbubint[icn]*divv*fac;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][icol]     += aux*funct[irn];
      estif[irow+1][icol+1] += aux*funct[irn];
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*
   Calculate velocity bubble part of matrix Nr(u):
 *----------------------------------------------------------------------*/
if (fdyn->nir != 0) /* evaluate for Newton iteraton */
{
/*----------------------------------------------------------------------*
    /
   |  v * u_bub * grad(u_old)     d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      aux = funct[irn]*vbubint[icn]*fac;
      estif[irow][icol]     += aux*vderxy[0][0];
      estif[irow][icol+1]   += aux*vderxy[0][1];
      estif[irow+1][icol]   += aux*vderxy[1][0];
      estif[irow+1][icol+1] += aux*vderxy[1][1];
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */

  if (fdyn->conte!=0)
  {
/*----------------------------------------------------------------------*
    /
   |  beta * v * u_old * div(u_bub)     d_omega
  /
 *----------------------------------------------------------------------*/
    if (fdyn->conte==1) beta = ONE;
    else beta = ONE/TWO;
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
        aux = beta*funct[irn]*fac;
        estif[irow][icol]     += aux*velint[0]*vbubderxy[0][icn];
        estif[irow][icol+1]   += aux*velint[0]*vbubderxy[1][icn];
        estif[irow+1][icol]   += aux*velint[1]*vbubderxy[0][icn];
        estif[irow+1][icol+1] += aux*velint[1]*vbubderxy[1][icn];
        irow += 2;
      } /* end loop over irn */
      icol += 2;
    } /* end loop over icn */
  }
} /* endif (fdyn->nir != 0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbkvv */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Kvp for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kvp is calculated.

NOTE: there's only one elestif
      --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at INT point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **pbubint   DOUBLE	   (i)    pressure bubble functions
\param***pbubderxy DOUBLE	   (i)    global deriv. of pre. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbkvp(DOUBLE         **estif,
		DOUBLE          *velint,
		DOUBLE         **vderxy,
		DOUBLE          *funct,
		DOUBLE         **derxy,
		DOUBLE         **pbubint,
		DOUBLE        ***pbubderxy,
		DOUBLE           fac,
		DOUBLE           visc,
		INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 *----------------------------------------------------------------------*/
INT     irow,icol,irn,posc;
DOUBLE  con,aux,aux0,aux1,beta,divv;

#ifdef DEBUG
dstrc_enter("f2_calbkvp");
#endif

/*--------------------------------------------------------------------- */
fdyn = alldyn[genprob.numff].fdyn;

con=fac*visc;

if (fdyn->vite==0)
{
/*----------------------------------------------------------------------*
   Calculate pressure bubble part of matrix K:
                  /
 sum over j=1,2  |  nue * grad(v) : grad(p_bub(j))   d_omega
                /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    irow=0;
    posc = icol + 2*iel;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][posc]   += con*(derxy[0][irn]*pbubderxy[0][0][icol] \
                                 +derxy[1][irn]*pbubderxy[0][1][icol]);
      estif[irow+1][posc] += con*(derxy[1][irn]*pbubderxy[1][1][icol] \
                                 +derxy[0][irn]*pbubderxy[1][0][icol]);
      irow += 2;
    } /* end loop over irn */
  } /* end loop over icn */
}
else
{
/*----------------------------------------------------------------------*
   Calculate pressure bubble part of matrix K:
                    /
   sum over j=1,2  |  2 * nue * eps(v) : eps(p_bub(j))   d_omega
                  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    irow=0;
    posc = icol + 2*iel;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][posc]   += con*(TWO*derxy[0][irn]*pbubderxy[0][0][icol] \
              +derxy[1][irn]*(pbubderxy[1][0][icol]+pbubderxy[0][1][icol]));
      estif[irow+1][posc] += con*(TWO*derxy[1][irn]*pbubderxy[1][1][icol] \
              +derxy[0][irn]*(pbubderxy[1][0][icol]+pbubderxy[0][1][icol]));
      irow += 2;
    } /* end loop over irn */
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*
   Calculate pressure bubble part of matrix Nc(u)
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
                  /
 sum over j=1,2  |  v * u_old * grad(p_bub(j))     d_omega
                /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  aux0 = (velint[0]*pbubderxy[0][0][icol] \
         +velint[1]*pbubderxy[1][0][icol])*fac;
  aux1 = (velint[0]*pbubderxy[0][1][icol] \
         +velint[1]*pbubderxy[1][1][icol])*fac;
  irow=0;
  posc = icol + 2*iel;
  for (irn=0;irn<iel;irn++)
  {
    estif[irow][posc]   += aux0*funct[irn];
    estif[irow+1][posc] += aux1*funct[irn];
    irow += 2;
  } /* end loop over irn */
} /* end loop over icn */

if (fdyn->conte!=0)
{
/*----------------------------------------------------------------------*
                   /
  sum over j=1,2  | beta * v * p_bub(j) * div(u_old)   d_omega
                 /
 *----------------------------------------------------------------------*/
  if (fdyn->conte==1) beta = ONE;
  else beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1];
  for (icol=0;icol<iel;icol++)
  {
    irow=0;
    posc = icol + 2*iel;
    aux = beta*divv*fac;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][posc]   += aux*funct[irn]*pbubint[0][icol];
      estif[irow+1][posc] += aux*funct[irn]*pbubint[1][icol];
      irow += 2;
    } /* end loop over irn */
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*
   Calculate pressure bubble part of matrix Nr(u):
 *----------------------------------------------------------------------*/
if (fdyn->nir != 0) /* evaluate for Newton iteraton */
{
/*----------------------------------------------------------------------*
                   /
  sum over j=1,2  |  v * p_bub(j) * grad(u_old)     d_omega
                 /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    irow=0;
    posc = icol + 2*iel;
    aux0=(pbubint[0][icol]*vderxy[0][0]+pbubint[1][icol]*vderxy[0][1])*fac;
    aux1=(pbubint[1][icol]*vderxy[1][1]+pbubint[0][icol]*vderxy[1][0])*fac;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][posc]   += aux0*funct[irn];
      estif[irow+1][posc] += aux1*funct[irn];
      irow += 2;
    } /* end loop over irn */
  } /* end loop over icn */

  if (fdyn->conte!=0)
  {
/*----------------------------------------------------------------------*
                   /
  sum over j=1,2  |  beta * v * u_old * div(p_bub(j))     d_omega
                 /
 *----------------------------------------------------------------------*/
    if (fdyn->conte==1) beta = ONE;
    else beta = ONE/TWO;
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc = icol + 2*iel;
      aux=(pbubderxy[0][0][icol]+pbubderxy[1][1][icol])*fac*beta;
      for (irn=0;irn<iel;irn++)
      {
        estif[irow][posc]   += aux*velint[0]*funct[irn];
        estif[irow+1][posc] += aux*velint[1]*funct[irn];
        irow += 2;
      } /* end loop over irn */
    } /* end loop over icn */
  }
} /* endif (fdyn->nir != 0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbkvp */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Kpv for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kpv is calculated.

NOTE: there's only one elestif
      --> Kpv is stored in estif[(2*iel)..(3*iel-1)][0..(2*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **vbubderxy DOUBLE	   (i)    global deriv. of vel. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbkpv(DOUBLE         **estif,
		DOUBLE          *funct,
		DOUBLE         **vbubderxy,
		DOUBLE           fac,
		INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 *----------------------------------------------------------------------*/
INT     irow,icol,icn,icd,posr;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calbkpv");
#endif

/*----------------------------------------------------------------------*
   Calculate velocity bubble part of divergence matrix GT:
       /
      | - q * div(u_bub)    d_omega
     /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
  for (icd=0;icd<2;icd++)
  {
    aux=vbubderxy[icd][icn]*fac;
    for (irow=0;irow<iel;irow++)
    {
      posr = irow + 2*iel;
      estif[posr][icol] -= aux*funct[irow];
    } /* end loop over irn */
    icol++;
  }
} /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbkpv */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Kpp for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kpp is calculated.

NOTE: there's only one elestif
      --> Kpp is stored in estif[(2*iel)..(3*iel-1)][(2*iel)..(3*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param***pbubderxy DOUBLE	   (i)    global deriv. of pre. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbkpp(DOUBLE         **estif,
		DOUBLE          *funct,
		DOUBLE        ***pbubderxy,
		DOUBLE           fac,
		INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 *----------------------------------------------------------------------*/
INT     irow,icol,posr,posc;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calbkpp");
#endif

/*----------------------------------------------------------------------*
   Calculate pressure bubble part of divergence matrix GT:
                      /
   sum over j=1,2    | - q * div(p_bub(j))    d_omega
                    /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  posc = icol + 2*iel;
  aux=(pbubderxy[0][0][icol]+pbubderxy[1][1][icol])*fac;
  for (irow=0;irow<iel;irow++)
  {
    posr = irow + 2*iel;
    estif[posr][posc] -= aux*funct[irow];
  } /* end loop over irn */
} /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbkpp */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Mvv for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Mvv is calculated.

NOTE: there's only one elemass
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]

</pre>
\param **emass     DOUBLE	   (i/o)  element mass matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param  *vbubint   DOUBLE	   (i)    velocity bubble functions
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbmvv(DOUBLE         **emass,
	        DOUBLE          *funct,
	        DOUBLE          *vbubint,
	        DOUBLE           fac,
	        INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;
INT     nvdfe;             /* number of velocity dofs of actual element */
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calbmvv");
#endif

nvdfe = NUM_F2_VELDOF*iel;

/*----------------------------------------------------------------------*
   Calculate velocity bubble part of matrix M:
    /
   |  v * u_bub    d_omega
  /
 *----------------------------------------------------------------------*/
icn=-1;
for(icol=0;icol<nvdfe;icol+=2)
{
   icn++;
   irn=-1;
   for(irow=0;irow<nvdfe;irow+=2)
   {
      irn++;
      aux = vbubint[icn]*funct[irn]*fac;
      emass[irow][icol]     += aux;
      emass[irow+1][icol+1] += aux;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbmvv */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of matrix Mvp for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Mvp is calculated.

NOTE: there's only one elemass
      --> Mvp is stored in emass[0..(2*iel-1)][(2*iel)..(3*iel-1)]

</pre>
\param **emass     DOUBLE	   (i/o)  element mass matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **pbubint   DOUBLE	   (i)    pressure bubble functions
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbmvp(DOUBLE         **emass,
	        DOUBLE          *funct,
	        DOUBLE         **pbubint,
	        DOUBLE           fac,
	        INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow,icol,irn,posc;
INT     nvdfe;             /* number of velocity dofs of actual element */

#ifdef DEBUG
dstrc_enter("f2_calbmvp");
#endif

nvdfe = NUM_F2_VELDOF*iel;

/*----------------------------------------------------------------------*
   Calculate pressure bubble part of matrix M:
                    /
  sum over j=1,2   |  v * p_bub(j)    d_omega
                  /
 *----------------------------------------------------------------------*/
for(icol=0;icol<iel;icol++)
{
  irn=-1;
  posc = 2*iel+icol;
  for(irow=0;irow<nvdfe;irow+=2)
  {
    irn++;
    emass[irow][posc]   += pbubint[0][icol]*funct[irn]*fac;
    emass[irow+1][posc] += pbubint[1][icol]*funct[irn]*fac;
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbmvp */

#endif
