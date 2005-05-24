/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of submesh matrices for fluid2

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
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
\brief evaluate Galerkin part of submesh stiffness matrix SMK for fluid2

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh stiffness matrix SMK
is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void

------------------------------------------------------------------------*/
void f2_calsmk(FLUID_DYN_ML    *mlvar,
	       DOUBLE         **smestif,
	       DOUBLE          *velint,
	       DOUBLE         **vderxy,
	       DOUBLE          *smfunct,
	       DOUBLE         **smderxy,
	       DOUBLE           fac,
	       DOUBLE           visc,
	       INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT     irow,icol;
DOUBLE  con,aux,beta,divv;

#ifdef DEBUG
dstrc_enter("f2_calsmk");
#endif

/*-------------------------------------------------- subgrid viscosity */
fdyn = alldyn[genprob.numff].fdyn;

if (mlvar->smsgvi>0) con=fac*(visc+mlvar->smsgvisc);
else con=fac*visc;

/*----------------------------------------------------------------------*
   Calculate viscous part of matrix SMK (including subgrid viscosity):
    /
   |  (nue+nueT) * grad(w) : grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  for (irow=0;irow<smiel;irow++)
  {
    smestif[irow][icol] += con*(smderxy[0][irow]*smderxy[0][icol]\
                               +smderxy[1][irow]*smderxy[1][icol]);
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*
   Calculate convective part of matrix SMK:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  w * u_old[ls_u_old] * grad(bub)     d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  aux = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol])*fac;
  for (irow=0;irow<smiel;irow++)
  {
    smestif[irow][icol] += aux*smfunct[irow];
  } /* end loop over irow */
}/* end loop over icol */

if (fdyn->conte!=0)
{
/*----------------------------------------------------------------------*
    /
   | beta * w * bub * div(u_old[ls_u_old])   d_omega
  /
 *----------------------------------------------------------------------*/
  if (fdyn->conte==1) beta = ONE;
  else beta = ONE/TWO;
/*---------------------------------------------- divergence of velocity */
  divv= vderxy[0][0]+vderxy[1][1];

  for (icol=0;icol<smiel;icol++)
  {
    aux = beta*smfunct[icol]*divv*fac;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end loop over irow */
  } /* end loop over icol */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calsmk */

/*!---------------------------------------------------------------------
\brief evaluate Galerkin part of submesh mass matrix SMM for fluid2

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh mass matrix SMM is
calculated.

</pre>
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void

------------------------------------------------------------------------*/
void f2_calsmm(DOUBLE         **smemass,
	       DOUBLE          *smfunct,
	       DOUBLE           fac,
	       INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT     irow,icol;

#ifdef DEBUG
dstrc_enter("f2_calsmm");
#endif

/*----------------------------------------------------------------------*
   Calculate matrix SMM:
    /
   |  v * bub  d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  for (irow=0;irow<smiel;irow++)
  {
    smemass[irow][icol] += smfunct[irow]*smfunct[icol]*fac;
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calsmm */

/*!---------------------------------------------------------------------
\brief evaluate diffusive part of submesh stiffn. matrix SMK for fluid2

<pre>                                                       gravem 07/03

In this routine, the diffusive Galerkin part of the submesh stiffness
matrix SMK is calculated.

</pre>
\param **smiediff  DOUBLE	   (i/o)  sm ele stiffn. matrix (diff.)
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void

------------------------------------------------------------------------*/
void f2_calsmkd(DOUBLE         **smiediff,
        	DOUBLE         **smderxy,
	        DOUBLE           fac,
	        INT              smiel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT     irow,icol;

#ifdef DEBUG
dstrc_enter("f2_calsmkd");
#endif

/*----------------------------------------------------------------------*
   Calculate diffusive (viscous) part of matrix SMK (inc. subgrid visc.):
    /
   |  grad(w) : grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  for (irow=0;irow<smiel;irow++)
  {
    smiediff[irow][icol] += fac*(smderxy[0][irow]*smderxy[0][icol]\
                                +smderxy[1][irow]*smderxy[1][icol]);
  } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calsmkd */

/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kvv

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Kvv is calculated:

    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /

    /
   |  v * u_old * grad(u)     d_omega
  /

    /
   |  v * u * grad(u_old)     d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param  *velint    DOUBLE	   (i)    vel at int point
\param **vderxy    DOUBLE	   (i)    global vel derivatives
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   iel	   INT  	   (i)	  number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_lscalkvv(DOUBLE	      **estif,
	       DOUBLE	       *velint,
	       DOUBLE	      **vderxy,
	       DOUBLE	       *funct,
	       DOUBLE	      **derxy,
	       DOUBLE		fac,
	       DOUBLE		visc,
	       INT		iel)
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
dstrc_enter("f2_lscalkvv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*------------------------------------------------- subgrid viscosity? */
if (fdyn->sgvisc>0) con=fac*(visc+fdyn->sugrvisc);
else con=fac*visc;

if (fdyn->vite==0)
{
/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K (including subgrid visc.):
    /
   |  (nue+nueT) * grad(v) : grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][icol]     += con*derxy[0][irn]*derxy[0][icn];
      estif[irow][icol+1]   += con*derxy[1][irn]*derxy[0][icn];
      estif[irow+1][icol]   += con*derxy[0][irn]*derxy[1][icn];
      estif[irow+1][icol+1] += con*derxy[1][irn]*derxy[1][icn];
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */
}
else
{
/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K (including subgrid visc.):
    /
   |  2 * (nue+nueT) * eps(v) : eps(u)   d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      estif[irow][icol]     += con*(TWO*derxy[0][irn]*derxy[0][icn] \
                                    + derxy[1][irn]*derxy[1][icn]);
      estif[irow][icol+1]   += con*(    derxy[1][irn]*derxy[0][icn]);
      estif[irow+1][icol]   += con*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += con*(TWO*derxy[1][irn]*derxy[1][icn] \
                                    + derxy[0][irn]*derxy[0][icn]);
      irow += 2;
    } /* end loop over irn */
    icol += 2;
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nc(u)
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  v * u_old * grad(u)     d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      aux = (velint[0]*derxy[0][icn] \
            +velint[1]*derxy[1][icn])*funct[irn]*fac;
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
   | beta * v * u * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  if (fdyn->conte==1) beta = ONE;
  else beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1];
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    aux = beta*funct[icn]*divv*fac;
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
   Calculate full Galerkin part of matrix Nr(u):
 *----------------------------------------------------------------------*/
if (fdyn->nir != 0) /* evaluate for Newton iteraton */
{
/*----------------------------------------------------------------------*
    /
   |  v * u * grad(u_old)     d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      aux = funct[irn]*funct[icn]*fac;
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
   |  beta * v * u_old * div(u)     d_omega
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
        estif[irow][icol]     += aux*velint[0]*derxy[0][icn];
        estif[irow][icol+1]   += aux*velint[0]*derxy[1][icn];
        estif[irow+1][icol]   += aux*velint[1]*derxy[0][icn];
        estif[irow+1][icol+1] += aux*velint[1]*derxy[1][icn];
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
} /* end of f2_lscalkvv */

/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kvp

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Kvp is calculated:

    /
   |  - div(v) * p     d_omega
  /

    /
   | - q * div(u)      d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: there's only one elestif
      --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]
      --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)]

</pre>
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac       DOUBLE	   (i)    weighting factor
\param   iel       INT  	   (i)    number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_lscalkvp(DOUBLE         **estif,
	       DOUBLE	       *funct,
	       DOUBLE	      **derxy,
	       DOUBLE		fac,
	       INT		iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |
 *----------------------------------------------------------------------*/
INT     irow, icol,irn,ird;
INT     posc;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_lscalkvp");
#endif

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kvp:
    /
   |  - div(v) * p     d_omega
  /

   and matrix Kpv:
    /
   | - q * div(u)      d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  irow=-1;
  posc = icol + 2*iel;
  for (irn=0;irn<iel;irn++)
  {
     for(ird=0;ird<2;ird++)
     {
	aux = funct[icol]*derxy[ird][irn]*fac;
	irow++;
	estif[irow][posc] -= aux;
	estif[posc][irow] -= aux;
     } /* end loop over ird */
  } /* end loop over irn */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalkvp */

/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Mvv

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Mvv is calculated:

    /
   |  v * u    d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: there's only one elemass
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]

</pre>
\param **emass     DOUBLE	   (i/o)  ele mass matrix
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac       DOUBLE	   (i)    weighting factor
\param   iel       INT   	   (i)    number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_lscalmvv(DOUBLE         **emass,
	       DOUBLE	       *funct,
	       DOUBLE		fac,
	       INT		iel)
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
dstrc_enter("f2_lscalmvv");
#endif

nvdfe = NUM_F2_VELDOF*iel;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mvv:
    /
   |  v * u    d_omega
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
      aux = funct[icn]*funct[irn]*fac;
      emass[irow][icol]     += aux;
      emass[irow+1][icol+1] += aux;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalmvv */



#endif
