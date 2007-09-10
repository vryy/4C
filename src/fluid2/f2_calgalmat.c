/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
\brief evaluate galerkin part of Kvv

<pre>                                                         genk 04/02
                                           modified for ALE   genk 10/02

In this routine the galerkin part of matrix Kvv is calculated:

EULER/ALE:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /

EULER/ALE:
    /
   |  v * u_old * grad(u)     d_omega
  /

    /
   |  v * u * grad(u_old)     d_omega
  /

ALE:

ale-convective velocity is split into
 c = u - u_G
   --> the known nonlinear term known from the EULER-case
   --> a new term:

    /
 - |  v * u_G * grad(u)     d_omega
  /

   --> this is a linear term inside the domain and for explicit treatement
       of  a free surface
   --> for implicit treatement of the free surface this term is nonlinear
       so it depends on the nonlinear iteration scheme if this term has to
       be taken into account on the LHS!

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]

</pre>
\param  *ele       ELEMENT         (i)    actual element
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param  *velint    DOUBLE	   (i)    vel at INT point
\param  *gridvint  DOUBLE          (i)    gridvel at INT point
\param **vderxy    DOUBLE	   (i)    global vel derivatives
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   iel	   INT  	   (i)	  number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_calkvv( ELEMENT         *ele,
		DOUBLE         **estif,
		DOUBLE          *velint,
		DOUBLE          *gridvint,
		DOUBLE         **vderxy,
		DOUBLE          *funct,
		DOUBLE         **derxy,
		DOUBLE           fac,
		DOUBLE           visc,
		INT              iel
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;
DOUBLE  c,aux;

#ifdef DEBUG
dstrc_enter("f2_calkvv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
c=fac*visc;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      estif[irow][icol]     += c*(TWO*derxy[0][irn]*derxy[0][icn] \
                                    + derxy[1][irn]*derxy[1][icn]);
      estif[irow+1][icol]   += c*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += c*(TWO*derxy[1][irn]*derxy[1][icn] \
                                    + derxy[0][irn]*derxy[0][icn]);
      estif[irow][icol+1]   += c*(    derxy[1][irn]*derxy[0][icn]);
      irow += 2;
   } /* end loop over irn */
   icol += 2;
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nc(u):
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

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nr(u):
    /
   |  v * u * grad(u_old)     d_omega
  /
 *----------------------------------------------------------------------*/
if (fdyn->nir != 0) /* evaluate for Newton iteraton */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
         aux = funct[irn]*funct[icn]*fac;
         estif[irow][icol]     += aux*vderxy[0][0];
         estif[irow+1][icol]   += aux*vderxy[1][0];
         estif[irow+1][icol+1] += aux*vderxy[1][1];
         estif[irow][icol+1]   += aux*vderxy[0][1];
         irow += 2;
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* endif (fdyn->nir != 0) */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part due to split of ALE-convective velocity:
     /
  - |  v * u_G_old * grad(u)     d_omega
   /
   REMARK:
    - this term is linear inside the domain and for explicit free surface
    - for implicit free surface this term is nonlinear (Nc), since u_G is
      unknown at the free surface. So it depends on the nonlinear iteration
      scheme if this term has to be taken into account:
	 - Newton: YES
	 - fixpoint-like: YES
 *----------------------------------------------------------------------*/
if(ele->e.f2->is_ale !=0) /* evaluate only for ALE */
{
   dsassert(gridvint!=NULL,"no grid velocity calculated!\n");
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
         aux = (gridvint[0]*derxy[0][icn] \
               +gridvint[1]*derxy[1][icn])*funct[irn]*fac;
         estif[irow][icol]     -= aux;
         estif[irow+1][icol+1] -= aux;
         irow += 2;
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* endif (is_ale !=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calkvv */

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
\param **estif     DOUBLE     (i/o)  ele stiffness matrix
\param  *funct     DOUBLE     (i)    nat. shape funcs
\param **derxy     DOUBLE     (i)    global coord. deriv.
\param   fac       DOUBLE     (i)    weighting factor
\param   iel       INT        (i)    number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_calkvp(
               DOUBLE         **estif,
               DOUBLE          *funct,
               DOUBLE         **derxy,
               DOUBLE           fac,
               INT              iel
              )
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
dstrc_enter("f2_calkvp");
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
} /* end of f2_calkvp */
/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kvg

<pre>                                                         genk 01/03

In this routine the galerkin part of matrix Kvg is calculated:

ALE:

    /
 - |  v * u_G * grad(u)     d_omega
  /


NOTE: there's only one elestif
      --> Kvg is stored in estif[0..(2*iel-1)](3*iel)..(5*iel-1)]

</pre>
\param **estif     DOUBLE           (i/o)  ele stiffness matrix
\param **vderxy    DOUBLE           (i)    global vel derivatives
\param  *funct     DOUBLE           (i)    nat. shape funcs
\param   fac       DOUBLE           (i)    weighting factor
\param   iel       INT              (i)    number of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_calkvg(
               DOUBLE         **estif,
               DOUBLE         **vderxy,
               DOUBLE          *funct,
               DOUBLE           fac,
               INT              iel
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calkvg");
#endif


/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nr(u_G):
       /
 (-)  |  v * u_G * grad(u_old)     d_omega
     /
   REMARK:
   for implicit treatement this is a reactive nonlinear term
   (u_G is unknown in this case) which is only evlueated for
   Newton iteration
 *----------------------------------------------------------------------*/
if (fdyn->nir != 0) /* evaluate for Newton iteraton */
{
   icol=NUMDOF_FLUID2*iel;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
         aux = funct[irn]*funct[icn]*fac;
         estif[irow][icol]     -= aux*vderxy[0][0];
         estif[irow+1][icol]   -= aux*vderxy[1][0];
         estif[irow+1][icol+1] -= aux*vderxy[1][1];
         estif[irow][icol+1]   -= aux*vderxy[0][1];
         irow += 2;
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* endif (fdyn->nir != 0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calkvg */

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
void f2_calmvv(
               DOUBLE         **emass,
               DOUBLE          *funct,
               DOUBLE           fac,
               INT              iel
              )
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
dstrc_enter("f2_calmvv");
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
} /* end of f2_calmvv */

/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kgv and Kgg

<pre>                                                         genk 01/03

In this routine the galerkin part of matrix Kgv is calculated:

    /
   |  w * u    d_gamma_FS
  /

In this routine the galerkin part of matrix Kgg is calculated:

     /
  - |  w * u_G    d_gamma_FS
   /

NOTE: there's only one estif

</pre>
\param **estif     DOUBLE     (i/o)  ele stiffness matrix
\param  *funct     DOUBLE     (i)    nat. shape funcs at edge
\param   fac       DOUBLE     (i)    weighting factor
\param   iel       INT        (i)    number of nodes of act. element
\param  *iedgnod   INT        (i)    edge nodes
\param   ngnode    INT        (i)    number of nodes of act. edge
\return void

------------------------------------------------------------------------*/
void f2_calkgedge(
                  DOUBLE         **estif,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT             *iedgnod,
                  INT              iel,
                  INT              ngnode
                )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |
 *----------------------------------------------------------------------*/
INT     irow,icolv,icolg,irn,icn;
INT     nd;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calkgedge");
#endif

nd = NUMDOF_FLUID2*iel;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kgv and Kgg:
      /
     |  w * u    d_gamma_FS
    /

     /
  - |  w * u_G    d_gamma_FS
   /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icolv = NUM_F2_VELDOF*iedgnod[icn];
   icolg = icolv+nd;
   for(irn=0;irn<ngnode;irn++)
   {
      irow = NUM_F2_VELDOF*iedgnod[irn]+nd;
      aux = funct[icn]*funct[irn]*fac;
      estif[irow][icolv]     += aux;
      estif[irow+1][icolv+1] += aux;
      estif[irow][icolg]     -= aux;
      estif[irow+1][icolg+1] -= aux;
   } /* end loop over irow */
} /* end loop over icol */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calkgedge */


/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kgg

<pre>                                                         genk 01/03

NOTE: there's only one estif

</pre>
\param **estif     DOUBLE     (i/o)  ele stiffness matrix
\param  *funct     DOUBLE     (i)    nat. shape funcs at edge
\param   fac       DOUBLE     (i)    weighting factor
\param   iel       INT        (i)    number of nodes of act. element
\param  *iedgnod   INT        (i)    edge nodes
\param   ngnode    INT        (i)    number of nodes of act. edge
\return void

------------------------------------------------------------------------*/
void f2_calgfskgg(
                  DOUBLE         **estif,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **xjm,
                  DOUBLE           det,
                  DOUBLE           fac,
                  INT              iel

                )
{
INT inode, node_start;
INT nd;
INT icol,irow;
INT i,j,k,l,m;
DOUBLE bop[3][2*MAXNOD_F2];
DOUBLE xji[2][2];
DOUBLE d[3][3];
DOUBLE db[3];
DOUBLE e1,e2,e3;
DOUBLE dum;

#ifdef DEBUG
dstrc_enter("f2_calgfskgg");
#endif

nd = 2*iel;

/*---------------------------------------------------------- initialise */
for (i=0;i<iel*2;i++)
{
   bop[0][i]=ZERO;
   bop[1][i]=ZERO;
   bop[2][i]=ZERO;
}

/*---------------------------------------------- inverse of jacobian ---*/
xji[0][0] = xjm[1][1]/ det;
xji[0][1] =-xjm[0][1]/ det;
xji[1][0] =-xjm[1][0]/ det;
xji[1][1] = xjm[0][0]/ det;

/*----------------------------- get operator b of global derivatives ---*/
for (inode=0; inode<iel; inode++)
{
  node_start = inode*2;

  bop[0][node_start+0] += xji[0][0] * deriv[0][inode]
                       +  xji[0][1] * deriv[1][inode];
  bop[1][node_start+1] += xji[1][0] * deriv[0][inode]
                       +  xji[1][1] * deriv[1][inode];
  bop[2][node_start+1] = bop[0][node_start+0];
  bop[2][node_start+0] = bop[1][node_start+1];
} /* end of loop over nodes */

/*---------------------------------------- get ficticous material data */
e1=ONE;
e2=ZERO;
e3=e1/TWO;

d[0][0]=e1;
d[0][1]=e2;
d[0][2]=ZERO;
d[1][0]=e2;
d[1][1]=e1;
d[1][2]=ZERO;
d[2][0]=ZERO;
d[2][1]=ZERO;
d[2][2]=e3;

/*---------------------------------------- get element stiffness matrix */
icol = iel * NUMDOF_FLUID2;
for (j=0; j<nd; j++)
{
   for (k=0; k<3; k++)
   {
      db[k] = ZERO ;
      for (l=0; l<3; l++)
      {
         db[k] += d[k][l]*bop[l][j]*fac ;
      }
   }
   irow = iel * NUMDOF_FLUID2;
   for (i=0; i<nd; i+=2)
   {
      dum = ZERO ;
      for (m=0; m<3; m++)
      {
         dum += bop[m][i]*db[m] ;
      }
      estif[irow+i][icol] += dum ;
   }
   icol++;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calkgedge */

/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kgg

<pre>                                                         genk 01/03

NOTE: there's only one estif

</pre>
\param **estif     DOUBLE     (i/o)  ele stiffness matrix
\param  *funct     DOUBLE     (i)    nat. shape funcs at edge
\param   fac       DOUBLE     (i)    weighting factor
\param   iel       INT        (i)    number of nodes of act. element
\param  *iedgnod   INT        (i)    edge nodes
\param   ngnode    INT        (i)    number of nodes of act. edge
\return void

------------------------------------------------------------------------*/
void f2_calgfskgedge(
                     DOUBLE         **estif,
                     DOUBLE          *funct,
                     DOUBLE          *vnint,
                     DOUBLE           fac,
                     INT             *iedgnod,
                     INT              iel,
                     INT              ngnode

                  )
{
INT     nd;
INT     irow, icolv,icolg;
INT     irn,icn;
DOUBLE  aux;

#ifdef DEBUG
dstrc_enter("f2_calgfskgedge");
#endif

nd = NUMDOF_FLUID2*iel;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kgv and Kgg:
      /
     |  w * n * u    d_gamma_FS
    /

     /
  - |  w * n * u_G    d_gamma_FS
   /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icolv = NUM_F2_VELDOF*iedgnod[icn];
   icolg = icolv+nd;
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+1+NUM_F2_VELDOF*iedgnod[irn];
      aux = funct[icn]*funct[irn]*fac;
      estif[irow][icolv]     += aux*vnint[0];
      estif[irow][icolv+1]   += aux*vnint[1];
      estif[irow][icolg]     -= aux*vnint[0];
      estif[irow][icolg+1]   -= aux*vnint[1];
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgfskgedge */

#endif
/*! @} (documentation module close)*/
#endif
