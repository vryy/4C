/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "xfem_prototypes.h"
/*! 
\addtogroup XFEM 
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA            *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC      *fdyn;



/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kvv

<pre>                                                            irhan 05/04

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
\param   index	   INT  	   (i)	  index for local assembly
\param   DENS	   DOUBLE  	   (i)	  fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calkvv(
  ELEMENT         *ele,
  DOUBLE         **estif,   
  DOUBLE          *velint,
  DOUBLE          *gridvint,
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     c,aux;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calkvv");
#endif		
/*----------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  c = fac*visc;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /
 *----------------------------------------------------------------------*/  
  for (icn=0; icn<TWO*iel; icn++) 
  {
    icol = index[icn];
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];      
      estif[irow  ][icol  ] += c*(TWO*derxy[0][irn]*derxy[0][icn] +
                                      derxy[1][irn]*derxy[1][icn]); 	     
      estif[irow+1][icol  ] += c*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += c*(TWO*derxy[1][irn]*derxy[1][icn] +
                                      derxy[0][irn]*derxy[0][icn]);
      estif[irow  ][icol+1] += c*(    derxy[1][irn]*derxy[0][icn]);
    }
  }

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nc(u):
    /
   |  v * u_old * grad(u)     d_omega
  /
 *----------------------------------------------------------------------*/
  if(fdyn->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];  
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][icn] +
               velint[1]*derxy[1][icn])*funct[irn]*fac*DENS;
        estif[irow  ][icol  ] += aux;
        estif[irow+1][icol+1] += aux;
      }
    }
  }
  
/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nr(u):
    /
   |  v * u * grad(u_old)     d_omega
  /
 *----------------------------------------------------------------------*/
  if (fdyn->nir!=0) /* evaluate for Newton iteration */
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];  
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = funct[irn]*funct[icn]*fac;
        estif[irow  ][icol  ] += aux*vderxy[0][0]*DENS;
        estif[irow+1][icol  ] += aux*vderxy[1][0]*DENS;
        estif[irow+1][icol+1] += aux*vderxy[1][1]*DENS;
        estif[irow  ][icol+1] += aux*vderxy[0][1]*DENS;
      }
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calkvv */



/*!--------------------------------------------------------------------- 
\brief evaluate galerkin part of Kvp

  
<pre>                                                            irhan 05/04

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
\param   index	   INT  	   (i)	  index for local assembly
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calkvp(
  DOUBLE         **estif,   
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  INT              iel,
  INT             *index     
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
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,ird;  
  INT        posc;
  DOUBLE     aux;

#ifdef DEBUG 
  dstrc_enter("xfem_f2_calkvp");
#endif		
/*----------------------------------------------------------------------*/

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
  for (icol=0; icol<iel; icol++)
  {
    posc = icol + 2*iel;
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];
      for(ird=0; ird<2; ird++)
      {      
	aux = -funct[icol]*derxy[ird][irn]*fac;
	estif[irow][posc] += aux;
	estif[posc][irow] += aux;
        irow++;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calkvp */



/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Mvv

<pre>                                                            irhan 05/04

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
\param   index	   INT  	   (i)	  index for local assembly
\param   DENS	   DOUBLE  	   (i)	  fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calmvv(
  DOUBLE         **emass,  
  DOUBLE          *funct, 
  DOUBLE           fac,   
  INT              iel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;  
  DOUBLE     aux;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calmvv");
#endif		
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mvv:
    /
   |  v * u    d_omega
  /
 *----------------------------------------------------------------------*/
  for (icn=0; icn<TWO*iel; icn++) 
  {
    icol = index[icn];
    for (irn=0; irn<TWO*iel; irn++)
    {
      irow = index[irn];      
      aux = funct[icn]*funct[irn]*fac;
      emass[irow  ][icol  ] += aux*DENS;
      emass[irow+1][icol+1] += aux*DENS;   
    }
  }

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of xfem_f2_calmvv */
/*! @} (documentation module close)*/
#endif
