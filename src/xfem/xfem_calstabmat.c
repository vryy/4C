/*!----------------------------------------------------------------------
\file
\brief stabilisation part of element stiffness matrix for fluid2

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
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
\brief evaluate stabilisaton part of Kvv

<pre>                                                            irhan 05/04
					    
In this routine the stabilisation part of matrix Kvv is calculated:

EULER/ALE:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /

EULER:

    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
  
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /  

EULER:
    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /

EULER/ALE:
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v))  div(eps(u))d_omega
  /  

EULER/ALE:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /
  
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]
      
NOTE: there are two velocities at integration point
      EULER: velint  = vel2int = u_old

NOTE: for EULER-case grid-velocity is not used      

</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *vel2int   DOUBLE          (i)     vel. at integr. point
\param  *gridvint  DOUBLE          (i)     grid vel at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)	   weighting factor	   
\param   visc	   DOUBLE	   (i)	   fluid viscosity
\param   iel	   INT		   (i)	   num. of nodes in ele
\param   ihoel	   INT		   (i)	   flag for higer ord. ele
\param   index	   INT  	   (i)	   index for local assembly
\param   DENS	   DOUBLE  	   (i)	   fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabkvv(			      
  ELEMENT         *ele,    
  DOUBLE         **estif,  
  DOUBLE          *velint,
  DOUBLE          *vel2int, 
  DOUBLE          *gridvint,
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |  
 *----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     tauc;
  DOUBLE     c,cc;
  DOUBLE     aux,auxr,auxc;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkvv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  tauc  = fdyn->tau[2];

/*----------------------------------------------------------------------*
   Calculate continuity stabilisation part:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
  if (gls->icont!=0)
  {
    c = fac*tauc;
    for(icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      for(irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        estif[irow  ][icol  ] += derxy[0][icn]*derxy[0][irn]*c;
        estif[irow+1][icol  ] += derxy[0][icn]*derxy[1][irn]*c;
        estif[irow  ][icol+1] += derxy[1][icn]*derxy[0][irn]*c;
        estif[irow+1][icol+1] += derxy[1][icn]*derxy[1][irn]*c;
      }
    }
  }
  
  c = fac*taumu;
  cc = c;
/*------------------------------ calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nc(u):
EULER:   
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega   
  /
 *----------------------------------------------------------------------*/    
      for (icn=0; icn<TWO*iel; icn++)
      {
        auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*cc;
        icol = index[icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*auxc;
          estif[irow  ][icol  ] += aux*DENS*DENS;
          estif[irow+1][icol+1] += aux*DENS*DENS;
        }
      }
    
/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nr(u):
EULER:   
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
    if (fdyn->nir!=0) /* evaluate for Newton iteration */
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        auxc = funct[icn]*cc;
        icol = index[icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];          
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*auxc;
          estif[irow  ][icol  ] += aux*vderxy[0][0]*DENS*DENS;
          estif[irow+1][icol  ] += aux*vderxy[1][0]*DENS*DENS;
          estif[irow  ][icol+1] += aux*vderxy[0][1]*DENS*DENS;
          estif[irow+1][icol+1] += aux*vderxy[1][1]*DENS*DENS;
        }
      }
    }

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part for higher order elements:
EULER:
    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
    if (ihoel!=0)
    {
      cc = c*visc;
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];        
        auxc = derxy2[0][icn] + derxy2[1][icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                    
          aux = (vel2int[0]*derxy[0][irn] + vel2int[1]*derxy[1][irn])*cc;
          estif[irow  ][icol  ] -= aux*(derxy2[0][icn] + auxc)*DENS;
          estif[irow+1][icol  ] -= aux*derxy2[2][icn]*DENS;;
          estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc)*DENS;;
          estif[irow  ][icol+1] -= aux*derxy2[2][icn]*DENS;;
        }
      }
    }
  }
  
/*-------------------------------- calculate viscous stabilisation part */
  if (ihoel!=0 && gls->ivisc!=0)
  {   
    switch (gls->ivisc) /* choose stabilisation type --> sign */
    {
        case 1: /* GLS- */
          sign = ONE;
          break;
        case 2: /* GLS+ */
          sign = -ONE;
          break;
        default:
          dserror("viscous stabilisation parameter unknown: IVISC");
    }

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part for higher order elements:
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) * div(eps(u))d_omega
  /
 *----------------------------------------------------------------------*/   
    cc = fac * taump * visc*visc * sign;
    
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];              
      auxc = derxy2[0][icn] + derxy2[1][icn];
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];        
        auxr = derxy2[0][irn] + derxy2[1][irn];
        aux  = auxc*auxr;
        estif[irow  ][icol  ] +=  (derxy2[0][icn]*(derxy2[0][irn]+auxr) +
                                   derxy2[2][icn]* derxy2[2][irn] +
                                   derxy2[0][irn]* auxc + aux)*cc;
        estif[irow+1][icol  ] +=  (derxy2[0][icn]* derxy2[2][irn] +
                                   derxy2[2][icn]*(derxy2[1][irn]+auxr) +
                                   derxy2[2][irn]* auxc)*cc;  
        estif[irow+1][icol+1] +=  (derxy2[2][icn]* derxy2[2][irn] +
                                   derxy2[1][icn]*(derxy2[1][irn]+auxr) +
                                   derxy2[1][irn]* auxc + aux)*cc;
        estif[irow  ][icol+1] +=  (derxy2[2][icn]*(derxy2[0][irn]+auxr) +
                                   derxy2[1][icn]* derxy2[2][irn] +
                                   derxy2[2][irn]* auxc)*cc;
      }
    }

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
    cc = fac * taump * visc * sign;
    
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];              
        aux = velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn];
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                  
          auxr = derxy2[0][irn] + derxy2[1][irn];
          estif[irow  ][icol  ] -= (derxy2[0][irn]+auxr)*aux*cc*DENS;
          estif[irow+1][icol  ] -=  derxy2[2][irn]*aux*cc*DENS;
          estif[irow  ][icol+1] -=  derxy2[2][irn]*aux*cc*DENS;
          estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*aux*cc*DENS;
        }
      }

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
    if (fdyn->nir!=0)  /* evaluate for Newton iteration */
    {
      for (icn=0; icn<TWO*iel; icn++)
      {
        icol = index[icn];              
        aux = funct[icn]*cc;
        for (irn=0; irn<TWO*iel; irn++)
        {
          irow = index[irn];                  
          auxr = derxy2[0][irn]*derxy2[1][irn];
          estif[irow  ][icol  ] -=  (derxy2[0][irn]*vderxy[0][0] +
                                     derxy2[2][irn]*vderxy[1][0] +
                                     auxr*vderxy[0][0])*aux*DENS;
          estif[irow+1][icol  ] -=  (derxy2[2][irn]*vderxy[0][0] +
                                     derxy2[1][irn]*vderxy[1][0] +
                                     auxr*vderxy[1][0])*aux*DENS;
          estif[irow  ][icol+1] -=  (derxy2[0][irn]*vderxy[0][1] +
                                     derxy2[2][irn]*vderxy[1][1] +
                                     auxr*vderxy[0][1])*aux*DENS;
          estif[irow+1][icol+1] -=  (derxy2[2][irn]*vderxy[0][1] +
                                     derxy2[1][irn]*vderxy[1][1] +
                                     auxr*vderxy[1][1])*aux*DENS;
        }
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of f2_calstabkvv */



/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvp

<pre>

In this routine the stabilisation part of matrix Kvv is calculated:

EULER:
    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /

EULER/ALE:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
 NOTE: there's only one elestif 				    
       --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]

NOTE: if the function is called from f2_calint, we have EULER case
      and so velint = u_old
      if the function is called from f2_calinta, we have ALE case
      and so velint = c (ale-convective velocity)
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\param   index	   INT  	   (i)	   index for local assembly
\param   DENS	   DOUBLE  	   (i)	   fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabkvp(
  ELEMENT         *ele,    
  DOUBLE         **estif, 
  DOUBLE          *velint,
  DOUBLE          *funct, 
  DOUBLE         **derxy, 
  DOUBLE         **derxy2,
  DOUBLE           fac,   
  DOUBLE           visc,  
  INT              iel,   
  INT              ihoel,
  INT             *index,
  DOUBLE           DENS
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
  INT        irow,icol,irn,posc;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     c;
  DOUBLE     aux;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkvp");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];

  c = fac * taumu;
/*------------------------------ calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
/*----------------------------------------------------------------------*
   Calculate advection stabilisation:
EULER:
    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /
 *----------------------------------------------------------------------*/    
    for (icol=0; icol<iel; icol++)
    {
      posc=2*iel+icol;
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*c;
        estif[irow  ][posc] += derxy[0][icol]*aux*DENS;
        estif[irow+1][posc] += derxy[1][icol]*aux*DENS;
      }
    }
  }

  /*-------------------------------- calculate viscous stabilisation part */
  if (gls->ivisc!=0 && ihoel!=0)
  {
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/    
    switch (gls->ivisc) /* choose stabilisation type --> sign */
    {
        case 1: /* GLS- */
          sign = ONE;
          break;
        case 2: /* GLS+ */
          sign = -ONE;
          break;
        default:
          dserror("viscous stabilisation parameter unknown: IVISC");
    }
    c = fac * taump * visc * sign;
    
    for (icol=0; icol<iel; icol++)
    {
      posc=2*iel+icol;
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow=index[irn];
        aux = derxy2[0][irn] + derxy2[1][irn];
        estif[irow  ][posc] -=  (derxy2[0][irn]*derxy[0][icol] +
                                 derxy2[2][irn]*derxy[1][icol] +
                                 aux*derxy[0][icol])*c;
        estif[irow+1][posc] -=  (derxy2[2][irn]*derxy[0][icol] +
                                 derxy2[1][irn]*derxy[1][icol] +
                                 aux*derxy[1][icol])*c;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkvp */



/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mvv

<pre>                                                            irhan 05/04

In this routine the stabilisation part of matrix Mvv is calculated:

EULER:
    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /

EULER/ALE:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /  
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]

NOTE: if the function is called from f2_calint, we have EULER case
      and so velint = u_old
      if the function is called from f2_calinta, we have ALE case
      and so velint = c (ale-convective velocity)
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\param   index	   INT  	   (i)	   index for local assembly
\param   DENS	   DOUBLE  	   (i)	   fluid density

\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabmvv(
  ELEMENT         *ele,     
  DOUBLE         **emass,  
  DOUBLE          *velint, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
/*----------------------------------------------------------------------*/
  INT        irow,icol,irn,icn;
  DOUBLE     taumu;
  DOUBLE     taump;
  DOUBLE     c,cc;
  DOUBLE     aux,auxc;
  DOUBLE     sign;

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/  
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabmvv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  gls  = ele->e.f2->stabi.gls;

  /* set stabilisation parameter */
  taumu = fdyn->tau[0];
  taump = fdyn->tau[1];
  
  c = fac * taumu;
  cc = c;

/*------------------------------ calculate advection stabilisation part */
  if (gls->iadvec!=0)
  {
/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
EULER:
    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
 *----------------------------------------------------------------------*/    
    for (icn=0; icn<TWO*iel; icn++)
    {
      auxc = funct[icn]*cc;
      icol = index[icn];
      for (irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
        emass[irow  ][icol  ] += aux*DENS*DENS;
        emass[irow+1][icol+1] += aux*DENS*DENS;
      }
    }
  }

/*-------------------------------- calculate viscous stabilisation part */
  if (gls->ivisc!=0 && ihoel!=0)
  {
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /
 *----------------------------------------------------------------------*/    
    switch (gls->ivisc) /* choose stabilisation type --> sign */
    {
        case 1: /* GLS- */
          sign = ONE;
          break;
        case 2: /* GLS+ */
          sign = -ONE;
          break;
        default:
          dserror("viscous stabilisation parameter unknown: IVISC");
    }
    c = fac * taump * visc * sign;
    
    for (icn=0; icn<TWO*iel; icn++)
    {      
      aux = funct[icn]*c;
      icol = index[icn];
      for(irn=0; irn<TWO*iel; irn++)
      {
        irow = index[irn];
        emass[irow  ][icol  ] -= (TWO*derxy2[0][irn] + derxy2[1][irn])*aux;
        emass[irow+1][icol  ] -=      derxy2[2][irn]*aux*DENS;;
        emass[irow+1][icol+1] -= (TWO*derxy2[1][irn] + derxy2[0][irn])*aux;
        emass[irow  ][icol+1] -=      derxy2[2][irn]*aux*DENS;;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
 return;
} /* end of xfem_f2_calstabmvv */



/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpv

<pre>                                                            irhan 05/04


In this routine the stabilisation part of matrix Kpv is calculated:

EULER/ALE:
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
  
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /  

EULER/ALE:
    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				    
      --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param  *ele       ELEMENT         (i)     actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *gridvint  DOUBLE          (i)     grid-vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\param   index	   INT  	   (i)	  index for local assembly
\param   DENS	   DOUBLE  	   (i)	  fluid density
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabkpv(
  ELEMENT         *ele,
  DOUBLE         **estif,   
  DOUBLE          *velint,
  DOUBLE          *gridvint, 
  DOUBLE         **vderxy, 
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
  DOUBLE         **derxy2, 
  DOUBLE           fac,    
  DOUBLE           visc,   
  INT              iel,    
  INT              ihoel,
  INT             *index,
  DOUBLE           DENS
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
  INT        irow,icol,icn,posr;
  DOUBLE     c;
  DOUBLE     aux;
  DOUBLE     taump;

#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkpv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump;
/*----------------------------------------------------------------------*
   Calculate stabilisation part Nc(u):
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
  /* evaluate for Newton- and fixed-point-like-iteration */
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*c;
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] -= derxy[0][irow]*aux*DENS;
        estif[posr][icol+1] -= derxy[1][irow]*aux*DENS;
      }
    }

  /*----------------------------------------------------------------------*
   Calculate stabilisation part Nr(u):
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
  if (fdyn->nir!=0) /* evaluate for Newton iteration */
  {
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = funct[icn]*c;
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] -= aux*(derxy[0][irow]*vderxy[0][0] +
                                    derxy[1][irow]*vderxy[1][0])*DENS;
        estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1] +
                                    derxy[1][irow]*vderxy[1][1])*DENS;
      }
    }
  }

/*----------------------------------------------------------------------*
   Calculate stabilisation part for higher order elements:
    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/
  if (ihoel!=0)
  {
    c = c * visc;
    
    for (icn=0; icn<TWO*iel; icn++)
    {
      icol = index[icn];
      aux = derxy2[0][icn] + derxy2[1][icn];
      for (irow=0; irow<iel; irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol  ] += ((derxy2[0][icn]+aux)*derxy[0][irow] +
                                 derxy2[2][icn]     *derxy[1][irow])*c;
        estif[posr][icol+1] +=  (derxy2[2][icn]     *derxy[0][irow] +
                                (derxy2[1][icn]+aux)*derxy[1][irow])*c;
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkpv */



/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mpv

<pre>                                                            irhan 05/04

In this routine the stabilisation part of matrix Mpv is calculated:

    /
   |  - tau_mp * grad(q) * u d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elemass  				    
      --> Mpv is stored in emass[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT		   (i)	   num. of nodes in ele
\param   index	   INT  	   (i)	  index for local assembly
\param   DENS	   DOUBLE  	   (i)	  fluid density

\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabmpv(
  DOUBLE         **emass,   
  DOUBLE          *funct,  
  DOUBLE         **derxy,  
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
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
  INT        irow,icol,icn,posr;
  DOUBLE     c;
  DOUBLE     taump;
  DOUBLE     auxc;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabmpv");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;
  
  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Mpv:
    /
   |  - tau_mp * grad(q) * u d_omega
  /
 *----------------------------------------------------------------------*/
  for (icn=0; icn<TWO*iel; icn++)
  {
    icol = index[icn];
    auxc = funct[icn]*c;
    for (irow=0; irow<iel; irow++)
    {
      posr = irow + 2*iel;
      emass[posr][icol  ] -= derxy[0][irow]*auxc*DENS;
      emass[posr][icol+1] -= derxy[1][irow]*auxc*DENS;
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabmpv */



/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpp

<pre>                                                            irhan 05/04

In this routine the stabilisation part of matrix Kpp is calculated:

    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Kpp is stored in				     
	      estif[((2*iel)..(3*iel-1)][((2*iel)..(3*iel-1)] 
      
</pre>
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT  	   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void xfem_f2_calstabkpp(
  DOUBLE         **estif,   
  DOUBLE         **derxy,  
  DOUBLE           fac,    
  INT              iel
  )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |
/*----------------------------------------------------------------------*/
  INT        irow,icol,posc,posr;
  DOUBLE     c;
  DOUBLE     taump;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calstabkpp");
#endif
/*----------------------------------------------------------------------*/

  /* initialise */
  fdyn = alldyn[genprob.numff].fdyn;

  /* set stabilisation parameter */
  taump = fdyn->tau[1];
  
  c = fac * taump;
/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Kpp:
    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0; icol<iel; icol++)
  {
    posc = icol + 2*iel;
    for (irow=0; irow<iel; irow++)
    {
      posr = irow + 2*iel;
      estif[posr][posc] -= (derxy[0][irow]*derxy[0][icol] +
                            derxy[1][irow]*derxy[1][icol])*c;
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_f2_calstabkpp */
/*! @} (documentation module close)*/	    
#endif
