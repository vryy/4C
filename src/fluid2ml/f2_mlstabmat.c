/*!----------------------------------------------------------------------
\file
\brief evaluate stabilization part of submesh matrices for fluid2

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm stiffness matrix SMK for fluid2

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh stiffness matrix 
SMK is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabsmk(FLUID_DYN_CALC  *dynvar,
	           FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smestif,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol;
DOUBLE tau;
DOUBLE con,ccon,beta,divv,cb,cbb,ccb;
DOUBLE aux,auxc;
DOUBLE sign;

#ifdef DEBUG 
dstrc_enter("f2_calstabsmk");
#endif

/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;
con = fac*tau;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate convective stabilization part:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  tau * u_old * grad(w) * u_old * grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  auxc = (velint[0]*smderxy[0][icol] + velint[1]*smderxy[1][icol])*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
    smestif[irow][icol] += aux;
  } /* end of loop over irow */
} /* end of loop over icol */

if (dynvar->conte!=0)  
{
  cb  = con*beta;
  cbb = cb*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  beta * tau * u_old * grad(w) * bub * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
      smestif[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * w * div(u_old) * u_old * grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  cb = cb*sign;
  for (icol=0;icol<smiel;icol++)
  {
    aux = (velint[0]*smderxy[0][icol] + velint[1]*smderxy[1][icol])*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau * w * div(u_old) * bub * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    aux = smfunct[icol]*divv*divv*cbb;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */
} 

/*----------------------------------------------------------------------*
   Calculate convective stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con*visc;
/*----------------------------------------------------------------------*
    /
   |  - tau * u_old * grad(w) * nue * delta(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (smderxy2[0][icol] + smderxy2[1][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
      smestif[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
  
  if (dynvar->conte!=0)  
  {
    ccb  = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * nue * w * div(u_old) * delta(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      aux = (smderxy2[0][icol] + smderxy2[1][icol])*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
  	smestif[irow][icol] -= smfunct[irow]*aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  } 

/*----------------------------------------------------------------------*
   Calculate viscous stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/   
  ccon = con * visc * visc * sign;
/*----------------------------------------------------------------------*
    /
   |  +/- tau  * nue * delta(w) * nue * delta(bub) d_omega
  /
 *----------------------------------------------------------------------*/   
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (smderxy2[0][icol] + smderxy2[1][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smestif[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

/*----------------------------------------------------------------------*
    /
   |  -/+ tau  * nue * delta(w) * u_old * grad(bub) d_omega
  /
 *----------------------------------------------------------------------*/
  ccon = con * visc * sign;
  
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (velint[0]*smderxy[0][icol] + velint[1]*smderxy[1][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smestif[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

  if (dynvar->conte!=0)  
  {
    ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau * nue * delta(w) * bub * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      auxc = smfunct[icol]*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
    	aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
    	smestif[irow][icol] -= aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  }
} /* end of stabilization for higher-order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabsmk */

/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of submesh mass matrix SMM for fluid2

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh mass matrix SMM 
is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabsmm(FLUID_DYN_CALC  *dynvar,
	           FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smemass,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol;
DOUBLE tau;
DOUBLE con,ccon,beta,divv,cb;
DOUBLE aux,auxc;
DOUBLE sign;

#ifdef DEBUG 
dstrc_enter("f2_calstabsmm");
#endif

/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate temporal stabilization part:
 *----------------------------------------------------------------------*/
if (ABS(mlvar->smstado)<3)
{ 
  if (mlvar->smstado==-1) con = fac*tau;
  else                    con = fac*tau*sign; 
/*---------------------------------------------------------------------*
    /
   |  +/- tau * (1/(theta*dt)) * w * bub d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*con*(ONE/dynvar->thsl);
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end loop over irow */
  } /* end loop over icol */

/*---------------------------------------------------------------------*
    /
   |  +/- tau * w * u_old * grad(bub) d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (velint[0]*smderxy[0][icol] + velint[1]*smderxy[1][icol])*con;
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end loop over irow */
  } /* end loop over icol */

  if (dynvar->conte!=0)  
  {
    cb = con*beta;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * w * bub * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      auxc = smfunct[icol]*divv*cb;
      for (irow=0;irow<smiel;irow++)
      {
        aux = smfunct[irow]*auxc;
        smemass[irow][icol] += aux;
      } /* end loop over irow */
    } /* end loop over icol */
  }  

/*----------------------------------------------------------------------*
   Calculate temporal stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
  if (ihoelsm!=0)
  {
    ccon = con * visc;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * w * nue * delta(bub)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {	 
      aux = (smderxy2[0][icol] + smderxy2[1][icol])*ccon;
      for(irow=0;irow<smiel;irow++)
      {
        smemass[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
    } /* end loop over icol */
  } /* end of temporal stabilization for higher-order elements */
} /* end of temporal stabilization */   

/*----------------------------------------------------------------------*
   Calculate convective stabilization part:
 *----------------------------------------------------------------------*/
con = fac*tau;
/*----------------------------------------------------------------------*
    /
   |  tau * u_old * grad(w) * bub d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  auxc = smfunct[icol]*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
    smemass[irow][icol]	+= aux;
  } /* end loop over irow */
} /* end loop over icol */
    
if (dynvar->conte!=0)  
{
  cb = con*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * w * div(u_old) * bub  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
}  

/*----------------------------------------------------------------------*
   Calculate viscous stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con * visc * sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * nue * delta(w) * bub  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {	 
    auxc = smfunct[icol]*ccon;
    for(irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smemass[irow][icol] -= aux;
    } /* end loop over irow */
  } /* end loop over icol */
} /* end of viscous stabilization for higher-order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabsmm */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_c * div(v) * div(u)   d_omega
  /

    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
  
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /  

    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /

    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) d_omega
  /  

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /

  
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)	   weighting factor	   
\param   visc	   DOUBLE	   (i)	   fluid viscosity
\param   iel	   INT		   (i)	   num. of nodes in ele
\param   ihoel	   INT		   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabkvv(ELEMENT         *ele,    
		   FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **estif,  
		   DOUBLE	   *velint, 
		   DOUBLE	  **vderxy, 
		   DOUBLE	   *funct,  
		   DOUBLE	  **derxy,  
		   DOUBLE	  **derxy2, 
		   DOUBLE	    fac,    
		   DOUBLE	    visc,   
		   INT  	    iel,    
                   INT  	    ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |  
/*----------------------------------------------------------------------*/
INT    irow,icol,irn,icn,ird;
DOUBLE taumu;
DOUBLE taump;
DOUBLE tauc;
DOUBLE con,ccon,beta,divv,cb,cbb,ccb;
DOUBLE aux,auxr,auxc;
DOUBLE signvi,signre;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabkvv");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];
tauc  = dynvar->tau[2];

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

switch (ele->e.f2->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ONE;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f2->ivisc) */

/*----------------------------------------------------------------------*
   Calculate continuity stabilisation part:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
if (ele->e.f2->icont!=0)
{
   con = fac*tauc;

   icol=0;
   for(icn=0;icn<iel;icn++)
   {
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
         estif[irow][icol]     += derxy[0][irn]*derxy[0][icn]*con;
         estif[irow][icol+1]   += derxy[0][irn]*derxy[1][icn]*con;
         estif[irow+1][icol]   += derxy[1][irn]*derxy[0][icn]*con;
         estif[irow+1][icol+1] += derxy[1][irn]*derxy[1][icn]*con;
         irow += 2;
      } /* end of loop over irn */
      icol += 2;
   } /* end of loop over icn */
} /* endif (ele->e.f2->icont!=0) */

con = fac*taumu;

if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
   Calculate convection stabilisation part Nc(u):
 *----------------------------------------------------------------------*/
  if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*con;
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
        estif[irow][icol]     += aux;
        estif[irow+1][icol+1] += aux;
        irow += 2;
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */

    if (dynvar->conte!=0)  
    {
      cb  = con*beta;
      cbb = cb*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  beta * tau_mu * u_old * grad(v) * u * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        auxc = funct[icn]*divv*cb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
          estif[irow][icol]     += aux;
          estif[irow+1][icol+1] += aux;
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      cb = cb*signre;
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*divv*cb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          estif[irow][icol]     += aux*funct[irn];
          estif[irow+1][icol+1] += aux*funct[irn];
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau_mu * v * div(u_old) * u * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        aux = funct[icn]*divv*divv*cbb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          estif[irow][icol]     += aux*funct[irn];
          estif[irow+1][icol+1] += aux*funct[irn];
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
    }
  } /* endif (dynvar->nic!=0) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part Nr(u):
 *----------------------------------------------------------------------*/
  if (dynvar->nir!=0)  /* evaluate for Newton iteraton */
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      auxc = funct[icn]*con;
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
        estif[irow][icol]     += aux*vderxy[0][0];
        estif[irow][icol+1]   += aux*vderxy[0][1];
        estif[irow+1][icol]   += aux*vderxy[1][0];
        estif[irow+1][icol+1] += aux*vderxy[1][1];
        irow += 2;
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */

    if (dynvar->conte!=0)  
    {
      cb  = con*beta;
      cbb = cb*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  beta * tau_mu * u_old * grad(v) * u_old * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*cb;
          estif[irow][icol]     += aux*velint[0]*derxy[0][icn];
          estif[irow][icol+1]   += aux*velint[0]*derxy[1][icn];
          estif[irow+1][icol]   += aux*velint[1]*derxy[0][icn];
          estif[irow+1][icol+1] += aux*velint[1]*derxy[1][icn];
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      cb = cb*signre;
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        auxc = funct[icn]*divv*cb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = funct[irn]*auxc;
          estif[irow][icol]     += aux*vderxy[0][0];
          estif[irow][icol+1]   += aux*vderxy[0][1];
          estif[irow+1][icol]   += aux*vderxy[1][0];
          estif[irow+1][icol+1] += aux*vderxy[1][1];
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau_mu * v * div(u_old) * u_old * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = funct[irn]*divv*cbb;
          estif[irow][icol]     += aux*velint[0]*derxy[0][icn];
          estif[irow][icol+1]   += aux*velint[0]*derxy[1][icn];
          estif[irow+1][icol]   += aux*velint[1]*derxy[0][icn];
          estif[irow+1][icol+1] += aux*velint[1]*derxy[1][icn];
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
    }
  } /* endif (dynvar->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part for higher order elements:
 *----------------------------------------------------------------------*/
  if (ihoel!=0)
  {
    ccon = con*visc;
    
    if (dynvar->vite==0)
    {
/*----------------------------------------------------------------------*
    /
   |  - tau_mu * nue * u_old * grad(v) * delta(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        auxc = (derxy2[0][icn] + derxy2[1][icn])*ccon;
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
          estif[irow][icol]     -= aux;
          estif[irow+1][icol+1] -= aux;
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
      
      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * nue * v * div(u_old) * delta(u)   d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          auxc = (derxy2[0][icn] + derxy2[1][icn])*divv*ccb;
          irow=0;
          for (irn=0;irn<iel;irn++)
          {
            aux = funct[irn]*auxc;
            estif[irow][icol]     -= aux;
            estif[irow+1][icol+1] -= aux;
            irow += 2;
          } /* end of loop over irn */
          icol += 2;
        } /* end of loop over icn */
      }	
    }
    else
    { 
/*----------------------------------------------------------------------*
    /
   |  - tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        auxc = derxy2[0][icn] + derxy2[1][icn];
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*ccon;
          estif[irow][icol]     -= aux*(derxy2[0][icn] + auxc);
          estif[irow][icol+1]   -= aux* derxy2[2][icn];
          estif[irow+1][icol]   -= aux* derxy2[2][icn];
          estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
          irow += 2;
        } /* end of loop over irn */
        icol += 2;
      } /* end of loop over icn */
      
      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * 2 * nue * v * div(u_old) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          auxc = derxy2[0][icn] + derxy2[1][icn];
          irow=0;
          for (irn=0;irn<iel;irn++)
          {
            aux = funct[irn]*divv*ccb;
            estif[irow][icol]     -= aux*(derxy2[0][icn] + auxc);
            estif[irow][icol+1]   -= aux* derxy2[2][icn];
            estif[irow+1][icol]   -= aux* derxy2[2][icn];
            estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
            irow += 2;
          } /* end of loop over irn */
          icol += 2;
        } /* end of loop over icn */
      }	
    }      
  } /* end of stabilisation for higher order elements */
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part for higher order elements:
 *----------------------------------------------------------------------*/   
if (ihoel!=0 && ele->e.f2->ivisc>0)
{   
  ccon = fac * taump * visc * visc * signvi;
  
  if (dynvar->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  +/- tau_mp  * nue**2 * delta(v) * delta(u) d_omega
  /
 *----------------------------------------------------------------------*/   
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      auxc = (derxy2[0][icn] + derxy2[1][icn])*ccon;
      for (irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn])*auxc;
	estif[irow][icol]     += aux;
	estif[irow+1][icol+1] += aux;
	irow += 2;				     
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * nue * delta(v) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
    ccon = fac * taump * visc * signvi;
   
    if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
    {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  aux = (derxy2[0][irn] + derxy2[1][irn])*auxc;
	  estif[irow][icol]	-= aux;
	  estif[irow+1][icol+1] -= aux;
	  irow += 2;
	} /* end of loop over irn */
	icol += 2;
      } /* end of loop over icn */

      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * 2 * nue * div(eps(v)) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  auxc = funct[icn]*divv*ccb;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = (derxy2[0][irn] + derxy2[1][irn])*auxc;
	    estif[irow][icol]	  -= aux;
	    estif[irow+1][icol+1] -= aux;
	    irow += 2;
	  } /* end of loop over irn */
	  icol += 2;
        } /* end of loop over icn */
      }
    } /* endif (dynvar->nic!=0) */
   
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * nue * delta(v) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
    if (dynvar->nir!=0)  /* evaluate for Newton iteraton */
    {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = funct[icn]*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  aux = (derxy2[0][irn] + derxy2[1][irn])*auxc;
	  estif[irow][icol]	-= vderxy[0][0]*aux;
	  estif[irow][icol+1]	-= vderxy[0][1]*aux;
	  estif[irow+1][icol]	-= vderxy[1][0]*aux;
	  estif[irow+1][icol+1] -= vderxy[1][1]*aux;
	  irow += 2;
	} /* end of loop over irn */
	icol += 2;
      } /* end of loop over icn */

      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * nue * delta(v) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = (derxy2[0][irn] + derxy2[1][irn])*ccb;
	    estif[irow][icol]	  -= velint[0]*derxy[0][icn]*aux;
	    estif[irow][icol+1]	  -= velint[0]*derxy[1][icn]*aux;
	    estif[irow+1][icol]	  -= velint[1]*derxy[0][icn]*aux;
	    estif[irow+1][icol+1] -= velint[1]*derxy[1][icn]*aux;
	    irow += 2;
	  } /* end of loop over irn */
	  icol += 2;
        } /* end of loop over icn */
      }
    } /* endif (dynvar->nir!=0) */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/   
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      auxc = derxy2[0][icn] + derxy2[1][icn];
      for (irn=0;irn<iel;irn++)
      {
        auxr = derxy2[0][irn] + derxy2[1][irn];
	aux  = auxc*auxr;
	estif[irow][icol]     +=  (derxy2[0][icn]*(derxy2[0][irn]+auxr) \
				 + derxy2[2][icn]* derxy2[2][irn]	\
	        		 + derxy2[0][irn]* auxc + aux)*ccon	;
	estif[irow][icol+1]   +=  (derxy2[2][icn]*(derxy2[0][irn]+auxr) \
				 + derxy2[1][icn]* derxy2[2][irn]	\
	        		 + derxy2[2][irn]* auxc)*ccon		;
	estif[irow+1][icol]   +=  (derxy2[0][icn]* derxy2[2][irn]	\
				 + derxy2[2][icn]*(derxy2[1][irn]+auxr) \
	        		 + derxy2[2][irn]* auxc)*ccon		;  
	estif[irow+1][icol+1] +=  (derxy2[2][icn]* derxy2[2][irn]	\
				 + derxy2[1][icn]*(derxy2[1][irn]+auxr) \
	        		 + derxy2[1][irn]* auxc + aux)*ccon	;
	irow += 2;				     
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
    ccon = fac * taump * visc * signvi;
   
    if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
    {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  auxr = derxy2[0][irn] + derxy2[1][irn];
	  estif[irow][icol]	-= (derxy2[0][irn]+auxr)*auxc;
	  estif[irow][icol+1]	-=  derxy2[2][irn]*auxc;
	  estif[irow+1][icol]	-=  derxy2[2][irn]*auxc;
	  estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*auxc;
	  irow += 2;
	} /* end of loop over irn */
	icol += 2;
      } /* end of loop over icn */

      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * 2 * nue * div(eps(v)) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  auxc = funct[icn]*divv*ccb;
	  for (irn=0;irn<iel;irn++)
	  {
	    auxr = derxy2[0][irn] + derxy2[1][irn];
	    estif[irow][icol]	  -= (derxy2[0][irn]+auxr)*auxc;
	    estif[irow][icol+1]   -=  derxy2[2][irn]*auxc;
	    estif[irow+1][icol]   -=  derxy2[2][irn]*auxc;
	    estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*auxc;
	    irow += 2;
	  } /* end of loop over irn */
	  icol += 2;
        } /* end of loop over icn */
      }
    } /* endif (dynvar->nic!=0) */
   
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
    if (dynvar->nir!=0)  /* evaluate for Newton iteraton */
    {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = funct[icn]*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  auxr = derxy2[0][irn] + derxy2[1][irn];
	  estif[irow][icol]	-=  (derxy2[0][irn]*vderxy[0][0]   \
	        		   + derxy2[2][irn]*vderxy[1][0]   \
	  			   + auxr*vderxy[0][0])*auxc	   ;
	  estif[irow][icol+1]	-=  (derxy2[0][irn]*vderxy[0][1]   \
	        		   + derxy2[2][irn]*vderxy[1][1]   \
	  			   + auxr*vderxy[0][1])*auxc	   ;
	  estif[irow+1][icol]	-=  (derxy2[2][irn]*vderxy[0][0]   \
	        		   + derxy2[1][irn]*vderxy[1][0]   \
	  			   + auxr*vderxy[1][0])*auxc	   ;
	  estif[irow+1][icol+1] -=  (derxy2[2][irn]*vderxy[0][1]   \
	        		   + derxy2[1][irn]*vderxy[1][1]   \
	  			   + auxr*vderxy[1][1])*auxc	   ;
	  irow += 2;
	} /* end of loop over irn */
	icol += 2;
      } /* end of loop over icn */

      if (dynvar->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * 2 * nue * div(eps(v)) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = derxy2[0][irn] + derxy2[1][irn];
	    estif[irow][icol]	  -=  ((derxy2[0][irn]+aux)*velint[0] 
	        	  	      + derxy2[2][irn]*velint[1])      
				      *derxy[0][icn]*ccb              ;
	    estif[irow][icol+1]	  -=  ((derxy2[0][irn]+aux)*velint[0] 
	        	  	      + derxy2[2][irn]*velint[1])      
				      *derxy[1][icn]*ccb              ;
	    estif[irow+1][icol]	  -=  ((derxy2[1][irn]+aux)*velint[1] 
	        	  	      + derxy2[2][irn]*velint[0])      
				      *derxy[0][icn]*ccb              ;
	    estif[irow+1][icol+1] -=  ((derxy2[1][irn]+aux)*velint[1] 
	        	  	      + derxy2[2][irn]*velint[0])      
				      *derxy[1][icn]*ccb              ;
	    irow += 2;
	  } /* end of loop over irn */
	  icol += 2;
        } /* end of loop over icn */
      }
    } /* endif (dynvar->nir!=0) */
  }
} /* end of viscous stabilisation for higher order elments */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabkvv */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvp

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /

    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
 NOTE: there's only one elestif 				    
       --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabkvp(ELEMENT         *ele,    
		   FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **estif, 
		   DOUBLE	   *velint,
		   DOUBLE	  **vderxy,
		   DOUBLE	   *funct, 
		   DOUBLE	  **derxy, 
		   DOUBLE	  **derxy2,
		   DOUBLE	    fac,   
		   DOUBLE	    visc,  
		   INT  	    iel,   
		   INT  	    ihoel)
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
INT    irow,icol,irn,ird,posc;
DOUBLE taumu;
DOUBLE taump;
DOUBLE con,beta,divv,cb;
DOUBLE aux;
DOUBLE signvi,signre;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabkvp");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

con = fac * taumu;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

switch (ele->e.f2->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ONE;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f2->ivisc) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation:
 *----------------------------------------------------------------------*/
if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    irow=0;
    posc=2*iel+icol;
    for (irn=0;irn<iel;irn++)
    {
      aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*con;
      estif[irow][posc]   += derxy[0][icol]*aux;
      estif[irow+1][posc] += derxy[1][icol]*aux;
      irow += 2;
    } /* end loop over irn */
  } /* end loop over icol */
    
  if (dynvar->conte!=0)  
  {
    cb = con*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc=2*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
    	aux = funct[irn]*divv*cb;
        estif[irow][posc]   += derxy[0][icol]*aux;
        estif[irow+1][posc] += derxy[1][icol]*aux;
    	irow += 2;
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */
  }  
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (ele->e.f2->ivisc>0 && ihoel!=0)
{
  con = fac * taump * visc * signvi;
  
  if (dynvar->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * nue * delta(v) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc=2*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn])*con;
	estif[irow][posc]   -= aux*derxy[0][icol];
	estif[irow+1][posc] -= aux*derxy[1][icol];
	irow += 2;			       
      } /* end loop over irn */
    } /* end loop over icol */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc=2*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
        aux = derxy2[0][irn] + derxy2[1][irn];
	estif[irow][posc]   -=  (derxy2[0][irn]*derxy[0][icol] \
			       + derxy2[2][irn]*derxy[1][icol] \
	        	       + aux*derxy[0][icol])*con       ;
	estif[irow+1][posc] -=  (derxy2[2][irn]*derxy[0][icol] \
			       + derxy2[1][irn]*derxy[1][icol] \
	        	       + aux*derxy[1][icol])*con       ;
	irow += 2;			       
      } /* end loop over irn */
    } /* end loop over icol */
  }
} /*------------ end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabkvp */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mvv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Mvv is calculated:

    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
  
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /  
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabmvv(ELEMENT         *ele,     
		   FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **emass,  
		   DOUBLE	   *velint, 
		   DOUBLE	  **vderxy,
    		   DOUBLE	   *funct,  
		   DOUBLE	  **derxy,  
		   DOUBLE	  **derxy2, 
		   DOUBLE	    fac,    
		   DOUBLE	    visc,   
		   INT  	    iel,    
		   INT  	    ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
/*----------------------------------------------------------------------*/
INT    irow,icol,irn,icn,ird;
DOUBLE taumu;
DOUBLE taump;
DOUBLE con,beta,divv,cb;
DOUBLE aux,auxc;
DOUBLE signvi,signre;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabmvv");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

con = fac * taumu;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

switch (ele->e.f2->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ONE;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f2->ivisc) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
 *----------------------------------------------------------------------*/
if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    auxc = funct[icn]*con;
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
      emass[irow][icol]     += aux;
      emass[irow+1][icol+1] += aux;
      irow += 2;
    } /* end loop over irn */
    icol += 2;      
  } /* end loop over icn */
    
  if (dynvar->conte!=0)  
  {
    cb = con*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      auxc = funct[icn]*divv*cb;
      for (irn=0;irn<iel;irn++)
      {
    	aux = funct[irn]*auxc;
        emass[irow][icol]     += aux;
        emass[irow+1][icol+1] += aux;
    	irow += 2;
      } /* end of loop over irn */
      icol += 2;
    } /* end of loop over icn */
  }  
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (ele->e.f2->ivisc>0 && ihoel!=0)
{
  con = fac * taump * visc * signvi;
  
  if (dynvar->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * nue * delta(v) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {      
      auxc = funct[icn]*con;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn])*auxc;
	emass[irow][icol]     -= aux;
	emass[irow+1][icol+1] -= aux;
	irow += 2;
      } /* end loop over irn */
      icol += 2;
    } /* end loop over icn */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {      
      aux = funct[icn]*con;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
        emass[irow][icol]     -= (TWO*derxy2[0][irn] + derxy2[1][irn])*aux;
	emass[irow][icol+1]   -=      derxy2[2][irn]*aux;
	emass[irow+1][icol]   -=      derxy2[2][irn]*aux;
	emass[irow+1][icol+1] -= (TWO*derxy2[1][irn] + derxy2[0][irn])*aux;
	irow += 2;
      } /* end loop over irn */
      icol += 2;
    } /* end loop over icn */
  }  
} /* end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabmvv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kpv is calculated:

    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
  
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /  

    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				    
      --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)      
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabkpv(FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **estif,   
		   DOUBLE	   *velint, 
		   DOUBLE	  **vderxy, 
		   DOUBLE	   *funct,  
		   DOUBLE	  **derxy,  
		   DOUBLE	  **derxy2, 
		   DOUBLE	    fac,    
		   DOUBLE	    visc,   
		   INT  	    iel,    
		   INT  	    ihoel)
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
INT    irow,icol,irn,ird,icn,posr;
DOUBLE con,beta,divv,cb;
DOUBLE aux;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabkpv");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

con = fac * taump;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nc(u):
 *----------------------------------------------------------------------*/
if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
{
/*----------------------------------------------------------------------*
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*con;
    for (irow=0;irow<iel;irow++)
    {
      posr = irow + 2*iel;
      estif[posr][icol]   -= derxy[0][irow]*aux;
      estif[posr][icol+1] -= derxy[1][irow]*aux;
    } /* end loop over irow */
    icol += 2;
  } /* end loop over icn */
  
  if (dynvar->conte!=0)  
  {
    cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau_mp * grad(q) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = funct[icn]*divv*cb;
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol]   -= derxy[0][irow]*aux;
        estif[posr][icol+1] -= derxy[1][irow]*aux;
      } /* end loop over irow */
      icol += 2;
    } /* end loop over icn */
  }
} /* endif (dynvar->nic!=0) */ 

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nr(u):
 *----------------------------------------------------------------------*/
if (dynvar->nir!=0) /* evaluate for Newton iteration */
{
/*----------------------------------------------------------------------*
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    aux = funct[icn]*con;
    for (irow=0;irow<iel;irow++)
    {
      posr = irow + 2*iel;
      estif[posr][icol]   -= aux*(derxy[0][irow]*vderxy[0][0]	\
        			+ derxy[1][irow]*vderxy[1][0])  ;
      estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1]	\
        			+ derxy[1][irow]*vderxy[1][1])  ;
    } /* end loop over irow */
    icol += 2;
  } /* end loop over icn */
  
  if (dynvar->conte!=0)  
  {
    cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau_mp * grad(q) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 2*iel;
        aux = (velint[0]*derxy[0][irow] + velint[1]*derxy[1][irow])*cb;
	estif[posr][icol]   -= aux*derxy[0][icn];
	estif[posr][icol+1] -= aux*derxy[1][icn];
      } /* end loop over irow */
      icol += 2;
    } /* end loop over icn */
  }
} /* endif (dynvar->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
  con = con * visc;
  
  if (dynvar->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mp * nue *grad(q) * delta(u) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = (derxy2[0][icn] + derxy2[1][icn])*con;
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol]   += aux*derxy[0][irow];
        estif[posr][icol+1] += aux*derxy[1][irow];
      } /* end loop over irow */
      icol += 2;
    } /* end loop over icn */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = derxy2[0][icn] + derxy2[1][icn];
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 2*iel;
        estif[posr][icol]   += ((derxy2[0][icn]+aux)*derxy[0][irow]   \
			       + derxy2[2][icn]     *derxy[1][irow])*con;
        estif[posr][icol+1] +=  (derxy2[2][icn]     *derxy[0][irow]   \
			       +(derxy2[1][icn]+aux)*derxy[1][irow])*con;
      } /* end loop over irow */
      icol += 2;
    } /* end loop over icn */
  }
} /* end of stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabkpv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpp

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kpp is calculated:

    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Kpp is stored in				     
	      estif[((2*iel)..(3*iel-1)][((2*iel)..(3*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT  	   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabkpp(FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **estif,   
		   DOUBLE	  **derxy,  
		   DOUBLE	    fac,    
		   INT  	    iel)
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
INT    irow,icol,posc,posr;
DOUBLE con;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabkpp");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

con = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Kpp:
    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
   posc = icol + 2*iel;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 2*iel;
      estif[posr][posc] -= (derxy[0][irow]*derxy[0][icol]      \
                           +derxy[1][irow]*derxy[1][icol])*con ;
   } /* end of loop over irow */
} /* end of loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabkpp */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mpv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Mpv is calculated:

    /
   |  - tau_mp * grad(q) * u d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elemass  				    
      --> Mpv is stored in emass[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT		   (i)	   num. of nodes in ele

\return void                                                                       

------------------------------------------------------------------------*/
void f2_lscalstabmpv(FLUID_DYN_CALC  *dynvar,
		   DOUBLE	  **emass,   
		   DOUBLE	   *funct,  
		   DOUBLE	  **derxy,  
		   DOUBLE	    fac,    
		   INT  	    iel)
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
INT    irow,icol,icn,posr;
DOUBLE con;
DOUBLE taump;
DOUBLE auxc;

#ifdef DEBUG 
dstrc_enter("f2_lscalstabmpv");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

con = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Mpv:
    /
   |  - tau_mp * grad(q) * u d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   auxc = funct[icn]*con;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 2*iel;
      emass[posr][icol]   -= derxy[0][irow]*auxc;
      emass[posr][icol+1] -= derxy[1][irow]*auxc;
   } /* end loop over irow */
   icol +=2;
} /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_lscalstabmpv */



#endif
