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

#endif
