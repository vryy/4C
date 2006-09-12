/*!----------------------------------------------------------------------
\file
\brief evaluate "Time" force vectors for submesh element for fluid2

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
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
\brief evaluate Galerkin part of submesh "Time" force vector for fluid2

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh "Time" force vector
on the rhs is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc 	   DOUBLE	   (i)    physical viscosity
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void

------------------------------------------------------------------------*/
void f2_calsmft(FLUID_DYN_ML	*mlvar,
		DOUBLE         **smetfor,
		DOUBLE  	*velintn,
		DOUBLE  	*velintnt,
		DOUBLE  	*velintnc,
		DOUBLE         **vderxyn,
		DOUBLE         **vderxync,
		DOUBLE         **vderxynv,
		DOUBLE         **vderxy2n,
		DOUBLE          *pderxyn,
		DOUBLE  	*smfunct,
		DOUBLE         **smderxy,
		DOUBLE  	 fac,
		DOUBLE  	 visc,
		INT		 smiel,
		INT		 iel,
		INT		 ihoelsm,
                INT		 ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT    irow,icol,icn;
DOUBLE con,ccon,beta,divvc,cb,facsr,facpr;
DOUBLE aux;

#ifdef DEBUG
dstrc_enter("f2_calsmft");
#endif

fdyn = alldyn[genprob.numff].fdyn;

facsr = fac*fdyn->thsr;
facpr = fac*fdyn->thpr;
con   = facsr*visc;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0:
  beta = ZERO;
break;
case 1:
  beta = ONE;
  divvc= vderxync[0][0]+vderxync[1][1];
break;
case 2:
  beta = ONE/TWO;
  divvc= vderxync[0][0]+vderxync[1][1];
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

/*----------------------------------------------------------------------*
   Calculate temporal forces (u_transient = u_n or ls_u_n or sm_u_n):
 *----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*
    /
   |  w * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
icol=mlvar->nelbub-2;
for (icn=0;icn<2;icn++)
{
  aux = velintnt[icn]*fac;
  for (irow=0;irow<smiel;irow++)
  {
    smetfor[irow][icol] += smfunct[irow]*aux;
  } /* end loop over irow */
  icol++;
} /* end loop over icol */

if (fdyn->thsr!=ZERO)
{
/*----------------------------------------------------------------------*
   Calculate convective forces (u_convective = u_n or ls_u_n):
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * w * u_convective * grad(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    aux = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1])*facsr;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (fdyn->conte!=0)
  {
    cb = beta*divvc*facsr;
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * beta * v * u_n * div(u_convective)  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = velintn[icn]*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

/*----------------------------------------------------------------------*
   Calculate viscous forces (separate large- and small-scale part of u_n):
 *----------------------------------------------------------------------*/
  if (mlvar->quastabub==0)
  {
    ccon = facsr*(visc + mlvar->smsgvisc);
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * grad(w) * (nue+nueT) * grad(sm_u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= (smderxy[0][irow]*vderxynv[0][icn]\
                               +smderxy[1][irow]*vderxynv[1][icn])*ccon;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

  if (ihoel!=0 && ihoelsm!=0)
  {
/*----------------------------------------------------------------------*
    /
   |  ((1-theta)*dt) * w * nue * delta(ls_u_n)  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = (vderxy2n[icn][0]+vderxy2n[icn][1])*con;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] += smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }
}

/*----------------------------------------------------------------------*
   Calculate pressure forces:
 *----------------------------------------------------------------------*/
if (fdyn->thpr!=ZERO && fdyn->iprerhs>0)
{
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * w * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    aux = pderxyn[icn]*facpr;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calsmft */

/*!---------------------------------------------------------------------
\brief evaluate stabilization part of sm "Time" force vector for fluid2

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh "Time" force
vector on the rhs is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc 	   DOUBLE	   (i)    physical viscosity
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void

------------------------------------------------------------------------*/
void f2_calstabsmft(FLUID_DYN_ML    *mlvar,
		    DOUBLE         **smetfor,
		    DOUBLE  	    *velintn,
		    DOUBLE  	    *velintnt,
		    DOUBLE  	    *velintnc,
		    DOUBLE         **vderxyn,
		    DOUBLE         **vderxync,
		    DOUBLE         **vderxynv,
		    DOUBLE         **vderxy2n,
		    DOUBLE          *pderxyn,
	 	    DOUBLE  	    *smfunct,
		    DOUBLE         **smderxy,
		    DOUBLE         **smderxy2,
		    DOUBLE  	     fac,
		    DOUBLE  	     visc,
		    INT		     smiel,
		    INT		     iel,
                    INT		     ihoelsm,
                    INT		     ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT    irow,icol,icn;
DOUBLE con,ccon,beta,divvc,cb,cbb;
DOUBLE aux,auxc;
DOUBLE sign,tau;

#ifdef DEBUG
dstrc_enter("f2_calstabsmft");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0:
  beta = ZERO;
break;
case 1:
  beta = ONE;
  divvc= vderxync[0][0]+vderxync[1][1];
break;
case 2:
  beta = ONE/TWO;
  divvc= vderxync[0][0]+vderxync[1][1];
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate temporal stabilization:
 *----------------------------------------------------------------------*/
if (ABS(mlvar->smstado)<3)
{
  if (mlvar->smstado==-1) con = fac*tau;
  else  		  con = fac*tau*sign;
/*---------------------------------------------------------------------*
    /
   |  +/-  tau * (1/(theta*dt)) * w * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    aux = velintnt[icn]*con*(ONE/fdyn->thsl);
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] += smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  ccon = con*(fdyn->thsr/fdyn->thsl);
/*---------------------------------------------------------------------*
    /
   |  -/+  tau * ((1-theta)*dt/(theta*dt)) * w * u_convective * grad(u_n)
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    aux = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (fdyn->conte!=0)
  {
    cb = beta*divvc*ccon;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt/(theta*dt)) * beta * w * u_n * div(u_conv.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = velintn[icn]*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

  if (ihoel!=0 && ihoelsm!=0)
  {
    ccon = ccon*visc;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt/(theta*dt)) * w * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = (vderxy2n[icn][0] + vderxy2n[icn][1])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] += smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

  if (fdyn->thpr!=ZERO && fdyn->iprerhs>0)
  {
    ccon = con*(fdyn->thpr/fdyn->thpl);
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt/(theta*dt)) * w * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = pderxyn[icn]*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    }
  }
}

/*----------------------------------------------------------------------*
   Calculate convective stabilization:
 *----------------------------------------------------------------------*/
con = fac*tau;
/*----------------------------------------------------------------------*
    /
   |  tau * u_convective * grad(w) * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
icol=mlvar->nelbub-2;
for (icn=0;icn<2;icn++)
{
  auxc = velintnt[icn]*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow])*auxc;
    smetfor[irow][icol] += aux;
  } /* end loop over irow */
  icol++;
} /* end loop over icol */

if (fdyn->conte!=0)
{
  cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * beta * w * div(u_convective) * u_transient  d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    aux = velintnt[icn]*divvc*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] += smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
}

if (fdyn->thsr!=ZERO)
{
  ccon = con*fdyn->thsr;
/*----------------------------------------------------------------------*
    /
   |  - tau * ((1-theta)*dt) * u_con. * grad(w) * u_con. * grad(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    auxc = (velintnc[0]*vderxyn[icn][0] + velintnc[1]*vderxyn[icn][1])*ccon;
     for (irow=0;irow<smiel;irow++)
    {
      aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (fdyn->conte!=0)
  {
    cb  = ccon*beta;
    cbb = cb*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau * ((1-th)*dt) * u_con. * grad(w) * u_n * div(u_con.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      auxc = velintn[icn]*divvc*cb;
       for (irow=0;irow<smiel;irow++)
      {
        aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow])*auxc;
        smetfor[irow][icol] -= aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */

    cb = cb*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau * ((1-th)*dt) * w * div(u_con.) * u_con. * grad(u_n)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = (velintnc[0]*vderxyn[icn][0] + velintnc[1]*vderxyn[icn][1])\
            *divvc*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= aux*smfunct[irow];
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */

/*----------------------------------------------------------------------*
    /
   |  -/+ beta * beta * tau * ((1-th)*dt) * v * div(u_con.) * u_n * div(u_con.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = velintn[icn]*divvc*divvc*cbb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= aux*smfunct[irow];
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

  if (ihoel!=0 && ihoelsm!=0)
  {
    ccon = con*fdyn->thsr*visc;
/*----------------------------------------------------------------------*
    /
   |  tau * ((1-theta)*dt) * u_con. * grad(w) * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      auxc = (vderxy2n[icn][0]+vderxy2n[icn][1])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow])*auxc;
        smetfor[irow][icol] += aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */

    if (fdyn->conte!=0)
    {
      cb = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * ((1-th)*dt) * w * div(u_con.) * nue * delta(u_n)
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-2;
      for (icn=0;icn<2;icn++)
      {
        aux = (vderxy2n[icn][0]+vderxy2n[icn][1])*divvc*cb;
        for (irow=0;irow<smiel;irow++)
        {
          smetfor[irow][icol] += aux*smfunct[irow];
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }
  }
}

if (fdyn->thpr!=ZERO && fdyn->iprerhs>0)
{
  ccon = con*fdyn->thpr;
/*----------------------------------------------------------------------*
    /
   |  - tau * ((1-theta)*dt) * u_convective * grad(w) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    auxc = pderxyn[icn]*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  }

  if (fdyn->conte!=0)
  {
    cb = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt) * beta * w * div(u_old) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      aux = pderxyn[icn]*divvc*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    }
  }
}

/*----------------------------------------------------------------------*
   Calculate viscous stabilization for higher-order elements:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * nue * delta(w) * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-2;
  for (icn=0;icn<2;icn++)
  {
    auxc = velintnt[icn]*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (fdyn->thsr!=ZERO)
  {
    ccon = con*fdyn->thsr*visc*sign;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * nue * delta(w) * u_convective * grad(u_n)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      auxc = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
        smetfor[irow][icol] += aux;
      }  /* end loop over irow */
      icol++;
    } /* end loop over icol */

    if (fdyn->conte!=0)
    {
      cb = beta*ccon;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * beta * nue * delta(w) *  * u_n * div(u_conv.)
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-2;
      for (icn=0;icn<2;icn++)
      {
        auxc = velintn[icn]*divvc*cb;
        for (irow=0;irow<smiel;irow++)
        {
         aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
         smetfor[irow][icol] += aux;
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }

    if (ihoel!=0)
    {
      ccon = ccon*visc;
/*---------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt) * nue * delta(w) * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-2;
      for (icn=0;icn<2;icn++)
      {
        auxc = (vderxy2n[icn][0] + vderxy2n[icn][1])*ccon;
        for (irow=0;irow<smiel;irow++)
        {
          aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
          smetfor[irow][icol] -= aux;
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }
  }

  if (fdyn->thpr!=ZERO && fdyn->iprerhs>0)
  {
    ccon = con*fdyn->thpr*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * nue * delta(w) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-2;
    for (icn=0;icn<2;icn++)
    {
      auxc = pderxyn[icn]*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
        smetfor[irow][icol] += aux;
      } /* end loop over irow */
      icol++;
    }
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabsmft */

/*!---------------------------------------------------------------------
\brief galerkin part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

                   /
   + (1-THETA)*dt |  v * b_old   d_omega
                 /

               /
   + THETA*dt |  v * b   d_omega
             /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'


</pre>
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param	 *funct       DOUBLE	      (i)    nat. shape functions
\param   *edeadn      DOUBLE          (i)    ele dead load at n
\param   *edeadng     DOUBLE          (i)    ele dead load at n+1
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_lscalgalexfv(DOUBLE          *eforce,
		     DOUBLE          *funct,
                     DOUBLE          *edeadn,
		     DOUBLE          *edeadng,
		     DOUBLE           fac,
		     INT              iel)
{
DOUBLE  facsl, facsr;
INT     inode,irow,isd;

fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------- set some factors */
facsl = fac*fdyn->thsl;
facsr = fac*fdyn->thsr;

/*----------------------------------------------------------------------*
   Calculate galerkin part of external forces:

                   /
   + (1-THETA)*dt |  v * b_old   d_omega
                 /

               /
   + THETA*dt |  v * b   d_omega
             /

 *----------------------------------------------------------------------*/
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*(edeadn[isd]*facsr+edeadng[isd]*facsl);
   } /* end of loop over isd */
} /* end of loop over inode */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalgalexfv */

/*!---------------------------------------------------------------------
\brief galerkin part of iteration forces for vel dofs

<pre>                                                         genk 04/02

In this routine the galerkin part of the iteration forces for vel dofs
is calculated:

                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
                 /


see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

</pre>
\param  *eforce    DOUBLE	   (i/o)  element force vector
\param  *covint    DOUBLE	   (i)	  conv. vels at int. point
\param  *velint    DOUBLE	   (i)    vel. at integr. point
\param **vderxy    DOUBLE	   (i)    vel. der. at integr. point
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac 	   DOUBLE	   (i)    weighting factor
\param	 iel	   INT		   (i)	  num. of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f2_lscalgalifv(DOUBLE          *eforce,
		  DOUBLE          *covint,
		  DOUBLE          *velint,
		  DOUBLE         **vderxy,
		  DOUBLE          *funct,
		  DOUBLE           fac,
		  INT              iel)
{
INT    inode,isd;
INT    irow;
DOUBLE facsl,betsl,beta,divv;

#ifdef DEBUG
dstrc_enter("f2_lscalgalifv");
#endif


fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------- set some factors */
facsl = fac * fdyn->thsl * fdyn->sigma;

switch (fdyn->conte) /* determine parameter beta of convective term */
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
} /* end switch (fdyn->conte) */

/*----------------------------------------------------------------------*
   Calculate convective forces of iteration force vector:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
                   /
   (+/-) THETA*dt |  v * u_old * grad(u_old)  d_omega
    |            /
    |-> signs due to nonlin. iteration scheme (fdyn->sigma)
 *----------------------------------------------------------------------*/
irow = -1;
for (inode=0;inode<iel;inode++)
{
  for (isd=0;isd<2;isd++)
  {
    irow++;
    eforce[irow] += funct[inode]*covint[isd]*facsl;
  } /* end loop over isd */
} /* end loop over inode */

if (fdyn->conte!=0)
{
  betsl = facsl * beta;
/*----------------------------------------------------------------------*
                   /
   (+/-) THETA*dt |  v * u_old * div(u_old)  d_omega
    |            /
    |-> signs due to nonlin. iteration scheme (fdyn->sigma)
 *----------------------------------------------------------------------*/
  irow = -1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      eforce[irow] += funct[inode]*velint[isd]*divv*betsl;
    } /* end loop over isd */
  } /* end loop over inode */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalgalifv */

/*!---------------------------------------------------------------------
\brief galerkin part of time forces for vel dofs

<pre>                                                         genk 04/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

      /
   + |  v * u     d_omega
    /

                      /
   (-) (1-THETA)*dt  |  v * u * grad(u)     d_omega
                    /

                      /
   (-) (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                    /

                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
   in ONESTEP methods: velint  = vel2int = U(n)
   in TWOSTEP methods: velint  = U(n+gamma)
   in TWOSTEP methods: vel2int = U(n)

</pre>
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param   *vel2int     DOUBLE	      (i)    vel. at integr. point
\param   *covint      DOUBLE	      (i)    conv. vel. at integr. p.
\param	 *funct       DOUBLE	      (i)    nat. shape functions
\param	**derxy       DOUBLE	      (i)    global derivatives
\param	**vderxy      DOUBLE	      (i/    global vel. deriv.
\param	  preint      DOUBLE	      (i)    pres. at integr. point
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_lscalgaltfv(DOUBLE          *eforce,
		  DOUBLE          *velint,
		  DOUBLE          *vel2int,
		  DOUBLE          *covint,
		  DOUBLE          *funct,
		  DOUBLE         **derxy,
		  DOUBLE         **vderxy,
		  DOUBLE           preint,
		  DOUBLE           visc,
		  DOUBLE           fac,
		  INT              iel)
{
INT    j,irow,isd,inode;
DOUBLE con,beta,divv,betsr;
DOUBLE aux;
DOUBLE facsr;
DOUBLE facpr;
DOUBLE fact[2];

#ifdef DEBUG
dstrc_enter("f2_lscalgaltfv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------- set some factors */
facsr = fac * fdyn->thsr;
facpr = fac * fdyn->thpr;
con   = facsr * visc;

switch (fdyn->conte) /* determine parameter beta of convective term */
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
} /* end switch (fdyn->conte) */

/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
   + |  v * u     d_omega
    /
 *----------------------------------------------------------------------*/
fact[0] = vel2int[0]*fac;
fact[1] = vel2int[1]*fac;
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*fact[isd];
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate convective forces of time force vector:
                    /
   - (1-THETA)*dt  |  v * u * grad(u)     d_omega
                  /
 *----------------------------------------------------------------------*/
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] -= funct[inode]*covint[isd]*facsr;
   } /* end of loop over isd */
} /* end of loop over irwo */

if (fdyn->conte!=0)
{
  betsr = facsr * beta;
/*----------------------------------------------------------------------*
                   /
   - (1-THETA)*dt |  v * u * div(u)  d_omega
                 /
 *----------------------------------------------------------------------*/
  irow = -1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      eforce[irow] += funct[inode]*velint[isd]*divv*betsr;
    } /* end loop over isd */
  } /* end loop over inode */
}

/*----------------------------------------------------------------------*
   Calculate viscous forces of time force vector:
 *----------------------------------------------------------------------*/
if (fdyn->vite==0)
{
/*----------------------------------------------------------------------*
                    /
   - (1-THETA)*dt  |  nue * grad(v) : grad(u)  d_omega
                  /
 *----------------------------------------------------------------------*/
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      eforce[irow] -= (derxy[0][inode]*vderxy[0][isd] \
                     + derxy[1][inode]*vderxy[1][isd])*con;
    } /* end of loop over isd */
  } /* end of loop over inode */
}
else
{
/*----------------------------------------------------------------------*
                    /
   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                  /
 *----------------------------------------------------------------------*/
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      for (j=0;j<2;j++)
      {
        eforce[irow] -= derxy[j][inode]*(vderxy[isd][j]+vderxy[j][isd])*con;
      } /* end of loop over j */
    } /* end of loop over isd */
  } /* end of loop over inode */
}

/*----------------------------------------------------------------------*
   Calculate pressure forces of time force vector:
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /
 *----------------------------------------------------------------------*/
if (fdyn->iprerhs>0)
{
   aux = preint * facpr;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      for (isd=0;isd<2;isd++)
      {
         irow++;
         eforce[irow] += derxy[isd][inode]*aux;
      } /* end of loop over isd */
   } /* end of loop over inode */
} /*endif (fdyn->iprerhs>0) */

/*----------------------------------------------------------------------*
   Calculate external forces of time force vector
                         (due to surface tension):
                    /
   + (1-THETA)*dt  |  v * h  d_gamma
                  /
 *----------------------------------------------------------------------*/
/* NOT IMPLEMENTED YET!!!!!!!!!!!!!!!
   maybe its better to do the integration over the element edge in
   a seperate loop (different gauss-points!!!!)
*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalgaltfv */

/*!---------------------------------------------------------------------
\brief galerkin part of time forces for pre dofs

<pre>                                                         genk 04/02

In this routine the galerkin part of the time forces for pre dofs
is calculated:

                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
    there's only one full element force vector
    for pre-dofs the pointer eforce points to the entry
    eforce[2*iel]

</pre>
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *funct       DOUBLE	      (i)    nat. shape functions
\param  **vderxy      DOUBLE	      (i)    global vel. deriv.
\param    fac         DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_lscalgaltfp(DOUBLE	  *eforce,
		  DOUBLE	  *funct,
		  DOUBLE	 **vderxy,
		  DOUBLE	   fac,
		  INT		   iel)
{
INT      inode;
DOUBLE   aux;
DOUBLE   facsr;

#ifdef DEBUG
dstrc_enter("f2_lscalgaltfp");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------- set some factors */
facsr = fac * fdyn->thsr;

/*----------------------------------------------------------------------*
   Calculate continuity forces of time force vector:
                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /
 *----------------------------------------------------------------------*/
aux = facsr * (vderxy[0][0] + vderxy[1][1]);
for(inode=0;inode<iel;inode++)
{
   eforce[inode] += funct[inode]*aux;
} /* end of loop over inode */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lscalgaltfp */

#endif
