/*!----------------------------------------------------------------------
\file
\brief evaluate "VMM" force vectors for fluid2 element

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
\brief evaluate Galerkin part of submesh "VMM" force vector for fluid2

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh "VMM" force vector on
the rhs is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smevfor   DOUBLE	   (i/o)  submesh element vmm force vec.
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **derxy2    DOUBLE	   (i)    2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void

------------------------------------------------------------------------*/
void f2_calsmfv(FLUID_DYN_ML	*mlvar,
		DOUBLE         **smevfor,
		DOUBLE  	*velint,
		DOUBLE         **vderxy,
		DOUBLE  	*smfunct,
		DOUBLE  	*funct,
		DOUBLE         **derxy,
		DOUBLE         **derxy2,
		DOUBLE  	 fac,
		DOUBLE  	 visc,
		INT		 smiel,
		INT		 iel,
                INT		 ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 *----------------------------------------------------------------------*/
INT    irow,icol,icn,isd;
DOUBLE con,beta,divv,cb,facsl,facpl;
DOUBLE aux;

#ifdef DEBUG
dstrc_enter("f2_calsmfv");
#endif

fdyn = alldyn[genprob.numff].fdyn;

facsl = fac*fdyn->thsl;
facpl = fac*fdyn->thpl;
con   = facsl*visc;

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
   Calculate temporal forces for velocity bubble function:
 *----------------------------------------------------------------------*/
if (fdyn->nis==0 && mlvar->transterm<2)
{
/*---------------------------------------------------------------------*
    /
   |  -  w * ls_shape_function d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    aux = funct[icol]*fac;
    for (irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
  } /* end loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate convective forces for velocity bubble function:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - (theta*dt) * w * u_old * grad(ls_shape_function) d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  aux = (velint[0]*derxy[0][icol] + velint[1]*derxy[1][icol])*facsl;
  for (irow=0;irow<smiel;irow++)
  {
    smevfor[irow][icol] -= smfunct[irow]*aux;
  } /* end loop over irow */
} /* end loop over icol */

if (fdyn->conte!=0)
{
  cb = beta*divv*facsl;
/*----------------------------------------------------------------------*
    /
   |  - (theta*dt) * beta * w * ls_shape_function * div(u_old)  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    aux = funct[icol]*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] -= smfunct[irow]*aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate viscous forces for velocity bubble function:
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
/*----------------------------------------------------------------------*
    /
   |  (theta*dt) * w * nue * delta(ls_shape_function)  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    aux = (derxy2[0][icol] + derxy2[1][icol])*con;
    for(irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] += smfunct[irow]*aux;
    } /* end loop over irow */
  } /* end loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate forces for pressure bubble function:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - (theta*dt) * w * d/dxi(ls_shape_function) d_omega
  /
 *----------------------------------------------------------------------*/
icol=iel;
for (icn=0;icn<iel;icn++)
{
  for (isd=0;isd<2;isd++)
  {
    aux = derxy[isd][icn]*facpl;
    for (irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  }
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calsmfv */

/*!---------------------------------------------------------------------
\brief evaluate stabilization part of sm "VMM" force vector for fluid2

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh "VMM" force
vector on the rhs is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smevfor   DOUBLE	   (i/o)  submesh element vmm force vec.
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param **smderxy2  DOUBLE	   (i)    sm 2nd glo. der. of shape fun.
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **derxy2    DOUBLE	   (i)    2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void

------------------------------------------------------------------------*/
void f2_calstabsmfv(FLUID_DYN_ML    *mlvar,
		    DOUBLE         **smevfor,
		    DOUBLE  	    *velint,
		    DOUBLE         **vderxy,
		    DOUBLE  	    *smfunct,
		    DOUBLE         **smderxy,
		    DOUBLE         **smderxy2,
		    DOUBLE  	    *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
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
INT    irow,icol,icn,isd;
DOUBLE con,ccon,beta,divv,cb,ccb,cbb;
DOUBLE aux,auxc;
DOUBLE sign,tau;

#ifdef DEBUG
dstrc_enter("f2_calstabsmfv");
#endif

fdyn = alldyn[genprob.numff].fdyn;

/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;
con = fac*tau;

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

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate temporal stabilization:
 *----------------------------------------------------------------------*/
if (fdyn->nis==0)
{
  if (mlvar->quastabub==0 && ABS(mlvar->smstado)<3)
  {
    if (mlvar->smstado==-1) con = fac*tau;
    else                    con = fac*tau*sign;

    if (mlvar->transterm<2)
    {
/*---------------------------------------------------------------------*
    /
   |  -/+  tau * (1/(theta*dt)) * w * ls_shape_function d_omega
  /
 *----------------------------------------------------------------------*/
      for (icol=0;icol<iel;icol++)
      {
        aux = funct[icol]*con*(ONE/fdyn->thsl);
        for (irow=0;irow<smiel;irow++)
        {
          smevfor[irow][icol] -= smfunct[irow]*aux;
        } /* end loop over irow */
      } /* end loop over icol */
    }

/*---------------------------------------------------------------------*
    /
   |  -/+  tau * w * u_old * grad(ls_shape_function) d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      aux = (velint[0]*derxy[0][icol] + velint[1]*derxy[1][icol])*con;
      for (irow=0;irow<smiel;irow++)
      {
    	smevfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
    } /* end loop over icol */

    if (fdyn->conte!=0)
    {
      cb = con*beta;
/*---------------------------------------------------------------------*
    /
   |  -/+  tau * w * ls_shape_function * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
      for (icol=0;icol<iel;icol++)
      {
        aux = funct[icol]*divv*cb;
        for (irow=0;irow<smiel;irow++)
        {
    	  smevfor[irow][icol] -= smfunct[irow]*aux;
        } /* end loop over irow */
      } /* end loop over icol */
    }

    if (ihoel!=0)
    {
      ccon = con*visc;
/*---------------------------------------------------------------------*
    /
   |  +/-  tau * w * nue * delta(ls_shape_function) d_omega
  /
 *----------------------------------------------------------------------*/
      for (icol=0;icol<iel;icol++)
      {
        aux = (derxy2[0][icol] + derxy2[1][icol])*ccon;
        for (irow=0;irow<smiel;irow++)
        {
    	  smevfor[irow][icol] += smfunct[irow]*aux;
        } /* end loop over irow */
      } /* end loop over icol */
    }

/*----------------------------------------------------------------------*
    /
   |  -/+ tau * w * d/dxi(ls_shape_function) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=iel;
    for (icn=0;icn<iel;icn++)
    {
      for (isd=0;isd<2;isd++)
      {
        aux = derxy[isd][icn]*con;
        for (irow=0;irow<smiel;irow++)
        {
          smevfor[irow][icol] -= smfunct[irow]*aux;
        } /* end loop over irow */
        icol++;
      }
    } /* end loop over icol */
  }

  if (mlvar->transterm<2)
  {
    con = fac*tau;
/*----------------------------------------------------------------------*
   Calculate convective stabilization for temporal forces of vel. bub.:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - tau * u_old * grad(w) * ls_shape_function d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      auxc = funct[icol]*con;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (velint[0]*smderxy[0][irow]+velint[1]*smderxy[1][irow])*auxc;
        smevfor[irow][icol] -= aux;
      } /* end loop over irow */
    } /* end loop over icol */

    if (fdyn->conte!=0)
    {
      cb = con*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * w * div(u_old) * ls_shape_function  d_omega
  /
 *----------------------------------------------------------------------*/
      for (icol=0;icol<iel;icol++)
      {
        aux = funct[icol]*divv*cb;
        for (irow=0;irow<smiel;irow++)
        {
          smevfor[irow][icol] -= smfunct[irow]*aux;
        } /* end of loop over irow */
      } /* end of loop over icol */
    }

/*----------------------------------------------------------------------*
   Calculate viscous stabilization for temporal forces of vel. bubble:
 *----------------------------------------------------------------------*/
    if (ihoel!=0)
    {
      ccon = con*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- nue * delta(w) * ls_shape_function  d_omega
  /
 *----------------------------------------------------------------------*/
      for (icol=0;icol<iel;icol++)
      {
        auxc = funct[icol]*ccon;
        for(irow=0;irow<smiel;irow++)
        {
          aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
          smevfor[irow][icol] += aux;
        } /* end loop over irow */
      } /* end loop over icol */
    } /* end of viscous stabilization for higher order elements */
  } /* end of stabilization for temporal forces */
}

con = fac*tau*fdyn->thsl;

/*----------------------------------------------------------------------*
   Calculate convective stabilisation for convective forces of vel. bub.:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - tau * (theta*dt) * u_old * grad(w) * u_old * grad(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  auxc = (velint[0]*derxy[0][icol] + velint[1]*derxy[1][icol])*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
    smevfor[irow][icol] -= aux;
  } /* end of loop over irow */
} /* end of loop over icol */

if (fdyn->conte!=0)
{
  cb  = con*beta;
  cbb = cb*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau * (th*dt) * u_old * grad(w) * ls_s_f * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    auxc = funct[icol]*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
      smevfor[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau * (th*dt) * w * div(u_old) * u_old * grad(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
  cb = cb*sign;
  for (icol=0;icol<iel;icol++)
  {
    aux = (velint[0]*derxy[0][icol] + velint[1]*derxy[1][icol])*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] -= aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */

/*----------------------------------------------------------------------*
    /
   |  -/+ beta * beta * tau * (th*dt) * w * div(u_old) * bub * div(u_old)
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    aux = funct[icol]*divv*divv*cbb;
    for (irow=0;irow<smiel;irow++)
    {
      smevfor[irow][icol] -= aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate convective stabilisation for viscous forces of vel. bub.:
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
  ccon = con*visc;
/*----------------------------------------------------------------------*
    /
   |  tau * (theta*dt) * u_old * grad(w) * nue * delta(ls_s_f)  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    auxc = (derxy2[0][icol] + derxy2[1][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
      smevfor[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

  if (fdyn->conte!=0)
  {
    ccb  = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * (theta*dt) * nue * w * div(u_old) * delta(ls_s_f)
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      aux = (derxy2[0][icol] + derxy2[1][icol])*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
  	smevfor[irow][icol] += smfunct[irow]*aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  }
}

/*----------------------------------------------------------------------*
   Calculate viscous stabilization for higher-order ele. for vel. bub.:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con * visc * visc * sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * (theta*dt) * nue * delta(w) * nue * delta(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<iel;icol++)
  {
    auxc = (derxy2[0][icol] + derxy2[1][icol])*ccon;
    for (irow=0;irow<iel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smevfor[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

/*----------------------------------------------------------------------*
    /
   |  +/- tau * (theta*dt) * nue * delta(w) * u_old * grad(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
  ccon = con * visc * sign;

  for (icol=0;icol<iel;icol++)
  {
    auxc = (velint[0]*derxy[0][icol] + velint[1]*derxy[1][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
      smevfor[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

  if (fdyn->conte!=0)
  {
    ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * (theta*dt) * nue * delta(w) * ls_s_f * div(u_old)
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      auxc = funct[icol]*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
    	aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
    	smevfor[irow][icol] += aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  }
} /* end of viscous stabilization for higher-order elements for vel. bub. */

/*----------------------------------------------------------------------*
   Calculate stabilization for forces of pressure bubble function:
 *----------------------------------------------------------------------*/
con = fac*tau*fdyn->thpl;
/*----------------------------------------------------------------------*
    /
   |  - tau * (theta*dt) * u_old * grad(w) * d/dxi(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
icol=iel;
for (icn=0;icn<iel;icn++)
{
  for (isd=0;isd<2;isd++)
  {
    auxc = derxy[isd][icn]*con;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow] + velint[1]*smderxy[1][irow])*auxc;
      smevfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  }
} /* end loop over icol */

if (fdyn->conte!=0)
{
  cb = con*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * (theta*dt) * beta * w * div(u_old) * d/dxi(ls_s_f) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=iel;
  for (icn=0;icn<iel;icn++)
  {
    for (isd=0;isd<2;isd++)
    {
      aux = derxy[isd][icn]*divv*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smevfor[irow][icol] -= smfunct[irow]*aux;
      } /* end of loop over irow */
      icol++;
    }
  } /* end of loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate viscous stabilization for higher-order ele. for pre. bub.:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * (theta*dt) * nue * delta(w) * d/dxi(ls_s_f)  d_omega
  /
 *----------------------------------------------------------------------*/
  icol=iel;
  for (icn=0;icn<iel;icn++)
  {
    for (isd=0;isd<2;isd++)
    {
      auxc = derxy[isd][icn]*ccon;
      for(irow=0;irow<smiel;irow++)
      {
        aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
        smevfor[irow][icol] += aux;
      } /* end loop over irow */
      icol++;
    }
  } /* end loop over icol */
} /* end of viscous stabilization for higher-order elements for pre. bub. */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabsmfv */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of large-scale rhs (vel. dofs) for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the large-scale rhs for the
velocity dofs is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **eforce   DOUBLE	   (i/o)  element iterative rhs
\param  *velint    DOUBLE	   (i)    velocity at int point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param  *smfint    DOUBLE	   (i)    rhs bubble functions
\param **smfderxy  DOUBLE	   (i)    global deriv. of rhs bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   visc      DOUBLE	   (i)    fluid viscosity
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbfv(FLUID_DYN_ML    *mlvar,
	       DOUBLE          *eforce,
	       DOUBLE          *velint,
	       DOUBLE         **vderxy,
	       DOUBLE          *funct,
	       DOUBLE         **derxy,
	       DOUBLE          *smfint,
	       DOUBLE         **smfderxy,
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
INT     irow,inode,isd,j;
DOUBLE  con,facsl,aux,auxc,beta,divv;

#ifdef DEBUG
dstrc_enter("f2_calbfv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------------------------- */
facsl=fac*fdyn->thsl;
con=facsl*visc;

/*-----------------------------------------------------------------------
   Calculate bubble part of temporal forces:
    /
   |  - v * f_bub   d_omega
   /
----------------------------------------------------------------------- */
if ((fdyn->nis==0 && mlvar->transterm==0) ||
    (fdyn->nis==0 && mlvar->transterm==2))
{
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    aux=funct[inode]*fac;
    for (isd=0;isd<2;isd++)
    {
      irow++;
      eforce[irow] -= aux*smfint[isd];
    } /* end loop over irn */
  } /* end loop over icn */
}

/*-----------------------------------------------------------------------
   Calculate bubble part of convective forces (convective part):
----------------------------------------------------------------------- */
/*-----------------------------------------------------------------------
    /
   |  - v * u_old * grad(f_bub)   d_omega
   /
----------------------------------------------------------------------- */
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    auxc=funct[inode]*facsl;
    for (isd=0;isd<2;isd++)
    {
      aux = (velint[0]*smfderxy[0][isd]+velint[1]*smfderxy[1][isd])*auxc;
      irow++;
      eforce[irow] -= aux;
    } /* end loop over irn */
  } /* end loop over icn */

  if (fdyn->conte!=0)
  {
/*----------------------------------------------------------------------*
    /
   |  - beta * v * f_bub * div(u_old)   d_omega
   /
*----------------------------------------------------------------------*/
    if (fdyn->conte==1) beta = ONE;
    else beta = ONE/TWO;
    divv= vderxy[0][0]+vderxy[1][1];
    irow=-1;
    for (inode=0;inode<iel;inode++)
    {
      aux=funct[inode]*divv*beta*facsl;
      for (isd=0;isd<2;isd++)
      {
        irow++;
        eforce[irow] -= aux*smfint[isd];
      } /* end loop over irn */
    } /* end loop over icn */
  }

/*-----------------------------------------------------------------------
   Calculate bubble part of convective forces (reactive part):
----------------------------------------------------------------------- */
if(fdyn->nir != 0)
{
/*-----------------------------------------------------------------------
    /
   |  - v * f_bub * grad(u_old)   d_omega
   /
----------------------------------------------------------------------- */
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    auxc=funct[inode]*facsl;
    for (isd=0;isd<2;isd++)
    {
      aux = (smfint[0]*vderxy[isd][0]+smfint[1]*vderxy[isd][1])*auxc;
      irow++;
      eforce[irow] -= aux;
    } /* end loop over irn */
  } /* end loop over icn */

  if (fdyn->conte!=0)
  {
/*----------------------------------------------------------------------*
    /
   |  - beta * v * u_old * div(f_bub)   d_omega
   /
*----------------------------------------------------------------------*/
    if (fdyn->conte==1) beta = ONE;
    else beta = ONE/TWO;
    irow=-1;
    auxc= smfderxy[0][0]+smfderxy[1][1];
    for (inode=0;inode<iel;inode++)
    {
      aux=funct[inode]*beta*facsl*auxc;
      for (isd=0;isd<2;isd++)
      {
        irow++;
        eforce[irow] -= aux*velint[isd];
      } /* end loop over irn */
    } /* end loop over icn */
  }
}

/*-----------------------------------------------------------------------
   Calculate bubble part of viscous forces:
----------------------------------------------------------------------- */
if (fdyn->vite==0)
{
/*-----------------------------------------------------------------------
    /
   |  - nue * grad(v) : grad(f_bub)  d_omega
   /
----------------------------------------------------------------------- */
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      eforce[irow] -= (derxy[0][inode]*smfderxy[isd][0] \
                       +derxy[1][inode]*smfderxy[isd][1])*con;
    } /* end loop over irn */
  } /* end loop over icn */
}
else
{
/*-----------------------------------------------------------------------
    /
   |  - 2 * nue * eps(v) : eps(f_bub)  d_omega
   /
----------------------------------------------------------------------- */
  irow=-1;
  for (inode=0;inode<iel;inode++)
  {
    for (isd=0;isd<2;isd++)
    {
      irow++;
      for (j=0;j<2;j++)
      {
        eforce[irow] -= derxy[j][inode]*(smfderxy[isd][j] \
	                                 +smfderxy[j][isd])*con;
      }
    } /* end loop over irn */
  } /* end loop over icn */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbfv */

/*!---------------------------------------------------------------------
\brief evaluate bubble part of large-scale rhs (pre. dofs) for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the large-scale rhs for the
pressure dofs is calculated.

</pre>
\param **eforce   DOUBLE	   (i/o)  element iterative rhs
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **smfderxy  DOUBLE	   (i)    global deriv. of rhs bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor
\param   iel	   INT  	   (i)	  number of element nodes
\return void

------------------------------------------------------------------------*/
void f2_calbfp(DOUBLE          *eforce,
	       DOUBLE          *funct,
	       DOUBLE         **smfderxy,
	       DOUBLE           fac,
	       INT              iel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 *----------------------------------------------------------------------*/
INT     inode,posr;
DOUBLE  facsl,aux;

#ifdef DEBUG
dstrc_enter("f2_calbfp");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------------------------- */
facsl=fac*fdyn->thsl;

/*-----------------------------------------------------------------------
   Calculate bubble part of continuity forces:
     /
    |  +/(-) q * div(f_bub)   d_omega
   /
----------------------------------------------------------------------- */
aux=(smfderxy[0][0]+smfderxy[1][1])*facsl;
for (inode=0;inode<iel;inode++)
{
  posr = 2*iel + inode;
  eforce[posr] += aux*funct[inode];
} /* end loop over icn */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calbfp */

#endif
