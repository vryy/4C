/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
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

static FLUID_DYNAMIC   *fdyn;
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
void f3_calgalexfv(
                  DOUBLE          *eforce,
		  DOUBLE          *funct,
                  DOUBLE          *edeadn,
		  DOUBLE          *edeadng,
		  DOUBLE           fac,
		  INT              iel
              )
{
DOUBLE  facsl, facsr;
INT     inode,irow,isd;

#ifdef DEBUG
dstrc_enter("f3_calgalexfv");
#endif

fdyn    = alldyn[genprob.numff].fdyn;
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
   for (isd=0;isd<3;isd++)
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
} /* end of f3_calgalexfv */

/*!---------------------------------------------------------------------
\brief stabilisation part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER:
                     /
   + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
                   /

ALE:
                     /
   + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
                   /

EULER/ALE:
                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b^  d_omega
                     /

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'


</pre>
\param   *gls         STAB_PAR_GLS    (i)    stabilisation
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param  **derxy2,     DOUBLE	      (i)    2nd. global derivatives
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  ihoel       INT	      (i)    flag for higer ord. ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void

------------------------------------------------------------------------*/
void f3_calstabexfv(
                     STAB_PAR_GLS    *gls,
                     ELEMENT         *ele,
                     DOUBLE          *eforce,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE          *edead,
                     DOUBLE          *velint,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              iel,
                     INT              ihoel,
                     INT              flag
                   )
{

INT    irow,inode,isd;
DOUBLE sign=ONE;
DOUBLE aux;
DOUBLE taumu,taump;
DOUBLE c,fvts;
DOUBLE fact[3];

#ifdef DEBUG
dstrc_enter("f3_calstabexfv");
#endif

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;

dsassert(ele->e.f3->stab_type == stab_gls,
         "routine with no or wrong stabilisation called");

/*--------------------------------------------------- set some factors */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = fdyn->thsr;
break;
case 1: /* evaluation at n+1 */
   c = fdyn->thsl;
break;
default:
   c = 0.0;
   dserror("value of flag not valid!!!\n");
}

/*----------------------------------------------------------------------*
   Calculate external/convective stab-forces of time force vector:
EULER:
                     /
   + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
                   /
ALE:
                     /
   + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
                   /
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
   fact[0] = edead[0]*fac*taumu*c;
   fact[1] = edead[1]*fac*taumu*c;
   fact[2] = edead[2]*fac*taumu*c;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1]
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
         eforce[irow] += aux*fact[isd];
      } /* end loop over isd */
   } /* end loop over inode */
} /* endif (ele->e.f3->iadvec!=0) */

/*----------------------------------------------------------------------*
   Calculate external/viscous stab-forces of time force vector:
                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b  d_omega
                     /
 *----------------------------------------------------------------------*/
if (gls->ivisc!=0 && ihoel!=0)
{
   if (gls->ivisc==2) sign*=-ONE; /* GLS+ stabilisation */

   fvts = fac*visc*taump*sign*c;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   -= ( derxy2[0][inode]*edead[0]  \
                        + derxy2[3][inode]*edead[1]  \
                        + derxy2[4][inode]*edead[2] + aux*edead[0])*fvts;
      eforce[irow+1] -= ( derxy2[3][inode]*edead[0]  \
                        + derxy2[1][inode]*edead[1]  \
                        + derxy2[5][inode]*edead[2] + aux*edead[1])*fvts;
      eforce[irow+2] -= ( derxy2[4][inode]*edead[0]  \
                        + derxy2[5][inode]*edead[1]  \
                        + derxy2[2][inode]*edead[2] + aux*edead[2])*fvts;
      irow += 3;
   } /* end loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_stabexfv */

/*!---------------------------------------------------------------------
\brief stabilisation part of external forces for pre dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for pre dofs
is calculated:

                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b^  d_omega
                    /

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
    there's only one full element force vector
    for pre-dofs the pointer eforce points to the entry
    eforce[2*iel]

</pre>
\param   *gls      STAB_PAR_GLS    (i)    stabilisation
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void

------------------------------------------------------------------------*/
void f3_calstabexfp(
                     STAB_PAR_GLS    *gls,
                     DOUBLE          *eforce,
                     DOUBLE         **derxy,
                     DOUBLE          *edead,
                     DOUBLE           fac,
                     INT              iel,
                     INT              flag
                   )
{
INT    inode;
DOUBLE c;
DOUBLE taump;
DOUBLE fact[3];

#ifdef DEBUG
dstrc_enter("f3_calstabexfp");
#endif

if (gls->ipres==0) goto end; /* no pressure stabilisation */
fdyn    = alldyn[genprob.numff].fdyn;

/*--------------------------------------------------- set some factors */
taump = fdyn->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = fdyn->thpr;
break;
case 1: /* evaluation at n+1 */
   c = fdyn->thpl;
break;
default:
   c = 0.0;
   dserror("value of flag not valid!!!\n");
}

/*----------------------------------------------------------------------*
   Calculate inertia/pressure stab forces of time force vector:
                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b  d_omega
                    /
 *----------------------------------------------------------------------*/
fact[0] = edead[0]*taump*fac*c;
fact[1] = edead[1]*taump*fac*c;
fact[2] = edead[2]*taump*fac*c;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]
                   + derxy[2][inode]*fact[2]);
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_stabexfp */


/*----------------------------------------------------------------------*/
void f3_calneumann(ELEMENT         *ele,
                   DOUBLE          *eforce,
                   DOUBLE         **xyze,
                   DOUBLE          *funct,
                   DOUBLE         **deriv,
                   DOUBLE         **xjm,
                   DOUBLE          *edeadng
  )
{
  INT       i;
  INT       iel;
  INT       imyrank=0;
  DIS_TYP   typ;
  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA   *data;

  INT       nir,nis;    /*  number of integration nodesin r,s  direction */
  INT       nil;
  INT       ihoel=0;    /* flag for higher order elements               */
  INT       icode=2;    /* flag for eveluation of shape functions -> f3_tri       */
  INT       intc;

  INT       foundsurf;
  INT       ngsurf;
  INT       surf;
  INT       ngline;
  INT       line;
  INT       ngnode;
  GSURF    *gsurf[6];
  NEUM_CONDITION *surfneum[6];
  INT       iedgnod[MAXNOD_F3];
  const INT numdf  = 3; /* dof per node                    */
  INT ili;


  ili = ele->Id_loc;

#ifdef DEBUG
  dstrc_enter("f3_calneumann");
#endif

  /*--------------------------------------------------- initialisation */
  iel = ele->numnp;  /*  number of nodes to me -> 4 (Tet 4) */
  typ = ele->distyp;

  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;
  /*-------------------------------- check for presence of freesurface  */
  foundsurf=0;
  /*----------------------------------- number of lines to this element */
  ngsurf=ele->g.gvol->ngsurf;

  /*-------- loop over lines, check for freesurface conditions on lines */

  for (i=0; i<ngsurf; i++)
  {
    DOUBLE pupu;
    pupu = fdyn->step;
    gsurf[i] = ele->g.gvol->gsurf[i];
    surfneum[i] = gsurf[i]->neum;
#if 0 /* to be checked */
    /*    if (gsurf[i]->proc!=imyrank) surfneum[i]=NULL; */
#endif
    if (surfneum[i]==NULL) continue;
    foundsurf++;
  }

  if (foundsurf==0) goto end;

  /*----- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: /* --> tri - element */
    typ=quad4;
    icode   = 2;
    ihoel   = 1;
    nir = ele->e.f3->nGP[0];
    nis = ele->e.f3->nGP[1];  /* number Gps per element */
/*    dserror("no neumann bc on triangular elements yet"); */
    break;
  case tet4: /* --> tri - element */
    typ=tri3;
    icode   = 3; /* siehe f3_tri, for tri3, icode should be > 1 */
    ihoel   = 1;
    nil = 3;   /* 1 or 3 Gps per element */
/*    dserror("no neumann bc on triangular elements yet"); */
    break;
  default:
    dserror("type %d unknown", typ);
  }


  /*-------------------------------------- set number of gauss points */
  /* We don't expect different numbers of gauss points for different
   * directions. But in case, let's use the larger one on all lines. */
/*  nil = MAX(nir,nis); */

  /*---------------------------------- loop over lines at free surface */

  if (typ==tri3)
  {
    for (surf=0; surf<ngsurf; surf++)
    {
      INT lr;
      INT ls;
      if (surfneum[surf]==NULL) continue;

      /*--------------------------------- check number of nodes on line */
      ngnode = gsurf[surf]->ngnode;

      /*------------------------------------------------ get edge nodes */
      f3_iedg(iedgnod,ele,surf);

      /*------------------------------ integration loop on actual gline */
      for (lr=0;lr<nil;lr++)
      {
/*      for (ls=0;ls<nil;ls++)
        {  */
	DOUBLE    e1;
	DOUBLE    e2;
	DOUBLE    facr;
	DOUBLE    facs;
	DOUBLE    fac;
	DOUBLE    det;
	INT       inode;
	INT       actcurve;
	DOUBLE    acttimefac;
	DOUBLE    acttimefacn;

	/*---------- get values of  shape functions and their derivatives */
/*      e1   = data->txgr[lr][0];
	e2   = data->txgs[ls][0];
	facr = data->twgt[lr][0]; */
/*  Reduced integrated */
/*      e1=1.0/3.0;
	e2=1.0/3.0;
	facr=1.0/2.0;  */
/* Three Gauss points */
	if (lr==0)
	{
	  e1=0.5;
	  e2=0.0;
	  facr=1.0/6.0;
	}
	if (lr==1)
	{
	  e1=0.5;
	  e2=0.5;
	  facr=1.0/6.0;
	}
	if (lr==2)
	{
	  e1=0.0;
	  e2=0.5;
	  facr=1.0/6.0;
	}
/* -----------------------------------------*/

	facs = ONE;
	f3_tri(funct,deriv,NULL,e1,e2,typ,icode);
/*      f2_degrectri(funct,deriv,e1,typ,1); -> only for rectangles
 */
	/*------------------------------- compute jacobian determinant */
	f3_edgejaco1(xyze,deriv,xjm,&det,iedgnod,ngnode,ele);
/*      f3_edgejaco(xyze,deriv,xjm,&det,iedgnod,ngnode,ele); */
	fac = det*facr*facs;

	actcurve = surfneum[surf]->curve-1;

	if (actcurve<0)
	{
	  acttimefac = 1.;
	  acttimefacn= 1.;
	}
	else
	{
	  DOUBLE acttime;
	  acttime = fdyn->acttime;
	  dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
	  acttime = fdyn->acttime-fdyn->dta;
	  dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
	}
	for (i=0;i<numdf;i++)
	{
	  if (surfneum[surf]->neum_onoff.a.iv[i]==0)
	  {
	    edeadng[i] = 0.;
	  }
	  if (surfneum[surf]->neum_onoff.a.iv[i]!=0)
	  {
	    edeadng[i] = surfneum[surf]->neum_val.a.dv[i]*acttimefac;
	  }
	}

	/*-------------------------------------------------------------*/
	/*                                uniform prescribed line load */
	/*-------------------------------------------------------------*/

	/*
	  Calculate external forces:

	  /
          + facs |  v * h   d_gamma
	  /
	*/
	for (inode=0;inode<ngnode;inode++)
	{
	  INT irow, isd;

	  /* matrices already permuted */
	  irow = NUMDOF_FLUID3*iedgnod[inode];
	  for (isd=0;isd<numdf;isd++)
	  {
	    /* NOTE
	     * Right now neumann bc are constant along a line. */
	    eforce[irow] += funct[inode]*(edeadng[isd]*fdyn->thsl)*fac;
	    irow++;
	  }
	}
      }
      /* To avoid double calculation. Won't work in parallel. (?) */
      ele->g.gvol->gsurf[surf]->neum=NULL;
    }
  }
  else
  {
    for (surf=0; surf<ngsurf; surf++)
    {
      INT lr;
      INT ls;
      if (surfneum[surf]==NULL) continue;

      /*--------------------------------- check number of nodes on line */
      ngnode = gsurf[surf]->ngnode;

      /*------------------------------------------------ get edge nodes */
      f3_iedg(iedgnod,ele,surf);

      /*------------------------------ integration loop on actual gline */

      for (lr=0;lr<nir;lr++)
      {
	for (ls=0;ls<nis;ls++)
        {
	  DOUBLE    e1;
	  DOUBLE    e2;
	  DOUBLE    facr;
	  DOUBLE    facs;
	  DOUBLE    fac;
	  DOUBLE    det;
	  INT       inode;
	  INT       actcurve;
	  DOUBLE    acttimefac;
	  DOUBLE    acttimefacn;

/*---------- get values of  shape functions and their derivatives */
	  nil = MAX(nir,nis);
	  e1   = data->qxg[lr][nil-1];
	  e2   = data->qxg[ls][nil-1];
	  facr = data->qwgt[lr][nil-1];
	  facs = data->qwgt[ls][nil-1];
	  f3_rec(funct,deriv,NULL,e1,e2,typ,icode);


/*------------------------------- compute jacobian determinant */
	  f3_edgejaco(xyze,deriv,xjm,&det,iedgnod,ngnode,ele);
	  fac = det*facr*facs;

	  actcurve = surfneum[surf]->curve-1;

	  if (actcurve<0)
	  {
	    acttimefac = 1.;
	    acttimefacn= 1.;
	  }
	  else
	  {
	    DOUBLE acttime;
	    acttime = fdyn->acttime;
	    dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
	    acttime = fdyn->acttime-fdyn->dta;
	    dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
	  }
	  for (i=0;i<numdf;i++)
	  {
	    if (surfneum[surf]->neum_onoff.a.iv[i]==0)
	    {
	      edeadng[i] = 0.;
	    }
	    if (surfneum[surf]->neum_onoff.a.iv[i]!=0)
	    {
	      edeadng[i] = surfneum[surf]->neum_val.a.dv[i]*acttimefac;
	    }
	  }

	  /*-------------------------------------------------------------*/
	  /*                                uniform prescribed line load */
	  /*-------------------------------------------------------------*/

	  /*
	    Calculate external forces:

	    /
	    + facs |  v * h   d_gamma
	    /
	  */
	  for (inode=0;inode<ngnode;inode++)
	  {
	    INT irow, isd;

	    /* matrices already permuted */
	    irow = NUMDOF_FLUID3*iedgnod[inode];
	    for (isd=0;isd<numdf;isd++)
	    {
	      /* NOTE
	       * Right now neumann bc are constant along a line. */
	      eforce[irow] += funct[inode]*(edeadng[isd]*fdyn->thsl)*fac;
	      irow++;
	    }
	  }
	}
	/* To avoid double calculation. Won't work in parallel. (?) */
	ele->g.gvol->gsurf[surf]->neum=NULL;
      }
    }
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


#endif
