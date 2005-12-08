/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
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
\param   *eforce     DOUBLE         (i/o)  element force vector
\param	 *funct      DOUBLE         (i)    nat. shape functions
\param   *edeadn     DOUBLE         (i)    ele dead load at n
\param   *edeadn     DOUBLE         (i)    ele dead load at n+1
\param	  fac        DOUBLE         (i)    weighting factor
\param	  iel        INT            (i)    num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_calgalexfv(
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
dstrc_enter("f2_calgalexfv");
#endif

/*--------------------------------------------------- set some factors */
fdyn  = alldyn[genprob.numff].fdyn;

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
} /* end of f2_calgalexfv */

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

NOTE: for EULER
      velint = U

NOTE: for ALE
      velint = C (ale-convective velocity)

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
void f2_calstabexfv(
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
DOUBLE aux;
DOUBLE sign=ONE;
DOUBLE taumu,taump;
DOUBLE c,fvts;
DOUBLE fact[2];

#ifdef DEBUG
dstrc_enter("f2_calstabexfv");
#endif

/*---------------------------------------------------------- initialise */
fdyn  = alldyn[genprob.numff].fdyn;

dsassert(ele->e.f2->stab_type == stab_gls,
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
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
         irow++;
         eforce[irow] += aux*fact[isd];
      } /* end loop over isd */
   } /* end loop over inode */
} /* endif (ele->e.f2->iadvec!=0) */
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
      eforce[irow]   -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*edead[0] \
                            + derxy2[2][inode]*edead[1])*fvts;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*edead[1] \
                            + derxy2[2][inode]*edead[0])*fvts;
      irow += 2;
   } /* end loop over inode */
} /* endif (ele->e.f2->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_stabexfv */

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
\param   *gls         STAB_PAR_GLS    (i)    stabilisation
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void

------------------------------------------------------------------------*/
void f2_calstabexfp(
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
DOUBLE fact[2];

#ifdef DEBUG
dstrc_enter("f2_calstabexfp");
#endif

if (gls->ipres==0) goto end; /* no pressure stabilisation */

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
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_stabexfp */



/*----------------------------------------------------------------------*/
/*!
  \brief calculate Neumann BC

  \param  *ele       ELEMENT          (i)    actual element
  \param   imyrank   INT              (i)    proc number
  \param  *eforce    DOUBLE           (o)    element force vector
  \param **xyze      DOUBLE           (-)    nodal coordinates at n+theta
  \param  *funct     DOUBLE           (-)    natural shape functions
  \param **deriv     DOUBLE           (-)    deriv. of nat. shape funcs
  \param **xjm       DOUBLE           (-)    jacobian matrix
  \param   *edeadn   DOUBLE           (-)    line load at n
  \param   *edeadng  DOUBLE           (-)    line load at n+g

  \author u.kue (based on the work of genk)
  \date 04/05
 */
/*----------------------------------------------------------------------*/
void f2_calneumann(ELEMENT         *ele,
                   INT              imyrank,
                   DOUBLE          *eforce,
                   DOUBLE         **xyze,
                   DOUBLE          *funct,
                   DOUBLE         **deriv,
                   DOUBLE         **xjm,
                   DOUBLE          *edeadn,
                   DOUBLE          *edeadng
  )
{
  INT       i;
  INT       iel;
  DIS_TYP   typ;
  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA   *data;

  INT       nir,nis;    /* number of integration nodesin r,s direction  */
  INT       nil;
  INT       ihoel=0;    /* flag for higher order elements               */
  INT       icode=2;    /* flag for eveluation of shape functions       */
  INT       intc;

  INT       foundline;
  INT       ngline;
  INT       line;
  INT       ngnode;
  GLINE    *gline[4];
  NEUM_CONDITION *lineneum[4];
  INT       iedgnod[MAXNOD_F2];
  const INT numdf  = 2; /* dof per node                    */


#ifdef DEBUG
  dstrc_enter("f2_calneumann");
#endif

  /*--------------------------------------------------- initialisation */
  iel = ele->numnp;
  typ = ele->distyp;

  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*-------------------------------- check for presence of freesurface  */
  foundline=0;
  /*----------------------------------- number of lines to this element */
  ngline=ele->g.gsurf->ngline;
  /*-------- loop over lines, check for freesurface conditions on lines */
  for (i=0; i<ngline; i++)
  {
    gline[i] = ele->g.gsurf->gline[i];
    lineneum[i] = gline[i]->neum;
#if 0 /* to be checked */
    if (gline[i]->proc!=imyrank) lineneum[i]=NULL;
#endif
    if (lineneum[i]==NULL) continue;
    foundline++;
  }
  if (foundline==0) goto end;

  /*----- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case quad4: case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2->nGP[0];
    nis = ele->e.f2->nGP[1];
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tri3:    /* initialise integration */
    nir  = ele->e.f2->nGP[0];
    nis  = 1;
    intc = ele->e.f2->nGP[1];
    dserror("no neumann bc on triangular elements yet");
    break;
  default:
    dserror_args(__FILE__, __LINE__, "type %d unknown", typ);
  }


  /*-------------------------------------- set number of gauss points */
  /* We don't expect different numbers of gauss points for different
   * directions. But in case, let's use the larger one on all lines. */
  nil = MAX(nir,nis);

  /*---------------------------------- loop over lines at free surface */
  for (line=0; line<ngline; line++)
  {
    INT lr;
    if (lineneum[line]==NULL) continue;

    /*--------------------------------- check number of nodes on line */
    ngnode = gline[line]->ngnode;

    /*------------------------------------------------ get edge nodes */
    f2_iedg(iedgnod,ele,line,0);

    /*------------------------------ integration loop on actual gline */
    for (lr=0;lr<nil;lr++)
    {
      DOUBLE    e1;
      DOUBLE    facr;
      DOUBLE    fac;
      DOUBLE    det;
      INT       inode;
      INT       actcurve;
      DOUBLE    acttimefac;
      DOUBLE    acttimefacn;

      /*---------- get values of  shape functions and their derivatives */
      e1   = data->qxg[lr][nil-1];
      facr = data->qwgt[lr][nil-1];
      f2_degrectri(funct,deriv,e1,typ,1);

      /*------------------------------- compute jacobian determinant */
      f2_edgejaco(xyze,deriv,xjm,&det,ngnode,iedgnod);
      fac = det*facr;

      actcurve = lineneum[line]->curve-1;
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
        if (lineneum[line]->neum_onoff.a.iv[i]==0)
        {
          edeadn[i]  = 0.;
          edeadng[i] = 0.;
        }
        if (lineneum[line]->neum_onoff.a.iv[i]!=0)
        {
          edeadn[i]  = lineneum[line]->neum_val.a.dv[i]*acttimefacn;
          edeadng[i] = lineneum[line]->neum_val.a.dv[i]*acttimefac;
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
        irow = NUMDOF_FLUID2*iedgnod[inode];
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
    ele->g.gsurf->gline[line]->neum=NULL;
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


#endif
/*! @} (documentation module close)*/
