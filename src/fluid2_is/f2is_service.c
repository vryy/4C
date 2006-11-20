/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2_is element

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifdef D_FLUID2_IS

#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2_is.h"

static INT PREDOF = 2;
#ifdef D_FSI
static INT NUMDF = 3;
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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
\brief set all arrays for element calculation

<pre>                                                         genk 04/02

get the element velocities and the pressure at different times

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the
      transformation in every time step

NOTE: if there is no classic time rhs (as described in WAW) the array
         eveln is misused and does NOT contain the velocity at time (n)
	 but rather a linear combination of old velocities and
	 accelerations depending upon the time integration scheme!!!!!
</pre>
\param   *ele      ELEMENT          (i)  actual element
\param  **xyze     DOUBLE           (o)  nodal coordinates
\param  **eveln    DOUBLE           (o)  ele vels at time n
\param  **evelng   DOUBLE           (o)  ele vels at time n+g
\param  **evhist   DOUBLE           (o)  history vector
\param   *epren    DOUBLE           (o)  ele pres at time n
\param   *edeadn   DOUBLE           (o)  ele dead load at n (selfweight)
\param   *edeadng  DOUBLE           (o)  ele dead load at n+g (selfweight)
\param   *ipos                      (i)  node array positions
\param   *hasext   INT              (o)  flag for external loads
\return void

------------------------------------------------------------------------*/
void f2is_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *epren,
  DOUBLE          *edeadn,
  DOUBLE          *edeadng,
  ARRAY_POSITION *ipos,
  INT             *hasext
  )
{
  INT    i;           /* simply some counters                             */
  INT    actcurve;    /* actual time curve                                */
  INT    nodesol, nodehist; /* position flags for sol_increment           */
  INT    ndsolold;          /* position flags for sol_increment           */
  DOUBLE acttimefac;  /* time factor from actual curve                    */
  DOUBLE acttimefacn; /* time factor at time (n)                          */
  DOUBLE acttime;
  NODE  *actnode;     /* actual node                                      */
  GSURF *actgsurf;

#ifdef DEBUG
  dstrc_enter("f2is_calset");
#endif

/*------------------------------------------------------- initialise ---*/
  fdyn = alldyn[genprob.numff].fdyn;
  nodesol  = ipos->velnp;
  nodehist = ipos->hist;
  ndsolold = ipos->veln;

/*-------------------------------------------- set element coordinates -*/
  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  /* -> implicit time integration method ---------*/

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[nodesol][0];
    evelng[1][i]=actnode->sol_increment.a.da[nodesol][1];

/*--------------------------------------- set vel. histories at (n) ---*/
    evhist[0][i] = actnode->sol_increment.a.da[nodehist][0];
    evhist[1][i] = actnode->sol_increment.a.da[nodehist][1];
/*------------------------------------ set element velocities at (n) ---*/
    eveln[0][i]=actnode->sol_increment.a.da[ndsolold][0];
    eveln[1][i]=actnode->sol_increment.a.da[ndsolold][1];
/*-------------------------------------- set supported pressures (n+1) */
  }

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];
    if (actnode->sol_increment.sdim <= PREDOF)
      break;
/*---------------------------------------------- set pressures (n+1) ---*/
    epren[i]   =actnode->sol_increment.a.da[nodesol][PREDOF];
  }

  /* check for dead load */
  actgsurf = ele->g.gsurf;
  if (actgsurf->neum!=NULL)
  {
    actcurve = actgsurf->neum->curve-1;
    if (actcurve<0)
    {
      acttimefac =ONE;
      acttimefacn=ONE;
    }
    else
    {
      acttime = fdyn->acttime;
      dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
      acttime = fdyn->acttime-fdyn->dta;
      dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
    }
    for (i=0;i<2;i++)
    {
      if (actgsurf->neum->neum_onoff.a.iv[i]==0)
      {
	edeadn[i]  = ZERO;
	edeadng[i] = ZERO;
      }
      if (actgsurf->neum->neum_type==neum_dead  &&
          actgsurf->neum->neum_onoff.a.iv[i]!=0)
      {
	edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*acttimefacn;
	edeadng[i] = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
	(*hasext)++;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



/*!---------------------------------------------------------------------
\brief set all arrays for element calculation for ALE

<pre>                                                         genk 10/02

get the element velocities and the pressure at different times

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the
      transformation in every time step

</pre>

\param   *ele       ELEMENT         (i)    actual element
\param   *ele       ELEMENT	    (i)    actual element
\param  **xyz0      DOUBLE          (o)    nodal coordinates at initial time
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\param  **eveln     DOUBLE          (o)    ele vels at time n
\param  **evelng    DOUBLE          (o)    ele vels at time n+g
\param  **evhist    DOUBLE          (o)    history vector
\param  **ealecovn  DOUBLE          (o)    ALE-convective vels at time n
\param  **ealecovng DOUBLE          (o)    ALE-convective vels at time n+g
\param  **egridv    DOUBLE          (o)    element grid velocity
\param   *epren     DOUBLE          (o)    ele pres at time n
\param   *edeadn    DOUBLE          (o)    ele dead load at n (selfweight)
\param   *edeadng   DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *ekappan   DOUBLE          (o)    nodal curvature at n
\param   *ekappang  DOUBLE          (o)    nodal curvature at n+g
\param   *ephin     DOUBLE          (o)    nodal height function at n
\param   *ephing    DOUBLE          (o)    nodal height function at n+g
\param   *ipos                      (i)    node array positions
\param   *hasext    INT             (o)    flag for external loads
\param    is_relax  INT             (i)    flag, if it's for relax.-param
\return void

------------------------------------------------------------------------*/
void f2is_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **ealecovn,
  DOUBLE         **ealecovng,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadn,
  DOUBLE          *edeadng,
  DOUBLE          *ekappan,
  DOUBLE          *ekappang,
  DOUBLE          *ephin,
  DOUBLE          *ephing,
  DOUBLE         **evnng,
  DOUBLE         **evnn,
  ARRAY_POSITION *ipos,
  INT             *hasext,
  INT              is_relax
  )
{
#ifdef D_FSI
  INT    i;           /* simply some counters                             */
  INT    actcurve;    /* actual time curve                                */
  INT    actmat;
  DOUBLE acttimefac;  /* time factor from actual curve                    */
  DOUBLE acttimefacn; /* time factor at time (n)                          */
  DOUBLE acttime;
  DOUBLE dens;
  NODE  *actfnode;    /* actual fluid node                                */
  GSURF *actgsurf;    /* actual gsurf                                     */

#ifdef DEBUG
  dstrc_enter("f2is_calseta");
#endif

  fdyn  = alldyn[genprob.numff].fdyn;

/*-------------------------------------------- set element coordinates */
  if (is_relax)
    f2_alecoor_sd(ele,xyze);
  else
    f2_alecoor(ele,xyze);

/*----------------------------------------------------------------------*
  | position of the different solutions:                                 |
  | node->sol_incement: solution history used for calculations           |
  |       sol_increment.a.da[0][i]: solution at (n-1)                    |
  |       sol_increment.a.da[1][i]: solution at (n)                      |
  |       sol_increment.a.da[2][i]: solution at (n+g)                    |
  |       sol_increment.a.da[3][i]: solution at (n+1)                    |
  |       sol_increment.a.da[4][i]: grid velocity                        |
  |       sol_increment.a.da[5][i]: convective velocity at (n)           |
  |       sol_increment.a.da[6][i]: convective velocity at (n+1)         |
  *----------------------------------------------------------------------*/

/* -> implicit time integration method ----------*/

  for (i=0;i<ele->numnp;i++)
  {
    actfnode=ele->node[i];
/*------------------------------------ set element velocities (n+gamma) */
    evelng[0][i]   =actfnode->sol_increment.a.da[ipos->velnp][0];
    evelng[1][i]   =actfnode->sol_increment.a.da[ipos->velnp][1];
    ealecovng[0][i]=actfnode->sol_increment.a.da[ipos->convnp][0];
    ealecovng[1][i]=actfnode->sol_increment.a.da[ipos->convnp][1];
    egridv[0][i]   =actfnode->sol_increment.a.da[ipos->gridv][0];
    egridv[1][i]   =actfnode->sol_increment.a.da[ipos->gridv][1];

/*--------------------------------------- set vel. histories at (n) ---*/
    evhist[0][i] = actfnode->sol_increment.a.da[ipos->hist][0];
    evhist[1][i] = actfnode->sol_increment.a.da[ipos->hist][1];
/*------------------------------------ set element velocities at (n) ---*/
    eveln[0][i]=actfnode->sol_increment.a.da[ipos->veln][0];
    eveln[1][i]=actfnode->sol_increment.a.da[ipos->veln][1];

    /* noch einmal kurz zurueck zur time rhs */
    ealecovn[0][i]=actfnode->sol_increment.a.da[ipos->convn][0];
    ealecovn[1][i]=actfnode->sol_increment.a.da[ipos->convn][1];
  }

  for(i=0;i<ele->numnp;i++)
  {
    actfnode=ele->node[i];
    if (actfnode->sol_increment.sdim <= PREDOF)
      break;
/*---------------------------------------------- set pressures (n+1) ---*/
    epren[i]   =actfnode->sol_increment.a.da[ipos->velnp][PREDOF];
  }

/*----------------------------------------------- check for dead load */
  actgsurf = ele->g.gsurf;
  if (actgsurf->neum!=NULL)
  {
    if (actgsurf->neum->neum_type==neum_LAS)
    {
      actcurve = actgsurf->neum->curve-1;
      if (actcurve<0)
	dserror("No Time curve given for neum_LAS!\n");
      acttime = fdyn->acttime;
      dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
      acttime=fdyn->acttime-fdyn->dta;
      dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
      edeadn[0]  = actgsurf->neum->neum_val.a.dv[0]*acttimefacn;
      edeadng[0] = actgsurf->neum->neum_val.a.dv[0]*acttimefac;
      edeadn[1]  = actgsurf->neum->neum_val.a.dv[1];
      edeadng[1] = actgsurf->neum->neum_val.a.dv[1];
      (*hasext)++;
    }
    else
    {
      actmat=ele->mat-1;
      dens = mat[actmat].m.fluid->density;
      actcurve = actgsurf->neum->curve-1;
      if (actcurve<0) acttimefac=ONE;
      else  dyn_facfromcurve(actcurve,fdyn->acttime,&acttimefac) ;
      for (i=0;i<2;i++)
      {
	actcurve = actgsurf->neum->curve-1;
	if (actcurve<0)
	{
	  acttimefac =ONE;
	  acttimefacn=ONE;
	}
	else
	{
	  acttime = fdyn->acttime;
	  dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
	  acttime = fdyn->acttime-fdyn->dta;
	  dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
	}
	for (i=0;i<2;i++)
	{
	  if (actgsurf->neum->neum_onoff.a.iv[i]==0)
	  {
	    edeadn[i]  = ZERO;
	    edeadng[i] = ZERO;
	  }
	  if (actgsurf->neum->neum_type==neum_dead  &&
	      actgsurf->neum->neum_onoff.a.iv[i]!=0)
	  {
	    edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*acttimefacn;
	    edeadng[i] = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
	    (*hasext)++;
	  }
	}
      }
    }
  }

/*------------------------------------------ curvature at free surface */
  if (ele->e.f2is->fs_on>0 && fdyn->surftens!=0)
  {
    if (fdyn->fsstnif!=0)
    {
      for (i=0;i<ele->numnp;i++)
	ekappan[i]=ele->e.f2is->kappa_ND.a.da[i][0];
    }
    if (fdyn->fsstnii!=0)
    {
      for (i=0;i<ele->numnp;i++)
	ekappang[i]=ele->e.f2is->kappa_ND.a.da[i][1];
    }
  }

/*------------------------------------ height function at free surface */
  if (ele->e.f2is->fs_on==5)
  {
    for (i=0;i<ele->numnp;i++)
    {
      actfnode=ele->node[i];
      if(actfnode->xfs==NULL) continue;
#if 0
      ephing[i]=actfnode->sol_increment.a.da[ipos->velnp][3];
#endif
      ephing[i]=actfnode->xfs[1];
      ephin[i] =actfnode->sol_increment.a.da[ipos->velnp][3];
    }
  }

/*------------------------------------------- normal at free surface */
  if (ele->e.f2is->fs_on)
  {
    for(i=0;i<ele->numnp;i++)
    {
      actfnode=ele->node[i];
      if (actfnode->actn==NULL) continue;
      evnng[0][i]=actfnode->actn[0];
      evnng[1][i]=actfnode->actn[1];
      if (actfnode->oldn==NULL) continue;
      evnn[0][i]=actfnode->oldn[0];
      evnn[1][i]=actfnode->oldn[1];
    }
  }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

#else
  dserror("FSI-functions not compiled in!\n");
#endif
  return;
}


#endif

