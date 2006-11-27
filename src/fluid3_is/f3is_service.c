/*!----------------------------------------------------------------------
\file
\brief service routines for fluid3 element

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>


 *----------------------------------------------------------------------*/
#ifdef D_FLUID3_IS
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3_is.h"

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

static INT PREDOF = 3;
static FLUID_DYNAMIC   *fdyn;
/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         genk 05/02

get the element velocities and the pressure at different times

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the
      transformation in every time step

</pre>
\param   *ele      ELEMENT          (i)    actual element
\param  **xyze     DOUBLE           (o)    nodal coordinates
\param  **ehist    DOUBLE           (o)    ele history data
\param  **evelng   DOUBLE           (o)    ele vels at time n+g
\param   *epren    DOUBLE           (o)    ele pres at time n
\param   *edeadng  DOUBLE           (o)    ele dead load at n+g (selfweight)
\param   *ipos                      (i)    node array positions
\param   *hasext   INT              (o)    flag for external loads
\return void

------------------------------------------------------------------------*/
void f3is_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **ehist,
  DOUBLE         **evelng,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos,
  INT             *hasext
  )
{
  INT    i;
  INT    actcurve;    /* actual time curve                                */
  DOUBLE acttime;
  DOUBLE acttimefac;  /* time factor from actual curve                    */
  DOUBLE acttimefacn; /* time factor at time (n)                          */
  NODE  *actnode;     /* actual node                                      */
  GVOL  *actgvol;

#ifdef DEBUG
  dstrc_enter("f3is_calset");
#endif

  fdyn    = alldyn[genprob.numff].fdyn;

/*-------------------------------------------- set element coordinates */
  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
    xyze[2][i]=ele->node[i]->x[2];
  }


/*---------------------------------------------------------------------*
  | position of the different solutions:                                |
  | node->sol_incement: solution history used for calculations          |
  |       sol_increment[ipos][i]: solution at some time level           |
  |       ipos-flags:  ipos.velnp ... solution at time (n+1)            |
  |                    ipos.veln  ... solution at time (n)              |
  |                    ipos.velnm ... solution at time (n-1)            |
  |                    ipos.gridv ... mesh velocity in actual time step |
  |                    ipos.convn ... convective velocity at time (n)   |
  |                    ipos.convnp... convective velocity at time (n+1) |
  |                    ipos.hist  ... sol. data for rhs                 |
  *---------------------------------------------------------------------*/


  for(i=0;i<ele->numnp;i++) /* loop nodes */
  {
    actnode=ele->node[i];
/*------------------------------------- set element velocities (n+1) ---*/
    evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
    evelng[2][i]=actnode->sol_increment.a.da[ipos->velnp][2];
/*---------------------------------------- set vel. histories at (n) ---*/
    ehist[0][i] = actnode->sol_increment.a.da[ipos->hist][0];
    ehist[1][i] = actnode->sol_increment.a.da[ipos->hist][1];
    ehist[2][i] = actnode->sol_increment.a.da[ipos->hist][2];
  }

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];
    if (actnode->sol_increment.sdim <= PREDOF)
      break;
/*---------------------------------------------- set pressures (n+1) ---*/
    epren[i]   =actnode->sol_increment.a.da[ipos->velnp][PREDOF];
  }

/*------------------------------------------------ check for dead load */
  actgvol = ele->g.gvol;
  if (actgvol->neum!=NULL)
  {
    actcurve = actgvol->neum->curve-1;
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
    for (i=0;i<3;i++)
    {
      if (actgvol->neum->neum_onoff.a.iv[i]==0)
      {
	edeadng[i] = ZERO;
      }
      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
	edeadng[i] = actgvol->neum->neum_val.a.dv[i]*acttimefac;
	(*hasext)++;
      }
    }
  }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         genk 02/04

get the element velocities and the pressure at different times

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the
      transformation in every time step

</pre>
\param   *ele       ELEMENT         (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\param  **ehist     DOUBLE          (o)    ele history data
\param  **evelng    DOUBLE          (o)    ele vels at time n+g
\param  **ealecovng DOUBLE          (o)    ALE-convective vels at time n+g
\param  **egridv    DOUBLE          (o)    element grid velocity
\param   *epren     DOUBLE          (o)    ele pres at time n
\param   *edeadng   DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *ipos                      (i)    node array positions
\param   *hasext    INT             (o)    flag for external loads
\return void

------------------------------------------------------------------------*/
void f3is_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **ehist,
  DOUBLE         **evelng,
  DOUBLE         **ealecovng,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos,
  INT             *hasext
  )
{
  INT    i;
  INT    actcurve;    /* actual time curve                                */
  DOUBLE acttimefac;  /* time factor from actual curve                    */
  DOUBLE acttimefacn; /* time factor at time (n)                          */
  DOUBLE acttime;
  NODE  *actnode;     /* actual node                                      */
  GVOL  *actgvol;

#ifdef DEBUG
  dstrc_enter("f3is_calseta");
#endif

  fdyn    = alldyn[genprob.numff].fdyn;

/*-------------------------------------------- set element coordinates */
  f3_alecoor(ele,xyze);


 /*---------------------------------------------------------------------*
  | position of the different solutions:                                |
  | node->sol_incement: solution history used for calculations          |
  |       sol_increment[0][i]: solution at (n-1)                        |
  |	  sol_increment[1][i]: solution at (n)                          |
  |	  sol_increment[2][i]: solution at (n+g)                        |
  |	  sol_increment[3][i]: solution at (n+1)                        |
  *---------------------------------------------------------------------*/


  for(i=0;i<ele->numnp;i++) /* loop nodes */
  {
    actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]   =actnode->sol_increment.a.da[ipos->velnp][0];
    evelng[1][i]   =actnode->sol_increment.a.da[ipos->velnp][1];
    evelng[2][i]   =actnode->sol_increment.a.da[ipos->velnp][2];
    ealecovng[0][i]=actnode->sol_increment.a.da[ipos->convnp][0];
    ealecovng[1][i]=actnode->sol_increment.a.da[ipos->convnp][1];
    ealecovng[2][i]=actnode->sol_increment.a.da[ipos->convnp][2];
    egridv[0][i]   =actnode->sol_increment.a.da[ipos->gridv][0];
    egridv[1][i]   =actnode->sol_increment.a.da[ipos->gridv][1];
    egridv[2][i]   =actnode->sol_increment.a.da[ipos->gridv][2];
    ehist[0][i] = actnode->sol_increment.a.da[ipos->hist][0];
    ehist[1][i] = actnode->sol_increment.a.da[ipos->hist][1];
    ehist[2][i] = actnode->sol_increment.a.da[ipos->hist][2];
  }

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];
    if (actnode->sol_increment.sdim <= PREDOF)
      break;
/*---------------------------------------------- set pressures (n+1) ---*/
    epren[i]   =actnode->sol_increment.a.da[ipos->velnp][PREDOF];
  }

/*------------------------------------------------ check for dead load */
  actgvol = ele->g.gvol;
  if (actgvol->neum!=NULL)
  {
    if (actgvol->neum->neum_type==neum_LAS)
    {
      actcurve = actgvol->neum->curve-1;
      if (actcurve<0)
	dserror("No Time curve given for neum_LAS!\n");
      acttime = fdyn->acttime;
      dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
      acttime=fdyn->acttime-fdyn->dta;
      dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
      edeadng[0] = actgvol->neum->neum_val.a.dv[0]*acttimefac;
      edeadng[1] = ZERO;
      edeadng[2] = actgvol->neum->neum_val.a.dv[2];
      (*hasext)++;
    }
    else
    {
      actcurve = actgvol->neum->curve-1;
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
      for (i=0;i<3;i++)
      {
	if (actgvol->neum->neum_onoff.a.iv[i]==0)
	{
	  edeadng[i] = ZERO;
	}
	if (actgvol->neum->neum_type==neum_dead  &&
	    actgvol->neum->neum_onoff.a.iv[i]!=0)
	{
	  edeadng[i] = actgvol->neum->neum_val.a.dv[i]*acttimefac;
	  (*hasext)++;
	}
      }
    }
  }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}

#endif
