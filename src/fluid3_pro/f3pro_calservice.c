/*!----------------------------------------------------------------------
\file
\brief service routines for fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID3_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3_PRO
#include "../headers/standardtypes.h"
#include "fluid3pro_prototypes.h"
#include "fluid3pro.h"
#include "../fluid3/fluid3_prototypes.h"

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

/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         basol 11/02

get the element velocities and the pressure at different times

</pre>
\param   *elevel      ELEMENT	   (i)    actual element for velocity
\param   *elepre      ELEMENT	   (i)    actual element for pressure
\param   **xyze       DOUBLE	   (o)    coordinates of the element
\param   **eveln      DOUBLE	   (o)    ele velocities at time n
\param   *epren       DOUBLE	   (o)    ele pressures at time n
\param   *ipos                     (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f3pro_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **evelnm,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *eprenm,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos
  )
{
  INT i, velnm, veln, velnp, hist;
  INT numpdof=0;
  NODE *actnode;                /* actual node for element */
  GVOL *actgvol;

#ifdef DEBUG
  dstrc_enter("f3pro_calset");
#endif

  velnm = ipos->velnm;
  veln  = ipos->veln;
  velnp = ipos->velnp;
  hist  = ipos->hist;

  /*--------------------------------------------- set element coordinates */
  for (i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
    xyze[2][i]=ele->node[i]->x[2];
  }

  /*---------------------------------------------------------------------*
   | position of the different solutions:                                |
   | node->sol_incement: solution history used for calculations          |
   |       sol_increment[ipos][i]: solution at some time level           |
   |       ipos->velnp .. solution at time n+1                           |
   |       ipos->veln  .. solution at time n                             |
   |       ipos->velnm .. solution at time n-1                           |
   *---------------------------------------------------------------------*/

  /* -> computation of time forces -------------------
     -> velocities and pressure at (n) are needed ----*/

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];

    /*----------------------------------- set element velocities at (n-1) */
    evelnm[0][i]=actnode->sol_increment.a.da[velnm][0];
    evelnm[1][i]=actnode->sol_increment.a.da[velnm][1];
    evelnm[2][i]=actnode->sol_increment.a.da[velnm][2];

    /*------------------------------------- set element velocities at (n) */
    eveln[0][i]=actnode->sol_increment.a.da[veln][0];
    eveln[1][i]=actnode->sol_increment.a.da[veln][1];
    eveln[2][i]=actnode->sol_increment.a.da[veln][2];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];
    evelng[2][i]=actnode->sol_increment.a.da[velnp][2];

    /*--------------------------------------- set vel. histories at (n) ---*/
    evhist[0][i] = actnode->sol_increment.a.da[hist][0];
    evhist[1][i] = actnode->sol_increment.a.da[hist][1];
    evhist[2][i] = actnode->sol_increment.a.da[hist][2];
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*---------------------------------------------- set pressures (n+1) ---*/
  if ((numpdof==-1) || (numpdof==-2))
  {
    ELEMENT* pele;
    pele = ele->e.f3pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i]  = pele->node[i]->sol_increment.a.da[1][0];
      eprenm[i] = pele->node[i]->sol_increment.a.da[2][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f3pro->press[i];
      eprenm[i]  = ele->e.f3pro->pressm[i];
    }
  }

  edeadng[0] = 0;
  edeadng[1] = 0;
  edeadng[2] = 0;

  /*------------------------------------------------ check for dead load */
  actgvol = ele->g.gvol;
  if (actgvol->neum!=NULL)
  {
    INT    actcurve;    /* actual time curve                                */
    DOUBLE acttime;
    DOUBLE acttimefac;  /* time factor from actual curve                    */
    DOUBLE acttimefacn; /* time factor at time (n)                          */

    actcurve = actgvol->neum->curve-1;
    if (actcurve<0)
    {
      acttimefac  = 1.;
      acttimefacn = 1.;
    }
    else
    {
      FLUID_DYNAMIC *fdyn;
      fdyn = alldyn[genprob.numff].fdyn;
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
      if (actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
	if (actgvol->neum->neum_type!=neum_dead)
	  dserror("unsupported neumann type %d",actgvol->neum->neum_type);
	edeadng[i] = actgvol->neum->neum_val.a.dv[i]*acttimefac;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of f3pro_calset */


/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         basol 11/02

get the element velocities and the pressure at different times

</pre>
\param   *elevel      ELEMENT	   (i)    actual element for velocity
\param   *elepre      ELEMENT	   (i)    actual element for pressure
\param   **xyze       DOUBLE	   (o)    coordinates of the element
\param   **eveln      DOUBLE	   (o)    ele velocities at time n
\param   *epren       DOUBLE	   (o)    ele pressures at time n
\param   *ipos                     (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f3pro_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **evelnm,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *eprenm,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos
  )
{
  INT i, velnm, veln, velnp, hist;
  INT numpdof=0;
  NODE *actnode;                /* actual node for element */
  GVOL *actgvol;

#ifdef DEBUG
  dstrc_enter("f3pro_calseta");
#endif

  velnm = ipos->velnm;
  veln  = ipos->veln;
  velnp = ipos->velnp;
  hist  = ipos->hist;

  /*--------------------------------------------- set element coordinates */
  f3_alecoor(ele,xyze);


  /*---------------------------------------------------------------------*
   | position of the different solutions:                                |
   | node->sol_incement: solution history used for calculations          |
   |       sol_increment[ipos][i]: solution at some time level           |
   |       ipos->velnp .. solution at time n+1                           |
   |       ipos->veln  .. solution at time n                             |
   |       ipos->velnm .. solution at time n-1                           |
   *---------------------------------------------------------------------*/

  /* -> computation of time forces -------------------
     -> velocities and pressure at (n) are needed ----*/

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];

    /*----------------------------------- set element velocities at (n-1) */
    evelnm[0][i]=actnode->sol_increment.a.da[velnm][0];
    evelnm[1][i]=actnode->sol_increment.a.da[velnm][1];
    evelnm[2][i]=actnode->sol_increment.a.da[velnm][2];

    /*------------------------------------- set element velocities at (n) */
    eveln[0][i]=actnode->sol_increment.a.da[veln][0];
    eveln[1][i]=actnode->sol_increment.a.da[veln][1];
    eveln[2][i]=actnode->sol_increment.a.da[veln][2];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];
    evelng[2][i]=actnode->sol_increment.a.da[velnp][2];

    /*--------------------------------------- set vel. histories at (n) ---*/
    evhist[0][i] = actnode->sol_increment.a.da[hist][0];
    evhist[1][i] = actnode->sol_increment.a.da[hist][1];
    evhist[2][i] = actnode->sol_increment.a.da[hist][2];

/*     ealecovng[0][i]=actnode->sol_increment.a.da[ipos->convnp][0]; */
/*     ealecovng[1][i]=actnode->sol_increment.a.da[ipos->convnp][1]; */
/*     ealecovng[2][i]=actnode->sol_increment.a.da[ipos->convnp][2]; */

    egridv[0][i]   =actnode->sol_increment.a.da[ipos->gridv][0];
    egridv[1][i]   =actnode->sol_increment.a.da[ipos->gridv][1];
    egridv[2][i]   =actnode->sol_increment.a.da[ipos->gridv][2];
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*---------------------------------------------- set pressures (n+1) ---*/
  if ((numpdof==-1) || (numpdof==-2))
  {
    ELEMENT* pele;
    pele = ele->e.f3pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i]  = pele->node[i]->sol_increment.a.da[1][0];
      eprenm[i] = pele->node[i]->sol_increment.a.da[2][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f3pro->press[i];
      eprenm[i]  = ele->e.f3pro->pressm[i];
    }
  }

  edeadng[0] = 0;
  edeadng[1] = 0;
  edeadng[2] = 0;

  /*------------------------------------------------ check for dead load */
  actgvol = ele->g.gvol;
  if (actgvol->neum!=NULL)
  {
    INT    actcurve;    /* actual time curve                                */
    DOUBLE acttime;
    DOUBLE acttimefac;  /* time factor from actual curve                    */
    DOUBLE acttimefacn; /* time factor at time (n)                          */

    actcurve = actgvol->neum->curve-1;
    if (actcurve<0)
    {
      acttimefac  = 1.;
      acttimefacn = 1.;
    }
    else
    {
      FLUID_DYNAMIC *fdyn;
      fdyn = alldyn[genprob.numff].fdyn;
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
      if (actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
	if (actgvol->neum->neum_type!=neum_dead)
	  dserror("unsupported neumann type %d",actgvol->neum->neum_type);
	edeadng[i] = actgvol->neum->neum_val.a.dv[i]*acttimefac;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of f3pro_calset */



/*----------------------------------------------------------------------*/
/*!
  \brief calculate shape functions and derivatives of the pressure

  \param pfunct       (o) shape function array
  \param pderiv       (o) shape function derivative array
  \param r            (i) natural coordinate
  \param s            (i) natural coordinate
  \param t            (i) natural coordinate
  \param dm           (i) discontinuous mode
  \param numpdof      (i) number of pressure dofs

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void f3pro_phex(DOUBLE* pfunct, DOUBLE** pderiv, DOUBLE r, DOUBLE s, DOUBLE t, DISMODE dm, INT* numpdof)
{
#ifdef DEBUG
  dstrc_enter("f3pro_prec");
#endif

  switch (dm)
  {
  case dm_q2pm1:
    /* discontinuous pressure on quadratic velocity */
    /* 27/4 */

    *numpdof = 4;

    pfunct[0]=1;
    pfunct[1]=r;
    pfunct[2]=s;
    pfunct[3]=t;

    pderiv[0][0]= 0;
    pderiv[1][0]= 0;
    pderiv[2][0]= 0;

    pderiv[0][1]= 1;
    pderiv[1][1]= 0;
    pderiv[2][1]= 0;

    pderiv[0][2]= 0;
    pderiv[1][2]= 1;
    pderiv[2][2]= 0;

    pderiv[0][3]= 0;
    pderiv[1][3]= 0;
    pderiv[2][3]= 1;
    break;
  case dm_q1p0:
    /* constant pressure */
    *numpdof = 1;

    pfunct[0]=1;

    pderiv[0][0]= 0;
    pderiv[1][0]= 0;
    pderiv[2][0]= 0;

    break;
  default:
    dserror("unsupported discretisation mode %d", dm);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
