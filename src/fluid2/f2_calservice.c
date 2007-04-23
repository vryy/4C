/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.def
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
void f2_calset(
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
dstrc_enter("f2_calset");
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
for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */
   evelng[0][i]=actnode->sol_increment.a.da[nodesol][0];
   evelng[1][i]=actnode->sol_increment.a.da[nodesol][1];
/*---------------------------------------------- set pressures (n+1) ---*/
   epren[i]   =actnode->sol_increment.a.da[nodesol][PREDOF];

/*if (fdyn->turbu > 2)*/ /* for laminar calculations only*/
/*{*/
/*--------------------------------------- set vel. histories at (n) ---*/
   evhist[0][i] = actnode->sol_increment.a.da[nodehist][0];
   evhist[1][i] = actnode->sol_increment.a.da[nodehist][1];
/*}*/
/*------------------------------------ set element velocities at (n) ---*/
   eveln[0][i]=actnode->sol_increment.a.da[ndsolold][0];
   eveln[1][i]=actnode->sol_increment.a.da[ndsolold][1];
/*   Kleine Ersetzung, fuer das Oseen-Problem!!! */
/*   evelng[0][i]=actnode->sol.a.da[0][0] * exp(-8.0*PI*PI*fdyn->acttime*0.01);
   evelng[1][i]=actnode->sol.a.da[0][1] * exp(-8.0*PI*PI*fdyn->acttime*0.01);
*/
/*-------------------------------------- set supported pressures (n+1) */
} /* end of loop over nodes of element */



/*------------------------------------------------ check for dead load */
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

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calset */



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
void f2_calseta(
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
dstrc_enter("f2_calseta");
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

for(i=0;i<ele->numnp;i++) /* loop nodes of element */
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
/*---------------------------------------------- set pressures (n+1) ---*/
   epren[i]   =actfnode->sol_increment.a.da[ipos->velnp][PREDOF];

   /* noch einmal kurz zurueck zur time rhs */
      ealecovn[0][i]=actfnode->sol_increment.a.da[ipos->convn][0];
      ealecovn[1][i]=actfnode->sol_increment.a.da[ipos->convn][1];
} /* end of loop over nodes of element */


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
if (ele->e.f2->fs_on>0 && fdyn->surftens!=0)
{
   if (fdyn->fsstnif!=0)
   {
      for (i=0;i<ele->numnp;i++)
         ekappan[i]=ele->e.f2->kappa_ND.a.da[i][0];
   }
   if (fdyn->fsstnii!=0)
   {
      for (i=0;i<ele->numnp;i++)
         ekappang[i]=ele->e.f2->kappa_ND.a.da[i][1];
   }
}

/*------------------------------------ height function at free surface */
if (ele->e.f2->fs_on==5)
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
if (ele->e.f2->fs_on)
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
} /* end of f2_calseta */


/*!---------------------------------------------------------------------
\brief set element coordinates during ALE calculations

<pre>                                                         genk 03/02

   nodal coordinates of actual element are evaluated. since nodes at the
   free surface have changing coordinates during the nonlin. iteration
   one has to treat them separately.
   Correct integration over moving domains has to be performed over the
   newest spatial configuration (in any case of theta)! Refere to me
   (Christiane Foerster) for details on geometric conservation and time
   integration on deforming domains!

</pre>

\warning Coordinates at new time level n+1 are returned! chfoe 03/05

\param   *ele       ELEMENT         (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\return void

------------------------------------------------------------------------*/
void f2_alecoor(
                  ELEMENT         *ele,
                  DOUBLE         **xyze
	       )
{
#ifdef D_FSI
INT i,j;
DOUBLE xy;
NODE  *actfnode;    /* actual fluid node                                */
NODE  *actanode;    /* actual ale node                                  */
GNODE *actfgnode;   /* actual fluid gnode                               */


#ifdef DEBUG
dstrc_enter("f2_alecoor");
#endif

/*-------------------------------------------- set element coordinates */
for(i=0;i<ele->numnp;i++)
{
   actfnode = ele->node[i];
   actfgnode = actfnode->gnode;
   actanode = actfgnode->mfcpnode[genprob.numaf];
   dsassert(actanode!=NULL,"ALE node not found (NULL-pointer)!\n");

   if(actfnode->xfs==NULL) /* no free surface */
   {
      for (j=0;j<2;j++)
      {
         xy     = actfnode->x[j];
         xyze[j][i] = xy + actanode->sol_mf.a.da[1][j];
      }
   }
   else /* free surface */
      for (j=0;j<2;j++)
         xyze[j][i] = actfnode->xfs[j];
}
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#else
dserror("FSI-functions not compiled in!\n");
#endif

return;
} /* end of f2_alecoor */


/*!---------------------------------------------------------------------
\brief set element coordinates during ALE calculations for relaxation
parameter of steepest descent method

<pre>                                                         chfoe 08/03


</pre>

\param   *ele       ELEMENT         (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\return void

------------------------------------------------------------------------*/
void f2_alecoor_sd(
                     ELEMENT         *ele,
                     DOUBLE         **xyze
	          )
{
#ifdef D_FSI
INT i,j;
DOUBLE theta;
NODE  *actfnode;    /* actual fluid node                                */
NODE  *actanode;    /* actual ale node                                  */
GNODE *actfgnode;   /* actual fluid gnode                               */


#ifdef DEBUG
dstrc_enter("f2_alecoor_sd");
#endif

fdyn  = alldyn[genprob.numff].fdyn;

theta = fdyn->theta;

/*-------------------------------------------- set element coordinates */

/* computation of relaxation parameter always goes from configuration 0
   to the one indicated by the actual incremental residuum vector g_i.
   See diss. Mok for details p. 119 ff */

for(i=0;i<ele->numnp;i++)
{
   actfnode = ele->node[i];
   actfgnode = actfnode->gnode;
   actanode = actfgnode->mfcpnode[genprob.numaf];
   /* if(actfnode->numdf==NUMDF) */
   {
      for (j=0;j<2;j++)
      {
         /* xyze[j][i] =  theta * ( actanode->sol_mf.a.da[2][j] ); */
         xyze[j][i] = ( actanode->sol_mf.a.da[2][j] );
      }
   }
#if 0
   else /* node on implicit free surface */
   {
   dserror("implicit free surface with steepest descent method not yet implementd");
   }
#endif
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#else
dserror("FSI-functions not compiled in!\n");
#endif

return;
} /* end of f2_alecoor_sd */


/*!---------------------------------------------------------------------
\brief routine interpolate fluid2 vector nodal quantaties
       at integration point

<pre>                                                         genk 04/02

</pre>

\param   *vecint   DOUBLE        (o)   values at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evec     DOUBLE        (i)   values at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_veci(
               DOUBLE  *vecint,
               DOUBLE  *funct,
               DOUBLE **evec,
               INT      iel
	    )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_veci");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   vecint[i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      vecint[i] += funct[j]*evec[i][j];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_veci */

/*!---------------------------------------------------------------------
\brief routine to interpolate nodal values at integration point

<pre>                                                         genk 04/02

</pre>
\param   *vecint   DOUBLE        (o)   values at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evec     DOUBLE        (i)   values at element nodes
\param    ngnode   INT           (i)   number of nodes on this edge
\param   *iedgnod  INT           (i)   edge node numbers
\return void

------------------------------------------------------------------------*/
void f2_edgeveci(
                  DOUBLE  *vecint,
                  DOUBLE  *funct,
                  DOUBLE **evec,
                  INT      ngnode,
                  INT     *iedgnod
               )
{
INT     i,j;
INT     node;

#ifdef DEBUG
dstrc_enter("f2_edgeveci");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   vecint[i]=ZERO;
   for (j=0;j<ngnode;j++) /* loop over all nodes j of the edge */
   {
      node=iedgnod[j];
      vecint[i] += funct[j]*evec[i][node];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_edgeveci */

/*!---------------------------------------------------------------------
\brief routine to calculate scalar values at integration point

<pre>                                                         genk 04/02

</pre>

\param  *funct     DOUBLE        (i)   shape functions
\param  *escal     DOUBLE        (i)   values at element nodes
\param   iel       INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
DOUBLE f2_scali(
               DOUBLE  *funct,
               DOUBLE  *escal,
               INT      iel
	     )
{
INT     j;
DOUBLE  scalint=ZERO;

#ifdef DEBUG
dstrc_enter("f2_scali");
#endif

for (j=0;j<iel;j++)
{
   scalint += funct[j] * escal[j];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return (scalint);
} /* end of f2_scali */

/*!---------------------------------------------------------------------
\brief routine to calculate curvature at integration point

<pre>                                                         genk 02/03

</pre>

\param  *funct     DOUBLE        (i)   shape functions
\param  *escal     DOUBLE        (i)   nodal kappa
\param   iedgnod   INT           (i)   local node numbers of actual line
\param   ngnode	   INT           (i)   number of nodes on actual line
\return DOUBLE

------------------------------------------------------------------------*/
DOUBLE f2_edgescali(
                    DOUBLE  *funct,
                    DOUBLE  *escal,
                    INT     *iedgnod,
                    INT      ngnode
	           )
{
INT     j;
DOUBLE  scalint= ZERO;

#ifdef DEBUG
dstrc_enter("f2_edgescali");
#endif

for (j=0;j<ngnode;j++)
{
   scalint += funct[j] * escal[iedgnod[j]];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return (scalint);
} /* end of f2_edgescali */



/*!---------------------------------------------------------------------
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the derivatives of the velocity w.r.t x/y are calculated

</pre>

\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_vder(
               DOUBLE **vderxy,
               DOUBLE **derxy,
               DOUBLE **evel,
               INT      iel
	    )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_vder");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   vderxy[0][i]=ZERO;
   vderxy[1][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_vder */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the 2nd derivatives of the velocity
w.r.t x/y are calculated

</pre>

\param  **vderxy2  DOUBLE        (o)   2nd velocity derivativs
\param  **derxy2   DOUBLE        (i)   2nd global derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_vder2(
               DOUBLE **vderxy2,
               DOUBLE **derxy2,
               DOUBLE **evel,
               INT      iel
             )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_vder2");
#endif

for (i=0;i<3;i++)
{
   vderxy2[0][i]=ZERO;
   vderxy2[1][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_vder2 */

/*!---------------------------------------------------------------------
\brief routine to calculate pressure derivatives at integration point

<pre>                                                         genk 04/02

In this routine derivatives of the pressure w.r.t x/y are calculated

</pre>

\param   *pderxy   DOUBLE        (o)   pressure derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param   *epre     DOUBLE        (i)   pressure at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_pder(
               DOUBLE  *pderxy,
               DOUBLE **derxy,
               DOUBLE  *epre,
               INT      iel
	    )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_pder");
#endif

for (i=0;i<2;i++) /* loop over directions i */
{
   pderxy[i] =  ZERO;
   for (j=0;j<iel;j++) /* loop over nodes j of the element */
   {
      pderxy[i] += derxy[i][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_pder */

/*!---------------------------------------------------------------------
\brief routine to calculate height function derivatives at integration
       point

<pre>                                                         genk 05/03

In this routine derivative of the height function w.r.t x is calculated

</pre>
\param  **derxy      DOUBLE        (i)   globael derivatives
\param   *ephi       DOUBLE        (i)   height funct at edge nodes
\param    iel        INT           (i)   number of nodes in this element
\param  *iedgnod     INT           (i)    edge node numbers
\return DOUBLE phiderx

------------------------------------------------------------------------*/
DOUBLE f2_phider(
                  DOUBLE **derxy,
                  DOUBLE  *ephi,
                  INT      ngnode,
                  INT     *iedgnod
	        )
{
INT    j;
INT    node;
DOUBLE phiderx=ZERO;

#ifdef DEBUG
dstrc_enter("f2_phider");
#endif

for (j=0;j<ngnode;j++)
{
   node=iedgnod[j];
   phiderx += derxy[0][j]*ephi[node];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return(phiderx);
} /* end of f2_phider */

/*!---------------------------------------------------------------------
\brief convective velocities

<pre>                                                         genk 04/02

in this routine the convective velocity is calculated at the
integration point:
 u * grad(u)
 e.g. 2D: COVx = Ux*Ux,x + Uy*Ux,y

</pre>

\param  **vderxy   DOUBLE        (i)   velocity derivativs
\param   *velint   DOUBLE        (i)   velocity at integration point
\param   *covint   DOUBLE        (o)   convective velocity at INT point
\return void

------------------------------------------------------------------------*/
void f2_covi(
               DOUBLE **vderxy,
               DOUBLE  *velint,
               DOUBLE  *covint
	    )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_covi");
#endif

for (i=0;i<2;i++)
{
   covint[i]=ZERO;
   for (j=0;j<2;j++)
   {
      covint[i] +=velint[j]*vderxy[i][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_covi */

/*!---------------------------------------------------------------------
\brief permutation of element force vector

<pre>                                                         genk 04/02

routine to rearrange the entries of the element force vector
this is necessary since we would like to use the existing assembly
routines for the RHS
hence a splitting of vel- and pre dof in the element force vector
is not possible any more!!!!


</pre>
\param   *eforce  DOUBLE        (i/o) element force vector
\param  **tmp     DOUBLE        (i)   working array
\param    iel     DOUBLE        (i)   number of nodes in this ele
\return void

------------------------------------------------------------------------*/
void f2_permeforce(
                     DOUBLE   *eforce,
                     DOUBLE  **tmp,
                     INT       iel
	          )
{
INT i,irow;
INT nvdof;      /* number of vel dofs                                   */
INT totdof;     /* total number of dofs                                 */

#ifdef DEBUG
dstrc_enter("f2_permeforce");
#endif

/*----------------------------------------------------- set some values */
nvdof  = NUM_F2_VELDOF*iel;
totdof = (NUM_F2_VELDOF+1)*iel;

/*-------------------------------------------------- rearrange vel-dofs */
irow = 0;
for (i=0;i<nvdof;i+=2)
{
   tmp[irow][0]   = eforce[i];
   tmp[irow+1][0] = eforce[i+1];
   irow += 3;
}

/*-------------------------------------------------- rearrange pre-dofs */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   tmp[irow][0] = eforce[i];
   irow += 3;
}

/*------------------------------------------------- copy back to eforce */
for (i=0;i<totdof;i++)
{
   eforce[i] = tmp[i][0];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_permeforce */

/*!---------------------------------------------------------------------
\brief permutation of element force vector for implicit free surface

<pre>                                                         genk 02/03

routine to rearrange the entries of the element force vector
this is necessary since we would like to use the existing assembly
routines for the RHS
hence a splitting of vel- and pre dof in the element force vector
is not possible any more!!!!


</pre>

\param   *eforce   DOUBLE        (i/o) element force vector
\param  **tmp      DOUBLE        (i)   working array
\param   *ele      ELEMENT       (i)   actual element
\return void

------------------------------------------------------------------------*/
void f2_permeforce_ifs(
                        DOUBLE   *eforce,
                        DOUBLE  **tmp,
                        ELEMENT  *ele
	              )
{
INT i,irow;
INT iel;
INT nvdof;      /* number of vel dofs                                   */
INT nvpdof;
INT rdist;
INT posr;

#ifdef DEBUG
dstrc_enter("f2_permeforce");
#endif

/*----------------------------------------------------- set some values */
iel    = ele->numnp;
nvdof  = NUM_F2_VELDOF*iel;
nvpdof = NUMDOF_FLUID2*iel;

irow=0;
rdist=nvdof;
for(i=0;i<iel;i++)
{
   posr=NUM_F2_VELDOF*i;
   switch(ele->node[i]->numdf)
   {
   case 3:
      tmp[irow][0]   = eforce[posr];
      tmp[irow+1][0] = eforce[posr+1];
      tmp[irow+2][0] = eforce[posr+rdist];
      irow+=3;
      rdist--;
   break;
   case 4:
      tmp[irow][0]   = eforce[posr];
      tmp[irow+1][0] = eforce[posr+1];
      tmp[irow+2][0] = eforce[posr+rdist];
      tmp[irow+3][0] = eforce[posr+rdist+iel];
      irow+=4;
      rdist--;
   break;
   case 5:
      tmp[irow][0]   = eforce[posr];
      tmp[irow+1][0] = eforce[posr+1];
      tmp[irow+2][0] = eforce[posr+rdist];
      tmp[irow+3][0] = eforce[posr+nvpdof];
      tmp[irow+4][0] = eforce[posr+nvpdof+1];
      irow+=5;
      rdist--;
   break;
   }
}

/*------------------------------------------------- copy back to eforce */
for (i=0;i<irow;i++)
{
   eforce[i] = tmp[i][0];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_permeforce_ifs */


/*!---------------------------------------------------------------------
\brief permutation of element stiffness matrix

<pre>                                                         genk 04/02

routine to add galerkin and stabilisation parts of the elment
stiffness matrix and to rearrange its entries!
this is necessary since we would like to use the existing assembly
routines for the stiffness matrix
hence a splitting of vel- and pre dofs is not possible any more!!!!

</pre>
\param  **estif   DOUBLE         (i/o) ele stiffnes matrix
\param  **emass   DOUBLE         (i)   ele mass matrix
\param  **tmp     DOUBLE         (-)   working array
\param	 *ele     ELEMENT        (i)   actual element
\return void

------------------------------------------------------------------------*/
void f2_permestif(
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE         **tmp,
                  ELEMENT         *ele
                 )
{
INT    i,j,icol,irow;     /* simply some counters                       */
INT    iel;               /* number of element nodes                    */
INT    nvdof;             /* number of vel dofs                         */
INT    npdof;             /* number of pre dofs                         */
INT    totdof;            /* total number of dofs                       */
DOUBLE thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG
dstrc_enter("f2_permestif");
#endif

/*----------------------------------------------------- set some values */
fdyn   = alldyn[genprob.numff].fdyn;
iel    = ele->numnp;
nvdof  = NUM_F2_VELDOF*iel;
npdof  = iel;
totdof = NUMDOF_FLUID2*iel;
thsl   = fdyn->thsl;

/*--------------------------------------------- copy estif to tmp-array *
                           and mutlitply stiffniss matrix with THETA*DT */
for (i=0;i<totdof;i++)
{
   for (j=0;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
   } /* end of loop over j */
} /* end of loop over i */
/*------------------------------- add mass matrix for instationary case */
if (fdyn->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (fdyn->nis==0) */

/*------------------------------------------------------- rearrange Kvv */
irow = 0;
for (i=0;i<nvdof;i+=2)
{
   icol = 0;
   for (j=0;j<nvdof;j+=2)
   {
      estif[irow][icol]     = tmp[i][j];
      estif[irow+1][icol]   = tmp[i+1][j];
      estif[irow][icol+1]   = tmp[i][j+1];
      estif[irow+1][icol+1] = tmp[i+1][j+1];
      icol += 3;
   } /* end of loop over j */
   irow += 3;
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kvp */
irow = 0;
for (i=0;i<nvdof;i+=2)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      icol += 3;
   } /* end of loop over j */
   irow += 3;
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kpv */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=2)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      icol += 3;
   } /* end of loop over j */
   irow += 3;
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kpp */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 3;
   } /* end of loop over j */
   irow += 3;
} /* end of loop over i */

/*---------------------------------------------------------------------*/


#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_permestif */

/*!---------------------------------------------------------------------
\brief permutation of element stiffness matrix for implict free surface
       (via local lagrange & height function)

<pre>                                                         genk 01/03
                                modified for height function  genk 05/03

routine to add galerkin and stabilisation parts of the elment
stiffness matrix and to rearrange its entries!
this is necessary since we would like to use the existing assembly
routines for the stiffness matrix
hence a splitting of vel- pre-  and grid/height function dofs is
not possible any more!!!!

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array
\param	  iel	  INT		 (i)   number of nodes in ele
\return void

------------------------------------------------------------------------*/
void f2_permestif_ifs(
                        DOUBLE         **estif,
                        DOUBLE         **emass,
                        DOUBLE         **tmp,
                        ELEMENT         *ele
                     )
{
INT    i,j,icol,irow;     /* simply some counters                       */
INT    nvdof;             /* number of vel dofs                         */
INT    rdist,cdist;
INT    totdof;            /* total number of dofs                       */
INT    iel;               /* actuel number of element nodes             */
INT    nvpdof;            /* number of vel+pres dofs                    */
INT    posc,posr;         /* positions in matrix                        */
DOUBLE thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG
dstrc_enter("f2_permestif_ifs");
#endif

/*----------------------------------------------------- set some values */
fdyn   = alldyn[genprob.numff].fdyn;
iel    = ele->numnp;
nvdof  = NUM_F2_VELDOF*iel;
nvpdof = NUMDOF_FLUID2*iel;
thsl   = fdyn->thsl;

switch (ele->e.f2->fs_on)
{
case 2: case 6:
   totdof = (NUMDOF_FLUID2+NUM_F2_VELDOF)*iel;
break;
case 5:
   totdof = (NUMDOF_FLUID2+1)*iel;
break;
default:
   dserror("FLUID2 element parameter fs_on out of range!\n");
}

/*--------------------------------------------- copy estif to tmp-array *
                           and mutlitply stiffness matrix with THETA*DT
                           and add mass matrix                          */
for (i=0;i<totdof;i++)
{
   for (j=0;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl + emass[i][j];
   } /* end of loop over j */
} /* end of loop over i */

/*--------------------------------------------------------------------- */
icol=0;
cdist=nvdof;
for(i=0;i<iel;i++)
{
   irow=0;
   posc=NUM_F2_VELDOF*i;
   rdist=nvdof;
   switch(ele->node[i]->numdf)
   {
   case 3:
      for (j=0;j<iel;j++)
      {
         posr=NUM_F2_VELDOF*j;
         switch(ele->node[j]->numdf)
         {
         case 3:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            irow += 3;
            rdist--;
         break;
         case 4:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+3][icol]   = tmp[posr+rdist+iel][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist+iel][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist+iel][posc+cdist];
            irow += 4;
            rdist--;
         break;
         case 5:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+3][icol]   = tmp[posr+nvpdof][posc];
            estif[irow+3][icol+1] = tmp[posr+nvpdof][posc+1];
            estif[irow+3][icol+2] = tmp[posr+nvpdof][posc+cdist];
            estif[irow+4][icol]   = tmp[posr+nvpdof+1][posc];
            estif[irow+4][icol+1] = tmp[posr+nvpdof+1][posc+1];
            estif[irow+4][icol+2] = tmp[posr+nvpdof+1][posc+cdist];
            irow += 5;
            rdist--;
         break;
         default:
            dserror("numdf invalid!\n");
         }
      }
      icol+=3;
   break;
   case 4:
      for (j=0;j<iel;j++)
      {
         posr=NUM_F2_VELDOF*j;
         switch(ele->node[j]->numdf)
         {
         case 3:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow][icol+3]   = tmp[posr][posc+cdist+iel];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist+iel];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+2][icol+3] = tmp[posr+rdist][posc+cdist+iel];
            irow += 3;
            rdist--;
         break;
         case 4:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow][icol+3]   = tmp[posr][posc+cdist+iel];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist+iel];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+2][icol+3] = tmp[posr+rdist][posc+cdist+iel];
            estif[irow+3][icol]   = tmp[posr+rdist+iel][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist+iel][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist+iel][posc+cdist];
            estif[irow+3][icol+3] = tmp[posr+rdist+iel][posc+cdist+iel];
            irow += 4;
            rdist--;
         break;
         default:
            dserror("numdf invalid!\n");
         }
      }
      icol+=4;
   break;
   case 5:
      for (j=0;j<iel;j++)
      {
         posr=NUM_F2_VELDOF*j;
         switch(ele->node[j]->numdf)
         {
         case 3:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow][icol+3]   = tmp[posr][posc+nvpdof];
            estif[irow][icol+4]   = tmp[posr][posc+nvpdof+1];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+3] = tmp[posr+1][posc+nvpdof];
            estif[irow+1][icol+4] = tmp[posr+1][posc+nvdof+1];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+2][icol+3] = tmp[posr+rdist][posc+nvpdof];
            estif[irow+2][icol+4] = tmp[posr+rdist][posc+nvpdof+1];
            irow += 3;
            rdist--;
         break;
         case 5:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+cdist];
            estif[irow][icol+3]   = tmp[posr][posc+nvpdof];
            estif[irow][icol+4]   = tmp[posr][posc+nvpdof+1];
            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+3] = tmp[posr+1][posc+nvpdof];
            estif[irow+1][icol+4] = tmp[posr+1][posc+nvpdof+1];
            estif[irow+2][icol]   = tmp[posr+rdist][posc];
            estif[irow+2][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+2][icol+2] = tmp[posr+rdist][posc+cdist];
            estif[irow+2][icol+3] = tmp[posr+rdist][posc+nvpdof];
            estif[irow+2][icol+4] = tmp[posr+rdist][posc+nvpdof+1];
            estif[irow+3][icol]   = tmp[posr+nvpdof][posc];
            estif[irow+3][icol+1] = tmp[posr+nvpdof][posc+1];
            estif[irow+3][icol+2] = tmp[posr+nvpdof][posc+cdist];
            estif[irow+3][icol+3] = tmp[posr+nvpdof][posc+nvpdof];
            estif[irow+3][icol+4] = tmp[posr+nvpdof][posc+nvpdof+1];
            estif[irow+4][icol]   = tmp[posr+nvpdof+1][posc];
            estif[irow+4][icol+1] = tmp[posr+nvpdof+1][posc+1];
            estif[irow+4][icol+2] = tmp[posr+nvpdof+1][posc+cdist];
            estif[irow+4][icol+3] = tmp[posr+nvpdof+1][posc+nvpdof];
            estif[irow+4][icol+4] = tmp[posr+nvpdof+1][posc+nvpdof+1];
            irow += 5;
            rdist--;
            break;
            default:
            dserror("numdf invalid!\n");
         }
      }
      icol+=5;
   break;
   default:
      dserror("numdf invalid!\n");
   }
   cdist--;
}

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_permestif_ifs */

/*!---------------------------------------------------------------------
\brief get edgnodes for element line

<pre>                                                         genk 01/03


</pre>

\param   iegnod    INT        (o)   edge nodes
\param  *ele       ELEMENT    (i)   actual element
\param   line      DOUBLE     (i)   actual line number
\param   init      INT        (i)   flag
\return void

------------------------------------------------------------------------*/
void f2_iedg(
               INT     *iegnod,
               ELEMENT *ele,
               INT      line,
               INT      init
	     )
{
INT i;
static INT iegq[4][4][2];
static INT iegt[4][4][2];

#ifdef DEBUG
dstrc_enter("f2_iedg");
#endif

/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          *
 * meaning of indices iegq[i][j][k]:
   i: node on edge
   j: line number of element
   k: index for linear of quadratic element
 *---------------------------------------------------------------------*/
if (init==1)
{
   /*-------------------------------------------- egde nodes for quad4 */
   iegq[0][0][0] = 0;
   iegq[1][0][0] = 1;
   iegq[0][1][0] = 1;
   iegq[1][1][0] = 2;
   iegq[0][2][0] = 2;
   iegq[1][2][0] = 3;
   iegq[0][3][0] = 3;
   iegq[1][3][0] = 0;
   /*----------------------------------- egde nodes for quad8 and quad9 */
   iegq[0][0][1] = 0;
   iegq[1][0][1] = 4;
   iegq[2][0][1] = 1;
   iegq[0][1][1] = 1;
   iegq[1][1][1] = 5;
   iegq[2][1][1] = 2;
   iegq[0][2][1] = 2;
   iegq[1][2][1] = 6;
   iegq[2][2][1] = 3;
   iegq[0][3][1] = 3;
   iegq[1][3][1] = 7;
   iegq[2][3][1] = 0;
   /*---------------------------------------------- egde nodes for tri3 */
   iegt[0][0][0] = 0;
   iegt[1][0][0] = 1;
   iegt[0][1][0] = 1;
   iegt[1][1][0] = 2;
   iegt[0][2][0] = 2;
   iegt[1][2][0] = 0;
   /*---------------------------------------------- egde nodes for tri6 */
   iegt[0][0][1] = 0;
   iegt[1][0][1] = 3;
   iegt[2][0][1] = 1;
   iegt[0][1][1] = 1;
   iegt[1][1][1] = 4;
   iegt[2][1][1] = 2;
   iegt[0][2][1] = 2;
   iegt[1][2][1] = 5;
   iegt[2][2][1] = 0;
}

/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
else if (init==0)
{
   switch(ele->distyp)
   {
   case quad4:
      for(i=0;i<2;i++) iegnod[i] = iegq[i][line][0];
   break;
   case quad8: case quad9:
      for(i=0;i<3;i++) iegnod[i] = iegq[i][line][1];
   break;
   case tri3:
      for(i=0;i<2;i++) iegnod[i] = iegt[i][line][0];
   break;
   case tri6:
      for(i=0;i<3;i++) iegnod[i] = iegt[i][line][1];
      dserror("iegnode for tri6 not tested yet\n");
   break;
   default:
      dserror("distyp unknown\n");
   } /*end switch(ele->distyp) */
}
else
   dserror("parameter 'init' out of range\n");

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f2_iedg */

#endif
/*! @} (documentation module close)*/
