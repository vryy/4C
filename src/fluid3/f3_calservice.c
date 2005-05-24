/*!----------------------------------------------------------------------
\file
\brief service routines for fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>


 *----------------------------------------------------------------------*/
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
\param   *edeadn   DOUBLE           (o)    ele dead load at n (selfweight)
\param   *edeadng  DOUBLE           (o)    ele dead load at n+g (selfweight)
\param   *ipos                      (i)    node array positions
\param   *hasext   INT              (o)    flag for external loads
\return void

------------------------------------------------------------------------*/
void f3_calset(
               ELEMENT         *ele,
               DOUBLE         **xyze,
               DOUBLE         **ehist,
               DOUBLE         **evelng,
               DOUBLE          *epren,
               DOUBLE          *edeadn,
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
dstrc_enter("f3_calset");
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
/*----------------------------------- set element velocities (n+gamma) */
   evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
   evelng[2][i]=actnode->sol_increment.a.da[ipos->velnp][2];
   ehist[0][i] = actnode->sol_increment.a.da[ipos->hist][0];
   ehist[1][i] = actnode->sol_increment.a.da[ipos->hist][1];
   ehist[2][i] = actnode->sol_increment.a.da[ipos->hist][2];
} /* end of loop over nodes */

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
         edeadn[i]  = ZERO;
         edeadng[i] = ZERO;
      }
      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i]  = actgvol->neum->neum_val.a.dv[i]*acttimefacn;
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
} /* end of f3_calset */

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
\param  **ealecovn  DOUBLE          (o)    ALE-convective vels at time n
\param  **ealecovng DOUBLE          (o)    ALE-convective vels at time n+g
\param  **egridv    DOUBLE          (o)    element grid velocity
\param   *epren     DOUBLE          (o)    ele pres at time n
\param   *edeadn    DOUBLE          (o)    ele dead load at n (selfweight)
\param   *edeadng   DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *ipos                      (i)    node array positions
\param   *hasext    INT             (o)    flag for external loads
\return void

------------------------------------------------------------------------*/
void f3_calseta(
                  ELEMENT         *ele,
                  DOUBLE         **xyze,
                  DOUBLE         **ehist,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovn,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE          *epren,
                  DOUBLE          *edeadn,
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
dstrc_enter("f3_calseta");
#endif

fdyn    = alldyn[genprob.numff].fdyn;

/*-------------------------------------------- set element coordinates */
f3_alecoor(ele,xyze);


/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
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
} /* end of loop over nodes */

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
      edeadn[0]  = actgvol->neum->neum_val.a.dv[0]*acttimefacn;
      edeadng[0] = actgvol->neum->neum_val.a.dv[0]*acttimefac;
      edeadn[1]  = ZERO;
      edeadng[1] = ZERO;
      edeadn[2]  = actgvol->neum->neum_val.a.dv[2];
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
            edeadn[i]  = ZERO;
            edeadng[i] = ZERO;
         }
         if (actgvol->neum->neum_type==neum_dead  &&
             actgvol->neum->neum_onoff.a.iv[i]!=0)
         {
            edeadn[i]  = actgvol->neum->neum_val.a.dv[i]*acttimefacn;
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
} /* end of f3_calseta */

/*!---------------------------------------------------------------------
\brief set element coordinates during ALE calculations

<pre>                                                         genk 02/04

   nodal coordinates of actual element are evaluated. since nodes at the
   free surface have changing coordinates during the nonlin. iteration
   one has to treat them separately.

</pre>
\param   *ele       ELEMENT         (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\return void

------------------------------------------------------------------------*/
void f3_alecoor(
                  ELEMENT         *ele,
                  DOUBLE         **xyze
	       )
{
#ifdef D_FSI
INT i,j;
DOUBLE omt,theta;
DOUBLE xy,xyn,xyng;
DOUBLE dt;
NODE  *actfnode;    /* actual fluid node                                */
NODE  *actanode;    /* actual ale node                                  */
GNODE *actfgnode;   /* actual fluid gnode                               */


#ifdef DEBUG
dstrc_enter("f3_alecoor");
#endif

fdyn    = alldyn[genprob.numff].fdyn;
theta = fdyn->theta;
omt   = ONE-theta;
dt = fdyn->dta;

/*-------------------------------------------- set element coordinates */
for(i=0;i<ele->numnp;i++)
{
   actfnode = ele->node[i];
   actfgnode = actfnode->gnode;
   actanode = actfgnode->mfcpnode[genprob.numaf];
   if(actfnode->xfs==NULL) /* no free surface */
   {
      for (j=0;j<3;j++)
      {
         xy     = actfnode->x[j];
         xyng   = xy + actanode->sol_mf.a.da[1][j];
         xyn    = xy + actanode->sol_mf.a.da[0][j];
         xyze[j][i] =  theta*(xyng)+omt*(xyn);
      }
   }
   else /* free surface */
   {
      for (j=0;j<3;j++)
      {
         xy     = actfnode->x[j];
         xyn    = xy + actanode->sol_mf.a.da[0][j];
         xyng   = actfnode->xfs[j];
         xyze[j][i] =  theta*(xyng)+omt*(xyn);
      }
   }
}

#else
dserror("FSI-functions not compiled in!\n");
#endif
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_alecoor */

/*!---------------------------------------------------------------------
\brief routine to calculate velocities at integration point

<pre>                                                         genk 05/02

</pre>
\param   *vecint   DOUBLE        (o)   vector at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **ecel     DOUBLE        (i)   vector at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f3_veci(
             DOUBLE  *vecint,
             DOUBLE  *funct,
             DOUBLE **evec,
             INT      iel
            )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_veci");
#endif

for (i=0;i<3;i++)
{
   vecint[i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vecint[i] += funct[j]*evec[i][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_veli */

/*!---------------------------------------------------------------------
\brief routine to calculate velocities at integration point on edge

<pre>                                                         genk 04/04

</pre>
\param   *velint   DOUBLE        (o)   velocities at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    ngnode   INT           (i)   number of nodes on this edge
\param   *iedgnod  INT           (i)   edge node numbers
\return void

------------------------------------------------------------------------*/
void f3_edgeveli(
                  DOUBLE  *velint,
                  DOUBLE  *funct,
                  DOUBLE **evel,
                  INT      ngnode,
                  INT     *iedgnod
               )
{
INT     i,j;
INT     node;

#ifdef DEBUG
dstrc_enter("f3_edgeveli");
#endif

for (i=0;i<3;i++) /* loop directions i */
{
   velint[i]=ZERO;
   for (j=0;j<ngnode;j++) /* loop over all nodes j of the edge */
   {
      node=iedgnod[j];
      velint[i] += funct[j]*evel[i][node];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_edgeveli */

/*!---------------------------------------------------------------------
\brief routine to calculate pressure at integration point

<pre>                                                         genk 05/02

</pre>
\param  *funct     DOUBLE        (i)   shape functions
\param  *epre      DOUBLE        (i)   pressure at element nodes
\param   iel       INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
DOUBLE f3_scali(
               DOUBLE  *funct,
               DOUBLE  *epre,
               INT      iel
            )
{
INT     j;
DOUBLE scalint=ZERO;

#ifdef DEBUG
dstrc_enter("f3_scali");
#endif

for (j=0;j<iel;j++)
{
   scalint += funct[j] * epre[j];
} /* end of loop over j */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return (scalint);
} /* end of f3_prei */

/*!---------------------------------------------------------------------
\brief routine to calculate height function at integration point

<pre>                                                         genk 05/03

</pre>
\param  *funct     DOUBLE        (i)   shape functions
\param  *ephi      DOUBLE        (i)   nodal phi
\param   iedgnod   INT           (i)   local node numbers of actual line
\param   ngnode	   INT           (i)   number of nodes on actual line
\return double

------------------------------------------------------------------------*/
DOUBLE f3_phii(
                  DOUBLE  *funct,
                  DOUBLE  *ephi,
                  INT     *iedgnod,
                  INT      ngnode
               )
{
INT     j;
DOUBLE  phiint = ZERO;

#ifdef DEBUG
dstrc_enter("f3_phii");
#endif

for (j=0;j<ngnode;j++)
{
   phiint += funct[j] * ephi[iedgnod[j]];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return (phiint);
} /* end of f3_phii */


/*!---------------------------------------------------------------------
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 05/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
vderxy[0][2] = Ux,z

</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f3_vder(
               DOUBLE **vderxy,
               DOUBLE **derxy,
               DOUBLE **evel,
               INT      iel
            )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_vder");
#endif

for (i=0;i<3;i++)
{
   vderxy[0][i]=ZERO;
   vderxy[1][i]=ZERO;
   vderxy[2][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
      vderxy[2][i] += derxy[i][j]*evel[2][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_vder */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the 2nd derivatives of the velocity
w.r.t x/y/z are calculated
   vderxy2[0][0] = Ux,xx
   vderxy2[0][3] = Ux,xy
   vderxy2[1][4] = Uy,xz
   vderxy2[2][5] = Uz,yz

</pre>
\param  **vderxy2  DOUBLE        (o)   2nd velocity derivativs
\param  **derxy2   DOUBLE        (i)   2nd global derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f3_vder2(
               DOUBLE **vderxy2,
               DOUBLE **derxy2,
               DOUBLE **evel,
               INT      iel
               )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_vder2");
#endif

for (i=0;i<6;i++)
{
   vderxy2[0][i]=ZERO;
   vderxy2[1][i]=ZERO;
   vderxy2[2][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
      vderxy2[2][i] += derxy2[i][j]*evel[2][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_vder2 */

/*!---------------------------------------------------------------------
\brief routine to calculate pressure derivatives at integration point

<pre>                                                         genk 05/02

In this routine derivatives of the pressure w.r.t x/y/z are calculated

</pre>
\param   *pderxy   DOUBLE        (o)   pressure derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param   *epre     DOUBLE        (i)   pressure at element nodes
\param    iel      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f3_pder(
               DOUBLE  *pderxy,
               DOUBLE **derxy,
               DOUBLE  *epre,
               INT      iel
            )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_pder");
#endif

for (i=0;i<3;i++)
{
   pderxy[i] =  ZERO;
   for (j=0;j<iel;j++)
   {
      pderxy[i] += derxy[i][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_pder */

/*!---------------------------------------------------------------------
\brief routine to calculate height function derivatives at integration
       point

<pre>                                                         genk 04/04

In this routine derivative of the height function w.r.t x is calculated

</pre>
\param  **derxy      DOUBLE        (i)   globael derivatives
\param   *ephi       DOUBLE        (i)   height funct at edge nodes
\param    iel        INT           (i)   number of nodes in this element
\param  *iedgnod     INT           (i)    edge node numbers
\return void

------------------------------------------------------------------------*/
void f3_phider(
                  DOUBLE  *phiderxy,
                  DOUBLE **derxy,
                  DOUBLE  *ephi,
                  INT      ngnode,
                  INT     *iedgnod
	        )
{
INT    i,j;
INT    node;

#ifdef DEBUG
dstrc_enter("f3_phider");
#endif

phiderxy[0]=phiderxy[1]=ZERO;

for (i=0;i<2;i++)
{
   for (j=0;j<ngnode;j++)
   {
      node=iedgnod[j];
      phiderxy[i] += derxy[0][j]*ephi[node];
   }
}
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_phider */

/*!---------------------------------------------------------------------
\brief convective velocities

<pre>                                                         genk 05/02

in this routine the convective velocity is calculated at the
integration point:
 u * grad(u)
 e.g. 3D: COVx = Ux*Ux,x + Uy*Ux,y + Uz*Ux,z

</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param   *velint   DOUBLE        (i)   velocity at integration point
\param   *covint   DOUBLE        (i)   convective velocity at INT point
\return void

------------------------------------------------------------------------*/
void f3_covi(
               DOUBLE **vderxy,
               DOUBLE  *velint,
               DOUBLE  *covint
            )
{
INT     i,j;
#ifdef DEBUG
dstrc_enter("f3_covi");
#endif

for (i=0;i<3;i++)
{
   covint[i]=ZERO;
   for (j=0;j<3;j++)
   {
      covint[i] +=velint[j]*vderxy[i][j];
   } /* end of loop over j */
} /* end of loop over i */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_covi */

/*!---------------------------------------------------------------------
\brief permutation of element force vector

<pre>                                                         genk 05/02

routine to rearrange the entries of the element force vector
this is necessary since we would like to use the existing assembly
routines for the RHS
hence a splitting of vel- and pre dof in the element force vector
is not possible any more!!!!


</pre>
\param   *eforce   DOUBLE        (i/o) element force vector
\param  **tmp      DOUBLE        (i)   working array
\param    iel      DOUBLE        (i)   number of nodes in this ele
\return void

------------------------------------------------------------------------*/
void f3_permeforce(
                  DOUBLE    *eforce,
                  DOUBLE   **tmp,
                  INT        iel
                  )
{
INT i,irow;
INT nvdof;      /* number of vel dofs                                   */
INT totdof;     /* total number of dofs                                 */

#ifdef DEBUG
dstrc_enter("f3_permeforce");
#endif

nvdof  = NUM_F3_VELDOF*iel;
totdof = (NUM_F3_VELDOF+1)*iel;

/*---------------------------------------------------- compute vel-dofs */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   tmp[irow][0]   = eforce[i];
   tmp[irow+1][0] = eforce[i+1];
   tmp[irow+2][0] = eforce[i+2];
   irow += 4;
} /* end of loop over i */

/*---------------------------------------------------- compute pre-dofs */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   tmp[irow][0] = eforce[i];
   irow += 4;
} /* end of loop over i */

/*------------------------------------------------- copy back to eforce */
for (i=0;i<totdof;i++)
{
   eforce[i] = tmp[i][0];
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_permeforce */

/*!---------------------------------------------------------------------
\brief permutation of element force vector for implicit free surface

<pre>                                                         genk 02/04

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
void f3_permeforce_ifs(
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
dstrc_enter("f3_permeforce_ifs");
#endif

/*----------------------------------------------------- set some values */
iel    = ele->numnp;
nvdof  = NUM_F3_VELDOF*iel;
nvpdof = NUMDOF_FLUID3*iel;

irow=0;
rdist=nvdof;
for(i=0;i<iel;i++)
{
   posr=NUM_F3_VELDOF*i;
   switch(ele->node[i]->numdf)
   {
   case 4:
      tmp[irow][0]   = eforce[posr];
      tmp[irow+1][0] = eforce[posr+1];
      tmp[irow+2][0] = eforce[posr+2];
      tmp[irow+3][0] = eforce[posr+rdist];
      irow+=4;
      rdist-=2;
   break;
   case 7:
      tmp[irow][0]   = eforce[posr];
      tmp[irow+1][0] = eforce[posr+1];
      tmp[irow+2][0] = eforce[posr+2];
      tmp[irow+3][0] = eforce[posr+rdist];
      tmp[irow+4][0] = eforce[posr+nvpdof];
      tmp[irow+5][0] = eforce[posr+nvpdof+1];
      tmp[irow+6][0] = eforce[posr+nvpdof+2];
      irow+=7;
      rdist-=2;
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
} /* end of f3_permeforce_ifs */


/*!---------------------------------------------------------------------
\brief permutation of element stiffness matrix

<pre>                                                         genk 05/02

routine to add galerkin and stabilisation parts of the elment
stiffness matrix and to rearrange its entries!
this is necessary since we would like to use the existing assembly
routines for the stiffness matrix
hence a splitting of vel- and pre dofs is not possible any more!!!!

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array
\param	  iel	  INT		 (i)   number of nodes in ele
\return void

------------------------------------------------------------------------*/
void f3_permestif(
		   DOUBLE         **estif,
		   DOUBLE         **emass,
		   DOUBLE         **tmp,
		   INT              iel
	          )
{
INT i,j,icol,irow;          /* simply some counters  	        	*/
INT nvdof;                  /* number of vel dofs                       */
INT npdof;                  /* number of pre dofs                       */
INT totdof;                 /* total number of dofs                     */
DOUBLE thsl;	            /* factor for LHS (THETA*DT)                */

#ifdef DEBUG
dstrc_enter("f3_permestif");
#endif

fdyn    = alldyn[genprob.numff].fdyn;

nvdof  = NUM_F3_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F3_VELDOF+1)*iel;
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

/*--------------------------------------------------------- compute Kvv */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]     = tmp[i][j];
      estif[irow+1][icol]   = tmp[i+1][j];
      estif[irow+2][icol]   = tmp[i+2][j];
      estif[irow][icol+1]   = tmp[i][j+1];
      estif[irow+1][icol+1] = tmp[i+1][j+1];
      estif[irow+2][icol+1] = tmp[i+2][j+1];
      estif[irow][icol+2]   = tmp[i][j+2];
      estif[irow+1][icol+2] = tmp[i+1][j+2];
      estif[irow+2][icol+2] = tmp[i+2][j+2];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kvp */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      estif[irow+2][icol] = tmp[i+2][j];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpv */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      estif[irow][icol+2] = tmp[i][j+2];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpp */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_permestif */

/*!---------------------------------------------------------------------
\brief permutation of element stiffness matrix for implict free surface
       (via local lagrange)

<pre>                                                         genk 01/03

routine to add galerkin and stabilisation parts of the elment
stiffness matrix and to rearrange its entries!
this is necessary since we would like to use the existing assembly
routines for the stiffness matrix
hence a splitting of vel- pre-  and grid/height function dofs is
not possible any more!!!!

</pre>
\param  **estif   DOUBLE            (i/o) ele stiffnes matrix
\param  **emass   DOUBLE            (i)   ele mass matrix
\param  **tmp     DOUBLE            (-)   working array
\param	 *ele     ELEMENT           (i)   actual element
\return void

------------------------------------------------------------------------*/
void f3_permestif_ifs(
                        DOUBLE         **estif,
                        DOUBLE         **emass,
                        DOUBLE         **tmp,
                        ELEMENT         *ele
                     )
{
INT    i,j,icol,irow;     /* simply some counters                       */
INT    nvdof;             /* number of vel dofs                         */
INT    rdist,cdist;
INT    totdof = 0;        /* total number of dofs                       */
INT    iel;               /* actuel number of element nodes             */
INT    nvpdof;            /* number of vel+pres dofs                    */
INT    posc,posr;         /* positions in matrix                        */
DOUBLE thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG
dstrc_enter("f3_permestif_ifs");
#endif

/*----------------------------------------------------- set some values */
fdyn    = alldyn[genprob.numff].fdyn;
iel    = ele->numnp;
nvdof  = NUM_F3_VELDOF*iel;
nvpdof = NUMDOF_FLUID3*iel;
thsl   = fdyn->thsl;

switch (ele->e.f3->fs_on)
{
case 2:
   totdof = (NUMDOF_FLUID3+NUM_F3_VELDOF)*iel;
break;
default:
   dserror("FLUID3 element parameter fs_on out of range!\n");
}

/*--------------------------------------------- copy estif to tmp-array *
                           and mutlitply stiffniss matrix with THETA*DT
                           and add mass matrix                          */
for (i=0;i<totdof;i++)
{
   for (j=0;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl + emass[i][j];
      estif[i][j]=ZERO;
   } /* end of loop over j */
} /* end of loop over i */

/*--------------------------------------------------------------------- */
icol=0;
cdist=nvdof;
for(i=0;i<iel;i++)
{
   irow=0;
   posc=NUM_F3_VELDOF*i;
   rdist=nvdof;
   switch(ele->node[i]->numdf)
   {
   case 4:
      for (j=0;j<iel;j++)
      {
         posr=NUM_F3_VELDOF*j;
         switch(ele->node[j]->numdf)
         {
         case 4:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+2];
            estif[irow][icol+3]   = tmp[posr][posc+cdist];

            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+2];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist];

            estif[irow+2][icol]   = tmp[posr+2][posc];
            estif[irow+2][icol+1] = tmp[posr+2][posc+1];
            estif[irow+2][icol+2] = tmp[posr+2][posc+2];
            estif[irow+2][icol+3] = tmp[posr+2][posc+cdist];

            estif[irow+3][icol]   = tmp[posr+rdist][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist][posc+2];
            estif[irow+3][icol+3] = tmp[posr+rdist][posc+cdist];
            irow += 4;
            rdist-=2;
         break;
         case 7:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+2];
            estif[irow][icol+3]   = tmp[posr][posc+cdist];

            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+2];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist];

            estif[irow+2][icol]   = tmp[posr+2][posc];
            estif[irow+2][icol+1] = tmp[posr+2][posc+1];
            estif[irow+2][icol+2] = tmp[posr+2][posc+2];
            estif[irow+2][icol+3] = tmp[posr+2][posc+cdist];

            estif[irow+3][icol]   = tmp[posr+rdist][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist][posc+2];
            estif[irow+3][icol+3] = tmp[posr+rdist][posc+cdist];

            estif[irow+4][icol]   = tmp[posr+nvpdof][posc];
            estif[irow+4][icol+1] = tmp[posr+nvpdof][posc+1];
            estif[irow+4][icol+2] = tmp[posr+nvpdof][posc+2];
            estif[irow+4][icol+3] = tmp[posr+nvpdof][posc+cdist];

            estif[irow+5][icol]   = tmp[posr+nvpdof+1][posc];
            estif[irow+5][icol+1] = tmp[posr+nvpdof+1][posc+1];
            estif[irow+5][icol+2] = tmp[posr+nvpdof+1][posc+2];
            estif[irow+5][icol+3] = tmp[posr+nvpdof+1][posc+cdist];

            estif[irow+6][icol]   = tmp[posr+nvpdof+2][posc];
            estif[irow+6][icol+1] = tmp[posr+nvpdof+2][posc+1];
            estif[irow+6][icol+2] = tmp[posr+nvpdof+2][posc+2];
            estif[irow+6][icol+3] = tmp[posr+nvpdof+2][posc+cdist];
            irow += 7;
            rdist-=2;
         break;
         default:
            dserror("numdf invalid!\n");
         }
      }
      icol+=4;
   break;
   case 7:
      for (j=0;j<iel;j++)
      {
         posr=NUM_F3_VELDOF*j;
         switch(ele->node[j]->numdf)
         {
         case 4:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+2];
            estif[irow][icol+3]   = tmp[posr][posc+cdist];
            estif[irow][icol+4]   = tmp[posr][posc+nvpdof];
            estif[irow][icol+5]   = tmp[posr][posc+nvpdof+1];
            estif[irow][icol+6]   = tmp[posr][posc+nvpdof+2];

            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+2];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+4] = tmp[posr+1][posc+nvpdof];
            estif[irow+1][icol+5] = tmp[posr+1][posc+nvpdof+1];
            estif[irow+1][icol+6] = tmp[posr+1][posc+nvpdof+2];

            estif[irow+2][icol]   = tmp[posr+2][posc];
            estif[irow+2][icol+1] = tmp[posr+2][posc+1];
            estif[irow+2][icol+2] = tmp[posr+2][posc+2];
            estif[irow+2][icol+3] = tmp[posr+2][posc+cdist];
            estif[irow+2][icol+4] = tmp[posr+2][posc+nvpdof];
            estif[irow+2][icol+5] = tmp[posr+2][posc+nvpdof+1];
            estif[irow+2][icol+6] = tmp[posr+2][posc+nvpdof+2];

            estif[irow+3][icol]   = tmp[posr+rdist][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist][posc+2];
            estif[irow+3][icol+3] = tmp[posr+rdist][posc+cdist];
            estif[irow+3][icol+4] = tmp[posr+rdist][posc+nvpdof];
            estif[irow+3][icol+5] = tmp[posr+rdist][posc+nvpdof+1];
            estif[irow+3][icol+6] = tmp[posr+rdist][posc+nvpdof+2];
            irow += 4;
            rdist-=2;
         break;
         case 7:
            estif[irow][icol]     = tmp[posr][posc];
            estif[irow][icol+1]   = tmp[posr][posc+1];
            estif[irow][icol+2]   = tmp[posr][posc+2];
            estif[irow][icol+3]   = tmp[posr][posc+cdist];
            estif[irow][icol+4]   = tmp[posr][posc+nvpdof];
            estif[irow][icol+5]   = tmp[posr][posc+nvpdof+1];
            estif[irow][icol+6]   = tmp[posr][posc+nvpdof+2];

            estif[irow+1][icol]   = tmp[posr+1][posc];
            estif[irow+1][icol+1] = tmp[posr+1][posc+1];
            estif[irow+1][icol+2] = tmp[posr+1][posc+2];
            estif[irow+1][icol+3] = tmp[posr+1][posc+cdist];
            estif[irow+1][icol+4] = tmp[posr+1][posc+nvpdof];
            estif[irow+1][icol+5] = tmp[posr+1][posc+nvpdof+1];
            estif[irow+1][icol+6] = tmp[posr+1][posc+nvpdof+2];

            estif[irow+2][icol]   = tmp[posr+2][posc];
            estif[irow+2][icol+1] = tmp[posr+2][posc+1];
            estif[irow+2][icol+2] = tmp[posr+2][posc+2];
            estif[irow+2][icol+3] = tmp[posr+2][posc+cdist];
            estif[irow+2][icol+4] = tmp[posr+2][posc+nvpdof];
            estif[irow+2][icol+5] = tmp[posr+2][posc+nvpdof+1];
            estif[irow+2][icol+6] = tmp[posr+2][posc+nvpdof+2];

            estif[irow+3][icol]   = tmp[posr+rdist][posc];
            estif[irow+3][icol+1] = tmp[posr+rdist][posc+1];
            estif[irow+3][icol+2] = tmp[posr+rdist][posc+2];
            estif[irow+3][icol+3] = tmp[posr+rdist][posc+cdist];
            estif[irow+3][icol+4] = tmp[posr+rdist][posc+nvpdof];
            estif[irow+3][icol+5] = tmp[posr+rdist][posc+nvpdof+1];
            estif[irow+3][icol+6] = tmp[posr+rdist][posc+nvpdof+2];

            estif[irow+4][icol]   = tmp[posr+nvpdof][posc];
            estif[irow+4][icol+1] = tmp[posr+nvpdof][posc+1];
            estif[irow+4][icol+2] = tmp[posr+nvpdof][posc+2];
            estif[irow+4][icol+3] = tmp[posr+nvpdof][posc+cdist];
            estif[irow+4][icol+4] = tmp[posr+nvpdof][posc+nvpdof];
            estif[irow+4][icol+5] = tmp[posr+nvpdof][posc+nvpdof+1];
            estif[irow+4][icol+6] = tmp[posr+nvpdof][posc+nvpdof+2];

            estif[irow+5][icol]   = tmp[posr+nvpdof+1][posc];
            estif[irow+5][icol+1] = tmp[posr+nvpdof+1][posc+1];
            estif[irow+5][icol+2] = tmp[posr+nvpdof+1][posc+2];
            estif[irow+5][icol+3] = tmp[posr+nvpdof+1][posc+cdist];
            estif[irow+5][icol+4] = tmp[posr+nvpdof+1][posc+nvpdof];
            estif[irow+5][icol+5] = tmp[posr+nvpdof+1][posc+nvpdof+1];
            estif[irow+5][icol+6] = tmp[posr+nvpdof+1][posc+nvpdof+2];

            estif[irow+6][icol]   = tmp[posr+nvpdof+2][posc];
            estif[irow+6][icol+1] = tmp[posr+nvpdof+2][posc+1];
            estif[irow+6][icol+2] = tmp[posr+nvpdof+2][posc+2];
            estif[irow+6][icol+3] = tmp[posr+nvpdof+2][posc+cdist];
            estif[irow+6][icol+4] = tmp[posr+nvpdof+2][posc+nvpdof];
            estif[irow+6][icol+5] = tmp[posr+nvpdof+2][posc+nvpdof+1];
            estif[irow+6][icol+6] = tmp[posr+nvpdof+2][posc+nvpdof+2];
            irow += 7;
            rdist-=2;
            break;
            default:
            dserror("numdf invalid!\n");
         }
      }
      icol+=7;
   break;
   default:
      dserror("numdf invalid!\n");
   }
   cdist-=2;
}

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_permestif_ifs */

/*!---------------------------------------------------------------------
  \brief get edgnodes for element line

  <pre>                                                         mn 03/04


  </pre>

  \param   iegnod    INT   	 (o)   edge nodes
  \param  *ele       ELEMENT	 (i)   actual element
  \param   line      DOUBLE	 (i)   actual line number
  \return void

  ----------------------------------------------------------------------*/
void f3_iedg(
    INT     *iegnod,
    ELEMENT *ele,
    INT      surf
    )
{

  INT i;
  const INT iegh[2][6][8] =
   {{{ 0, 1, 2, 3, 0, 0, 0, 0 },  /* surf 0                   */
     { 0, 1, 5, 4, 0, 0, 0, 0 },  /* surf 1                   */
     { 1, 2, 6, 5, 0, 0, 0, 0 },  /* surf 2 for hex 8 element */
     { 2, 3, 7, 6, 0, 0, 0, 0 },  /* surf 3                   */
     { 3, 0, 4, 7, 0, 0, 0, 0 },  /* surf 4                   */
     { 4, 5, 6, 7, 0, 0, 0, 0 }}, /* surf 5                   */
    {{ 0, 0, 0, 0, 0, 0, 0, 0 },   /* surf 0                    */
     { 0, 0, 0, 0, 0, 0, 0, 0 },   /* surf 1                    */
     { 0, 0, 0, 0, 0, 0, 0, 0 },   /* surf 2 for hex 20 element */
     { 0, 0, 0, 0, 0, 0, 0, 0 },   /* surf 3                    */
     { 0, 0, 0, 0, 0, 0, 0, 0 },   /* surf 4                    */
     { 0, 0, 0, 0, 0, 0, 0, 0 }}}; /* surf 5                    */


#ifdef DEBUG
  dstrc_enter("f3_iedg");
#endif

  switch(ele->distyp)
  {
    case hex8:
      for(i=0;i<4;i++) iegnod[i] = iegh[0][surf][i];
      break;
    case hex20:
      dserror("iedg for hex20 not yet implemented !!\n");
      for(i=0;i<8;i++) iegnod[i] = iegh[1][surf][i];
      break;
    case tet4:
      dserror("iedg for tet4 not yet implemented !!\n");
      break;
    case tet10:
      dserror("iedg for tet10 not yet implemented !!\n");
      break;
    default:
      dserror("distyp unknown\n");
  } /* end switch(ele->distyp) */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_iedg */



#endif
