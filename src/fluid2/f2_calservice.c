/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element 

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
static INT PREDOF = 2;
static INT NUMDF = 3;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!--------------------------------------------------------------------- 
\brief set all arrays for element calculation

<pre>                                                         genk 04/02

get the element velocities and the pressure at different times 

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the	 
      transformation in every time step 			 
				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)  actual element
\param  **xyze     DOUBLE          (o)  nodal coordinates
\param  **eveln    DOUBLE	   (o)  ele vels at time n
\param  **evelng   DOUBLE	   (o)  ele vels at time n+g
\param   *epren    DOUBLE	   (o)  ele pres at time n
\param   *edeadn   DOUBLE          (o)  ele dead load at n (selfweight)
\param   *edeadng  DOUBLE          (o)  ele dead load at n+g (selfweight)
\param   *hasext   INT             (o)  flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                DOUBLE         **xyze,
                DOUBLE         **eveln,    
	        DOUBLE         **evelng,   
	        DOUBLE          *epren,
		DOUBLE          *edeadn,
		DOUBLE          *edeadng,
		INT             *hasext
	      )
{
INT    i;           /* simply some counters                             */
INT    actmat  ;    /* material number of the element                   */
INT    actcurve;    /* actual time curve                                */
DOUBLE acttimefac;  /* time factor from actual curve                    */
DOUBLE dens;        /* density                                          */
NODE  *actnode;     /* actual node                                      */
GSURF *actgsurf;

#ifdef DEBUG 
dstrc_enter("f2_calset");
#endif

/*-------------------------------------------- set element coordinates */
for(i=0;i<ele->numnp;i++)
{
   xyze[0][i]=ele->node[i]->x[0];
   xyze[1][i]=ele->node[i]->x[1];
}

/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
 *---------------------------------------------------------------------*/


if(dynvar->isemim==0)  /* -> implicit time integration method ---------*/
{
   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
      actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */      
      evelng[0][i]=actnode->sol_increment.a.da[3][0];
      evelng[1][i]=actnode->sol_increment.a.da[3][1]; 
/*-------------------------------------- set supported pressures (n+1) */   
   } /* end of loop over nodes of element */
   
} /* endif (dynvar->isemim==0) */
else  /* -> semi-implicit time integration method 
       | for semi-impl. methods one needs extra suup. velocities-------*/
{   
   dserror("semi-implicit methods not checked yet!!!\n");
   if(dynvar->itwost==0)   /* -> semi-implicit one-step ---------------*/
   {
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*-------------------------------- set element velocities at time (n)  */ 	 
         evelng[0][i]=actnode->sol_increment.a.da[1][0];
	 evelng[1][i]=actnode->sol_increment.a.da[1][1];
      } /* end of loop over nodes of element */
   } /* endif (dynvar->itwost==0) */
   else  /* -> semi-implicit two-step ---------------------------------*/
   {  
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*--------------------------- set element velocities at time (n+gamma) */	 
         evelng[0][i]=actnode->sol_increment.a.da[2][0];
	 evelng[1][i]=actnode->sol_increment.a.da[2][1];
      } /* end of loop over nodes of element */
   } /* endif else */  
} /* endif else */

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*------------------------------------- set element velocities at (n) */	 
         eveln[0][i]=actnode->sol_increment.a.da[1][0];
	 eveln[1][i]=actnode->sol_increment.a.da[1][1];    
/*------------------------------------------------- set pressures (n) */   
         epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
      } /* end of loop over nodes of element */  
} /* endif (dynvar->nif!=0) */		       

/*------------------------------------------------ check for dead load */
actgsurf = ele->g.gsurf;
if (actgsurf->neum!=NULL)
{
   actmat=ele->mat-1;
   dens = mat[actmat].m.fluid->density;
   actcurve = actgsurf->neum->curve-1;
   if (actcurve<0) acttimefac=ONE;
   else  dyn_facfromcurve(actcurve,dynvar->acttime,&acttimefac) ;    
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
         edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
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
\param   *dynvar    FLUID_DYN_CALC  (i)
\param   *ele       ELEMENT	    (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\param  **eveln     DOUBLE	    (o)    ele vels at time n
\param  **evelng    DOUBLE	    (o)    ele vels at time n+g
\param  **ealecovn  DOUBLE          (o)    ALE-convective vels at time n
\param  **ealecovng DOUBLE          (o)    ALE-convective vels at time n+g
\param  **egridv    DOUBLE          (o)    element grid velocity
\param   *epren     DOUBLE	    (o)    ele pres at time n
\param   *edeadn    DOUBLE          (o)    ele dead load at n (selfweight)
\param   *edeadng   DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *ekappan   DOUBLE          (o)    nodal curvature at n
\param   *ekappang  DOUBLE          (o)    nodal curvature at n+g
\param   *hasext    INT             (o)    flag for external loads
\param    is_relax  INT             (i)    flag, if it's for relax.-param
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calseta( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                DOUBLE         **xyze,
                DOUBLE         **eveln,    
	        DOUBLE         **evelng,
		DOUBLE         **ealecovn,
		DOUBLE         **ealecovng,
		DOUBLE         **egridv, 
	        DOUBLE          *epren,
		DOUBLE          *edeadn,
		DOUBLE          *edeadng,
		DOUBLE          *ekappan,
		DOUBLE          *ekappang,
		INT             *hasext,
		INT              is_relax
	      )
{
INT    i;           /* simply some counters                             */
INT    actmat  ;    /* material number of the element                   */
INT    actcurve;    /* actual time curve                                */
DOUBLE acttimefac;  /* time factor from actual curve                    */
DOUBLE dens;        /* density                                          */
NODE  *actfnode;    /* actual fluid node                                */
GSURF *actgsurf;    /* actual gsurf                                     */

#ifdef DEBUG 
dstrc_enter("f2_calseta");
#endif

/*-------------------------------------------- set element coordinates */ 
if (is_relax)
{
   f2_alecoor_sd(dynvar,ele,xyze);
}
else
   f2_alecoor(dynvar,ele,xyze);

/*----------------------------------------------------------------------*
 | position of the different solutions:                                 |
 | node->sol_incement: solution history used for calculations	    	|
 |       sol_increment.a.da[0][i]: solution at (n-1)			|
 |       sol_increment.a.da[1][i]: solution at (n)			|
 |       sol_increment.a.da[2][i]: solution at (n+g)			|
 |       sol_increment.a.da[3][i]: solution at (n+1)			|
 |       sol_increment.a.da[4][i]: grid velocity			|
 |       sol_increment.a.da[5][i]: convective velocity at (n)   	|
 |       sol_increment.a.da[6][i]: convective velocity at (n+1) 	|
 *----------------------------------------------------------------------*/

if(dynvar->isemim==0)  /* -> implicit time integration method ----------*/
{
   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
      actfnode=ele->node[i];
/*------------------------------------ set element velocities (n+gamma) */ 
      evelng[0][i]   =actfnode->sol_increment.a.da[3][0];
      evelng[1][i]   =actfnode->sol_increment.a.da[3][1];
      ealecovng[0][i]=actfnode->sol_increment.a.da[6][0];
      ealecovng[1][i]=actfnode->sol_increment.a.da[6][1];
      egridv[0][i]   =actfnode->sol_increment.a.da[4][0];
      egridv[1][i]   =actfnode->sol_increment.a.da[4][1];     
   } /* end of loop over nodes of element */
   
} /* endif (dynvar->isemim==0) */
else  /* -> semi-implicit time integration method 
       | for semi-impl. methods one needs extra suup. velocities-------*/
{   
   dserror("semi-implicit method not implemented for ALE!!!\n");  
} /* endif else */

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actfnode=ele->node[i];
/*------------------------------------- set element velocities at (n) */	 
         eveln[0][i] =actfnode->sol_increment.a.da[1][0];
	 eveln[1][i] =actfnode->sol_increment.a.da[1][1];
	 ealecovn[0][i]=actfnode->sol_increment.a.da[5][0];
	 ealecovn[1][i]=actfnode->sol_increment.a.da[5][1];  
/*------------------------------------------------- set pressures (n) */   
         epren[i]    =actfnode->sol_increment.a.da[1][PREDOF];      
      } /* end of loop over nodes of element */  
} /* endif (dynvar->nif!=0) */		       

/*----------------------------------------------- check for dead load */
actgsurf = ele->g.gsurf;
if (actgsurf->neum!=NULL)
{
   actmat=ele->mat-1;
   dens = mat[actmat].m.fluid->density;
   actcurve = actgsurf->neum->curve-1;
   if (actcurve<0) acttimefac=ONE;
   else  dyn_facfromcurve(actcurve,dynvar->acttime,&acttimefac) ;    
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
         edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
	 edeadng[i] = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
	 (*hasext)++;
      }
   }
}

/*------------------------------------------ curvature at free surface */
if (ele->e.f2->fs_on>0 && dynvar->surftens!=0)
{
   if (dynvar->fsstnif!=0)
   {
      for (i=0;i<ele->numnp;i++)
      {
         ekappan[i]=ele->e.f2->kappa_ND.a.da[i][0];
      }
   }
   if (dynvar->fsstnii!=0)
   {
      for (i=0;i<ele->numnp;i++)
      {
        ekappang[i]=ele->e.f2->kappa_ND.a.da[i][1];
      }
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calseta */


/*!--------------------------------------------------------------------- 
\brief set element coordinates during ALE calculations

<pre>                                                         genk 03/02
			 
				      
</pre>
\param   *dynvar    FLUID_DYN_CALC  (i)
\param   *ele       ELEMENT	    (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\return void                                                                       

------------------------------------------------------------------------*/
void f2_alecoor( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                DOUBLE         **xyze
	       )
{
INT i,j;
DOUBLE omt,theta;
DOUBLE xy,xyn,xyng;
DOUBLE dt;
NODE  *actfnode;    /* actual fluid node                                */
NODE  *actanode;    /* actual ale node                                  */
GNODE *actfgnode;   /* actual fluid gnode                               */


#ifdef DEBUG 
dstrc_enter("f2_alecoor");
#endif

#ifdef D_FSI

omt = dynvar->omt;
theta = ONE-omt;
dt = dynvar->dta;

/*-------------------------------------------- set element coordinates */ 
for(i=0;i<ele->numnp;i++)
{
   actfnode = ele->node[i];
   actfgnode = actfnode->gnode;
   actanode = actfgnode->mfcpnode[genprob.numaf];
   if(actfnode->numdf==NUMDF)
   {
      for (j=0;j<2;j++)
      {
         xy     = actfnode->x[j];
         xyng   = xy + actanode->sol_mf.a.da[1][j];
         xyn    = xy + actanode->sol_mf.a.da[0][j];
         xyze[j][i] =  theta*(xyng)+omt*(xyn); 
      }
   }
   else /* node on implicit free surface */
   {
      for (j=0;j<2;j++)
      {
         xy    = actfnode->x[j];
	 xyn   = xy + actanode->sol_mf.a.da[0][j];
	 xyng  = xyn + actfnode->sol_increment.a.da[3][j+NUMDF]*dt;
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
} /* end of f2_alecoor */


/*!--------------------------------------------------------------------- 
\brief set element coordinates during ALE calculations for relaxation 
parameter of steepest descent method

<pre>                                                         ChFoe 08/03
			 
				      
</pre>
\param   *dynvar    FLUID_DYN_CALC  (i)
\param   *ele       ELEMENT	    (i)    actual element
\param  **xyze      DOUBLE          (o)    nodal coordinates at time n+theta
\return void                                                                       

------------------------------------------------------------------------*/
void f2_alecoor_sd( 
                   FLUID_DYN_CALC  *dynvar, 
   	           ELEMENT         *ele,     
                   DOUBLE         **xyze
	          )
{
INT i,j;
DOUBLE omt,theta;
DOUBLE xy;
NODE  *actfnode;    /* actual fluid node                                */
NODE  *actanode;    /* actual ale node                                  */
GNODE *actfgnode;   /* actual fluid gnode                               */


#ifdef DEBUG 
dstrc_enter("f2_alecoor_sd");
#endif

#ifdef D_FSI

omt = dynvar->omt;
theta = ONE-omt;

/*-------------------------------------------- set element coordinates */ 

/* computation of relaxation parameter always goes from configuration 0
   to the one indicated by the actual incremental residuum vector g_i.
   See diss. Mok for details p. 119 ff */

for(i=0;i<ele->numnp;i++)
{
   actfnode = ele->node[i];
   actfgnode = actfnode->gnode;
   actanode = actfgnode->mfcpnode[genprob.numaf];
   if(actfnode->numdf==NUMDF)
   {
      for (j=0;j<2;j++)
      {
         xyze[j][i] =  theta * ( actanode->sol_mf.a.da[2][j] ); 
      }
   }
   else /* node on implicit free surface */
   {
   dserror("implicit free surface with steepest descent method not yet implementd");
/*      for (j=0;j<2;j++)
      {
         xy    = actfnode->x[j];
	 xyn   = xy + actanode->sol_mf.a.da[0][j];
	 xyng  = xyn + actfnode->sol_increment.a.da[3][j+NUMDF]*dt;
         xyze[j][i] =  theta*(xyng)+omt*(xyn);
      } */
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
} /* end of f2_alecoor_sd */


/*!--------------------------------------------------------------------- 
\brief routine to calculate velocities at integration point

<pre>                                                         genk 04/02
				      
</pre>
\param   *velint   DOUBLE        (o)   velocities at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_veli(
             DOUBLE  *velint,     
             DOUBLE  *funct,    
	     DOUBLE **evel,     
	     INT      iel       
	    ) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f2_veli");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   velint[i]=ZERO; 
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      velint[i] += funct[j]*evel[i][j];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_veli */

/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure at integration point

<pre>                                                         genk 04/02
				      
</pre>
\param  *preint    DOUBLE        (o)   pressure at integration point
\param  *funct     DOUBLE        (i)   shape functions
\param  *epre      DOUBLE        (i)   pressure at element nodes
\param   iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_prei(
             DOUBLE  *preint,     
             DOUBLE  *funct,    
	     DOUBLE  *epre,     
	     INT      iel       
	    ) 
{
INT     j;

#ifdef DEBUG 
dstrc_enter("f2_prei");
#endif

*preint = ZERO;
for (j=0;j<iel;j++)
{
   *preint += funct[j] * epre[j];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_prei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate surface tension at integration point

<pre>                                                         genk 02/03
				      
</pre>
\param  *funct     DOUBLE        (i)   shape functions
\param  *ekappa    DOUBLE        (i)   nodal kappa
\param   iedgnod   INT           (i)   local node numbers of actual line
\param   ngnode	   INT           (i)   number of nodes on actual line
\return DOUBLE                                                                       

------------------------------------------------------------------------*/
DOUBLE f2_kappai(  
		 DOUBLE  *funct,    
	         DOUBLE  *ekappa,  
	         INT     *iedgnod,   
	         INT      ngnode       
	        ) 
{
INT     j;
DOUBLE  kappaint;

#ifdef DEBUG 
dstrc_enter("f2_kappai");
#endif

kappaint = ZERO;
for (j=0;j<ngnode;j++)
{
   kappaint += funct[j] * ekappa[iedgnod[j]];
}
  
/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return (kappaint); 
} /* end of f2_prei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
				      
</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
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
\param    iel	   INT           (i)   number of nodes in this element
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
\param    iel	   INT           (i)   number of nodes in this element
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
\param   *eforce   DOUBLE        (i/o) element force vector
\param  **tmp      DOUBLE        (i)   working array
\param    iel	   DOUBLE        (i)   number of nodes in this ele
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
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array		
\param	  iel	  INT		 (i)   number of nodes in ele
\param	 *dynvar  FLUID_DYN_CALC
\return void                                                                       

------------------------------------------------------------------------*/
void f2_permestif(                  
		   DOUBLE         **estif,   
		   DOUBLE         **emass, 
		   DOUBLE         **tmp,   
		   ELEMENT         *ele,   
		   FLUID_DYN_CALC  *dynvar		   		    
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
iel    = ele->numnp;
nvdof  = NUM_F2_VELDOF*iel;
npdof  = iel;
totdof = NUMDOF_FLUID2*iel;
thsl   = dynvar->thsl;

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
if (dynvar->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (dynvar->nis==0) */

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

<pre>                                                         genk 01/03

routine to add galerkin and stabilisation parts of the elment	   
stiffness matrix and to rearrange its entries!  		   
this is necessary since we would like to use the existing assembly 
routines for the stiffness matrix				   
hence a splitting of vel- pre-  and grid dofs is not possible any more!!!!

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array		
\param	  iel	  INT		 (i)   number of nodes in ele
\param	 *dynvar  FLUID_DYN_CALC
\return void                                                                       

------------------------------------------------------------------------*/
void f2_permestif_ifs(                  
		      DOUBLE         **estif,   
		      DOUBLE         **emass, 
		      DOUBLE         **tmp,   
		      ELEMENT         *ele,   
		      FLUID_DYN_CALC  *dynvar		   		    
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
iel    = ele->numnp;
nvdof  = NUM_F2_VELDOF*iel;
nvpdof = NUMDOF_FLUID2*iel;
totdof = (NUMDOF_FLUID2+NUM_F2_VELDOF)*iel;
thsl   = dynvar->thsl;

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
if (dynvar->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (dynvar->nis==0) */

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

\param   iegnod    INT   	 (o)   edge nodes
\param  *ele       ELEMENT	 (i)   actual element
\param   line      DOUBLE	 (i)   actual line number		
\param	 iinit	   INT		 (i)   flag
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
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/
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
