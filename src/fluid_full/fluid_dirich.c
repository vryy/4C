/*!----------------------------------------------------------------------
\file
\brief setting dirichlet conditions for fluid

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "fluid_prototypes.h"
#include "../solver/solver.h"
/*-------------------------------------------- prototypes for this file */
#ifdef D_FSI
static DOUBLE fluid_usd(DOUBLE r, DOUBLE alpha);
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;

static FLUID_DYNAMIC *fdyn;

/*!---------------------------------------------------------------------
\brief routine to initialise the dirichlet boundary conditions

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are initialised:
- real pressure from input-file is transformed to kinematic pressure:
  actgnode->dirich->dirich_val.a.dv[predof] /= dens

- if the initial field is a zero-field, then set the dirichlet values in
  the solution history of the node:
    'actnode->sol.a.da[0][j] = initval*acttimefac' (--> output)
    'actnode->sol_increment.a.da[ipos.veln][j] = initval*acttimefac'
                                         |        |         |
                                    time (n)      |         |
                                  initial value from input  |
                                               factor from timecurve (T=0.0)

</pre>
\param *actfield FIELD         (i)  actual field (fluid)

\return void

------------------------------------------------------------------------*/
void fluid_initdirich(  FIELD          *actfield )
{
INT        i,j;
INT        numnp_total;               /* total number of fluid nodes    */
INT        numele_total;              /* total number of fluid elements */
INT        predof;	              /* number of pressure dof	        */
INT        numdf;	              /* number of fluid dofs	        */
INT        actcurve;	              /* actual timecurve  	        */
INT        numveldof;
DOUBLE     dens;	              /* density			*/
DOUBLE     timefac[MAXTIMECURVE];     /* factors from time-curve        */
DOUBLE     T=0.0;	              /* starting time		        */
DOUBLE     acttimefac;                /* actual factor from timecurve   */
DOUBLE     initval;	              /* intial dirichlet value	        */
DOUBLE     funct_fac;	             /* factor from spatial function */
GNODE     *actgnode;	              /* actual GNODE		        */
NODE      *actnode;	              /* actual NODE		        */
ELEMENT   *actele;	              /* actual ELEMENT		        */
#ifdef D_FSI
NODE      *firstnode, *lastnode;
NODE      *actanode,*alnode,*afnode;
DOUBLE     ract,rtotal,r;
DOUBLE     rtotalv[3],ractv[3];
DOUBLE     alpha,length;
DOUBLE     uact,u10;
#endif

INT counter=0;

#ifdef DEBUG
dstrc_enter("fluid_initdirich");
#endif

fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
predof       = numdf-1;
numveldof    = numdf-1;

/*-------------------------- since different materials are not allowed
              one can work with the material parameters of any element */
actele = &(actfield->dis[0].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

/*------------------------------------------ check dirichlet conditions */
for (i=0;i<numnp_total;i++) /* loop all nodes */
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
      continue;
   if (actgnode->dirich->dirich_type==dirich_FSI)
      counter++;
   if (actgnode->dirich->dirich_type==dirich_none)
   {
      for (j=0;j<actnode->numdf;j++) /* loop all dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
            actcurve = actgnode->dirich->curve.a.iv[j];
	    if(actcurve>numcurve)
	       dserror("Load curve: actual curve > number defined curves\n");
      } /* end of loop over all dofs */
#ifdef D_FSI
      if (fdyn->freesurf==5 && actnode->xfs!=NULL)
      {
         if (actgnode->dirich->dirich_onoff.a.iv[numdf]!=0)
	    dserror("No DBC for height function dof allowed!\n");
      }
#endif
      /* transform real pressure from input to kinematic pressure ---*/
      if (actgnode->dirich->dirich_onoff.a.iv[predof]!=0)
          actgnode->dirich->dirich_val.a.dv[predof] /= dens;
   }
} /* end of loop over all nodes */

if (counter>0 && par.myrank==0)
{
printf("\n");
printf("          | FIELD FLUID     | number of nodes coupled with structure: %d \n",counter);
printf("\n");
}

/*---------- set dirichlet conditions at time (0) for zero intial field */
if (fdyn->init==0)
{
/*------------------------------------------ get values from time curve */
   for (actcurve=0;actcurve<numcurve;actcurve++)
   {
     dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
   }/* end loop over active timecurves */
   /*--------------------- before applying DBC transform to xyz* co-sys */
   locsys_trans_sol_dirich(actfield,0,0,0,0);
   locsys_trans_sol_dirich(actfield,0,1,1,0);
/*------------------------------------------------- loop over all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode  = &(actfield->dis[0].node[i]);
      actgnode = actnode->gnode;
      if (actgnode->dirich==NULL)
         continue;
      switch(actgnode->dirich->dirich_type)
      {
      case dirich_none:
         for (j=0;j<actnode->numdf;j++) /* loop all dofs */
         {
            if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
               continue;
            actcurve = actgnode->dirich->curve.a.iv[j]-1;
            if (actcurve<0)
               acttimefac = ONE;
            else
               acttimefac = timefac[actcurve];
            initval   = actgnode->dirich->dirich_val.a.dv[j];
            funct_fac = actgnode->d_funct[j];
            actnode->sol_increment.a.da[ipos.veln][j] = initval*acttimefac*funct_fac;
            actnode->sol.a.da[0][j] = initval*acttimefac*funct_fac;
         } /* end loop over dofs */
      break;
#ifdef D_FSI
      case dirich_FSI: /* FSI --> dirichvalues = grid velocity!!! */
	 for (j=0;j<numveldof;j++) /* loop vel-dofs */
	 {
	    initval = actnode->sol_increment.a.da[ipos.gridv][j];
	    actnode->sol_increment.a.da[ipos.veln][j] = initval;
	    actnode->sol.a.da[0][j] = initval;
	 }
      break;
      case dirich_FSI_pseudo: /* FSI --> dirichvalues = grid velocity!!! */
	 for (j=0;j<numveldof;j++) /* loop vel-dofs */
	 {
	    initval = actnode->sol_increment.a.da[ipos.gridv][j];
	    actnode->sol_increment.a.da[ipos.veln][j] = initval;
	    actnode->sol.a.da[0][j] = initval;
	 }
      break;
      case dirich_slip:
         /*------------------------------------ get values and pointers */
         dsassert(numveldof==2,"slip-dirich only for 2D up to now!\n");
         firstnode=actgnode->slipdirich->firstnode;
         lastnode =actgnode->slipdirich->lastnode;
         afnode = firstnode->gnode->mfcpnode[genprob.numaf];
         alnode = lastnode->gnode->mfcpnode[genprob.numaf];
         actanode=actgnode->mfcpnode[genprob.numaf];
         dsassert(actanode!=NULL,"cannot read from NULL-pointer!\n");
         dsassert(afnode!=NULL,"cannot read from NULL-pointer!\n");
         dsassert(alnode!=NULL,"cannot read from NULL-pointer!\n");
         length=actgnode->slipdirich->length;
         alpha =actgnode->slipdirich->alpha;
         actcurve = actgnode->slipdirich->curve-1;
         /*------------------------------------- evaluate time function */
         if (actcurve<0)
            acttimefac = ONE;
         else
            acttimefac = timefac[actcurve];
         initval = actgnode->slipdirich->dirich_val;

         /*------------------------ actual total length of slip BC line */
         for(j=0;j<3;j++)
            rtotalv[j]= (lastnode->x[j]+alnode->sol_mf.a.da[1][j])
                       -(firstnode->x[j]+afnode->sol_mf.a.da[1][j]);
         rtotal=sqrt(DSQR(rtotalv[0])+DSQR(rtotalv[1])+DSQR(rtotalv[2]));

         /*---------------------------- actual position on slip BC line */
         for(j=0;j<3;j++)
            ractv[j]= (actnode->x[j]+actanode->sol_mf.a.da[1][j])
                     -(firstnode->x[j]+afnode->sol_mf.a.da[1][j]);
         ract=sqrt(DSQR(ractv[0])+DSQR(ractv[1])+DSQR(ractv[2]));
         if (rtotal<length)
         {
            length=rtotal;
            printf("WARNING: SLIP BC LINE shorter than given length!\n");
         }
         if (ract<length) /* evaluate slip model function */
         {
            r=ract/length*TEN;
            uact = fluid_usd(r,alpha);
            u10  = fluid_usd(TEN,alpha);
            actnode->sol_increment.a.da[ipos.veln][1]=uact/u10*initval*acttimefac;
            actnode->sol_increment.a.da[ipos.veln][0]=ZERO;
            actnode->sol.a.da[0][1]=uact/u10*initval*acttimefac;
            actnode->sol.a.da[0][0]=ZERO;
         }
         else
         {
            actnode->sol_increment.a.da[ipos.veln][1]=initval*acttimefac;
            actnode->sol_increment.a.da[ipos.veln][0]=ZERO;
            actnode->sol.a.da[0][1]=uact/u10*initval*acttimefac;
            actnode->sol.a.da[0][0]=ZERO;
         }
      break;
#endif
      default:
         dserror("dirch_type unknown!\n");
      } /* end switch */
   } /*end loop over nodes */

   /*------------------------------------- transform back to XYZ co-sys */
   locsys_trans_sol_dirich(actfield,0,0,0,1);
   locsys_trans_sol_dirich(actfield,0,1,1,1);

} /* endif fdyn->init */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_initdirich*/


/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions at time <time>

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->acttime>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[pos][j] = initval*acttimefac'
                                  |         |         |
                            time (n+1)      |         |
                          initial value from input    |
                                         factor from timecurve
</pre>
\param *actfield FIELD		(i)  actual field (fluid)
\param  pos	 INT		(i)	position, where to write dbc

\return void

------------------------------------------------------------------------*/
void fluid_setdirich(   FIELD           *actfield,
			INT		 pos
	            )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        actcurve;	             /* actual timecurve		*/
INT        numveldof;
DOUBLE     timefac[MAXTIMECURVE];    /* factors from time-curve         */
DOUBLE     T;		             /* actual time		        */
DOUBLE     acttimefac;               /* actual factor from timecurve    */
DOUBLE     initval;	             /* intial dirichlet value	        */
DOUBLE     funct_fac;	             /* factor from spatial function */
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */
#ifdef D_FSI
NODE      *firstnode, *lastnode;
NODE      *actanode,*alnode,*afnode;
DOUBLE     ract,rtotal,r;
DOUBLE     rtotalv[3],ractv[3];
DOUBLE     alpha,length;
DOUBLE     uact,u10;
#endif

#ifdef DEBUG
dstrc_enter("fluid_setdirich");
#endif

/*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = fdyn->acttime;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
} /* end loop over active timecurves */

/*------------------------ dirichlet values are applied in the xyz* co-sys
   so transform nodal values with dirichlet conditions                   */
locsys_trans_sol_dirich(actfield,0,1,pos,0);

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
         actcurve = actgnode->dirich->curve.a.iv[j]-1;
         if (actcurve<0)
            acttimefac = ONE;
         else
            acttimefac = timefac[actcurve];
         initval   = actgnode->dirich->dirich_val.a.dv[j];
         funct_fac = actgnode->d_funct[j];
         actnode->sol_increment.a.da[pos][j] = initval*acttimefac*funct_fac;
      } /* end loop over dofs */
   break;
#ifdef D_FSI
   case dirich_FSI: /* dirichvalues = grid velocity!!! */
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[pos][j]
	=actnode->sol_increment.a.da[ipos.gridv][j];
   break;
   case dirich_FSI_pseudo: /* dirichvalues = grid velocity!!! */
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[pos][j]
	=actnode->sol_increment.a.da[ipos.gridv][j];
   break;
   case dirich_slip: /* slip boundary condition */
      /*--------------------------------------- get values and pointers */
      dsassert(numveldof==2,"slip-dirich only for 2D up to now!\n");
      firstnode=actgnode->slipdirich->firstnode;
      lastnode =actgnode->slipdirich->lastnode;
      afnode = firstnode->gnode->mfcpnode[genprob.numaf];
      alnode = lastnode->gnode->mfcpnode[genprob.numaf];
      actanode=actgnode->mfcpnode[genprob.numaf];
      dsassert(actanode!=NULL,"cannot read from NULL-pointer!\n");
      dsassert(afnode!=NULL,"cannot read from NULL-pointer!\n");
      dsassert(alnode!=NULL,"cannot read from NULL-pointer!\n");
      length=actgnode->slipdirich->length;
      alpha =actgnode->slipdirich->alpha;
      actcurve = actgnode->slipdirich->curve-1;
      /*---------------------------------------- evaluate time function */
#if 0
      if (actcurve<0)
         acttimefac = ONE;
      else
         acttimefac = timefac[actcurve];

      initval = actgnode->slipdirich->dirich_val;
#endif
      acttimefac=ONE;
      initval = lastnode->sol_increment.a.da[ipos.gridv][1];
      /*--------------------------- actual total length of slip BC line */
      for(j=0;j<3;j++)
         rtotalv[j]= (lastnode->x[j]+alnode->sol_mf.a.da[1][j])
                    -(firstnode->x[j]+afnode->sol_mf.a.da[1][j]);
      rtotal=sqrt(DSQR(rtotalv[0])+DSQR(rtotalv[1])+DSQR(rtotalv[2]));

      /*------------------------------- actual position on slip BC line */
      for(j=0;j<3;j++)
         ractv[j]= (actnode->x[j]+actanode->sol_mf.a.da[1][j])
                  -(firstnode->x[j]+afnode->sol_mf.a.da[1][j]);
      ract=sqrt(DSQR(ractv[0])+DSQR(ractv[1])+DSQR(ractv[2]));
      if (rtotal<length)
      {
         length=rtotal;
         printf("WARNING: SLIP BC LINE shorter than given length!\n");
      }
      if (ract<length) /* evaluate slip model function */
      {
         r=ract/length*TEN;
         uact = fluid_usd(r,alpha);
         u10  = fluid_usd(TEN,alpha);
         actnode->sol_increment.a.da[pos][1]=uact/u10*initval*acttimefac;
         actnode->sol_increment.a.da[pos][0]=ZERO;
      }
      else
      {
         actnode->sol_increment.a.da[pos][1]=initval*acttimefac;
         actnode->sol_increment.a.da[pos][0]=ZERO;
      }
   break;
#endif
   default:
      dserror("dirich_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*------------------------ dirichlet values are applied in the xyz* co-sys
   so transform back nodal values with dirichlet conditions             */
locsys_trans_sol_dirich(actfield,0,1,pos,1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich*/


/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions for a parabolic velocity
profile modified for projection algorithm

<pre>                                                         genk 04/02
                                                              basol 03/03

in this routine the dirichlet boundary conditions for fluid2pro
elements are set at time <T=fdyn->acttime>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac'
                                          |        |         |
                                     time (n+1)    |         |
                                   initial value from input  |
                                                factor from timecurve
</pre>
\param *actfield FIELD         (i)  actual field (fluid)

\return void

------------------------------------------------------------------------*/
void fluid_setdirich_parabolic(FIELD  *actfield)
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        actcurve;	             /* actual timecurve		*/
DOUBLE     timefac[MAXTIMECURVE];    /* factors from time-curve         */
DOUBLE     T;		            /* actual time		       */
DOUBLE     acttimefac;              /* actual factor from timecurve    */
DOUBLE     L=ONE;
DOUBLE     initval;	            /* intial dirichlet value	       */
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;	             /* actual NODE		        */

#ifdef DEBUG
dstrc_enter("fluid_setdirich_parabolic");
#endif
/*======================================================================*
* REMARK: Here the inflow profile is assumed to be set at point x=-0.5
* also the profile is formulated according to the parabolic formula below
* -(TWO/L)*(TWO/L)*actnode->x[1]*actnode->x[1]+ONE;
* a better way of doing this would be implementing the velocity profile into
* preprocessor code "GID"
*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = fdyn->acttime;
numdf        = fdyn->numdf;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
} /* end loop over active timecurves */

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   dsassert(actnode->locsysId==0,"locsys not implemented yet!\n");
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL) continue;
   for (j=0;j<numdf;j++) /* loop dofs */
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0) continue;
      actcurve = actgnode->dirich->curve.a.iv[j]-1;
      if (actcurve<0)
      {
         acttimefac = ONE;
      }else{
         acttimefac = timefac[actcurve];
      }
      if (actnode->x[0]==-0.5)
      {
         if (j==0) /*--the parabolic velocity profile is only for the x-dof--*/
	 {
            /*--parabolic profile is formulated as below---*/
            initval  = -(TWO/L)*(TWO/L)*actnode->x[1]*actnode->x[1]+ONE;
            /*--it is multiplied with a timefac------------*/
	    actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
         }else{
	    initval  = actgnode->dirich->dirich_val.a.dv[j];
            actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
	 }/*end of if (j==0)*/
      }else{
         initval  = actgnode->dirich->dirich_val.a.dv[j];
         actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
      }/*end of (actnode->x[0]==-0.5)*/
   } /* end loop over dofs */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich*/

/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions for the specific problem
"flow around a cylinder"

<pre>                                                         basol 05/03

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->acttime>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac'
                                          |        |         |
                                     time (n+1)    |         |
                                   initial value from input  |
                                                factor from timecurve
</pre>
\param *actfield FIELD         (i)  actual field (fluid)

\return void

------------------------------------------------------------------------*/
void fluid_setdirich_cyl(FIELD  *actfield)
{
INT        i,j;
INT        numnp_total;                /* total number of fluid nodes     */
INT        numele_total;               /* total number of fluid elements  */
INT        numdf;	               /* number of fluid dofs    	  */
INT        actcurve;	               /* actual timecurve		  */
DOUBLE   timefac[MAXTIMECURVE];      /* factors from time-curve         */
DOUBLE   T;  		             /* actual time		        */
DOUBLE   acttimefac;                 /* actual factor from timecurve    */
DOUBLE   initval;  	             /* intial dirichlet value	        */
/*DOUBLE   Um=1.5;                     /* maximum velocity                */
DOUBLE   Um=0.3;                     /* maximum velocity                */
DOUBLE   H=0.41;                     /* height of the channel           */
GNODE     *actgnode;	              /* actual GNODE		         */
NODE      *actnode;	              /* actual NODE		         */

#ifdef DEBUG
dstrc_enter("fluid_setdirich_cyl");
#endif
/*======================================================================*
* REMARK: Here the inflow profile is assumed to be set at point x=0.0
* also the profile is formulated according to the parabolic formula below
* FOUR*Um*actnode->x[1]*(H-actnode->x[1])/(H*H);
* see the DFG Benchmark paper for the flow around a cylinder for the proper
* description of the problem
* a better way of doing this would be implementing the velocity profile into
* preprocessor code "GID"
*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = fdyn->acttime;
numdf        = fdyn->numdf;
/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
} /* end loop over active timecurves */

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   dsassert(actnode->locsysId==0,"locsys not implemented yet!\n");
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   for (j=0;j<numdf;j++) /* loop dofs */
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
         continue;
      actcurve = actgnode->dirich->curve.a.iv[j]-1;
      if (actcurve<0)
         acttimefac = ONE;
      else
         acttimefac = timefac[actcurve];
      if (actnode->x[0]==0.0)
      {
         if (j==0)/*--the parabolic velocity profile is only for the x-dof--*/
	 {
            /*--parabolic profile is formulated as below---*/
	    initval  = FOUR*Um*actnode->x[1]*(H-actnode->x[1])/(H*H);
            /*--it is multiplied with a timefac------------*/
	    actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
         }else{
	    initval  = actgnode->dirich->dirich_val.a.dv[j];
            actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
	 }/*end of if (j==0)*/
      }else{
         initval  = actgnode->dirich->dirich_val.a.dv[j];
         actnode->sol_increment.a.da[ipos.velnp][j] = initval*acttimefac;
      }/*end of if (actnode->x[0]==0.0)*/
   } /* end loop over dofs */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_settdirich*/


/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions for fluid to determine
Relaxation parameter

<pre>                                                             cf 08/03

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set for the fluid solution needed to determine the Relaxation
parameter via steepest descent method.

the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[ipos.relax][j] = 0.0
                                 |          |
                fluid sol. for RelaxParam   |
					at Dirichlet boundaries
  AND:
    'actnode->sol_increment.a.da[ipos.relax][j] 
                                 |  = actnode->sol_increment.a.da[ipos.gridv][j]
                                 |                  |
                fluid sol. for RelaxParam           |
					at fsi coupling interface,
					      grid velocity
</pre>
\param *actfield FIELD         (i)  actual field (fluid)

\return void

------------------------------------------------------------------------*/
void fluid_setdirich_sd( FIELD *actfield )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        numveldof;
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG
dstrc_enter("fluid_setdirich_sd");
#endif

/*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   dsassert(actnode->locsysId==0,"locsys not implemented yet!\n");
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
         actnode->sol_increment.a.da[ipos.relax][j] = 0.0;
         actnode->sol_increment.a.da[ipos.convnp][j] = 0.0;
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[ipos.relax][j]
	=actnode->sol_increment.a.da[ipos.gridv][j];
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich_sd*/




/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions of acceleration

<pre>                                                             cf 09/03

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set for the fluid acceleration. These conditions depend
upon the dbc for the velocity set in the input file according to the
time stepping sceme.
The actual implementation serves the generalised-alpha scheme.

sol_increment[pos_to][i] = fac1 * sol_increment[pos_from1][i]
                         + fac2 * sol_increment[pos_from2][i]
                         + fac3 * sol_increment[pos_from3][i]


</pre>
\param *actfield	FIELD		(i)	actual field (fluid)
\param  pos_to		INT		(i)	pos in sol_increment to write to
\param  pos1_from	INT		(i)	1st pos to read from
\param  pos2_from	INT		(i)	2nd pos to read from
\param  pos3_from	INT		(i)	3rd pos to read from
\param	fac1		DOUBLE		(i)	factor of value at pos1_from
\param	fac2		DOUBLE		(i)	factor of value at pos2_from
\param	fac3		DOUBLE		(i)	factor of value at pos3_from

\return void

------------------------------------------------------------------------*/
void fluid_setdirich_acc(
                         FIELD		*actfield,
			 INT		 pos_to,
			 INT		 pos1_from,
			 INT		 pos2_from,
			 INT		 pos3_from,
			 DOUBLE		 fac1,
			 DOUBLE		 fac2,
			 DOUBLE		 fac3
	                )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        numveldof;
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG
dstrc_enter("fluid_setdirich_acc");
#endif

/*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   dsassert(actnode->locsysId==0,"locsys not implemented yet!\n");
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
	       actnode->sol_increment.a.da[pos_to][j]
	       = fac1 * actnode->sol_increment.a.da[pos1_from][j]
               + fac2 * actnode->sol_increment.a.da[pos2_from][j]
               + fac3 * actnode->sol_increment.a.da[pos3_from][j];
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
	 dserror("generalised alpha with FSI not yet implemented");
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich_acc*/



/*!--------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector

<pre>                                                         genk 04/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[velnp][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= estif[i][j] * dirich[j];
\endcode

\param  *actele    ELEMENT   (i)   actual element
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces
\param   readfrom  INT       (i)   position, where to read dbc from

\return void
------------------------------------------------------------------------*/
void fluid_caldirich(
                        ELEMENT         *actele,
		        DOUBLE          *dforces,
                        DOUBLE         **estif,
		        INT             *hasdirich,
			INT		 readfrom
                     )
{

INT         i,j;
INT         ilocsys;
INT         nrow;
INT         numdf;                      /* number of fluid dofs         */
INT         nd=0;
DOUBLE      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele  */
DOUBLE      val[MAXDOFPERNODE];
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */

#ifdef DEBUG
dstrc_enter("fluid_caldirich");
#endif

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;
   if (actgnode->dirich==NULL)
      continue;
   else
      *hasdirich=1;
      break;
} /* end loop over nodes */

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = ZERO;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                  written to the nodes (sol_increment[ipos.velnp][j]) */
nrow=0;
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   actnode  = actele->node[i];
   numdf    = actnode->numdf;
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
   {
      nrow+=numdf;
      continue;
   }
   ilocsys=actnode->locsysId-1;
   if (ilocsys>=0) /* local co-sys */
   {
      /*--------------------- transform values at node from XYZ to xyz* */
      for (j=0;j<numdf;j++)
         val[j]=actnode->sol_increment.a.da[readfrom][j];
      locsys_trans_nodval(actele,&(val[0]),numdf,ilocsys,0);
      for (j=0;j<numdf;j++)
      {
         dirich_onoff[nrow+j] = actgnode->dirich->dirich_onoff.a.iv[j];
         dirich[nrow+j] = val[j];
      }
   }
   else
   {
      for (j=0; j<numdf; j++) /* loop dofs */
      {
         dirich_onoff[nrow+j] = actgnode->dirich->dirich_onoff.a.iv[j];
         dirich[nrow+j] = actnode->sol_increment.a.da[readfrom][j];
      } /* end loop over dofs */
   }
   nrow+=numdf;
} /* end loop over nodes */

dsassert(nrow==nd,"failure during calculation of dirich forces\n");
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_caldirich*/

/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector for the
projection method (constant velocity profile)

<pre>                                                     basol 11/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[ipos.velnp][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= (emass[i][j]+dt*estif[i][j]) * dirich[j];
\endcode

\param  *actele    ELEMENT   (i)   actual element
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param **emass     DOUBLE    (i)   element mass matrix
\param   dt        DOUBLE    (i)   time increment
\param  theta      DOUBLE    (i)   variable for the time integration of viscousity matrix
                                   for the implicitly treated K, theta=1.0
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces
\return void
------------------------------------------------------------------------*/
void fluid_pm_caldirich(
                     ELEMENT   *actele,
		     DOUBLE   *dforces,
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,
		     DOUBLE   theta,
		     INT       *hasdirich
		    )
{

INT         i,j;
INT         numdf;                      /* number of fluid dofs         */
INT         nd=0;
DOUBLE    dirich[MAXDOFPERELE];      /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE];  /* dirichlet flags of act. ele  */
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */

#ifdef DEBUG
dstrc_enter("fluid_pm_caldirich");
#endif

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;
   if (actgnode->dirich==NULL)
      continue;
   else
      *hasdirich=1;
   break;
} /* end loop over nodes */

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                  written to the nodes (sol_increment[ipos.velnp][j]) */
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actnode->sol_increment.a.da[ipos.velnp][j];
   } /* end loop over dofs */
} /* end loop over nodes */

/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= (estif[i][j]*dt*theta+emass[i][j]) * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_caldirich*/

/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector for the
       projection method for a parabolic velocity profile

<pre>                                                     basol 11/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[ipos.velnp][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= (emass[i][j]+dt*estif[i][j]) * dirich[j];
\endcode

\param  *actele    ELEMENT   (i)   actual element
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param **emass     DOUBLE    (i)   element mass matrix
\param   dt        DOUBLE    (i)   time increment
\param  theta      DOUBLE    (i)   variable for the time integration of viscousity matrix
                                   for the implicitly treated K, theta=1.0
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces
\return void
------------------------------------------------------------------------*/
void fluid_pm_caldirich_parabolic(
                     ELEMENT   *actele,
		     DOUBLE   *dforces,
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,
                     DOUBLE   theta,
		     INT       *hasdirich
		    )
{

INT         i,j;
INT         numdf;                          /* number of fluid dofs         */
INT         nd=0;
DOUBLE    dirich[MAXDOFPERELE];           /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE];    /* dirichlet flags of act. ele  */
GNODE      *actgnode;	                  /* actual GNODE                 */
NODE       *actnode;	                  /* actual NODE                  */
DOUBLE     tol=1.0E-6;
DOUBLE     y_coor;
DOUBLE     L=ONE;
#ifdef DEBUG
dstrc_enter("fluid_pm_caldirich_parabolic");
#endif
/*======================================================================*
* REMARK: the profile is formulated according to the parabolic formula below
* -(TWO/L)*(TWO/L)*y_coor*y_coor+ONE;
* if the parabolic input profile is implemented into the preprocessor "GID"
* one doesn't need any further subroutines then "fluid_caldirich"
*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;
   if (actgnode->dirich==NULL)
      continue;
   else
      *hasdirich=1;
   break;
} /* end loop over nodes */

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                  written to the nodes (sol_increment[ipos.velnp][j]) */
/*----------------------------------------------------------------------*/
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      if (actgnode->dirich->dirich_val.a.dv[j] < tol)
      {
         dirich[i*numdf+j] = ZERO;
      }else{
          /*--the parabolic equation is formulated---*/
	  /*--according to the y-coordinate----------*/
          y_coor = actnode->x[1];
	  /*--parabolic value of the velocity--------*/
	  dirich[i*numdf+j] = -(TWO/L)*(TWO/L)*y_coor*y_coor+ONE;
      }/* end of if (actgnode->dirich->dirich_val.a.dv[j] < tol)*/
   } /* end loop over dofs */
} /* end loop over nodes */
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      /*--(M+dt*K)*dirich-----------------------------*/
      /*--K_eff*dirich--------------------------------*/
      dforces[i] -= (estif[i][j]*dt*theta+emass[i][j]) * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_pm_caldirich_parabolic*/

/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector for the
       projection method for a parabolic velocity profile & for the specific
       problem "flow around a cylinder"

<pre>                                                     basol 11/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[ipos.velnp][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= (emass[i][j]+dt*estif[i][j]) * dirich[j];
\endcode

\param  *actele    ELEMENT   (i)   actual element
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param **emass     DOUBLE    (i)   element mass matrix
\param   dt        DOUBLE    (i)   time increment
\param  theta      DOUBLE    (i)   variable for the time integration of viscousity matrix
                                   for the implicitly treated K, theta=1.0
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces
\return void
*----------------------------------------------------------------------*/
void fluid_pm_caldirich_cyl(
                     ELEMENT   *actele,
		     DOUBLE   *dforces,
                     DOUBLE   **estif,
		     DOUBLE   **emass,
		     DOUBLE   dt,
                     DOUBLE   theta,
		     INT       *hasdirich
		    )
{

INT         i,j;
INT         numdf;                        /* number of fluid dofs         */
INT         nd=0;
DOUBLE    dirich[MAXDOFPERELE];         /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE];  /* dirichlet flags of act. ele  */
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */
DOUBLE     tol=1.0E-6;
/*===============================================================*/
/*for detailed description of the problem see the DFG benchmark--*/
/*paper   ------------------------------------------------------ */
DOUBLE   Um=1.5;                       /* mean velocity                */
DOUBLE   H=0.41;                       /* height of the channel        */
/*===============================================================*/
#ifdef DEBUG
dstrc_enter("fluid_pm_caldirich_cyl");
#endif

/*======================================================================*
* REMARK: the profile is formulated according to the parabolic formula below
* FOUR*Um*actnode->x[1]*(H-actnode->x[1])/(H*H);
* if the parabolic input profile is implemented into the preprocessor "GID"
* one doesn't need any further subroutines then "fluid_caldirich"
*------------------------- check if there are any dirichlet conditions--*/
/*------------------------- for the nodes of this element----------------*/
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;
   if (actgnode->dirich==NULL)
      continue;
   else
      *hasdirich=1;
   break;
} /* end loop over nodes */

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                  written to the nodes (sol_increment[ipos.velnp][j]) */
/*----------------------------------------------------------------------*/
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      if (actgnode->dirich->dirich_val.a.dv[j] < tol)
      {
         dirich[i*numdf+j] = ZERO;
      }else{
         /*--parabolic value of the velocity-----------------*/
         dirich[i*numdf+j] = FOUR*Um*actnode->x[1]*(H-actnode->x[1])/(H*H);
      }/*end of if (actgnode->dirich->dirich_val.a.dv[j] < tol)*/
   } /* end loop over dofs */
} /* end loop over nodes */
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      /*--(M+dt*K)*dirich-----------------------------*/
      /*--K_eff*dirich--------------------------------*/
      dforces[i] -= (estif[i][j]*dt*theta+emass[i][j]) * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_pm_caldirich_parabolic*/
/*----------------------------------------------------------------------*
 | function for slip BC                                                 |
 *----------------------------------------------------------------------*/
#ifdef D_FSI
static DOUBLE fluid_usd(DOUBLE r, DOUBLE alpha)
{
DOUBLE twoa;
DOUBLE rpi;
DOUBLE ur;
DOUBLE fac1,fac2,fac3,fac4;

#ifdef DEBUG
dstrc_enter("fluid_usd");
#endif

twoa=TWO*alpha;
rpi = pow(r,PI/twoa);

if (FABS(rpi)<EPS8)
   dserror("cannot evaluate slip BC function at r=0.0");

fac1=twoa*rpi/(twoa-sin(twoa));
fac2=(-rpi+ONE+log(rpi))*sin(twoa);
fac3=twoa*DSQR(rpi-ONE);
fac4=ONE/(rpi+ONE);

ur=fac1*(fac2/fac3+fac4);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return(ur);
}
#endif
/*! @} (documentation module close)*/
#endif
