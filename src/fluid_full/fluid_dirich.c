#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#define ONE (1.0)
#define NUMTIMECURVE (5)
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 |  routine to initialise the dirichlet boundary conditions             |
 |                                                         genk 04/02   |
 *----------------------------------------------------------------------*/
void fluid_initdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn)
{
GNODE                *actgnode;
NODE                 *actnode;
ELEMENT              *actele;
int                   i,j;
int                   numnp_total;
int                   numele_total;
int                   predof;
int                   actmat;
double                dens;
int                   numdf;
int                   actcurve;
double                timefac[NUMTIMECURVE];
double                T=0.0;
double                acttimefac;
double                initval;

#ifdef DEBUG 
dstrc_enter("fluid_initdirich");
#endif  

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
predof       = numdf-1;


/*------------ transform real pressure from input to kinematic pressure */
for (i=0;i<numele_total;i++)
{
   actele = &(actfield->dis[0].element[i]);
   actmat = actele->mat-1;
   dens   = mat[actmat].m.fluid->density;
   for(j=0;j<actele->numnp;j++)
   {
      actgnode = actele->node[j]->gnode;
      if (actgnode->dirich==NULL)
         continue;
      if (actgnode->dirich->dirich_onoff.a.iv[predof]!=0)
         actgnode->dirich->dirich_val.a.dv[predof] /= dens;
   }
}

/*---------- set dirichlet conditions at time (0) for zero intial field */
if (fdyn->init==0)
{
/*------------------------------------------ get values from time curve */
   for (actcurve=0;actcurve<numcurve;actcurve++)
   {
     dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
   }   
/*------------------------------------------------- loop over all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode  = &(actfield->dis[0].node[i]); 
      actgnode = actnode->gnode;      
      if (actgnode->dirich==NULL)
         continue;
      for (j=0;j<numdf;j++)
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
         actcurve = actgnode->dirich->curve.a.iv[j]-1;
         if (actcurve<0)
            acttimefac = ONE;
         else
            acttimefac = timefac[actcurve];
         initval  = actgnode->dirich->dirich_val.a.dv[j];               
         actnode->sol_increment.a.da[1][j] = initval*acttimefac;
      }
   }   
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_initdirich*/

 /*----------------------------------------------------------------------*
 |  routine to set dirichlet boundary conditions on at time <time>       |
 |                                                           genk 04/02  |
 *----------------------------------------------------------------------*/
void fluid_setdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn)
{
GNODE                *actgnode;
NODE                 *actnode;
int                   i,j;
int                   numnp_total;
int                   numele_total;
int                   predof;
int                   numdf;
int                   actcurve;
double                timefac[NUMTIMECURVE];
double                T;
double                acttimefac;
double                initval;


#ifdef DEBUG 
dstrc_enter("fluid_setdirich");
#endif 

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
predof       = numdf-1; 
T            = fdyn->time;
numdf        = fdyn->numdf;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
}

/*------------------------------------------------- loop over all nodes */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]); 
   actgnode = actnode->gnode;      
   if (actgnode->dirich==NULL)
         continue;
   for (j=0;j<numdf;j++)
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
         continue;
      actcurve = actgnode->dirich->curve.a.iv[j]-1;
      if (actcurve<0)
         acttimefac = ONE;
      else
         acttimefac = timefac[actcurve];
      initval  = actgnode->dirich->dirich_val.a.dv[j];               
      actnode->sol_increment.a.da[3][j] = initval*acttimefac;
   }
}





/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_settdirich*/
