#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../headers/solution.h"
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
 |  routine to set dirichlet boundary conditions on at time <time>       |
 |                                                           genk 04/02  |
 *----------------------------------------------------------------------*/
void ale_setdirich(FIELD  *actfield, STRUCT_DYNAMIC *sdyn)
{
GNODE                *actgnode;
NODE                 *actnode;
int                   i,j;
int                   numnp_total;
int                   numele_total;
int                   actcurve;
double                timefac[NUMTIMECURVE];
double                T;
double                acttimefac;
double                initval;

double                cx,cy,win,wino,dd;

#ifdef DEBUG 
dstrc_enter("ale_setdirich");
#endif 

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = sdyn->time;

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
   for (j=0;j<actnode->numdf;j++)
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
         continue;
      actcurve = actgnode->dirich->curve.a.iv[j]-1;
      if (actcurve<0)
         acttimefac = ONE;
      else
         acttimefac = timefac[actcurve];
      initval  = actgnode->dirich->dirich_val.a.dv[j];               
/*=====================================================================*
 |    example: rotating hole (dirich_val.a.dv == 90)                   |
 |    sonst: Normalfall:                                               |
 |     actnode->sol.a.da[0][j] = initval*acttimefac;                   |
 *=====================================================================*/
      if (initval != 90)
        actnode->sol.a.da[0][j] = initval*acttimefac;
      else
      {
	cx = actnode->x[0]; 
	cy = actnode->x[1];
	win = (initval * acttimefac * 3.14159265359)/180.0;
	wino= atan(cy/cx);
	dd = sqrt(cx*cx+cy*cy);
	if(cx < 0.0) wino += 3.14159265359;
        if (j==0)
	  actnode->sol.a.da[0][j] = dd * cos(win+wino) - cx;
        else
	  actnode->sol.a.da[0][j] = dd * sin(win+wino) - cy;
      }
/*=====================================================================*/
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_settdirich*/





/*----------------------------------------------------------------------*
 |  routine to calculate the element dirichlet load vector              |
 |                                                  genk 05/02          |
 *----------------------------------------------------------------------*/
void ale_caldirich(
                     ELEMENT   *actele, 
		     double    *fullvec,
		     int        dim,
                     ARRAY     *estif_global
		    )     
{

int                   i,j;
int                   dof;
int                   numdf;
int                   nd=0;
double              **estif;
double                dirich[MAXDOFPERELE];
int                   dirich_onoff[MAXDOFPERELE];
double                dforces[MAXDOFPERELE];
GNODE                *actgnode;
NODE                 *actnode;
int                   lm[MAXDOFPERELE];

#ifdef DEBUG 
dstrc_enter("ale_caldirich");
#endif  
/*----------------------------------------------------------------------*/
estif  = estif_global->a.da;
/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
   dforces[i] = 0.0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                                 dirichlet values at (n) were already *
/*                                     written to the nodes (sol[0][j]) */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];   
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++)
   {
      lm[i*numdf+j] = actele->node[i]->dof[j];
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actnode->sol.a.da[0][j];
   }
}
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
/*-------- now assemble the vector dforces to the global vector fullvec */
for (i=0; i<nd; i++)
{
   if (lm[i] >= dim) continue;
   fullvec[lm[i]] += dforces[i];
}
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_caldirich*/ 
#endif
