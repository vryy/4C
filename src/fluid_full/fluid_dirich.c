/*!----------------------------------------------------------------------
\file
\brief setting dirichlet conditions for fluid

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
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
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;

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
    'actnode->sol_increment.a.da[1][j] = initval*acttimefac'
                                 |        |         |
                            time (n)      |         |               
                          initial value from input  |               
                                       factor from timecurve (T=0.0)

</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)    

\return void                                            

------------------------------------------------------------------------*/
void fluid_initdirich(  FIELD          *actfield, 
                        FLUID_DYNAMIC  *fdyn
		     )
{
int        i,j;
int        numnp_total;               /* total number of fluid nodes    */
int        numele_total;              /* total number of fluid elements */
int        predof;	              /* number of pressure dof	        */
int        numdf;	              /* number of fluid dofs	        */
int        actcurve;	              /* actual timecurve  	        */
int        numveldof;
double     dens;	              /* density			*/
double     timefac[MAXTIMECURVE];     /* factors from time-curve        */
double     T=0.0;	              /* starting time		        */
double     acttimefac;                /* actual factor from timecurve   */
double     initval;	              /* intial dirichlet value	        */
GNODE     *actgnode;	              /* actual GNODE		        */
NODE      *actnode;	              /* actual NODE		        */
ELEMENT   *actele;	              /* actual ELEMENT		        */

int counter=0;

#ifdef DEBUG 
dstrc_enter("fluid_initdirich");
#endif  

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
            initval  = actgnode->dirich->dirich_val.a.dv[j];               
            actnode->sol_increment.a.da[1][j] = initval*acttimefac;
	    actnode->sol.a.da[0][j] = initval*acttimefac;
         } /* end loop over dofs */
      break;
      case dirich_FSI: /* FSI --> dirichvalues = grid velocity!!! */
	 for (j=0;j<numveldof;j++) /* loop vel-dofs */
	 {
	    initval = actnode->sol_increment.a.da[4][j];  
	    actnode->sol_increment.a.da[1][j] = initval; 
	    actnode->sol.a.da[0][j] = initval;
	 }
      break;
      default:
         dserror("dirch_type unknown!\n");
      } /* end switch */
   } /*end loop over nodes */   
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
elements are set at time <T=fdyn->time>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[3][j] = initval*acttimefac'
                                 |        |         |
                            time (n+1)    |         |               
                          initial value from input  |               
                                       factor from timecurve
</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)  

\return void     

------------------------------------------------------------------------*/
void fluid_setdirich(   FIELD           *actfield, 
                        FLUID_DYNAMIC   *fdyn
	            )
{
int        i,j;
int        numnp_total;              /* total number of fluid nodes     */
int        numele_total;             /* total number of fluid elements  */
int        numdf;	             /* number of fluid dofs    	*/
int        actcurve;	             /* actual timecurve		*/
int        numveldof;
double     timefac[MAXTIMECURVE];    /* factors from time-curve         */
double     T;		             /* actual time		        */
double     acttimefac;               /* actual factor from timecurve    */
double     initval;	             /* intial dirichlet value	        */
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG 
dstrc_enter("fluid_setdirich");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = fdyn->time;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
} /* end loop over active timecurves */ 

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
         initval  = actgnode->dirich->dirich_val.a.dv[j];               
         actnode->sol_increment.a.da[3][j] = initval*acttimefac;	 
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */     
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[3][j]
	=actnode->sol_increment.a.da[4][j];
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
} /* end of fluid_settdirich*/

/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector

<pre>                                                         genk 04/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[3][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= estif[i][j] * dirich[j];
\endcode			     

\param  *actele    ELEMENT   (i)   actual element	  
\param  *dforces   double    (o)   dirichlet force vector
\param **estif     double    (i)   element stiffness matrix
\param  *hasdirich int       (o)   flag if s.th. was written to dforces

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_caldirich(
                        ELEMENT         *actele,  
		        double          *dforces, 
                        double         **estif,   
		        int             *hasdirich
		    )     
{

int         i,j;
int         nrow;
int         numdf;                      /* number of fluid dofs         */
int         nd=0;                      
double      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
int         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele  */ 
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
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                           written to the nodes (sol_increment[3][j]) */
nrow=0;
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];   
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[nrow+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[nrow+j] = actnode->sol_increment.a.da[3][j];
   } /* end loop over dofs */
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

#endif
/*! @} (documentation module close)*/
