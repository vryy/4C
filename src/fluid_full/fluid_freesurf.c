/*!----------------------------------------------------------------------
\file
\brief setting free surface conditions for fluid and ale

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
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
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;  
/*!--------------------------------------------------------------------- 
\brief create free surface condition

<pre>                                                         genk 01/03

in this function the DBCs at the free surface are created:
ALE    - free surface = DBC
FLUID  - no boundary condition at all

</pre>


\return void                                            

------------------------------------------------------------------------*/
void fluid_createfreesurf()
{
INT       i;                              /* simply some counters       */
INT       hasdirich,hascouple;
INT       hasfsi,hasneum;                 /* different flags            */
INT       hasfreesurf;
DLINE    *actdline;                       /* actual DLINE               */
DNODE    *actdnode;                       /* actual DNODE               */
FIELDTYP  fieldtyp;                       

#ifdef DEBUG 
dstrc_enter("fluid_createfreesurf");
#endif

#ifdef D_FSI

dsassert(genprob.numff>=0,"No fluid field in fluid function!\n");
dsassert(genprob.numaf>=0,"No ale field in fluid freesurf function!\n");

/*--------------------------------------------------------- loop dlines */
for (i=0; i<design->ndline; i++)
{ 
   hasdirich=0;
   hascouple=0;
   hasfsi   =0;
   hasneum  =0;
   hasfreesurf=0;
   actdline = &(design->dline[i]);
   /*--------------------------------------------- check for conditions */
   if (actdline->dirich!=NULL) hasdirich++;
   if (actdline->couple!=NULL) hascouple++;
   if (actdline->fsicouple!=NULL) hasfsi++;
   if (actdline->neum!=NULL) hasneum++;
   if (actdline->freesurf!=NULL) hasfreesurf++;
   if (hasfreesurf==0) continue;
   fieldtyp=actdline->freesurf->fieldtyp;
   
   switch (fieldtyp)
   {
   case fluid:
      /*------------------ just check if this is really a free surface! */
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DLINE\n");   
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DLINE\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DLINE\n");
   break;
   case ale:
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DLINE\n");   
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DLINE\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DLINE\n");
      /*-------- allocate space for a dirichlet condition in this dline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");  
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV"); 
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
      /*----------------------------------- initialise for free surface */
      actdline->dirich->dirich_onoff.a.iv[0]=1;
      actdline->dirich->dirich_onoff.a.iv[1]=1;
      actdline->dirich->dirich_type=dirich_freesurf;  
   break;
   case structure:
   break;
   default:
       dserror("fieldtyp unknown!");
   break;
   }
} /* end of loops over dlines */
/*--------------------------------------------------------- loop dlines */
for (i=0; i<design->ndnode; i++)
{ 
   hasdirich=0;
   hascouple=0;
   hasfsi   =0;
   hasneum  =0;
   hasfreesurf=0;
   actdnode = &(design->dnode[i]);
   /*--------------------------------------------- check for conditions */
   if (actdnode->dirich!=NULL) hasdirich++;
   if (actdnode->couple!=NULL) hascouple++;
   if (actdnode->fsicouple!=NULL) hasfsi++;
   if (actdnode->neum!=NULL) hasneum++;
   if (actdnode->freesurf!=NULL) hasfreesurf++;
   if (hasfreesurf==0) continue;
   fieldtyp=actdnode->freesurf->fieldtyp;
   
   switch (fieldtyp)
   {
   case fluid:
      /*------------------ just check if this is really a free surface! */
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DNODE\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DNODE\n");   
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DNODE\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DNODE\n");
   break;
   case ale:
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DNODE\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DNODE\n");   
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DNODE\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DNODE\n");
      /*-------- allocate space for a dirichlet condition in this dline */
      actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdnode->dirich) dserror("Allocation of memory failed");  
      amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdnode->dirich->dirich_onoff));
      amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV"); 
      amzero(&(actdnode->dirich->dirich_val));
      amzero(&(actdnode->dirich->curve));
      /*----------------------------------- initialise for free surface */
      actdnode->dirich->dirich_onoff.a.iv[0] = 1;
      actdnode->dirich->dirich_onoff.a.iv[1] = 1; 
      actdnode->dirich->dirich_type=dirich_freesurf;   
   break;
   case structure:
   break;
   default:
       dserror("fieldtyp unknown!");
   break;
   }
}

#else
dserror("FSI-functions not compiled in!\n");
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_createfreesurf*/

/*!--------------------------------------------------------------------- 
\brief dofs at free surface

<pre>                                                         genk 01/03

set the dof numbers at the free surface:
2D explicit: numdf=3
2D implicit: numdf=5

</pre>


\return void                                            

------------------------------------------------------------------------*/
void fluid_freesurf_setdofs()
{
#ifdef D_FSI
INT i,j;
FIELD *fluidfield;
FLUID_DYNAMIC *fdyn;
ELEMENT *actele;
NODE *actnode;
GNODE *actgnode;

#ifdef DEBUG 
dstrc_enter("fluid_freesurf_setdofs");
#endif

/*------------------------------------- leave if there's no fluid field */
dsassert(genprob.numff>=0,"No fluid field in fluid function!\n");

fluidfield = &(field[genprob.numff]);
fdyn = alldyn[genprob.numff].fdyn;

/*-------------------------------------------------------- free surface */
if (fdyn->freesurf>0)
{
   /*---------------------------------------------------- loop elements */
   for (i=0;i<fluidfield->dis[0].numele;i++)
   {
      actele = &(fluidfield->dis[0].element[i]);   
      for (j=0;j<actele->numnp;j++)
      {
         actnode=actele->node[j];
         actgnode=actnode->gnode;
         if (actgnode->freesurf!=NULL)
         {
            switch (actele->eltyp)
            {
            case el_fluid2:   
#ifdef D_FLUID2
               if(fdyn->freesurf==1)
	       actele->e.f2->fs_on=1;
	       else
	       {	 
                  actnode->numdf=5;
	          actele->e.f2->fs_on=2;
               }
#endif
	    break;
            case el_fluid3:
#ifdef D_FLUID3	
               dserror("implicit free surface not implemented yet for fluid3\n");
               if (fdyn->freesurf==1)
	       actele->e.f3->fs_on=1;
	       else 
	       {
	          actnode->numdf=7;
                  actele->e.f3->fs_on=2;
	       }	       
#endif
            break;
            default:
            dserror("no fluid element in fluid field!\n");
            }         
         }
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of fluid_createfreesurf*/
/*!--------------------------------------------------------------------- 
\brief set coordinates at free surface

<pre>                                                         genk 02/03

some examples require a specific shape of the intial fluid domain; the
modifactions are done in this function.

</pre>


\return void                                            

------------------------------------------------------------------------*/
void fluid_modcoor()
{
#ifdef D_FSI
INT     j;                              /* simply some counters       */
DOUBLE  a;
DOUBLE  y,x;
FIELD  *fluidfield;
FIELD  *alefield;
NODE   *actnode;
/*----------------------------------------- variables for solitary wave */
DOUBLE eta,H,d,fac,fac1,sech;

#ifdef DEBUG 
dstrc_enter("fluid_modcoor");
#endif


/*--------------------------- modify coordinates for free oszillation */ 
/*                            (Ramaswamy 1990)                        */
if (strncmp(allfiles.title[0],"freeosz",7)==0)
{
   dsassert(genprob.numff>=0,"No fluid field in freesurf fluid function!\n");
   dsassert(genprob.numaf>=0,"No ale field in fluid freesurf function!\n");
   
   if (par.myrank==0)
   {
      printf("\n");
      printf("Modifying coordinates at free surface ...\n");  
      printf("\n");
   }
   a=ONE/TEN/TEN;
   /*----------------------------------------------------- loop fields */
   fluidfield=&(field[genprob.numff]);
   alefield  =&(field[genprob.numaf]);
   /*---------------------------------------- loop nodes of fluidfield */
   for (j=0;j<fluidfield->dis[0].numnp;j++)
   {
      actnode  = &(fluidfield->dis[0].node[j]);
      x=actnode->x[0];
      y=actnode->x[1];
      /*----------------------------------------------- modification:
     	modify all y-coordinates; amplitude depends on old
     	y-coordinate: 0<=y<=1					       */ 
      actnode->x[1]=y-y*a*sin(PI*x);
   }   
   /*---------------------------------------- loop nodes of fluidfield */
   for (j=0;j<alefield->dis[0].numnp;j++)
   {
      actnode  = &(alefield->dis[0].node[j]);
      x=actnode->x[0];
      y=actnode->x[1];
      /*----------------------------------------------- modification:
     	modify all y-coordinates; amplitude depends on old
     	y-coordinate: 0<=y<=1					       */ 
      actnode->x[1]=y-y*a*sin(PI*x);
   } 
}

/*----------------------------- modify coordinates for solitary wave */
/*                              (Ramaswamy 1990)                     */
if (strncmp(allfiles.title[0],"solwave",7)==0)
{
   dsassert(genprob.numff>=0,"No fluid field in freesurf fluid function!\n");
   dsassert(genprob.numaf>=0,"No ale field in fluid freesurf function!\n");
   
   if (par.myrank==0)
   {
      printf("\n");
      printf("Modifying coordinates at free surface ...\n");  
      printf("\n");
   }
   
   fluidfield=&(field[genprob.numff]);
   alefield =&(field[genprob.numaf]);
   /*--------------------------------------------- set some constants */
   d	  = TEN;
   H	  = ONE/FIVE*d;
   fac1    = sqrt((THREE*H)/(FOUR*d*d*d));
   /*-------------------------------------- loop nodes of fluid field */
   for (j=0;j<fluidfield->dis[0].numnp;j++)
   {
      actnode=&(fluidfield->dis[0].node[j]);
      x   = actnode->x[0];
      y   = actnode->x[1];
      fac = fac1*x;
      /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 where xct=x!    */ 
      sech = cosh(fac);
      sech = ONE/sech/sech;
      eta  = d + H*sech;
      /*------------------------------------------ modify coordinates */
      y = y/d*eta;
      actnode->x[1] = y;  
   }
   /*---------------------------------------- loop nodes of ale field */
   for (j=0;j<alefield->dis[0].numnp;j++)
   {
      actnode=&(alefield->dis[0].node[j]);
      x   = actnode->x[0];
      y   = actnode->x[1];
      fac = fac1*x;
      /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 where xct=x!    */ 
      sech = cosh(fac);
      sech = ONE/sech/sech;
      eta  = d + H*sech;
      /*------------------------------------------ modify coordinates */
      y = y/d*eta;
      actnode->x[1] = y;       
   }
}

/*----------------------- modify coordinates for wavebreaking problem */
/*                        (Radovitzky & Ortiz 1998)                   */
if (strncmp(allfiles.title[0],"wavebreak",9)==0)
{
   dsassert(genprob.numff>=0,"No fluid field in freesurf fluid function!\n");
   dsassert(genprob.numaf>=0,"No ale field in fluid freesurf function!\n");
   
   if (par.myrank==0)
   {
      printf("\n");
      printf("Modifying coordinates at free surface ...\n");  
      printf("\n");
   }
   
   fluidfield=&(field[genprob.numff]);
   alefield =&(field[genprob.numaf]);
   /*--------------------------------------------- set some constants */
   d	  = TEN;
   H	  = FIVE;
   fac1    = sqrt((THREE*H)/(FOUR*d*d*d));
   /*-------------------------------------- loop nodes of fluid field */
   for (j=0;j<fluidfield->dis[0].numnp;j++)
   {
      actnode=&(fluidfield->dis[0].node[j]);
      x   = actnode->x[0];
      if (x > 50.0) continue;
      y   = actnode->x[1];
      fac = fac1*x;
      /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 where xct=x!    */ 
      sech = cosh(fac);
      sech = ONE/sech/sech;
      eta  = d + H*sech;
      /*------------------------------------------ modify coordinates */
      y = y/d*eta;
      actnode->x[1] = y;  
   }
   /*---------------------------------------- loop nodes of ale field */
   for (j=0;j<alefield->dis[0].numnp;j++)
   {
      actnode=&(alefield->dis[0].node[j]);
      x   = actnode->x[0];
      if (x > 50.0) continue;
      y   = actnode->x[1];
      fac = fac1*x;
      /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 where xct=x!    */ 
      sech = cosh(fac);
      sech = ONE/sech/sech;
      eta  = d + H*sech;
      /*------------------------------------------ modify coordinates */
      y = y/d*eta;
      actnode->x[1] = y;       
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of fluid_modcoor*/

#endif
/*! @} (documentation module close)*/
