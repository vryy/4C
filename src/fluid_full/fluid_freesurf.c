/*!----------------------------------------------------------------------
\file
\brief setting free surface conditions for fluid and ale

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
#include "../solver/solver.h"
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
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;
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
#ifdef D_FSI
INT	  i,j;				  /* simply some counters	*/
INT	  freesurf;
INT       dim;
INT       numdf;
INT	  hasdirich,hascouple;
INT	  hasfsi,hasneum;		  /* different flags		*/
INT       hasfreesurf;
DLINE    *actdline;                       /* actual DLINE               */
DNODE    *actdnode;                       /* actual DNODE               */
DSURF    *actdsurf;                       /* actual DSURF               */
FIELDTYP  fieldtyp;

#ifdef DEBUG
dstrc_enter("fluid_createfreesurf");
#endif

/* quit if there's no ALE field */
if (genprob.numaf==-1 || genprob.numfld==1) goto end;


dsassert(genprob.numff>=0,"No fluid field in fluid function!\n");
dsassert(genprob.numaf>=0,"No ale field in fluid freesurf function!\n");

freesurf = alldyn[genprob.numff].fdyn->freesurf;
numdf = alldyn[genprob.numff].fdyn->numdf;
dim = numdf-1;

if (freesurf==0) goto end;

/*--------------------------------------------------------- loop dsurfs */
for (i=0; i<design->ndsurf; i++)
{
   hasdirich=0;
   hascouple=0;
   hasfsi   =0;
   hasneum  =0;
   hasfreesurf=0;
   actdsurf = &(design->dsurf[i]);
   /*--------------------------------------------- check for conditions */
   if (actdsurf->dirich!=NULL)
   {
      hasdirich++;
      /*--------------- for implicit free surface copy dirich condition
                                                       to the grid dofs */
      if (freesurf==2 || freesurf==6)
      for (j=0;j<dim;j++)
      {
         dsassert(actdsurf->dirich->dirich_onoff.fdim>j+numdf,
            "index not allowed for array!\n");
         actdsurf->dirich->dirich_onoff.a.iv[j+numdf]
            =actdsurf->dirich->dirich_onoff.a.iv[j];
         dsassert(actdsurf->dirich->dirich_val.fdim>j+numdf,
            "index not allowed for array!\n");
         actdsurf->dirich->dirich_val.a.dv[j+numdf]
            =actdsurf->dirich->dirich_val.a.dv[j];
      }
   }
   if (actdsurf->couple!=NULL) hascouple++;
   if (actdsurf->fsicouple!=NULL) hasfsi++;
   if (actdsurf->neum!=NULL) hasneum++;
   if (actdsurf->freesurf!=NULL) hasfreesurf++;
   if (hasfreesurf==0) continue;
   fieldtyp=actdsurf->freesurf->fieldtyp;

   switch (fieldtyp)
   {
   case fluid:
      /*------------------ just check if this is really a free surface! */
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DSURF\n");
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DSURF\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DSURF\n");
   break;
   case ale:
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DSURF\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DSURF\n");
      dsassert(hasfsi==0,"fsi- and freesurface condition defined on same DSURF\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DSURF\n");
      /*-------- allocate space for a dirichlet condition in this dsurf */
      actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
      amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdsurf->dirich->dirich_onoff));
      amzero(&(actdsurf->dirich->dirich_val));
      amzero(&(actdsurf->dirich->curve));
      amzero(&(actdsurf->dirich->funct));
      /*----------------------------------- initialise for free surface */
      for (j=0;j<dim;j++)
         actdsurf->dirich->dirich_onoff.a.iv[j]=1;
      actdsurf->dirich->dirich_type=dirich_freesurf;
   break;
   case structure:
   break;
   default:
       dserror("fieldtyp unknown!");
   }
} /* end of loops over dsurfs */
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
   if (actdline->dirich!=NULL)
   {
      hasdirich++;
      /*--------------- for implicit free surface copy dirich condition
                                                       to the grid dofs */
      if (freesurf==2 || freesurf==6)
      for (j=0;j<dim;j++)
      {
         dsassert(actdline->dirich->dirich_onoff.fdim>j+numdf,
            "index not allowed for array!\n");
         actdline->dirich->dirich_onoff.a.iv[j+numdf]
            =actdline->dirich->dirich_onoff.a.iv[j];
         dsassert(actdline->dirich->dirich_val.fdim>j+numdf,
            "index not allowed for array!\n");
         actdline->dirich->dirich_val.a.dv[j+numdf]
            =actdline->dirich->dirich_val.a.dv[j];
      }
   }
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
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amzero(&(actdline->dirich->dirich_val));
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->curve));
      amdef("funct",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->funct));
      /*----------------------------------- initialise for free surface */
      for (j=0;j<dim;j++)
         actdline->dirich->dirich_onoff.a.iv[j]=1;
      actdline->dirich->dirich_type=dirich_freesurf;
   break;
   case structure:
   break;
   default:
       dserror("fieldtyp unknown!");
   break;
   }
} /* end of loops over dlines */
/*--------------------------------------------------------- loop dnodes */
for (i=0; i<design->ndnode; i++)
{
   hasdirich=0;
   hascouple=0;
   hasfsi   =0;
   hasneum  =0;
   hasfreesurf=0;
   actdnode = &(design->dnode[i]);
   /*--------------------------------------------- check for conditions */
   if (actdnode->dirich!=NULL)
   {
      hasdirich++;
      /*--------------- for implicit free surface copy dirich condition
                                                       to the grid dofs */
      if (freesurf==2 || freesurf==6)
      for (j=0;j<dim;j++)
      {
         dsassert(actdnode->dirich->dirich_onoff.fdim>j+numdf,
            "index not allowed for array!\n");
         actdnode->dirich->dirich_onoff.a.iv[j+numdf]
            =actdnode->dirich->dirich_onoff.a.iv[j];
         dsassert(actdnode->dirich->dirich_val.fdim>j+numdf,
            "index not allowed for array!\n");
         actdnode->dirich->dirich_val.a.dv[j+numdf]
            =actdnode->dirich->dirich_val.a.dv[j];
      }
   }
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
      amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdnode->dirich->dirich_onoff));
      amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amzero(&(actdnode->dirich->dirich_val));
      amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
      amzero(&(actdnode->dirich->curve));
      amdef("funct",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdnode->dirich->funct));
      /*----------------------------------- initialise for free surface */
      for (j=0;j<dim;j++)
         actdnode->dirich->dirich_onoff.a.iv[j]=1;
      actdnode->dirich->dirich_type=dirich_freesurf;
   break;
   case structure:
   break;
   default:
       dserror("fieldtyp unknown!");
   break;
   }
}
end:
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
2D local lagrange explicit:       numdf=3
2D local lagrange part. implicit: numdf=5
2D heightfunction separate:       numdf=3 + hfdof
</pre>


\return void

------------------------------------------------------------------------*/
void fluid_freesurf_setdofs()
{
#ifdef D_FSI
INT               i,j;
FIELD            *fluidfield;
FLUID_DYNAMIC    *fdyn;
ELEMENT          *actele;
NODE             *actnode;
GNODE            *actgnode;
DISCRET          *actdis;
INT               numdf;
INT               numnp;
INT               numele;
INT               counter=0;
INT               counter2=0;
INT               myrank;
INT               tmp;

#ifdef DEBUG
dstrc_enter("fluid_freesurf_setdofs");
#endif

/*------------------------------------- leave if there's no fluid field */
dsassert(genprob.numff>=0,"No fluid field in fluid function!\n");

fluidfield = &(field[genprob.numff]);
fdyn = alldyn[genprob.numff].fdyn;
numdf = fdyn->numdf;
actdis = &(fluidfield->dis[0]);
numnp  = actdis->numnp;
numele = actdis->numele;
myrank = par.myrank;


/*------------------------------------------------- plausability checks */
if (fdyn->freesurf==2 || fdyn->freesurf==6)
{
   tmp = MAXDOFPERNODE-5;
   if (numdf==3)
   {
      tmp = MAXDOFPERNODE-5;
      if (tmp<0)
         dserror("MAXDOF PER NODE < 5 !!!\n");
   }
   if (numdf==4)
   {
      tmp = MAXDOFPERNODE-7;
      if (tmp<0)
         dserror("MAXDOF PER NODE < 7 !!!\n");
   }
}

/*-------------------------------------- initialise heightfunction dofs */
for (i=0;i<numnp;i++)
{
   actnode = &(actdis->node[i]);
   actnode->hfdof=NULL;
}

/*-------------------------------------------------------- free surface */
if (fdyn->freesurf>0)
{
   /*---------------------------------------------------- loop elements */
   for (i=0;i<numele;i++)
   {
      actele = &(actdis->element[i]);
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
               if(fdyn->freesurf==1) /* loc. lag. expl. */
	          actele->e.f2->fs_on=1;
	       else if (fdyn->freesurf==2) /* loc. lag. part. impl. */
	       {
                  actnode->numdf=5;
	          actele->e.f2->fs_on=2;
               }
	       else if (fdyn->freesurf==3)
	       {
	          if (actnode->hfdof==NULL)
		  {
		     actnode->hfdof=(INT*)CCACALLOC(1,sizeof(INT));
                     actnode->hfdof[0]=counter;
		     actele->e.f2->fs_on=3;
		     counter++;
		     if (actnode->proc==myrank)
		     counter2++;
		  }
	       }
	       else if (fdyn->freesurf==5)
	       {
	          actnode->numdf=4;
		  actele->e.f2->fs_on=5;
	       }
               else if (fdyn->freesurf==6)
               {
                  actnode->numdf=5;
	          actele->e.f2->fs_on=6;
               }
	       else
	          dserror("Parameter freesurf out of range\n");
#endif
	    break;
            case el_fluid3:
#ifdef D_FLUID3
               if (fdyn->freesurf==1)
	          actele->e.f3->fs_on=1;
	       else if (fdyn->freesurf==2)
	       {
   	          actnode->numdf=7;
                  actele->e.f3->fs_on=2;
	       }
	       else if (fdyn->freesurf==3)
	       {
	          if (actnode->hfdof==NULL)
		  {
		     actnode->hfdof=(INT*)CCACALLOC(1,sizeof(INT));
                     actnode->hfdof[0]=counter;
		     actele->e.f3->fs_on=3;
		     counter++;
		     if (actnode->proc==myrank)
		     counter2++;
		  }
	       }
               else
	          dserror("Parameter freesurf out of range\n");
#endif
            break;

#ifdef D_FLUID3_F
            case el_fluid3_fast:
               if (fdyn->freesurf==1)
                  actele->e.f3->fs_on=1;
               else if (fdyn->freesurf==2)
               {
                  actnode->numdf=7;
                  actele->e.f3->fs_on=2;
               }
               else if (fdyn->freesurf==3)
               {
                  if (actnode->hfdof==NULL)
                  {
                     actnode->hfdof=(INT*)CCACALLOC(1,sizeof(INT));
                     actnode->hfdof[0]=counter;
                     actele->e.f3->fs_on=3;
                     counter++;
                     if (actnode->proc==myrank)
                     counter2++;
                  }
               }
               else
                  dserror("Parameter freesurf out of range\n");
            break;
#endif

            default:
               dserror("no fluid element in fluid field!\n");
            }
         }
      }
   }
   fdyn->hf_numeq_total=counter;
   fdyn->hf_numeq=counter2;
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
DOUBLE  y,x,z;
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
   /*------------------------------------------ loop nodes of alefield */
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

if (strncmp(allfiles.title[0],"f3_freeosz",10)==0)
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
      z=actnode->x[2];
      /*----------------------------------------------- modification:
     	modify all z-coordinates; amplitude depends on old
     	z-coordinate: 0<=z<=1					       */
      actnode->x[2]=z-z*a*sin(PI*x);
   }
   /*------------------------------------------ loop nodes of alefield */
   for (j=0;j<alefield->dis[0].numnp;j++)
   {
      actnode  = &(alefield->dis[0].node[j]);
      x=actnode->x[0];
      z=actnode->x[2];
      /*----------------------------------------------- modification:
     	modify all y-coordinates; amplitude depends on old
     	y-coordinate: 0<=z<=1					       */
      actnode->x[2]=z-z*a*sin(PI*x);
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
   H	  = ONE/FIVE*d;
   fac1    = sqrt((THREE*H)/(FOUR*d*d*d));
   /*-------------------------------------- loop nodes of fluid field */
   for (j=0;j<fluidfield->dis[0].numnp;j++)
   {
      actnode=&(fluidfield->dis[0].node[j]);
      x   = actnode->x[0];
      if (x > 80.0) continue;
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
      if (x > 80.0) continue;
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

/*!---------------------------------------------------------------------
\brief update coordinates at free surface

<pre>                                                         genk 01/04

   In this routine the position of the free surface is determined:
   flag = 0: predictor for next time step
   flag = 1: update during the nonlinear iteration

</pre>
\param *fluidfield FIELD          actual field               (i)
\param *fdyn	   FLUID_DYNAMIC                             (i)
\param  dt         DOUBLE         time increment             (i)
\param  flag       INT            evaluation flag            (i)

\return void

------------------------------------------------------------------------*/
void fluid_updfscoor(FIELD *fluidfield, FLUID_DYNAMIC *fdyn,
                     DOUBLE dt, INT flag)
{
#ifdef D_FSI
INT    j,k;
INT    numnp_total;
INT    numdf;
INT    dim;
INT    phipos;
DOUBLE pred;
NODE  *actfnode, *actanode;
GNODE *actfgnode;
static ARRAY   val1_a,val2_a;
static DOUBLE *val1,*val2;

#ifdef DEBUG
dstrc_enter("fluid_updfscoor");
#endif

numdf = fdyn->numdf;
dim   = numdf-1;
phipos= dim-1;
numnp_total = fluidfield->dis[0].numnp;

switch (flag)
{
case -1: /* predict coordinates at the very beginning */
   val1=amdef("val1",&val1_a,MAXDOFPERNODE,1,"DV");
   val2=amdef("val2",&val2_a,MAXDOFPERNODE,1,"DV");
   switch (fdyn->freesurf)
   {
   case 1: case 2: case 6:
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         for (k=0;k<dim;k++)
         {
            actfnode->xfs[k] = actfnode->x[k] + actanode->sol_mf.a.da[1][k]
	                      + actfnode->sol_increment.a.da[ipos.velnp][k]*dt;
         }
      }
   break;
   case 3: case 5:
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         dsassert(actfnode->locsysId==0,
          "no local co-system at free surface allowed!\n");
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         actfnode->xfs[phipos] = actfnode->x[phipos] + actanode->sol_mf.a.da[1][phipos]
	                       + actfnode->sol_increment.a.da[ipos.velnp][phipos]*dt;
      }
   break;
   default:
      dserror("don't know how to predict free surface!\n");
   }
break;
case 0: /* predict coordinates for the next step                         */
   switch (fdyn->freesurf)
   {
   case 1: /* local lagrange, explicit  */
      /*-------------------------------------- loop nodes of fluid field */
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         for (k=0;k<dim;k++)
         {
            actfnode->xfs[k] = actfnode->x[k] + actanode->sol_mf.a.da[1][k]
	                      + actfnode->sol_increment.a.da[ipos.velnp][k]*dt;
         }
      }
   break;
   case 2: case 6: /* local lagrange, part. implicit */
      /*-------------------------------------- loop nodes of fluid field */
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         for (k=0;k<dim;k++)
         {
            /*------------------- correct ALE solution at free surface */
            actanode->sol_mf.a.da[1][k] = actfnode->xfs[k]-actfnode->x[k];
            /*------- predict free surface position for next time step */
	    actfnode->xfs[k] += actfnode->sol_increment.a.da[ipos.velnp][k+numdf]*dt;

/*------------------------------------------------------------- old stuff
#if 0
             actfnode->xfs[k] = actfnode->x[k] + actanode->sol_mf.a.da[1][k]
	                     + actfnode->sol_increment.a.da[ipos.velnp][k+numdf]*dt;
#endif
#if 0
            actfnode->xfs[k] += actfnode->sol_increment.a.da[ipos.velnp][k+numdf]*dt;
#endif
*/
         }
      }
   break;
   case 3: case 5:/* height function separat & implicit */
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         dsassert(actfnode->locsysId==0,
          "no local co-system at free surface allowed!\n");
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         /*------------------------- correct ALE solution at free surface */
         actanode->sol_mf.a.da[1][phipos] = actfnode->xfs[phipos]-actfnode->x[phipos];
         /*------------- predict free surface position for next time step */
#if 0
         /*------------- precitor based on fluid velocity leads to local instabilities in the
            velocity field if the velocity at the free surface is tangential to the
            heightfunction */
	 actfnode->xfs[phipos] += actfnode->sol_increment.a.da[ipos.velnp][phipos]*dt;

         /* better: predictor based on phi, however this only reduces the effect */
         actfnode->xfs[phipos] += (actfnode->sol_increment.a.da[ipos.velnp][numdf]-
                                   actfnode->sol_increment.a.da[ipos.veln][numdf]);
#endif
         /* predictor of second order: */
         pred = (THREE/TWO*actfnode->sol_increment.a.da[ipos.velnp][numdf]
               - TWO*actfnode->sol_increment.a.da[ipos.veln][numdf]
               +ONE/TWO*actfnode->sol_increment.a.da[ipos.velnm][numdf]);
         actfnode->xfs[phipos] += pred;

      }
   break;
   default:
      dserror("don't know how to update free surface!\n");
   }
break;
case 1: /* update coordinates: correction during the fluid iteration */
   switch (fdyn->freesurf)
   {
   case 2: case 6: /* local lagrange, partitioned implicit */
      /*-------------------------------------- loop nodes of fluid field */
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         for (k=0;k<dim;k++)
         {
/*------------------------------------------------------------ old stuff
#if 0
            ---------- first order correction of free surface position
            actfnode->xfs[k] = actfnode->x[k] + actanode->sol_mf.a.da[0][k]
	                  + actfnode->sol_increment.a.da[ipos.velnp][k+numdf]*dt;
#endif
*/
             /*------ second order correction of free surface position */
             actfnode->xfs[k] = actfnode->x[k] + actanode->sol_mf.a.da[0][k]
	                + (actfnode->sol_increment.a.da[ipos.velnp][k+numdf]
	                +  actfnode->sol_increment.a.da[ipos.veln][k+numdf])/TWO*dt;
         }
      }
   break;
   case 3: case 5: /* height function separat & implicit */
      /*-------------------------------------- loop nodes of fluid field */
      for (j=0;j<numnp_total;j++)
      {
         actfnode=&(fluidfield->dis[0].node[j]);
         if (actfnode->xfs==NULL) continue;
         dsassert(actfnode->locsysId==0,
          "no local co-system at free surface allowed!\n");
	 actfnode->xfs[phipos] = actfnode->sol_increment.a.da[ipos.velnp][numdf];
      }
   break;
   default:
      dserror("don't know how to update free surface!\n");
   }
break;
default:
   dserror("flag out of range!\n");
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
