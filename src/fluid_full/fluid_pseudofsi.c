/*!---------------------------------------------------------------------
\file
\brief setting pseudo fsi-conditions for fluid and ale

<pre>
Maintainer: Steffen Genkinger
            genk@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

---------------------------------------------------------------------*/
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

in this function the DBCs at pseude fsi-interface is created:
FLUID  - dirichlet boundary condition

</pre>


\return void                                            

------------------------------------------------------------------------*/
void fluid_createpseudofsi()
{
INT       i;                              /* simply some counters       */
INT       hasdirich,hascouple;
INT       hasfsi,hasneum;                 /* different flags            */
INT       hasfreesurf;
DLINE    *actdline;                       /* actual DLINE               */
DNODE    *actdnode;                       /* actual DNODE               */
FIELDTYP  fieldtyp;                       

#ifdef DEBUG 
dstrc_enter("fluid_createpseudofsi");
#endif

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
   if (hasfsi==0) continue;
   if (actdline->fsicouple->fsi_typ!=fsi_pseudo) continue;
   fieldtyp=actdline->fsicouple->fieldtyp;   
   switch (fieldtyp)
   {
   case fluid:
      /*------------------ just check if this is really a free surface! */
      dsassert(hasdirich==0,"dirich- and freesurface condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and freesurface condition defined on same DLINE\n");   
      dsassert(hasfreesurf==0,"fsi- and freesurface condition defined on same DLINE\n");
      dsassert(hasneum==0,"neumann- and freesurface condition defined on same DLINE\n");
      /*----------- allocate space for a dirichlet condition in this dline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");  
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV"); 
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
      /*----------------------------------- initialise for fsi-coupling */
      actdline->dirich->dirich_onoff.a.iv[0] = 1;   
      actdline->dirich->dirich_onoff.a.iv[1] = 1;       
      actdline->dirich->dirich_type=dirich_FSI_pseudo;
   break;
   case ale:
      dserror("fieldtyp ale not allowed for pseudo fsi!");   
   break;
   case structure:
      dserror("fieldtyp structure not allowed for pseudo fsi!");   
   break;
   default:
       dserror("fieldtyp unknown!");
   break;
   }
} /* end of loops over dlines */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fsi_creatcoup*/
#endif
