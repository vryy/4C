/*!----------------------------------------------------------------------
\file
\brief postprocessing 

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD         *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of partitions, size numfld                                    |
 *----------------------------------------------------------------------*/
struct _PARTITION    *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
extern struct _PAR   par;                     
/*!---------------------------------------------------------------------                                         
\brief call of Visualisation tools

<pre>                                                         genk 07/02       

This routine checks the type of problem and based on the program  options
a visualisation tool is called.
At the moment implemented:
VISUAL2

</pre>  
\return void                                                                       

------------------------------------------------------------------------*/
void ntavisual()
{
int    i;        /* simply a counter                                    */
int    actnum;   /* field number to visualise                           */
FIELD *actfield; /* actual field                                        */

#ifdef DEBUG 
dstrc_enter("ntavisual");
#endif


/*------------------------do initial partitioning of nodes and elements */
part_fields();

#ifdef D_FLUID
/*---------------------------------- set dofs for implicit free surface */
if (genprob.numff>=0) fluid_freesurf_setdofs();
/*----------------------------- modify coordinates for special problems */
if (genprob.numff>=0) fluid_modcoor();
#endif

/*------------------------------------------------ assign dofs to nodes */
for(i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   if (actfield->ndis==1) assign_dof(actfield);
   if (actfield->ndis>1) assign_dof_ndis(actfield);
}
/*--------make the procs know their own nodes and elements a bit better */
part_assignfield();


switch(genprob.visual)
{
case 2: /* 2D - Problem: Visualisation with VISUAL2*/
/*----------------------------------------------------------------------*/
#ifdef VISUAL2_PACKAGE
/*----------------------------------------------------------------------*/
   printf("\n");
   printf("   Starting Visual2 filter ... \n");
   printf("   ----------------------------\n");      
/*---------------------------------------------- check number of fields */
   switch (genprob.numfld)
   {
   case 1: /* single field problem */
      actfield=&(field[0]);
      if (actfield->fieldtyp==fluid)
      {
         printf("   Visualisation of a single field problem: FLUID\n");
	 vis2caf(0,-1,-1); 
      }
      else if (actfield->fieldtyp==structure)
      {
         dserror("visualisation of structural problem not implemented yet!\n");
      }
      else if (actfield->fieldtyp==ale)
      {
         dserror("visualisation of ale problem not implemented yet!\n");
      }
   break;
   default: /* multi field problem */
      printf("\n");
      printf("   Visualisation of a multi field problem:\n");
      printf("   Number of fields: %d\n",genprob.numfld);
      printf("\n");
            
      if (genprob.numsf>=0)
      printf("   Actual number of STRUCTURE field:  %d\n",genprob.numsf); 
      if (genprob.numff>=0)
      printf("   Actual number of FLUID field:      %d\n",genprob.numff);
      if (genprob.numaf>=0)
      printf("   Actual number of ALE field:        %d\n",genprob.numaf);            
      printf("\n");
      printf("   Which field do you want to visualise?\n");
      scanf("%d",&actnum);
      if(actnum==genprob.numff) 
         vis2caf(genprob.numff,genprob.numaf,genprob.numsf);
      else
         dserror("Visualisation of ale/struct problem not implemented yet!\n");
/*      if(actnum==nums) vis2cas(nums);
      if(actnum==numa) vis2caa(numa);               */
   } /* end switch(numfld) */   
break;
/*----------------------------------------------------------------------*/
#else
   dserror("VISUAL2 package is not compiled in!!!\n");
#endif
/*----------------------------------------------------------------------*/   
case 3: /* 3D - Problem : Visualisation with VISUAL3 */
   dserror("VISUAL3 not implemented yet!\n");
break;

default:
   dserror("visualisation mode unknown!\n");
} /* switch(genprob.visual) */

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of ntavisual */
