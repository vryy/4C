/*!----------------------------------------------------------------------
\file
\brief postprocessing 

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
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
int    numfld;   /* number of fields                                    */
int    numf;     /* number of fluid field                               */
int    nums;     /* number of structural field                          */
int    numa;     /* number of ale field                                 */
int    actnum;   /* number of field, which will be visualised           */
FIELD *actfield; /* actual field                                        */

#ifdef DEBUG 
dstrc_enter("ntavisual");
#endif

/*------------------------do initial partitioning of nodes and elements */
part_fields();

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
   numfld = genprob.numfld;
   switch (numfld)
   {
   case 1: /* single field problem */
      actfield=&(field[0]);
      if (actfield->fieldtyp==fluid)
      {
         printf("   Visualisation of a single field problem: FLUID\n");
	 vis2caf(0); 
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
      dserror("Visualisation of multi field problem not implemented yet!\n");
/*      printf("Visualisation of a multi field problem:\n");
      printf("\n");
      printf("Number of fields: %d\n",numfld);
      printf("\n");
      for (i=0;i<numfld;i++)
      {
         actfield=&(field[i]);
         if (actfield->fieldtyp==fluid)
         break;
	 numf=i;
      }
      printf("Actual number of FLUID field:      %d\n",numf);
      for (i=0;i<numfld;i++)
      {
         actfield=&(field[i]);
         if (actfield->fieldtyp==strcuture)
         break;
	 nums=i;
      }
      printf("Actual number of STRUCTURE field:  %d\n",nums); 
      for (i=0;i<numfld;i++)
      {
         actfield=&(field[i]);
         if (actfield->fieldtyp==strcuture)
         break;
	 numa=i;
      }
      printf("Actual number of ALE field:        %d\n",numa);            
      printf("\n");
      printf("Which field do you want to visualise?\n);
      scanf(&actnum);
      if(actnum==numf) vis2caf(numf);
      if(actnum==nums) vis2cas(nums);
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
