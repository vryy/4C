/*!----------------------------------------------------------------------
\file
\brief postprocessing

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

-----------------------------------------------------------------------*/
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
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
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
INT    i;        /* simply a counter                                    */
#ifdef VISUAL2_PACKAGE
INT    actnum;   /* field number to visualise                           */
INT    screen;
INT    dummy;
#endif
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

/*---------------------------------- allocate fluid integration data ---*/
alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));
 
switch(genprob.visual)
{
case 2: /* 2D - Problem: Visualisation with VISUAL2*/
/*----------------------------------------------------------------------*/
#ifdef VISUAL2_PACKAGE
/*----------------------------------------------------------------------*/
   printf("\n");
   printf("     Starting Visual2 filter ... \n");
   printf("     ----------------------------\n\n");
   printf("     [ ] = DEFAULT = RETURN \n\n");
/*---------------------------------------------- check number of fields */
   switch (genprob.numfld)
   {
   case 1: /* single field problem */
      actfield=&(field[0]);
      if (actfield->fieldtyp==fluid)
      {
         printf("     Visualisation of a single field problem: FLUID\n");
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
      printf("     Visualisation of a multi field problem:\n");
      printf("     Number of fields: %d\n",genprob.numfld);
      printf("\n");

      input2:
      if (genprob.numsf>=0)
      printf("     Actual number of STRUCTURE field:   %d\n",genprob.numsf);
      if (genprob.numff>=0)
      printf("     Actual number of FLUID field:      [%d]\n",genprob.numff);
      if (genprob.numaf>=0)
      printf("     Actual number of ALE field:         %d\n",genprob.numaf);
      printf("\n");
      printf("     Which field do you want to visualise?\n");
      screen=getchar();
      switch(screen)
      {
      case 10: actnum=genprob.numff; break;
      case 48: actnum=0; dummy=getchar(); break;
      case 49: actnum=1; dummy=getchar(); break;
      case 50: actnum=2; dummy=getchar(); break;
      default:
         printf("\nTry again!\n");
         goto input2;
      }
      if(actnum==genprob.numff)
         vis2caf(genprob.numff,genprob.numaf,genprob.numsf);
      else
         dserror("Visualisation of ale/struct problem not implemented yet!\n");
   } /* end switch(numfld) */
break;
/*----------------------------------------------------------------------*/
#else
   dserror("VISUAL2 package is not compiled in!!!\n");
#endif
/*----------------------------------------------------------------------*/
case 3: /* 3D - Problem : Visualisation with VISUAL3 */
/*----------------------------------------------------------------------*/
#ifdef VISUAL3_PACKAGE
/*----------------------------------------------------------------------*/
   printf("\n");
   printf("     Starting Visual3 filter ... \n");
   printf("     ----------------------------\n\n");
   printf("     [ ] = DEFAULT = RETURN \n\n");

/*---------------------------------------------- check number of fields */
   switch (genprob.numfld)
   {
   case 1: /* single field problem */
      actfield=&(field[0]);
      if (actfield->fieldtyp==fluid)
      {
         printf("     Visualisation of a single field problem: FLUID\n");
	 vis3caf(0,-1,-1);
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
      printf("     Visualisation of a multi field problem:\n");
      printf("     Number of fields: %d\n",genprob.numfld);
      printf("\n");

      input1:
      if (genprob.numsf>=0)
      printf("     Actual number of STRUCTURE field:   %d\n",genprob.numsf);
      if (genprob.numff>=0)
      printf("     Actual number of FLUID field:      [%d]\n",genprob.numff);
      if (genprob.numaf>=0)
      printf("     Actual number of ALE field:         %d\n",genprob.numaf);
      printf("\n");
      printf("     Which field do you want to visualise?\n");
      screen=getchar();
      switch(screen)
      {
      case 10: actnum=genprob.numff; break;
      case 48: actnum=0; dummy=getchar(); break;
      case 49: actnum=1; dummy=getchar(); break;
      case 50: actnum=2; dummy=getchar(); break;
      default:
         printf("\nTry again!\n");
         goto input1;
      }
      if(actnum==genprob.numff)
         vis3caf(genprob.numff,genprob.numaf,genprob.numsf);
      else
         dserror("Visualisation of ale/struct problem not implemented yet!\n");
   } /* end switch(numfld) */
break;
#else
dserror("VISUAL3 package is not compiled in!!!\n");
#endif
default:
   dserror("visualisation mode unknown!\n");
} /* switch(genprob.visual) */

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of ntavisual */
