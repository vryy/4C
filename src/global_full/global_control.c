#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
struct _FIELD         *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
struct _GENPROB       genprob;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
struct _PARTITION    *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
struct _DESIGN       *design;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | struct _ALLDYNA       *alldyn;                                       |
 *----------------------------------------------------------------------*/
ALLDYNA             *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
struct _STATIC_VAR   *statvar;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
struct _MATERIAL     *mat;
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
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[])
{
int      ierr;
long int handle;
double   t0,ti,tc;
/*---------------------------------------- init devices, tracing, etc...*/
ntaini(argc,argv);
/*--------------------------------input phase, input of all information */
t0=ds_cputime();
ntainp();
ti=ds_cputime()-t0;
if (par.myrank==0)
{
 printf("\n");
 printf("Total CPU Time for INPUT:       %10.3#E sec \n",ti);
 printf("\n"); 
}
/*--------------------------------close the input file, delete the copy */
/*---------- You cannot read anything from input file beyond this point */
frend();
/*---------------------------------------- check for visualisation mode */
if (genprob.visual!=0) goto visualisation;
/*--------------------------------- set up field-specific communicators */
#ifdef PARALLEL 
create_communicators();
#endif
/*----------------------------------------------------calculation phase */
t0=ds_cputime();
ntacal();
tc=ds_cputime()-t0;
if (par.myrank==0)
{
 printf("\n");
 printf("Total CPU Time for CALCULATION: %10.3#E sec \n",tc);
 printf("\n"); 
}
goto endcal;
/*----------------------------------------------------------------------*/
visualisation:
if (genprob.visual==2 || genprob.visual==3)
ntavisual();
/*----------------------------------------------------------------------*/
endcal:
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of ntam */
    
