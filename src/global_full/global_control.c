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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of partitions, size numfld                                    |
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
/*----------------------------------------------------------------------*
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[])
{
int ierr;
int handle;
double t0,ti,tc;
/*---------------------------------------- init devices, tracing, etc...*/
ntaini(argc,argv);
/*--------------------------------input phase, input of all information */
t0=ds_cputime();
ntainp();
ti=ds_cputime()-t0;
/*------------------------ write output of all preprocessor information */
if (par.myrank==0 && genprob.restart==0)
{
   pss_write("input_file",allfiles.numrows,allfiles.numcol,sizeof(char),
                allfiles.input_file_hook,&handle,&ierr);
   if (ierr!=1) dserror("cannot write input_file to pss");
}
/*--------------------------------close the input file, delete the copy */
/*---------- You cannot read anything from input file beyond this point */
frend();
/*--------------------------------- set up field-specific communicators */
#ifdef PARALLEL 
create_communicators();
#endif
/*----------------------------------------- read information of restart */
if (genprob.restart!=0)
{
   /* nothing implemented yet */
}
/*----------------------------------------------------calculation phase */
t0=ds_cputime();
ntacal();
tc=ds_cputime();
#ifdef PARALLEL 
MPI_Barrier(MPI_COMM_WORLD);
#endif
if (par.myrank==0)
{
 printf("\n");
 printf("Total CPU Time for INPUT:       %10.3#E sec \n",ti);
 printf("Total CPU Time for CALCULATION: %10.3#E sec \n",tc);
 printf("\n"); 
}
 

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of ntam */
    
