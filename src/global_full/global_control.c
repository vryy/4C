#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[])
{
int ierr;
int handle;
/*---------------------------------------- init devices, tracing, etc...*/
ntaini(argc,argv);
/*--------------------------------input phase, input of all information */
ntainp();
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
ntacal();


 

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of ntam */
