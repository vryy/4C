#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 | open all files                                                       |
 |                                                                      |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntadev(int argc, char *argv[])
{
char  *charpointer;   
/*------------ cannot use dserror here, because .err files not yet open */
/*--------------------------------------------------cheque command line */
   if (argc <= 1)
   {
      printf("You forgot to give the input and output file names!\n");
      printf("Try again!\n");
#ifdef PARALLEL 
      MPI_Finalize();
#endif 
      exit(1);
   }
   else if (argc <= 2)
   {
      printf("You forgot to give the output file name!\n");
      printf("Try again!\n");
#ifdef PARALLEL 
      MPI_Finalize();
#endif 
      exit(1);
   }    
   else
   {
/*-----------------------------get input file name and open input file */
     if (par.myrank==0)
     {
        allfiles.inputfile_name=argv[1];
        if ( (allfiles.in_input=fopen(allfiles.inputfile_name,"r"))==NULL)
        {
           printf("Opening of input file %s failed\n",allfiles.inputfile_name);
           exit(1);
        }
        printf("input is read from         %s\n",allfiles.inputfile_name);
     }
/*-----------------------get output file kenner and open .out and .err */
     allfiles.outlenght=strlen(argv[2]);
     allfiles.outputfile_kenner=argv[2];
     allfiles.num_outputfiles=13;      
     if (allfiles.outlenght>=40)
     {
        if (par.myrank==0)
        {
           printf("Your outputfile kenner is too long!\n");
           fflush(stdout);
        }
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
     strcpy(allfiles.outputfile_name,allfiles.outputfile_kenner);
     charpointer=allfiles.outputfile_name+strlen(allfiles.outputfile_kenner);
     sprintf(charpointer,"%d",par.myrank);
     charpointer = allfiles.outputfile_name+strlen(allfiles.outputfile_name);


/*------------------------------------------------------open .out file */
if (par.myrank==0)
{     
     strncpy(charpointer,".out",4);
     if ( (allfiles.out_out=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .out failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("output is written to       %s\n",allfiles.outputfile_name);
}     
/*------------------------------------------------------open .err file */
     
     strncpy(charpointer,".err",4);
     if ( (allfiles.out_err=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .err failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    printf("errors are reported to     %s\n",allfiles.outputfile_name);
/*------------------------------------------------------open .pss file */
/* is opened on all procs */     
     strncpy(charpointer,".pss",4);
     if ( (allfiles.out_pss=fopen(allfiles.outputfile_name,"w+b"))==NULL)
     {
        printf("Opening of output file .pss failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    allfiles.pss_counter=0;
    printf("binary is written   to     %s\n",allfiles.outputfile_name);
/*-----------------------------------------------open .flavia.msh file */
if (par.myrank==0)
{     
     strncpy(charpointer,".flavia.msh",11);
     if ( (allfiles.gidmsh=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .out failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("output is written to       %s\n",allfiles.outputfile_name);
}     
/*-----------------------------------------------open .flavia.res file */
if (par.myrank==0)
{     
     strncpy(charpointer,".flavia.res",11);
     if ( (allfiles.gidres=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .out failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("output is written to       %s\n",allfiles.outputfile_name);
}     
/*--------------------------to open any other file just add block here */


   }

return;
}/* end of ntadev */
