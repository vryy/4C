#include "../headers/standardtypes.h"
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
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB  genprob;
/*----------------------------------------------------------------------*
 | open all files                                                       |
 |                                                                      |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntadev(int argc, char *argv[])
{
char  *charpointer; 
char  *resptr;  
int    isrestart;
int    length;
int    lasthandle,nrecords;
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
/*-------------------------------------------------- check for restart */
     genprob.restart=0;
     if (argc > 3)
     {
       resptr = argv[3];
       length = strlen(argv[3]);
       if (length == 7)
       if (strncmp("restart",resptr,length)==0) genprob.restart++;
     } 
/*----------------- in case of restart open the pss-file to apend mode */
if (genprob.restart)
{
     if ( (allfiles.out_pss=fopen(allfiles.outputfile_name,"a+b"))==NULL)
     {
        printf("Opening of restart file .pss failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    printf("binary is written   to     %s\n",allfiles.outputfile_name);
    printf("restart           from     %s\n",allfiles.outputfile_name);
    /*----------------------- set file pointer to the end of pss file */
    fseek(allfiles.out_pss,0,SEEK_END);
    /*---------------- write status report about pss-file to err file */
#ifdef DEBUG 
    pss_status_to_err();
#endif
    /*----------------------------------- get the last handle in file */
    pss_last_handle_in_file(&lasthandle,&nrecords);
    /*--- the initial value for the next record on file is therefore :*/
    allfiles.pss_counter = lasthandle + 1;
}
else/*------------------------------- no restart, open a new pss-file */
{
     if ( (allfiles.out_pss=fopen(allfiles.outputfile_name,"w+b"))==NULL)
     {
        printf("Opening of output file .pss failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    printf("binary is written   to     %s\n",allfiles.outputfile_name);
/* 
   processors must have unique handles
   process 0 takes handles 0..99999
   process 1 takes handles 100000..199999
   etc....
   if one processors handles are all used the process starts again 
   at myrank*100000 + nprocs*100000
*/   
    /* initial values: */
    allfiles.pss_counter = par.myrank*100000;
}
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
