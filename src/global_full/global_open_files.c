#include "../headers/standardtypes.h"
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
char  *vismode;
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
/*--------------------------------------- check for visualisation mode */ 
/* visulisation mode is called by:
   cca.exe input.dat output0 vis2 */
     genprob.visual=0;
     if ( argc > 3)
     {
       vismode=argv[3];
       length = strlen(argv[3]);
       if (length == 4)
       {
         if (strncmp("vis2",vismode,length)==0) genprob.visual=2;
         if (strncmp("VIS2",vismode,length)==0) genprob.visual=2;
         if (strncmp("vis3",vismode,length)==0) genprob.visual=3;
         if (strncmp("VIS3",vismode,length)==0) genprob.visual=3;
       }
     }
     if (genprob.visual) goto visualisation ;
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
/*------------------------------------------------------open .tur file */
if (par.myrank==0)
{     
     strncpy(charpointer,".tur",4);
     if ( (allfiles.out_tur=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .tur failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("turbulence is written to   %s\n",allfiles.outputfile_name);
}     
/*-----------------------------------------------------open .mon file */
if (par.myrank==0)
{     
     strncpy(charpointer,".mon",4);
     if ( (allfiles.out_mon=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of monitoring file .mon failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("monitoring is written to   %s\n",allfiles.outputfile_name);
} 

/*------------------------------------------------------open .pss file */
/* is opened on all procs */     
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
/* 
   old pss-file is opened read/write/append on allfiles.in_pss
   new pss-file is opened read/write        on allfiles.out_pss
*/   
if (genprob.restart)
{
     strncpy(charpointer,".pss",4);
     if ( (allfiles.in_pss=fopen(allfiles.outputfile_name,"a+b"))==NULL)
     {
        printf("Opening of restart file .pss failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    printf("restart           from     %s\n",allfiles.outputfile_name);
    /*----------------------- set file pointer to the end of pss file */
    fseek(allfiles.in_pss,0,SEEK_END);

     strncpy(charpointer,".respss",7);
     if ( (allfiles.out_pss=fopen(allfiles.outputfile_name,"w+b"))==NULL)
     {
        printf("Opening of restart file .respss failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    if (par.myrank==0)
    printf("binary is written   to     %s\n",allfiles.outputfile_name);
    /*----------------------- set file pointer to the end of pss file */
    fseek(allfiles.out_pss,0,SEEK_END);
    /*---------------- write status report about pss-file to err file */
}
else/*------------------------------- no restart, open a new pss-file */
{
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
    printf("binary is written   to     %s\n",allfiles.outputfile_name);
    allfiles.in_pss=NULL;
}
/*------------------------------------------------------open .plt file */
if (par.myrank==0)
{     
     strncpy(charpointer,".plt",4);
     if ( (allfiles.gnu=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .plt failed\n");
#ifdef PARALLEL 
        MPI_Finalize();
#endif 
        exit(1);
     }
    printf("output is written to       %s\n",allfiles.outputfile_name);
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

visualisation:
if (genprob.visual>0)
{
   /*------------------------------------------ check number of procs */
   if (par.nprocs-1) printf("Visualisation only with one proc!\n");
#ifdef PARALLEL 
       MPI_Finalize();
#endif
   /*------------------------------ open files for visual2 or visual3 */
   if (genprob.visual==2 || genprob.visual==3)
   {
/*---------------------------- at the moment the informations of the 
  input file is used for visualisation with Visual 2 / Visual 3       */
/*-----------------------------get input file name and open input file */
      allfiles.inputfile_name=argv[1];
      if ( (allfiles.in_input=fopen(allfiles.inputfile_name,"r"))==NULL)
      {
         printf("Opening of input file %s failed\n",allfiles.inputfile_name);
         exit(1);
      }
      printf("input is read from           %s\n",allfiles.inputfile_name);   
/*------------------------------------- kenner of pss file and open it */
     allfiles.vispsslength=strlen(argv[2]);
     allfiles.vispssfile_kenner=argv[2];    
     strcpy(allfiles.vispssfile_name,allfiles.vispssfile_kenner);
     charpointer=allfiles.vispssfile_name+strlen(allfiles.vispssfile_kenner);
     charpointer = allfiles.vispssfile_name+strlen(allfiles.vispssfile_name); 

     strncpy(charpointer,".pss",4);
     if ( (allfiles.in_pss=fopen(allfiles.vispssfile_name,"a+b"))==NULL)
     {
        printf("Opening of visual data file .pss failed\n");
     }
    printf("visual data are read from    %s\n",allfiles.vispssfile_name);
   }
/*------------------------------------------------------open .err file */     
     strcpy(allfiles.outputfile_name,allfiles.vispssfile_kenner);
     charpointer=allfiles.outputfile_name+strlen(allfiles.vispssfile_kenner);     
     strncpy(charpointer,".vis.err",8);
     if ( (allfiles.out_err=fopen(allfiles.outputfile_name,"w"))==NULL)
     {
        printf("Opening of output file .err failed\n");
        exit(1);
     }
    printf("errors are reported to       %s\n",allfiles.outputfile_name);   
}
/*----------------------------------------------- end of visualisation */
   }
return;
}/* end of ntadev */
