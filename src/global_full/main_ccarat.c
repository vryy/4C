/*!---------------------------------------------------------------------
\file
\brief main routine

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/compile_settings.h"
#ifndef DEFINE_STRING
#define DEFINE_STRING "\n\tunknown"
#endif
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
struct _PAR     par;
/*!---------------------------------------------------------------------

\brief main routine                                              

<pre>                                                        m.gee 8/00 
main is only printing the ccarat head and the finish
</pre>
\param argc     INT     (i)   number of arguments on command line including exe                                
\param argv     *char[] (i)   array of arguments from command line
\return void                                               

------------------------------------------------------------------------*/
INT main(INT argc, char *argv[])
{
static char release[13] = "01_20040623 ";
#ifdef PARALLEL 
static char *buff,*dbuff;
INT          buffsize=MPIBUFFSIZE;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
MPI_Comm_size(MPI_COMM_WORLD, &par.nprocs);
/*------------------------------------------------ attach buffer to mpi */
buff = (char*)malloc(buffsize);
if (!buff) 
{
   printf("Allocation of memory for mpi buffer failed");
   MPI_Finalize();
   exit(1);
}   
MPI_Buffer_attach(buff,buffsize);
#else
par.myrank=0;
par.nprocs=1;
#endif
if (par.myrank==0)
{
   printf("\n");
   printf("****************************************\n");
   printf("*                                      *\n");
   printf("*             C C A R A T              *\n");
   printf("*                                      *\n");
   printf("*                                      *\n");
#ifdef PARALLEL 
   printf("*           parallel version           *\n");
#else
   printf("*          sequential version          *\n");
#endif
   printf("*                                      *\n");
   printf("*         Release: %s        *\n",release);
   printf("*                                      *\n");
   printf("*       Institut fuer Baustatik        *\n");
   printf("*        Universitaet Stuttgart        *\n");
   printf("*                                      *\n");
   printf("*  Lehrstuhl fuer Numerische Mechanik  *\n");
   printf("*   Technische Universitaet Muenchen   *\n");
   printf("*                                      *\n");
   printf("*    (c) 2004 All Rights Reserved.     *\n");
   printf("*                                      *\n");
   printf("****************************************\n\n");
#ifdef PARALLEL 
   printf("number of processors: %d\n",par.nprocs);
#endif
}

if ((argc == 2) && (strcmp(argv[1], "-v") == 0)) {
  if (par.myrank==0) {
    printf("\nDefine flags used to build ccarat:\n%s\n\n", DEFINE_STRING);
  }
}
else {
/*----------------------------------------------- everything is in here */
ntam(argc,argv);
/*----------------------------------------------------------------------*/
}

#ifdef PARALLEL
MPI_Barrier(MPI_COMM_WORLD);
printf("processor %d finished normally\n",par.myrank);
MPI_Buffer_detach(&dbuff,&buffsize);
if (dbuff!=buff || buffsize != MPIBUFFSIZE)
dserror("Illegal modification of mpi buffer adress or size appeared"); 
free(dbuff);
MPI_Finalize();
#else
printf("processor %d finished normally\n",par.myrank);
#endif
return(0);
} /* end of main */
