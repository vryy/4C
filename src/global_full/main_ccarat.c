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
#include <../headers/compile_settings.h>

#ifdef TRAP_FE

#ifdef LINUX_MUENCH
/*
 * This is to get the GNU prototypes. Is there a better way to state
 * that we want to use them? */
#define __USE_GNU
#include <fenv.h>
#endif

#ifdef HPUX_MUENCH
#include <fenv.h>
#endif

#ifdef HPUXITA
#include <fenv.h>
#endif

#ifdef HPUX11
#include <fenv.h>
#endif

#endif /* TRAP_FE */

/* In case the settings header is brocken */
#ifndef COMPILE_SETTINGS_H
#define DEFINE_STRING "\n\tunknown"
#define CREATOR "unknown"
#define CREATION_DATE "unknown"
#define CONFIGURATION "unknown"
#endif

#define print_define(arg)  if (strstr(DEFINE_STRING, #arg )==NULL) printf("\n\t" #arg "=%d", arg);
#define print_define_dbl(arg)  if (strstr(DEFINE_STRING, #arg )==NULL) printf("\n\t" #arg "=%f", arg);

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
struct _PAR     par;

/* bla */
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
static char release[13] = "03_20050823*";
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
   printf("*    (c) 2005 All Rights Reserved.     *\n");
   printf("*                                      *\n");
   printf("****************************************\n\n");
#ifdef PARALLEL
   printf("number of processors: %d\n",par.nprocs);
#endif
}


if ((argc == 2) && (strcmp(argv[1], "-v") == 0)) {
  if (par.myrank==0) {
    printf("\nBuilt by %s on %s\n", CREATOR, CREATION_DATE);
    printf("using configuration: %s\n", CONFIGURATION);
    printf("\nDefine flags used to build ccarat:\n%s\n", DEFINE_STRING);
    printf("\nDefault values:\n");
    print_define(MAXNOD);
    print_define(MAXELE);
    print_define(MAXDOFPERNODE);
    print_define(MAXGAUSS);
    print_define(MAXFIELD);
    print_define(MAXTIMECURVE);
    print_define(MAXRECORDPERELE);
    print_define(MAXNUMMATRICES);
    print_define(MAXNOD_AXISHELL);
    print_define(MAXNOD_BEAM3);
    print_define(MAXNOD_BRICK1);
    print_define(NUMDOF_BRICK1);
    print_define(MAXQINTC);
    print_define(MAXQINTP);
    print_define(MAXTINTC);
    print_define(MAXTINTP);
    print_define(FLUID_NUM_LD);
    print_define(NUM_F2_VELDOF);
    print_define(NUMDOF_FLUID2);
    print_define(MAXNOD_F2);
    print_define(NUM_F3_VELDOF);
    print_define(MAXNOD_F3);
    print_define(MAXNOD_SHELL8);
    print_define(NUMDOF_SHELL8);
    print_define(MAXHYB_SHELL8);
    print_define(MAXNOD_SHELL9);
    print_define(MAXLAY_SHELL9);
    print_define(MAXKLAY_SHELL9);
    print_define(NUMDOF_SHELL9);
    print_define(MAXHYB_SHELL9);
    print_define_dbl(A3FAC_SHELL9);
    print_define(MAXNODESTRESS_SHELL9);
    print_define(MAXNOD_WALL1);
    printf("\n\n");
  }
}
else {
  /* Here we turn the NaN and inf numbers of. No need to calculate
   * those. If those appear the calculation needs much (!) more
   * time. Better stop immediately if some illegal operation occurs. */
#ifdef TRAP_FE

  /* Sadly, it seems the functions needed for this are different on
   * different maschines. */
#ifdef LINUX_MUENCH

  /* This is a GNU extention thus it's only available on linux. But
   * it's exactly what we want: SIGFPE just for the given
   * exceptions. We don't care about FE_INEXACT. (It happens all the
   * time.) */
  /* Over- and underflow seem to happen sometimes. Does it worry us?
   * Will that spoil the results? */
  /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
  feenableexcept(FE_INVALID | FE_DIVBYZERO);

  /* The hard GNU way. But it does too much. */
  /*fesetenv((fenv_t*)-2);*/

#endif

#ifdef HPUX_MUENCH
  /*
   * Don't ask me why they want this. The man page said it's needed on
   * itanium maschines. */
#pragma STDC FENV_ACCESS ON
  /*fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
  fesettrapenable(FE_INVALID | FE_DIVBYZERO);
#endif

#ifdef HPUXITA
#pragma STDC FENV_ACCESS ON
  fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
#endif

#ifdef HPUX11
#pragma STDC FENV_ACCESS ON
  fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
#endif

#endif /* TRAP_FE */

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
