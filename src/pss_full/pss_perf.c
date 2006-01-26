/*!---------------------------------------------------------------------
  \file
  \brief contains performance measuring routines

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

  ---------------------------------------------------------------------*/

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

/*---------------------------------------------------- local prototype */
double cputime_thread(void);



#include "../headers/standardtypes.h"
#include <sys/time.h>


#if defined(HPUX_MUENCH) || defined(LINUX_MUENCH)
#include <unistd.h>
#include <sys/times.h>
#endif

#if defined(SX6) || defined(SX8)
#include <unistd.h>
#include <sys/times.h>
#endif

#if defined(HPUX11) && defined(HPUX_GNU)
#include <unistd.h>
#include <sys/times.h>
#endif

#if defined(SUSE73)
#include <unistd.h>
#include <sys/times.h>
#endif


#ifdef PERF
static DOUBLE begtime[100];
static DOUBLE sumtime[100];

static INT    counter[100];
static char   name[100][25];
static INT    bezug[100];
#endif


#ifdef PERF

/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Initializes all counters with 0 and set the names and 'bezuege'.
  </pre>

  \return void

  ------------------------------------------------------------------------*/
void perf_init_all ()
{
  INT index;

  for (index=0; index<100; index++)
  {
    begtime[index] = 0.0;
    sumtime[index] = 0.0;
    counter[index] = 0;
    bezug[index]   = 0;
    sprintf(name[index],"perf slot %2i",index);
  }

  strcpy(name[0], "TOTAL");                     bezug[0]  = 0;
  strcpy(name[1], "INPUT");                     bezug[1]  = 0;
  strcpy(name[2], "CALCULATION");               bezug[2]  = 0;

  strcpy(name[3], "input control");             bezug[3]  = 1;
  strcpy(name[4], "input design");              bezug[4]  = 1;
  strcpy(name[5], "input design topology");     bezug[5]  = 1;
  strcpy(name[6], "input material");            bezug[6]  = 1;
  strcpy(name[7], "input fields");              bezug[7]  = 1;
  strcpy(name[8], "input detailed topology");   bezug[8]  = 1;
  strcpy(name[9], "input topology fe");         bezug[9]  = 1;
  strcpy(name[10],"input conditions");          bezug[10] = 1;
  strcpy(name[11],"input inherit");             bezug[11] = 1;

  strcpy(name[12],"part fields");               bezug[12] = 0;
  strcpy(name[13],"assign dofs");               bezug[13] = 0;
  strcpy(name[14],"part assign field");         bezug[14] = 0;
  strcpy(name[15],"mask global mat");           bezug[15] = 0;
  strcpy(name[16],"calc matrices");             bezug[16] = 0;
  strcpy(name[17],"assemble matrices");         bezug[17] = 0;
  strcpy(name[18],"assemble rhs");              bezug[18] = 0;
  strcpy(name[19],"exchange dofs");             bezug[19] = 0;
  strcpy(name[20],"solve system");              bezug[20] = 0;

  strcpy(name[40],"fsi_dyn");                   bezug[40] = 0;
  strcpy(name[41],"init");                      bezug[41] = 0;
  strcpy(name[42],"fsi_fluid");                 bezug[42] = 0;
  strcpy(name[43],"fsi_struct");                bezug[43] = 0;
  strcpy(name[44],"fsi_ale");                   bezug[44] = 0;
  strcpy(name[45],"comp relax");                bezug[45] = 0;
  strcpy(name[46],"fsi_fluid-init");            bezug[46] = 42;
  strcpy(name[47],"fsi_fluid-prep");            bezug[47] = 42;
  strcpy(name[48],"fsi_fluid-calc");            bezug[48] = 42;
  strcpy(name[49],"fsi_fluid-solv");            bezug[49] = 42;

  strcpy(name[50],"fsi_fluid-post");            bezug[50] = 42;
  strcpy(name[51],"fsi_fluid-coup");            bezug[51] = 42;
  strcpy(name[52],"fsi_fluid-fina");            bezug[52] = 42;
  strcpy(name[53],"allocate_solmf");            bezug[53] = 41;
  strcpy(name[54],"find coup nodes");           bezug[54] = 41;
  strcpy(name[55],"plausibility");              bezug[55] = 41;
  strcpy(name[56],"loop struct");               bezug[56] = 54;
  strcpy(name[57],"loop ale");                  bezug[57] = 54;
  strcpy(name[58],"");                          bezug[58] = 0;
  strcpy(name[59],"");                          bezug[59] = 0;

  strcpy(name[60],"");                          bezug[60] = 0;
  strcpy(name[61],"");                          bezug[61] = 0;
  strcpy(name[62],"");                          bezug[62] = 0;
  strcpy(name[63],"");                          bezug[63] = 0;
  strcpy(name[64],"");                          bezug[64] = 0;
  strcpy(name[65],"");                          bezug[65] = 0;
  strcpy(name[66],"");                          bezug[66] = 0;
  strcpy(name[67],"");                          bezug[67] = 0;
  strcpy(name[68],"");                          bezug[68] = 0;
  strcpy(name[69],"");                          bezug[69] = 0;

  strcpy(name[71],"out_results");               bezug[71] = 0;
  strcpy(name[72],"restart_write_bin");         bezug[72] = 0;
  strcpy(name[73],"out_gid");                   bezug[73] = 0;
  strcpy(name[74],"restart_write");             bezug[74] = 0;

  strcpy(name[80],"solver");                    bezug[80] = 0;
  strcpy(name[81],"ele_calc");                  bezug[81] = 0;

  strcpy(name[92],"amzero");                    bezug[92] = 0;
  strcpy(name[93],"amdel");                     bezug[93] = 0;
  strcpy(name[94],"amredef");                   bezug[94] = 0;
  strcpy(name[95],"amdef");                     bezug[95] = 0;

  /*----------------------------------------------------------------------*/
  return;
}
#endif /* ifdef PERF */



/*!---------------------------------------------------------------------
  \brief routine to get the exact time

  <pre>                                                        mn 01/04
  Gets the current system (cpu) time.
  </pre>
  \return void

  ------------------------------------------------------------------------*/
DOUBLE perf_time ()
{
  DOUBLE  ret=0;

#ifdef HPUX10
  struct timeval _tstart;
  DOUBLE t1;

  gettimeofday(&_tstart, NULL);
  t1 =  (double)_tstart.tv_sec + (double)_tstart.tv_usec/(1000*1000);
  ret = t1;
#endif

#if defined(HPUX11) && !defined(HPUX_GNU)
  DOUBLE time;

  time = cputime_thread( );
  ret = time;
#endif

#if defined(HPUX11) && defined(HPUX_GNU)
  clock_t c0;
  c0 = clock();
  ret = c0;
#endif


#if defined(HPUXITA) && !defined(HPUX_GNU)
  DOUBLE time;

  time = cputime_thread( );
  ret = time;
#endif

#if defined(HPUXITA) && defined(HPUX_GNU)
  clock_t c0;
  c0 = clock();
  ret = c0;
#endif


#ifdef SUSE73
  DOUBLE clk_tck;
  struct tms buf;

  times(&buf);
  clk_tck = (DOUBLE)sysconf(_SC_CLK_TCK);
  ret = (buf.tms_utime + buf.tms_stime)/clk_tck;
#if 0
  struct timeval _tstart;
  DOUBLE t1;

  gettimeofday(&_tstart, NULL);
  t1 =  (double)_tstart.tv_sec + (double)_tstart.tv_usec/(1000*1000);
  ret = t1;
#endif
#endif

#ifdef WIN
  ret = 0.0;
#endif


#if defined(SX6) || defined(SX8)
  /*ret = clock();*/
  DOUBLE clk_tck;
  struct tms buf;

  times(&buf);
  clk_tck = (DOUBLE)sysconf(_SC_CLK_TCK);
  ret = (buf.tms_utime + buf.tms_stime)/clk_tck;
#endif


#ifdef LINUX_MUENCH
  /*ret = clock();*/
  DOUBLE clk_tck;
  struct tms buf;

  times(&buf);
  clk_tck = (DOUBLE)sysconf(_SC_CLK_TCK);
  if (clk_tck > 0)
    ret = (buf.tms_utime + buf.tms_stime)/clk_tck;
#endif

#ifdef HPUX_MUENCH
  DOUBLE clk_tck;
  struct tms buf;

  times(&buf);
  clk_tck = (DOUBLE)sysconf(_SC_CLK_TCK);
  ret = (buf.tms_utime + buf.tms_stime)/clk_tck;
#endif

#ifdef BULL
  struct timeval _tstart;
  DOUBLE t1;

  gettimeofday(&_tstart, NULL);
  t1 =  (double)_tstart.tv_sec + (double)_tstart.tv_usec/(1000*1000);
  ret = t1/((DOUBLE)sysconf(_SC_CLK_TCK));
#endif

  /*----------------------------------------------------------------------*/
  return (ret);
}


#ifdef PERF
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Initializes one counter with 0.
  </pre>

  \param  index   INT   (i)   index of the counter to initialize
  \return void

  ------------------------------------------------------------------------*/
void perf_init (INT index)
{

  begtime[index] = 0.0;
  sumtime[index] = 0.0;
  counter[index] = 0;
  /*----------------------------------------------------------------------*/
  return;
}



/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  start of the region for one timer.
  </pre>

  \param  index   INT   (i)   index of the counter to use
  \return void

  ------------------------------------------------------------------------*/
void perf_begin (INT index)
{
  DOUBLE perf_time();

#if defined(SX6) || defined(SX8)
  char   name[28];
  sprintf(name, "%2i_%24s",index,name[index]);
  ftrace_region_begin(name);
#else
  begtime[index] = perf_time();
#endif

  /*----------------------------------------------------------------------*/
  return;
}
void perfbeginf (INT *index)
{
  DOUBLE perf_time();

  begtime[*index] = perf_time();

  /*----------------------------------------------------------------------*/
  return;
}
void perfbeginf_ (INT *index)
{
  DOUBLE perf_time();

  begtime[*index] = perf_time();

  /*----------------------------------------------------------------------*/
  return;
}


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  end of the region for one timer.
  </pre>

  \param  index   INT   (i)   index of the counter to use
  \return void


  ------------------------------------------------------------------------*/
void perf_end (INT index)
{
  DOUBLE perf_time();

  DOUBLE end_time, elapsed_time;

#if defined(SX6) || defined(SX8)
  char   name[28];
  sprintf(name, "%2i_%24s",index,name[index]);
  ftrace_region_end(name);
#else
  end_time = perf_time();

  elapsed_time = end_time - begtime[index];

#if defined(HPUX11) && defined(HPUX_GNU)
  elapsed_time = elapsed_time/CLOCKS_PER_SEC;
#endif

  sumtime[index] += elapsed_time;
  counter[index] += 1;
#endif

  /*----------------------------------------------------------------------*/
  return;
}
void perfendf (INT *index)
{
  DOUBLE perf_time();

  DOUBLE end_time, elapsed_time;

  end_time = perf_time();

  elapsed_time = end_time - begtime[*index];

  sumtime[*index] += elapsed_time;
  counter[*index] += 1;

  /*----------------------------------------------------------------------*/
  return;
}
void perfendf_ (INT *index)
{
  DOUBLE perf_time();

  DOUBLE end_time, elapsed_time;

  end_time = perf_time();

  elapsed_time = end_time - begtime[*index];

  sumtime[*index] += elapsed_time;
  counter[*index] += 1;

  /*----------------------------------------------------------------------*/
  return;
}



/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Print the results for all timers.
  - for the sequential case just print the timers to the screen and the
    error file
  - in the parallel case collect the timers on the first processor and
    print them there
  </pre>

  \param  out           (i)   file stream to write to
  \return void

  ------------------------------------------------------------------------*/
void perf_out ()
{
  INT             proc;
  INT             index;

  FILE*           screen = stdout;
  FILE*           err    = allfiles.out_err;

#ifdef PARALLEL
  struct _ARRAY   alltimer_s_a;
  struct _ARRAY   allcounter_s_a;
  DOUBLE        **alltimer_s;
  INT           **allcounter_s;
  struct _ARRAY   alltimer_r_a;
  struct _ARRAY   allcounter_r_a;
  DOUBLE        **alltimer_r;
  INT           **allcounter_r;
  INTRA          *actintra;

  DOUBLE          sum_timer;
  INT             sum_counter;
#endif



  /* for the sequential case just print the timers to the screen and the error file */
#ifndef PARALLEL


  /* print out the headings */
  fprintf(err, "\n%2s %25s: %7s %12s %5s%3s\n",
      "ID",
      "Functionname",
      "Count",
      "Time [sec]",
      "%",
      "of"
      );
  fprintf(err,
      "============================================================\n");

  fprintf(screen, "\n%2s %25s: %7s %12s %5s%3s\n",
      "ID",
      "Functionname",
      "Count",
      "Time [sec]",
      "%",
      "of"
      );
  fprintf(screen,
      "============================================================\n");


  /* loop all counters and print only those that were used */
  for (index=0; index<100; index++)
  {
    if (counter[index] != 0)
    {
      fprintf(err, "%2d %25s: %7d %12.6e %5.1f %2d \n",
          index,
          name[index],
          counter[index],
          sumtime[index],
          (sumtime[bezug[index]] != 0) ? (sumtime[index]/(sumtime[bezug[index]]))*100 : 0,
          bezug[index]
          );

      fprintf(screen, "%2d %25s: %7d %12.6e %5.1f %2d \n",
          index,
          name[index],
          counter[index],
          sumtime[index],
          (sumtime[bezug[index]] != 0) ? (sumtime[index]/(sumtime[bezug[index]]))*100 : 0,
          bezug[index]
          );
    } /* if (counter[index] != 0) */


    /* add separators at some points */
    if (index== 2 || index==11 || index==20 ||
        index==45 || index==52 || index==55 ||
        index==57 || index==74 || index==82 )
    {
      fprintf(err,
          "------------------------------------------------------------\n");
      fprintf(screen,
          "------------------------------------------------------------\n");
    }

  } /* for (index=0; index<100; index++) */


  /* print footer line */
  fprintf(err,
      "============================================================\n");
  fprintf(screen,
      "============================================================\n");




  /* in the parallel case collect the timers on the first processor and print them there */
#else

  actintra    = &(par.intra[0]);


  /* allocate memory for the timers and counters of all procs */
  alltimer_s   = amdef("alltimer_s",  &alltimer_s_a,  par.nprocs+1,100,"DA");
  alltimer_r   = amdef("alltimer_r",  &alltimer_r_a,  par.nprocs+1,100,"DA");
  allcounter_s = amdef("allcounter_s",&allcounter_s_a,par.nprocs+1,100,"IA");
  allcounter_r = amdef("allcounter_r",&allcounter_r_a,par.nprocs+1,100,"IA");

  amzero(&alltimer_s_a);
  amzero(&alltimer_r_a);
  amzero(&allcounter_s_a);
  amzero(&allcounter_r_a);


  /* copy the data into these arrays */
  for (index=0; index<100; index++)
  {
    alltimer_s[par.myrank][index]   = sumtime[index];
    allcounter_s[par.myrank][index] = counter[index];
  }


  /* allreduce the data */
  MPI_Allreduce(alltimer_s[0],  alltimer_r[0],  (par.nprocs+1)*100,
      MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
  MPI_Allreduce(allcounter_s[0],allcounter_r[0],(par.nprocs+1)*100,
      MPI_INT   ,MPI_SUM,actintra->MPI_INTRA_COMM);


  /* print only on the first processor */
  if (par.myrank==0)
  {
    /* print out the headings */
    /* first line */
    fprintf(err,   "\n                              ");
    fprintf(screen,"\n                              ");

    for (proc=0;proc<par.nprocs;proc++)
      fprintf(err, "     Processor %2i:   ", proc);

    fprintf(err,   "            Average:\n");
    fprintf(screen,"            Average:\n");

    /* second line */
    fprintf(err,   "                              ");
    fprintf(screen,"                              ");

    for (proc=0;proc<par.nprocs;proc++)
      fprintf(err, "  ------------------ ", proc);

    fprintf(err,   "  ---------------------------\n");
    fprintf(screen,"  ---------------------------\n");

    /* third line */
    fprintf(err,   "%2s %25s: ","ID","Functionname");
    fprintf(screen,"%2s %25s: ","ID","Functionname");

    for (proc=0;proc<par.nprocs;proc++)
      fprintf(err, "  Count   Time [sec] ", proc);

    fprintf(err,   "%7s %12s %5s%3s\n","Count","Time [sec]","%","of");
    fprintf(screen,"%7s %12s %5s%3s\n","Count","Time [sec]","%","of");

    for (proc=0;proc<par.nprocs;proc++)
      fprintf(err, "=====================");
    fprintf(err,
        "============================================================\n");
    fprintf(screen,
        "============================================================\n");


    /* loop all counters and print only those that were used */
    for (index=0; index<100; index++)
    {
      sum_timer   = 0.0;
      sum_counter = 0;

      if (allcounter_r[0][index] != 0)
      {
        fprintf(err,   "%2d %25s: ",index,name[index]);
        fprintf(screen,"%2d %25s: ",index,name[index]);

        for (proc=0;proc<par.nprocs;proc++)
        {
          fprintf(err, "%7d %12.6e ",
              allcounter_r[proc][index],
              alltimer_r[proc][index]
              );
          sum_timer   += alltimer_r[proc][index];
          sum_counter += allcounter_r[proc][index];
        }

        alltimer_r[par.nprocs][index]   = sum_timer/par.nprocs;
        allcounter_r[par.nprocs][index] = sum_counter/par.nprocs;

        fprintf(err,    "%7d %12.6e %5.1f %2d\n",
            allcounter_r[par.nprocs][index],
            alltimer_r[par.nprocs][index],
            (alltimer_r[par.nprocs][bezug[index]] != 0) ?
                (alltimer_r[par.nprocs][index]/(alltimer_r[par.nprocs][bezug[index]]))*100 : 0,
            bezug[index]
            );
        fprintf(screen, "%7d %12.6e %5.1f %2d\n",
            allcounter_r[par.nprocs][index],
            alltimer_r[par.nprocs][index],
            (alltimer_r[par.nprocs][bezug[index]] != 0) ?
                (alltimer_r[par.nprocs][index]/(alltimer_r[par.nprocs][bezug[index]]))*100 : 0,
            bezug[index]
            );
      } /* if (allcounter_r[0][index] != 0) */


      /* add separators at some points */
      if (index== 2 || index==11 || index==20 ||
          index==45 || index==52 || index==55 ||
          index==57 || index==74 || index==82 )
      {
        for (proc=0;proc<par.nprocs;proc++)
          fprintf(err, "---------------------");
        fprintf(err,
            "------------------------------------------------------------\n");
        fprintf(screen,
            "------------------------------------------------------------\n");
      }

    } /* for (index=0; index<100; index++) */


    /* print footer line */
    for (proc=0;proc<par.nprocs;proc++)
      fprintf(err, "=====================");
    fprintf(err,
        "============================================================\n");
    fprintf(screen,
        "============================================================\n");

    fflush(err);
    fflush(screen);

  } /* if (par.myrank==0) */


  MPI_Barrier(actintra->MPI_INTRA_COMM);

#endif


  /*----------------------------------------------------------------------*/
  return;
}
#endif /* ifdef PERF */




#if (defined(HPUX11) || defined(HPUXITA)) && !defined(HPUX_GNU)
/*
 * A lightweight and thread-safe procedure to return the thread time
 * on the HP V2250 with 1 microsec resolution.
 *
 * Note that the argument is NOT used.  This procedure is
 * generally intended to only return time intervals.  Its
 * usage is;
 *
 *   double begin_time , end_time , elapsed_time;
 *   double cputime_thread() , zero=0.0;
 *
 *   begin_time = cputime_thread( zero );
 *      < code to be timed >
 *   end_time = cputime_thread( zero );
 *   elapsed_time = end_time - begin_time;
 *
 */

#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#define HP_USER_CPU_TIME 1
#define HP_SYSTEM_CPU_TIME 2
#define HP_INTERRUPT_TIME 3
#define HP_MAX_THREADS 1024

/* global data */
static int first_call = 1;

double cputime_thread (void)
{


  /* Prototypes */
  extern int __lw_get_time_factors(double *,time_t *);
  extern void __lw_get_thread_times(int ,int64_t *,int64_t *);

  /* Declarations */
  static double base_time[HP_MAX_THREADS];
  static double sec_per_itick = -1.0;
  static double last_cpu_time[HP_MAX_THREADS];
  static double last_wall_time[HP_MAX_THREADS];

  int64_t cputime;
  int64_t walltime;
  int index;
  time_t boottime;
  int i;
  int return_value;
  double timex;

  /* first call initializations */
  /*    (1) get the tick time from the thread package
   *    (2) set the base time for all possible threads to zero
   *    (3) get the current times (processor and wall-clock)
   *    (4) establish the base time for this thread
   *    (5) return zero as the elapsed time
   *    (6) insure we won't initialize again.
   */

  if ( first_call ) {
    return_value = __lw_get_time_factors( &sec_per_itick , &boottime );
    if ( (return_value < 0) || (sec_per_itick < 0) ) {
      fprintf( stderr , "cputime_thread: __lw_get_time_factors() failed.\n" );
      return( (double) -1 );
    }
    printf("sec_per_itick = %g boottime = %d\n",sec_per_itick,boottime);

    for ( i=0 ; i<HP_MAX_THREADS ; i++ )
      base_time[i] = 0.0;

    cputime = 0;
    walltime = 0;
    __lw_get_thread_times( HP_USER_CPU_TIME , &cputime , &walltime );

#ifndef HPUX11
    index = pthread_self();
    if ( index < 0 ) index = 0;
#else
    index = 0;
#endif

    base_time[index] = ((double) cputime) * sec_per_itick;
    first_call = 0;
    return( (double) 0 );
  } /* first_call */


  /* usual entry - get thread time & subtract off base time */
  cputime = 0;
  walltime = 0;
  __lw_get_thread_times( HP_USER_CPU_TIME , &cputime , &walltime );

#ifndef HPUX11
  index = pthread_self();
  if ( index < 0 )  index = 0;
#else
    index = 0;
#endif

  return( ((double) cputime) * sec_per_itick - base_time[index] );

} /* cputime_thread */

#endif
