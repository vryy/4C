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
/*---------------------------------------------------- local prototype */
double cputime_thread();

#include "../headers/standardtypes.h"
#include <sys/time.h>


#if defined(HPUX_MUENCH) || defined(LINUX_MUENCH)
#include <unistd.h>
#include <sys/times.h>
#endif

#if defined(SX6)
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
#endif


/*!---------------------------------------------------------------------
  \brief routine to get the exact time

  <pre>                                                        mn 01/04
  Gets the current system time.
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
  DOUBLE cputime_thread() , zero=0.0;

  time = cputime_thread( zero );
  ret = time;
#endif

#if defined(HPUX11) && defined(HPUX_GNU)
  clock_t c0;
  c0 = clock();
  ret = c0;
#endif


#if defined(HPUXITA) && !defined(HPUX_GNU)
  DOUBLE time;
  DOUBLE cputime_thread() , zero=0.0;

  time = cputime_thread( zero );
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


#ifdef SX6
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
  Initializes all counters with 0.
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
  }
  /*----------------------------------------------------------------------*/
  return;
}


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Initializes one counters with 0.
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
  DOUBLE perf_cpu();
  DOUBLE perf_time();

  begtime[index] = perf_time();

  /*----------------------------------------------------------------------*/
  return;
}
void perfbeginf (INT *index)
{
  DOUBLE perf_cpu();
  DOUBLE perf_time();

  begtime[*index] = perf_time();

  /*----------------------------------------------------------------------*/
  return;
}
void perfbeginf_ (INT *index)
{
  DOUBLE perf_cpu();
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

  end_time = perf_time();

  elapsed_time = end_time - begtime[index];

#if defined(HPUX11) && defined(HPUX_GNU)
  elapsed_time = elapsed_time/CLOCKS_PER_SEC;
#endif

  sumtime[index] += elapsed_time;
  counter[index] += 1;

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
  Print the results for one timer.
  </pre>

  \param  out           (i)   file stream to write to
  \param  index   INT   (i)   index of the counter to use
  \param  string  char  (i)   name of this counter
  \param  bezug   INT   (i)   calculate percentage relative to this counter
  \param  ops     INT   (i)   number of FLOP in this counter region
  \return void

  ------------------------------------------------------------------------*/
void perf_print (FILE* out, INT index, char string[], INT bezug, INT ops)
{

  if (ops!=1)
  {
    fprintf(out, "%2d %25s: %7d %12.6e %5.1f %2d %8d %8.3f\n",
        index,
        string,
        counter[index],
        sumtime[index],
        (sumtime[bezug] != 0) ? (sumtime[index]/(sumtime[bezug]))*100 : 0,
        bezug,
        ops,
        ((sumtime[bezug] != 0) && (counter[index] != 0)) ? ops/(sumtime[index]/counter[index]*1.0e6) : 0
        );
  }
  else
  {
    fprintf(out, "%2d %25s: %7d %12.6e %5.1f %2d \n",
        index,
        string,
        counter[index],
        sumtime[index],
        (sumtime[bezug] != 0) ? (sumtime[index]/(sumtime[bezug]))*100 : 0,
        bezug
        );
  }

  /*----------------------------------------------------------------------*/
  return;
}


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Print the results for all timers.
  </pre>

  \param  out           (i)   file stream to write to
  \return void

  ------------------------------------------------------------------------*/
void perf_out (FILE* out)
{
  INT i;
  INT user_perf;

  /* --------------- print out time counters */
  fprintf(out, "%2s %25s: %7s %12s %5s%3s %8s %8s\n",
      "ID",
      "Functionname",
      "Count",
      "Time [sec]",
      "% ",
      " of",
      "FLOP",
      "MFLOPS"
      );

  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out,  1,"INPUT",                   0,1);
  perf_print(out,  2,"CALCULATION",             0,1);
  perf_print(out,  0,"TOTAL",                   0,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out,  3,"input control",           1,1);
  perf_print(out,  4,"input design",            1,1);
  perf_print(out,  5,"input design topology",   1,1);
  perf_print(out,  6,"input material",          1,1);
  perf_print(out,  7,"input fields",            1,1);
  perf_print(out,  8,"input detailed topology", 1,1);
  perf_print(out,  9,"input topology fe",       1,1);
  perf_print(out, 10,"input conditions",        1,1);
  perf_print(out, 11,"input inherit",           1,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 12,"part fields",             0,1);
  perf_print(out, 13,"assign dofs",             0,1);
  perf_print(out, 14,"part assign field",       0,1);
  perf_print(out, 15,"mask global mat",         0,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 16,"calc matrices",           0,1);
  perf_print(out, 17,"assemble matrices",       0,1);
  perf_print(out, 18,"assemble rhs vectors",    0,1);
  perf_print(out, 19,"exchange dofs",           0,1);
  perf_print(out, 20,"solve system",            0,1);
  if (counter[21]+counter[22]>0)
    fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  if (counter[21]>0)
    perf_print(out, 21,"permute fluid matrices",  2,1);
  if (counter[22]>0)
    perf_print(out, 22,"local co-ord. system",    2,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  user_perf = 0;
  for (i=30; i<40; ++i) {
    if (counter[i] > 0) {
      perf_print(out, i, "temp perf slot", 0, 1);
      user_perf = 1;
    }
  }
  if (user_perf)
    fprintf(out, "-----------------------------------------------------------------------------\n");
  perf_print(out, 40,"calc matrices fast", 0,1);
  perf_print(out, 41,"assemble_fast", 0,1);
  perf_print(out, 42,"init element", 0,1);
  perf_print(out, 43,"calele", 0,1);
  fprintf(out, "-----------------------------------------------------------------------------\n");
  perf_print(out, 44,"calset", 43,1);
  perf_print(out, 45,"elecord", 43,1);
  perf_print(out, 46,"elesize", 43,1);
  perf_print(out, 47,"calint", 43,1);
  perf_print(out, 48,"make_estif", 43,1);
  perf_print(out, 49,"caldirich", 43,1);
  perf_print(out, 50,"massrhs", 43,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 51,"functderiv", 47,1);
  perf_print(out, 52,"jacob", 47,1);
  perf_print(out, 53,"gder", 47,1);
  perf_print(out, 54,"gder2", 47,1);
  perf_print(out, 55,"veli", 47,1);
  perf_print(out, 56,"vder", 47,1);
  perf_print(out, 57,"elesize_2", 47,1);
  perf_print(out, 58,"gal_mat", 47,1);
  perf_print(out, 59,"stab_mat", 47,1);
  perf_print(out, 60,"iter_rhs", 47,1);
  perf_print(out, 61,"ext_rhs", 47,1);
  perf_print(out, 62,"time_rhs", 47,1);
  perf_print(out, 63,"ext_rhs", 47,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 64,"invupdate", 41,1);
  perf_print(out, 65,"invbindx", 41,1);
  perf_print(out, 66,"make lm", 41,1);
  perf_print(out, 67,"faddmsr", 41,1);
  perf_print(out, 68,"copy rhs", 41,1);
  perf_print(out, 69,"assemble rhs", 41,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 71,"out_results", 0,1);
  perf_print(out, 72,"restart_write_bin", 0,1);
  perf_print(out, 73,"out_gid", 0,1);
  perf_print(out, 74,"restart_write", 0,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");
  perf_print(out, 95,"amdef",                   0,1);
  perf_print(out, 94,"amredef",                 0,1);
  perf_print(out, 93,"amdel",                   0,1);
  perf_print(out, 92,"amzero",                  0,1);
  fprintf(out, "%77s\n","-----------------------------------------------------------------------------");

  /*----------------------------------------------------------------------*/
  return;
}
#endif




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
