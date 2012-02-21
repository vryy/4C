
#include <unistd.h>
#include <sys/times.h>

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_init_control.H"
#include "../drt_lib/global_inp_control2.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 *----------------------------------------------------------------------*/
struct _GENPROB       genprob;


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR     par;


static double cputime()
{
  double ret;
  double clk_tck;
  struct tms buf;

  times(&buf);
  clk_tck = (double)sysconf(_SC_CLK_TCK);
  ret = (buf.tms_utime + buf.tms_stime)/clk_tck;
  return ret;
}

void ntacal();

/*----------------------------------------------------------------------*
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(
    int                 argc,
    char               *argv[]
    )
{
  double   t0,ti,tc;

  ntaini_ccadiscret(argc,argv);

  /* input phase, input of all information */

  t0=cputime();

  ntainp_ccadiscret();

  ti=cputime()-t0;
  if (par.myrank==0)
  {
    printf("\nTotal CPU Time for INPUT:       %10.3E sec \n\n",ti);
  }

  /*--------------------------------------------------calculation phase */
  t0=cputime();

  ntacal();

  tc=cputime()-t0;
  if (par.myrank==0)
  {
    printf("\nTotal CPU Time for CALCULATION: %10.3E sec \n\n",tc);
  }
}

