
#include <unistd.h>
#include <sys/times.h>

#include "../drt_lib/standardtypes_cpp.H"
#include "global_init_control.H"
#include "global_inp_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"


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
  Teuchos::RCP<Epetra_Comm> gcomm = DRT::Problem::Instance()->GetNPGroup()->GlobalComm();

  double   t0,ti,tc;

  /// IO file names and kenners
  std::string inputfile_name;
  std::string outputfile_kenner;
  std::string restartfile_kenner;

  ntaini_ccadiscret(argc,argv,inputfile_name,outputfile_kenner,restartfile_kenner);

  /* input phase, input of all information */

  t0=cputime();

  ntainp_ccadiscret(inputfile_name,outputfile_kenner,restartfile_kenner);

  ti=cputime()-t0;
  if (gcomm->MyPID()==0)
  {
    printf("\nTotal CPU Time for INPUT:       %10.3E sec \n\n",ti);
  }

  /*--------------------------------------------------calculation phase */
  t0=cputime();

  ntacal();

  tc=cputime()-t0;
  if (gcomm->MyPID()==0)
  {
    printf("\nTotal CPU Time for CALCULATION: %10.3E sec \n\n",tc);
  }
}

