/*----------------------------------------------------------------------*/
/*! \file

\brief the basic baci routine

\level 0


*/
/*----------------------------------------------------------------------*/

#include "baci_comm_utils.H"
#include "baci_global_full_init_control.H"
#include "baci_global_full_inp_control.H"
#include "baci_io_pstream.H"
#include "baci_lib_globalproblem.H"



void ntacal();

/*----------------------------------------------------------------------*
 | main routine                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[])
{
  Teuchos::RCP<Epetra_Comm> gcomm = DRT::Problem::Instance()->GetCommunicators()->GlobalComm();

  double t0, ti, tc;

  /// IO file names and kenners
  std::string inputfile_name;
  std::string outputfile_kenner;
  std::string restartfile_kenner;

  ntaini_ccadiscret(argc, argv, inputfile_name, outputfile_kenner, restartfile_kenner);

  /* input phase, input of all information */

  t0 = DRT::Problem::Walltime();

  ntainp_ccadiscret(inputfile_name, outputfile_kenner, restartfile_kenner);

  ti = DRT::Problem::Walltime() - t0;
  if (gcomm->MyPID() == 0)
  {
    IO::cout << "\nTotal CPU Time for INPUT:       " << std::setw(10) << std::setprecision(3)
             << std::scientific << ti << " sec \n\n";
  }

  /*--------------------------------------------------calculation phase */
  t0 = DRT::Problem::Walltime();

  ntacal();

  tc = DRT::Problem::Walltime() - t0;
  if (gcomm->MyPID() == 0)
  {
    IO::cout << "\nTotal CPU Time for CALCULATION: " << std::setw(10) << std::setprecision(3)
             << std::scientific << tc << " sec \n\n";
  }
}
