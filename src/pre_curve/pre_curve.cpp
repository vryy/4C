
#include <string>
#include <fstream>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCP.hpp>

#include <cstdlib>

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_inputreader.H"

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


int main(int argc, char** argv)
{
  int myrank = 0;

#ifdef PARALLEL
  int nproc  = 1;
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif

  Teuchos::CommandLineProcessor CLP;
  CLP.setDocString("pre curve plot\n");

  std::string datfile;
  std::string file="xxx.plot";

  int curve = 1;

  double t1 = 0.;
  double t2 = 1.;
  double dt = 0.1;

  CLP.throwExceptions(false);

  CLP.setOption("datfile",&datfile,"dat file to read curve from");
  CLP.setOption("curve",&curve,"curve number to plot");
  CLP.setOption("file",&file,"plot file name");
  CLP.setOption("t1",&t1,"start time");
  CLP.setOption("t2",&t2,"end time");
  CLP.setOption("dt",&dt,"time step");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = CLP.parse(argc,argv);

  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
  {
    std::exit(0);
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
  {
    std::exit(1);
  }

  if (datfile=="")
    dserror("dat file name is required");

#ifdef PARALLEL
  Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(com);
#else
  Epetra_SerialComm* com = new Epetra_SerialComm();
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(com);
#endif

  // and now the actual reading
  DRT::INPUT::DatFileReader reader(datfile.c_str(), comm);
  reader.Activate();

  /*---------------------------------------------- input of time curves */
  DRT::UTILS::TimeCurveManager::Instance().ReadInput();

  if (myrank==0)
  {
    DRT::UTILS::TimeCurve& c = DRT::UTILS::TimeCurveManager::Instance().Curve(curve-1);

    std::ofstream f(file.c_str());
    f << "# count   time   f(t)   f'(t)   f''(t)\n";
    int i=0;
    for (double t=t1; t<=t2; t+=dt)
    {
      std::vector<double> values = c.FctDer(t,2);
      f << i << "   "
        << t << "   "
        << values[0] << "   "
        << values[1] << "   "
        << values[2] << "   "
        << "\n";
      i += 1;
    }
  }

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;
}
