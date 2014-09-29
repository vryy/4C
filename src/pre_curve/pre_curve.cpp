
#include <string>
#include <fstream>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCP.hpp>

#include <cstdlib>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_parobjectregister.H"

/*======================================================================*/
/*======================================================================*/
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

  int printparobjecttypes = 0;

  CLP.throwExceptions(false);

  CLP.setOption("datfile",&datfile,"dat file to read curve from");
  CLP.setOption("curve",&curve,"curve number to plot");
  CLP.setOption("file",&file,"plot file name");
  CLP.setOption("t1",&t1,"start time");
  CLP.setOption("t2",&t2,"end time");
  CLP.setOption("dt",&dt,"time step");
  CLP.setOption("printparobjecttypes",&printparobjecttypes,"print names of parobject types (registration hack)");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = CLP.parse(argc,argv);

  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
  {
    std::exit(0);
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
  {
    std::exit(1);
  }
  if (printparobjecttypes)
  {
    // hack so that the parobject types are registered
    PrintParObjectList();
    exit(0);
  }

  if (datfile=="")
    dserror("dat file name is required");

  // create "dummy" NP group which only sets the correct communicators
  COMM_UTILS::CreateComm(0,NULL);
  // create a problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(problem->GetNPGroup()->GlobalComm().get(), false);
  // and now the actual reading (no dump to the non-existing error file!)
  DRT::INPUT::DatFileReader reader(datfile.c_str(), comm, 0);

  /*---------------------------------------------- input of time curves */
  DRT::UTILS::TimeCurveManager manager;
  manager.ReadInput(reader);

  if (myrank==0)
  {
    DRT::UTILS::TimeCurve& c = manager.Curve(curve-1);

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
