/*----------------------------------------------------------------------*/
/*! \file

\brief general initialization routine of BACI

\maintainer Martin Kronbichler

\level 0
*/
/*----------------------------------------------------------------------*/

#include <string>
#include <sstream>

#include "global_init_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"

/*----------------------------------------------------------------------*/
/*
  Setup of input and output files. No actual read is performed here.
 */
/*----------------------------------------------------------------------*/
void ntaini_ccadiscret(int argc, char** argv, std::string& inputfile_name,
    std::string& outputfile_kenner, std::string& restartfile_kenner)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  int group = problem->GetNPGroup()->GroupId();
  int ngroups = problem->GetNPGroup()->NumGroups();
  NP_TYPE npType = problem->GetNPGroup()->NpType();
  int restartgroup = 0;
  int myrank = lcomm->MyPID();


  if (argc <= 1)
  {
    if (myrank == 0)
    {
      printf("You forgot to give the input and output file names!\n");
      printf("Try again!\n");
    }
    MPI_Finalize();
    exit(1);
  }
  else if (argc <= 2)
  {
    if (myrank == 0)
    {
      printf("You forgot to give the output file name!\n");
      printf("Try again!\n");
    }
    MPI_Finalize();
    exit(1);
  }


  // parse command line and separate input/output arguments
  std::vector<std::string> inout(0);
  for (int i = 1; i < argc; i++)
  {
    std::string temp = argv[i];
    if (temp.substr(0, 1) != "-") inout.push_back(temp);
  }

  // number of input/output arguments specified by the user
  int inoutargs = int(inout.size());

  std::stringstream infilename;
  std::stringstream outfilekenner;
  std::stringstream restartfilekenner;
  // set input file name in each group
  switch (npType)
  {
    case no_nested_parallelism:
      infilename << inout[0];
      outfilekenner << inout[1];
      restartgroup = 0;
      break;
    case every_group_read_dat_file:
    case copy_dat_file:
    {
      // TODO: to be removed after implementation has been done
      // if(npType == copy_dat_file)
      //  dserror("copyDatFile is not yet available for -nptype=");
      if (inoutargs > 4)
        dserror(
            "You specified too many arguments (%d). A maximum of four args is allowed", inoutargs);

      infilename << inout[0];
      // check whether outfilekenner includes a dash and in case separate the number at the end
      size_t pos = inout[1].rfind('-');
      if (pos != std::string::npos)
      {
        int number = atoi(inout[1].substr(pos + 1).c_str());
        inout[1] = inout[1].substr(0, pos);
        outfilekenner << inout[1] << "_group_" << group << "_"
                      << "-" << number;
      }
      else
      {
        outfilekenner << inout[1] << "_group_" << group << "_";
      }
      restartgroup = 0;
    }
    break;
    case separate_dat_files:
      if (inoutargs % ngroups != 0)
        dserror("Each group needs the same number of arguments for input/output.");
      inoutargs /= ngroups;
      infilename << inout[group * inoutargs];
      outfilekenner << inout[group * inoutargs + 1];
      restartgroup = group;
      break;
    default:
      dserror(
          "-nptype is not correct. Only copyDatFile, everyGroupReadDatFile and separateDatFiles "
          "are available");
      break;
  }


  // REMARK:
  // error files are opened by OpenErrorFile()
  // called in ntainp_ccadiscret()
  // inform user
  if (myrank == 0)
  {
    std::cout << "input is read from     " << infilename.str() << std::endl;
  }


  // default case is an identical restartfile_kenner and outputfile_kenner
  restartfilekenner << outfilekenner.str();
  for (int i = 2; i < inoutargs; i++)
  {
    std::string restart = inout[restartgroup * inoutargs + i];

    if (restart.substr(0, 8) == "restart=")
    {
      int r = atoi(restart.substr(8, std::string::npos).c_str());
      // tell the global problem about the restart step given in the command line
      problem->SetRestartStep(r);
    }
    else if (restart.substr(0, 12) == "restartfrom=")
    {
      restartfilekenner.str("");
      restartfilekenner << (restart.substr(12, std::string::npos).c_str());

      switch (npType)
      {
        case no_nested_parallelism:
        case separate_dat_files:
          // nothing to add to restartfilekenner
          break;
        case every_group_read_dat_file:
        case copy_dat_file:
        {
          // check whether restartfilekenner includes a dash and in case separate the number at the
          // end
          size_t pos = restartfilekenner.str().rfind('-');
          if (pos != std::string::npos)
          {
            int number = atoi(restartfilekenner.str().substr(pos + 1).c_str());
            std::string kenner = restartfilekenner.str().substr(0, pos);
            restartfilekenner.str("");
            restartfilekenner << kenner << "_group_" << group << "_"
                              << "-" << number;
          }
          else
          {
            restartfilekenner << "_group_" << group << "_";
          }
        }
        break;
        default:
          dserror(
              "-nptype is not correct. Only copyDatFile, everyGroupReadDatFile and "
              "separateDatFiles are available");
          break;
      }
    }
  }

  /// set IO file names and kenners
  inputfile_name = infilename.str();
  outputfile_kenner = outfilekenner.str();
  restartfile_kenner = restartfilekenner.str();

  return;
}
