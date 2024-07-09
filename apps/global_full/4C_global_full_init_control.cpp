/*----------------------------------------------------------------------*/
/*! \file

\brief general initialization routine of 4C


\level 0
*/
/*----------------------------------------------------------------------*/

#include "4C_global_full_init_control.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"

#include <sstream>
#include <string>

/*----------------------------------------------------------------------*/
/*
  Setup of input and output files. No actual read is performed here.
 */
/*----------------------------------------------------------------------*/
void ntaini_ccadiscret(int argc, char** argv, std::string& inputfile_name,
    std::string& outputfile_kenner, std::string& restartfile_kenner)
{
  using namespace FourC;

  Global::Problem* problem = Global::Problem::instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->get_communicators()->local_comm();
  int group = problem->get_communicators()->group_id();
  int ngroups = problem->get_communicators()->num_groups();
  Core::Communication::NestedParallelismType npType = problem->get_communicators()->np_type();
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
  auto inoutargs = int(inout.size());

  std::stringstream infilename;
  std::stringstream outfilekenner;
  std::stringstream restartfilekenner;
  // set input file name in each group
  switch (npType)
  {
    case Core::Communication::NestedParallelismType::no_nested_parallelism:
      infilename << inout[0];
      outfilekenner << inout[1];
      restartgroup = 0;
      break;
    case Core::Communication::NestedParallelismType::every_group_read_dat_file:
    {
      if (inoutargs > 4)
        FOUR_C_THROW(
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
    case Core::Communication::NestedParallelismType::separate_dat_files:
      if (inoutargs % ngroups != 0)
        FOUR_C_THROW("Each group needs the same number of arguments for input/output.");
      inoutargs /= ngroups;
      infilename << inout[group * inoutargs];
      outfilekenner << inout[group * inoutargs + 1];
      restartgroup = group;
      break;
    default:
      FOUR_C_THROW(
          "-nptype is not correct. Only everyGroupReadDatFile and separateDatFiles "
          "are available");
      break;
  }


  if (myrank == 0)
  {
    std::cout << "input is read from     " << infilename.str() << std::endl;
  }

  // bool parameter defining if input argument is given
  bool restartIsGiven = false;
  bool restartfromIsGiven = false;

  // default case is an identical restartfile_kenner and outputfile_kenner
  restartfilekenner << outfilekenner.str();
  for (int i = 2; i < inoutargs; i++)
  {
    std::string restart = inout[restartgroup * inoutargs + i];

    if (restart.substr(0, 8) == "restart=")
    {
      const std::string option = restart.substr(8, std::string::npos);
      int r;
      if (option.compare("last_possible") == 0)
      {
        r = -1;  // here we use a negative value to trigger the search in the control file in the
                 // later step. It does not mean a restart from a negative number is allowed from
                 // the user point of view.
      }
      else
      {
        r = atoi(option.c_str());
        if (r < 0) FOUR_C_THROW("Restart number must be a positive value");
      }
      // tell the global problem about the restart step given in the command line
      problem->set_restart_step(r);
      restartIsGiven = true;
    }
    else if (restart.substr(0, 12) == "restartfrom=")
    {
      restartfilekenner.str("");
      restartfilekenner << (restart.substr(12, std::string::npos).c_str());

      switch (npType)
      {
        case Core::Communication::NestedParallelismType::no_nested_parallelism:
        case Core::Communication::NestedParallelismType::separate_dat_files:
          // nothing to add to restartfilekenner
          break;
        case Core::Communication::NestedParallelismType::every_group_read_dat_file:
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
          FOUR_C_THROW(
              "-nptype is not correct. Only everyGroupReadDatFile and "
              "separateDatFiles are available");
          break;
      }

      restartfromIsGiven = true;
    }
  }

  // throw error in case restartfrom is given but no restart step is specified
  if (restartfromIsGiven && !restartIsGiven)
  {
    FOUR_C_THROW("You need to specify a restart step when using restartfrom.");
  }

  /// set IO file names and kenners
  inputfile_name = infilename.str();
  outputfile_kenner = outfilekenner.str();
  restartfile_kenner = restartfilekenner.str();
}
