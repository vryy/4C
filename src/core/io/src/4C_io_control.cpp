/*----------------------------------------------------------------------*/
/*! \file
 * \brief output control
\level 0
*/
/*----------------------------------------------------------------------*/


#include "4C_config.hpp"
#include "4C_config_revision.hpp"

#include "4C_io_control.hpp"

#include "4C_io_legacy_table.hpp"
#include "4C_io_pstream.hpp"

#include <Epetra_MpiComm.h>
#include <pwd.h>
#include <unistd.h>

#include <array>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::OutputControl::OutputControl(const Epetra_Comm& comm, std::string problemtype,
    const Core::FE::ShapeFunctionType type_of_spatial_approx, std::string inputfile,
    const std::string& outputname, const int ndim, const int restart_step, const int filesteps,
    const bool write_binary_output)
    : problemtype_(std::move(problemtype)),
      inputfile_(std::move(inputfile)),
      ndim_(ndim),
      filename_(outputname),
      restartname_(outputname),
      filesteps_(filesteps),
      restart_step_(restart_step),
      myrank_(comm.MyPID()),
      write_binary_output_(write_binary_output)
{
  if (restart_step)
  {
    if (myrank_ == 0)
    {
      int number = 0;
      size_t pos = RestartFinder(filename_);
      if (pos != std::string::npos)
      {
        number = atoi(filename_.substr(pos + 1).c_str());
        filename_ = filename_.substr(0, pos);
      }

      for (;;)
      {
        number += 1;
        std::stringstream name;
        name << filename_ << "-" << number << ".control";
        std::ifstream file(name.str().c_str());
        if (not file)
        {
          filename_ = name.str();
          filename_ = filename_.substr(0, filename_.length() - 8);
          std::cout << "restart with new output file: " << filename_ << "\n";
          break;
        }
      }
    }

    if (comm.NumProc() > 1)
    {
      int length = static_cast<int>(filename_.length());
      std::vector<int> name(filename_.begin(), filename_.end());
      int err = comm.Broadcast(&length, 1, 0);
      if (err) FOUR_C_THROW("communication error");
      name.resize(length);
      err = comm.Broadcast(name.data(), length, 0);
      if (err) FOUR_C_THROW("communication error");
      filename_.assign(name.begin(), name.end());
    }
  }

  std::stringstream name;
  name << filename_ << ".control";
  write_header(name.str(), type_of_spatial_approx);

  insert_restart_back_reference(restart_step, outputname);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::OutputControl::OutputControl(const Epetra_Comm& comm, std::string problemtype,
    const Core::FE::ShapeFunctionType type_of_spatial_approx, std::string inputfile,
    const std::string& restartname, std::string outputname, const int ndim, const int restart_step,
    const int filesteps, const bool write_binary_output, const bool adaptname)
    : problemtype_(std::move(problemtype)),
      inputfile_(std::move(inputfile)),
      ndim_(ndim),
      filename_(std::move(outputname)),
      restartname_(restartname),
      filesteps_(filesteps),
      restart_step_(restart_step),
      myrank_(comm.MyPID()),
      write_binary_output_(write_binary_output)
{
  if (restart_step)
  {
    if (myrank_ == 0 && adaptname)
    {
      // check whether filename_ includes a dash and in case separate the number at the end
      int number = 0;
      size_t pos = RestartFinder(filename_);
      if (pos != std::string::npos)
      {
        number = atoi(filename_.substr(pos + 1).c_str());
        filename_ = filename_.substr(0, pos);
      }

      // either add or increase the number in the end or just set the new name for the control file
      for (;;)
      {
        // if no number is found and the control file name does not yet exist -> create it
        if (number == 0)
        {
          std::stringstream name;
          name << filename_ << ".control";
          std::ifstream file(name.str().c_str());
          if (not file)
          {
            std::cout << "restart with new output file: " << filename_ << '\n';
            break;
          }
        }
        // a number was found or the file does already exist -> set number correctly and add it
        number += 1;
        std::stringstream name;
        name << filename_ << "-" << number << ".control";
        std::ifstream file(name.str().c_str());
        if (not file)
        {
          filename_ = name.str();
          filename_ = filename_.substr(0, filename_.length() - 8);
          std::cout << "restart with new output file: " << filename_ << '\n';
          break;
        }
      }
    }

    if (comm.NumProc() > 1)
    {
      int length = static_cast<int>(filename_.length());
      std::vector<int> name(filename_.begin(), filename_.end());
      int err = comm.Broadcast(&length, 1, 0);
      if (err) FOUR_C_THROW("communication error");
      name.resize(length);
      err = comm.Broadcast(name.data(), length, 0);
      if (err) FOUR_C_THROW("communication error");
      filename_.assign(name.begin(), name.end());
    }
  }

  if (write_binary_output_)
  {
    std::stringstream name;
    name << filename_ << ".control";

    write_header(name.str(), type_of_spatial_approx);

    insert_restart_back_reference(restart_step, restartname);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::OutputControl::OutputControl(const OutputControl& ocontrol, const char* new_prefix)
    : problemtype_(ocontrol.problemtype_),
      inputfile_(ocontrol.inputfile_),
      ndim_(ocontrol.ndim_),
      filename_(ocontrol.filename_),
      restartname_(ocontrol.restartname_),
      filesteps_(ocontrol.filesteps_),
      restart_step_(ocontrol.restart_step_),
      myrank_(ocontrol.myrank_),
      write_binary_output_(ocontrol.write_binary_output_)
{
  // replace file names if provided
  if (new_prefix)
  {
    // modify file name
    {
      std::string filename_path;
      std::string filename_suffix;
      size_t pos = filename_.rfind('/');

      if (pos != std::string::npos) filename_path = filename_.substr(0, pos + 1);

      filename_suffix = filename_.substr(pos + 1);
      filename_ = filename_path + new_prefix + filename_suffix;
    }

    // modify restart name
    {
      std::string restartname_path;
      std::string restartname_suffix;
      size_t pos = restartname_.rfind('/');

      if (pos != std::string::npos) restartname_path = restartname_.substr(0, pos + 1);

      restartname_suffix = restartname_.substr(pos + 1);
      restartname_ = restartname_path + new_prefix + restartname_suffix;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::OutputControl::OverwriteResultFile(const Core::FE::ShapeFunctionType& spatial_approx)
{
  std::stringstream name;
  name << filename_ << ".control";

  write_header(name.str(), spatial_approx);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::OutputControl::new_result_file(
    int numb_run, const Core::FE::ShapeFunctionType& spatial_approx)
{
  if (filename_.rfind("_run_") != std::string::npos)
  {
    size_t pos = filename_.rfind("_run_");
    if (pos == std::string::npos) FOUR_C_THROW("inconsistent file name");
    filename_ = filename_.substr(0, pos);
  }

  std::stringstream name;
  name << filename_ << "_run_" << numb_run;
  filename_ = name.str();
  name << ".control";


  write_header(name.str(), spatial_approx);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::OutputControl::new_result_file(const std::string& name_appendix, int numb_run,
    const Core::FE::ShapeFunctionType& spatial_approx)
{
  std::stringstream name;
  name << name_appendix;
  name << "_run_" << numb_run;

  new_result_file(name.str(), spatial_approx);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::OutputControl::new_result_file(
    std::string name, const Core::FE::ShapeFunctionType& spatial_approx)
{
  filename_ = name;
  name += ".control";

  write_header(name, spatial_approx);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::IO::OutputControl::write_header(
    const std::string& control_file_name, const Core::FE::ShapeFunctionType& spatial_approx)
{
  if (myrank_ == 0)
  {
    if (controlfile_.is_open()) controlfile_.close();

    controlfile_.open(control_file_name.c_str(), std::ios_base::out);
    if (not controlfile_)
      FOUR_C_THROW("Could not open control file '%s' for writing", control_file_name.c_str());

    time_t time_value;
    time_value = time(nullptr);

    std::array<char, 256> hostname;
    struct passwd* user_entry;
    user_entry = getpwuid(getuid());
    gethostname(hostname.data(), 256);

    controlfile_ << "# 4C output control file\n"
                 << "# created by " << user_entry->pw_name << " on " << hostname.data() << " at "
                 << ctime(&time_value) << "# using code version (git SHA1) "
                 << VersionControl::git_hash << " \n\n"
                 << "input_file = \"" << inputfile_ << "\"\n"
                 << "problem_type = \"" << problemtype_ << "\"\n"
                 << "spatial_approximation = \""
                 << Core::FE::ShapeFunctionTypeToString(spatial_approx) << "\"\n"
                 << "ndim = " << ndim_ << "\n"
                 << "\n";

    controlfile_ << std::flush;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::OutputControl::insert_restart_back_reference(
    int restart, const std::string& outputname)
{
  if (myrank_ != 0) return;

  // insert back reference
  if (restart)
  {
    size_t pos = outputname.rfind('/');
    controlfile_ << "restarted_run = \""
                 << ((pos != std::string::npos) ? outputname.substr(pos + 1) : outputname)
                 << "\"\n\n";

    controlfile_ << std::flush;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string Core::IO::OutputControl::FileNameOnlyPrefix() const
{
  std::string filenameonlyprefix = filename_;

  size_t pos = filename_.rfind('/');
  if (pos != std::string::npos)
  {
    filenameonlyprefix = filename_.substr(pos + 1);
  }

  return filenameonlyprefix;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string Core::IO::OutputControl::DirectoryName() const
{
  std::filesystem::path path(filename_);
  return path.parent_path();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::InputControl::InputControl(const std::string& filename, const bool serial)
    : filename_(filename)
{
  std::stringstream name;
  name << filename << ".control";

  if (!serial)
    parse_control_file(&table_, name.str().c_str(), MPI_COMM_WORLD);
  else
    parse_control_file_serial(&table_, name.str().c_str());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::IO::InputControl::InputControl(const std::string& filename, const Epetra_Comm& comm)
    : filename_(filename)
{
  std::stringstream name;
  name << filename << ".control";

  // works for parallel, as well as serial applications because we only have an Epetra_MpiComm
  const auto* epetrampicomm = dynamic_cast<const Epetra_MpiComm*>(&comm);
  if (!epetrampicomm)
    FOUR_C_THROW("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  else
  {
    const MPI_Comm lcomm = epetrampicomm->GetMpiComm();
    parse_control_file(&table_, name.str().c_str(), lcomm);
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::InputControl::~InputControl() { destroy_map(&table_); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
size_t Core::IO::RestartFinder(const std::string& filename)
{
  size_t pos;
  for (pos = filename.size(); pos > 0; --pos)
  {
    if (filename[pos - 1] == '-') return pos - 1;

    if (not std::isdigit(filename[pos - 1]) or filename[pos - 1] == '/') return std::string::npos;
  }
  return std::string::npos;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Core::IO::GetLastPossibleRestartStep(Core::IO::InputControl& inputcontrol)
{
  /* Go to the first symbol under the name "field" and get the
   * corresponding step. Note that it will find the last "field"
   * group starting from the end of the file and looking backwards. */

  SYMBOL* symbol = map_find_symbol(inputcontrol.ControlFile(), "field");
  if (symbol != nullptr && symbol_is_map(symbol))
  {
    MAP* map;
    symbol_get_map(symbol, &map);
    return map_read_int(map, "step");
  }

  FOUR_C_THROW(
      "No restart entry in symbol table. "
      "Control file corrupt?\n\nLooking for control file at: %s",
      inputcontrol.FileName().c_str());

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
