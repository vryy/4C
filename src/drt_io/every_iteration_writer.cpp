/*----------------------------------------------------------------------------*/
/*! \file
\brief Write output for each Newton step during one load step in an extra output file.


\level 3
*/
/*----------------------------------------------------------------------------*/

#include "every_iteration_writer.H"
#include "io.H"
#include "io_control.H"
#include "io_pstream.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <boost/filesystem.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
IO::EveryIterationWriter::EveryIterationWriter()
    : isinit_(false),
      issetup_(false),
      isnewton_initialized_(false),
      myrank_(-1),
      run_number_(-1),
      write_only_this_step_(-1),
      write_owner_each_newton_iteration_(false),
      base_filename_(),
      parent_writer_(NULL),
      interface_(NULL),
      every_iter_writer_(Teuchos::null)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::Init(const IO::DiscretizationWriter* parent_writer,
    IO::EveryIterationWriterInterface* interface, const Teuchos::ParameterList& params)
{
  issetup_ = false;

  parent_writer_ = parent_writer;
  interface_ = interface;


  myrank_ = parent_writer_->Output()->MyRank();

  run_number_ = params.get<int>("RUN_NUMBER");
  write_only_this_step_ = params.get<int>("STEP_NP_NUMBER");
  write_owner_each_newton_iteration_ =
      DRT::INPUT::IntegralValue<bool>(params, "WRITE_OWNER_EACH_NEWTON_ITER");

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::Setup()
{
  ThrowIfNotInitialized(__LINE__);

  /* Remove the restart counter from the folder name. Note that the restart
   * counter stays a part of the final file name of the corresponding step. */
  const std::string filename_without_restart =
      RemoveRestartCountFromFileName(ParentWriter().Output()->FileNameOnlyPrefix());

  const std::string dir_name(filename_without_restart + "_every_iter");

  std::string file_dir_path = ExtractPath(ParentWriter().Output()->FileName());
  file_dir_path += dir_name;
  CreateDirectory(file_dir_path);

  // modify the prefix and copy possible restart number
  std::string prefix;
  prefix = dir_name + "/" + CreateRunDirectory(file_dir_path);

  Teuchos::RCP<IO::OutputControl> control_iteration = Teuchos::null;
  control_iteration = Teuchos::rcp(new IO::OutputControl(*ParentWriter().Output(), prefix.c_str()));

  // adjust steps per file
  AdjustStepsPerFile(*control_iteration);

  // create new output writer object
  every_iter_writer_ = Teuchos::rcp(
      new IO::DiscretizationWriter(ParentWriter(), control_iteration, IO::CopyType::shape));

  // save base file name
  base_filename_ = every_iter_writer_->Output()->FileName();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::EveryIterationWriter::CreateRunDirectory(const std::string& file_dir_path) const
{
  if (run_number_ < 0) return "";

  std::ostringstream run_dir;
  run_dir << "run_" << run_number_;

  std::string full_dir_path = file_dir_path + "/" + run_dir.str();
  CreateDirectory(full_dir_path);

  return run_dir.str() + "/";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::CreateDirectory(const std::string& dir_path) const
{
  IO::CreateDirectory(dir_path, myrank_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::CreateDirectory(const std::string& dir_path, const int myrank)
{
  if (myrank != 0) return;

  boost::filesystem::path dir(dir_path);
  if (not boost::filesystem::is_directory(dir))
    if (not boost::filesystem::create_directory(dir))
      dserror("The directory \"%s\" could not be created!", dir_path.c_str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::ExtractPath(const std::string& full_filename)
{
  std::string filename_path;
  size_t pos = full_filename.rfind('/');

  if (pos != std::string::npos) filename_path = full_filename.substr(0, pos + 1);

  return filename_path;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::ExtractFileName(const std::string& full_filename)
{
  std::string filenameonly = full_filename;

  size_t pos = full_filename.rfind('/');
  if (pos != std::string::npos)
  {
    filenameonly = full_filename.substr(pos + 1);
  }

  return filenameonly;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::RemoveRestartCountFromFileName(const std::string& filename)
{
  const int restart_count = DRT::Problem::Instance()->Restart();
  if (restart_count == 0) return filename;

  size_t pos = filename.rfind('-');
  if (pos == std::string::npos) return filename;

  // potential restart counter
  const char rc = filename.at(pos + 1);

  // check if it is really an integer
  const int irc = rc - '0';
  if (irc >= 0 and irc <= 9)
    return filename.substr(0, pos);
  else
    return filename;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool IO::EveryIterationWriter::WriteThisStep() const
{
  return (write_only_this_step_ < 0 or write_only_this_step_ == interface_->GetStepNp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::AdjustStepsPerFile(IO::OutputControl& control) const
{
  int new_file_steps = control.FileSteps() * MAX_NUMBER_LINE_SEARCH_ITERATIONS_;
  if (new_file_steps > std::numeric_limits<int>::max())
    new_file_steps = std::numeric_limits<int>::max();

  control.SetFileSteps(new_file_steps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::PrintPath2Screen(const std::string& path) const
{
  IO::cout << "every iteration output path: \"" << path << "\"" << IO::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::InitNewtonIteration()
{
  ThrowIfNotSetup(__LINE__);

  if (not WriteThisStep()) return;

  const int curr_step = interface_->GetStepNp();

  // create new result and mesh files for each step
  std::ostringstream result_name;
  result_name << base_filename_ << "_step_" << curr_step;

  PrintPath2Screen(result_name.str());

  every_iter_writer_->NewResultFile(result_name.str());
  every_iter_writer_->WriteMesh(0, 0.0);

  every_iter_writer_->NewStep(0, 0.0);

  interface_->PrepareOutput();
  interface_->OutputDebugState(*every_iter_writer_, true);

  isnewton_initialized_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::AddNewtonIteration(const int newton_iteration)
{
  ThrowIfNotSetup(__LINE__);

  if (not WriteThisStep()) return;

  if (not isnewton_initialized_) dserror("Call InitNewtonIteration() first!");

  const int counter = MAX_NUMBER_LINE_SEARCH_ITERATIONS_ * newton_iteration;
  every_iter_writer_->WriteMesh(counter, counter);
  every_iter_writer_->NewStep(counter, counter);

  interface_->PrepareOutput();
  interface_->OutputDebugState(*every_iter_writer_, write_owner_each_newton_iteration_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::AddLineSearchIteration(
    const int newton_iteration, const int linesearch_iteration)
{
  ThrowIfNotSetup(__LINE__);

  if (not WriteThisStep()) return;

  if (not isnewton_initialized_) dserror("Call InitNewtonIteration() first!");

  if (linesearch_iteration >= static_cast<int>(MAX_NUMBER_LINE_SEARCH_ITERATIONS_))
    dserror(
        "The EveryIterationWriter does not support more than %d line search"
        " steps. If this number is exceeded, the counters will get mixed up.",
        MAX_NUMBER_LINE_SEARCH_ITERATIONS_ - 1);

  const int counter = MAX_NUMBER_LINE_SEARCH_ITERATIONS_ * newton_iteration + linesearch_iteration;
  every_iter_writer_->WriteMesh(counter, counter);
  every_iter_writer_->NewStep(counter, counter);

  interface_->PrepareOutput();
  interface_->OutputDebugState(*every_iter_writer_, write_owner_each_newton_iteration_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int IO::CountLinesInFile(const std::string& filepath)
{
  std::ifstream myfile(filepath);

  if (not myfile.is_open()) return -1;

  // new lines will be skipped unless we stop it from happening:
  myfile.unsetf(std::ios_base::skipws);

  // count the newlines with an algorithm specialized for counting:
  int line_count =
      std::count(std::istream_iterator<char>(myfile), std::istream_iterator<char>(), '\n');

  myfile.close();

  return line_count;
}
