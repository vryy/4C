/*----------------------------------------------------------------------------*/
/*!
\file every_iteration_writer.cpp

\brief Write output for each Newton step during one load step in
       an extra output file.

\maintainer Michael Hiermeier

\date Aug 10, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "every_iteration_writer.H"
#include "io.H"
#include "io_control.H"
#include "io_pstream.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <boost/filesystem.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
IO::EveryIterationWriter::EveryIterationWriter()
    : isinit_( false ),
      issetup_(false),
      isnewton_initialized_(false),
      myrank_(-1),
      run_number_(-1),
      write_only_this_step_(-1),
      write_owner_each_newton_iteration_(false),
      base_filename_(),
      parent_writer_( NULL ),
      interface_( NULL ),
      every_iter_writer_( Teuchos::null )
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::Init(
    const IO::DiscretizationWriter* parent_writer,
    IO::EveryIterationWriterInterface* interface,
    const Teuchos::ParameterList& params )
{
  issetup_ = false;

  parent_writer_ = parent_writer;
  interface_ = interface;

  myrank_ = parent_writer_->Output()->MyRank();

  run_number_ = params.get<int>("RUN_NUMBER");
  write_only_this_step_ = params.get<int>("STEP_NP_NUMBER");
  write_owner_each_newton_iteration_ = DRT::INPUT::IntegralValue<bool>(
      params, "WRITE_OWNER_EACH_NEWTON_ITER" );

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::Setup()
{
  ThrowIfNotInitialized( __LINE__ );

  const std::string dir_name( ParentWriter().Output()->FileNameOnlyPrefix() +
      "_every_iter");

  std::string file_dir_path = ExtractPath( ParentWriter().Output()->FileName() );
  file_dir_path += dir_name;
  CreateDirectory( file_dir_path );

  std::string restart_dir_path = ExtractPath( ParentWriter().Output()->RestartName() );
  restart_dir_path += dir_name;
  CreateDirectory( restart_dir_path );

  // modify the prefix and copy possible restart number
  std::string prefix;
  prefix = dir_name + "/" + CreateRunDirectory( file_dir_path );

  Teuchos::RCP<IO::OutputControl> control_iteration = Teuchos::null;
  control_iteration = Teuchos::rcp( new IO::OutputControl(
      *ParentWriter().Output(), prefix.c_str() ) );

  // adjust steps per file
  AdjustStepsPerFile( *control_iteration );

  // create new output writer object
  every_iter_writer_ = Teuchos::rcp( new IO::DiscretizationWriter(
      ParentWriter(), control_iteration, IO::CopyType::shape ) );

  // save base file name
  base_filename_ = every_iter_writer_->Output()->FileName();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::EveryIterationWriter::CreateRunDirectory(
    const std::string& file_dir_path ) const
{
  if ( run_number_ < 0 )
    return "";

  std::ostringstream run_dir;
  run_dir << "run_" << run_number_;

  std::string full_dir_path = file_dir_path + "/" + run_dir.str();
  CreateDirectory( full_dir_path );

  return run_dir.str() + "/";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::CreateDirectory( const std::string& dir_path ) const
{
  if ( myrank_ != 0)
    return;

  boost::filesystem::path dir( dir_path );
  if ( not boost::filesystem::is_directory( dir ) )
    if ( not boost::filesystem::create_directory( dir ) )
      dserror( "The directory \"%s\" could not be created!", dir_path.c_str() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string IO::EveryIterationWriter::ExtractPath(
    const std::string& full_filename ) const
{
  std::string filename_path;
  size_t pos = full_filename.rfind('/');

  if (pos!=std::string::npos)
    filename_path = full_filename.substr(0,pos+1);

  return filename_path;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool IO::EveryIterationWriter::WriteThisStep() const
{
  return ( write_only_this_step_ < 0 or
      write_only_this_step_ == interface_->GetStepNp() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::AdjustStepsPerFile( IO::OutputControl& control ) const
{
  int new_file_steps = control.FileSteps() * MAX_NUMBER_LINE_SEARCH_ITERATIONS_;
  if ( new_file_steps > std::numeric_limits<int>::max() )
    new_file_steps = std::numeric_limits<int>::max();

  control.SetFileSteps( new_file_steps );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::PrintPath2Screen( const std::string& path ) const
{
  IO::cout << "every iteration output path: \"" << path << "\"" << IO::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::InitNewtonIteration()
{
  ThrowIfNotSetup( __LINE__ );

  if ( not WriteThisStep() )
    return;

  const int curr_step = interface_->GetStepNp();

  // create new result and mesh files for each step
  std::ostringstream result_name;
  result_name << base_filename_ << "_step_" << curr_step;

  PrintPath2Screen( result_name.str() );

  every_iter_writer_->NewResultFile( result_name.str() );
  every_iter_writer_->WriteMesh( 0, 0.0 );

  every_iter_writer_->NewStep( 0, 0.0 );

  interface_->PrepareOutput();
  interface_->OutputState( *every_iter_writer_, true );

  isnewton_initialized_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void IO::EveryIterationWriter::AddNewtonIteration( const int newton_iteration )
{
  ThrowIfNotSetup( __LINE__ );

  if ( not WriteThisStep() )
    return;

  if ( not isnewton_initialized_ )
    dserror( "Call InitNewtonIteration() first!" );

  const int counter = MAX_NUMBER_LINE_SEARCH_ITERATIONS_ * newton_iteration;
  every_iter_writer_->WriteMesh( counter, counter );
  every_iter_writer_->NewStep( counter, counter );

  interface_->PrepareOutput();
  interface_->OutputState( *every_iter_writer_, write_owner_each_newton_iteration_ );
}
