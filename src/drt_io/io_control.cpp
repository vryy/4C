
#ifdef CCADISCRET

#include "io_control.H"

#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#ifndef WIN_MUENCH
#include <pwd.h>
#endif

#include <vector>
#include <iostream>
#include <sstream>

#include "../drt_lib/drt_dserror.H"

extern "C" {

#include "compile_settings.h"

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::OutputControl::OutputControl(const Epetra_Comm& comm,
                                 std::string problemtype,
                                 std::string inputfile,
                                 std::string outputname,
                                 int ndim,
                                 int restart,
                                 int filesteps)
  : problemtype_(problemtype),
    inputfile_(inputfile),
    ndim_(ndim),
    filename_(outputname),
    restartname_(outputname),
    filesteps_(filesteps)
{
  if (restart)
  {
    if (comm.MyPID()==0)
    {
      int number = 0;
      unsigned pos = filename_.rfind('-');
      if (pos!=string::npos)
      {
        number = atoi(filename_.substr(pos+1).c_str());
        filename_ = filename_.substr(0,pos);
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
          filename_ = filename_.substr(0,filename_.length()-8);
          std::cout << "restart with new output file: "
                    << filename_
                    << "\n";
          break;
        }
      }
    }

    if (comm.NumProc()>1)
    {
      int length = filename_.length();
      std::vector<int> name(filename_.begin(),filename_.end());
      int err = comm.Broadcast(&length, 1, 0);
      if (err)
        dserror("communication error");
      name.resize(length);
      err = comm.Broadcast(&name[0], length, 0);
      if (err)
        dserror("communication error");
      filename_.assign(name.begin(),name.end());
    }
  }

  if (comm.MyPID()==0)
  {
    std::stringstream name;
    name << filename_ << ".control";
    controlfile_.open(name.str().c_str(),std::ios_base::out);
    if (not controlfile_)
      dserror("could not open control file '%s' for writing", name.str().c_str());

    time_t time_value;
    time_value = time(NULL);

    char hostname[31];
    struct passwd *user_entry;
#ifndef WIN_MUENCH
    user_entry = getpwuid(getuid());
    gethostname(hostname, 30);
#else
    strcpy(hostname, "unknown host");
#endif

    controlfile_ << "# baci output control file\n"
                 << "# created by "
#if !defined(WIN_MUENCH) && !defined(HPUX_GNU)
                 << user_entry->pw_name
#else
                 << "unknown"
#endif
                 << " on " << hostname << " at " << ctime(&time_value)
                 << "# using code revision " CHANGEDREVISION " from " CHANGEDDATE " \n\n"
                 << "input_file = \"" << inputfile << "\"\n"
                 << "problem_type = \"" << problemtype << "\"\n"
                 << "ndim = " << ndim << "\n"
                 << "\n";

    // insert back reference
    if (restart)
    {
      unsigned pos = outputname.rfind('/');
      controlfile_ << "restarted_run = \""
                   << ((pos!=string::npos) ? outputname.substr(pos+1) : outputname)
                   << "\"\n\n";
    }

    controlfile_ << std::flush;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::OutputControl::OutputControl(const Epetra_Comm& comm,
                                 std::string problemtype,
                                 std::string inputfile,
                                 std::string restartname,
                                 std::string outputname,
                                 int ndim,
                                 int restart,
                                 int filesteps)
  : problemtype_(problemtype),
    inputfile_(inputfile),
    ndim_(ndim),
    filename_(outputname),
    restartname_(restartname),
    filesteps_(filesteps)
{
  if (restart)
  {
    if (comm.NumProc()>1)
    {
      int length = filename_.length();
      std::vector<int> name(filename_.begin(),filename_.end());
      int err = comm.Broadcast(&length, 1, 0);
      if (err)
        dserror("communication error");
      name.resize(length);
      err = comm.Broadcast(&name[0], length, 0);
      if (err)
        dserror("communication error");
      filename_.assign(name.begin(),name.end());
    }
  }

  if (comm.MyPID()==0)
  {
    std::stringstream name;
    name << filename_ << ".control";
    controlfile_.open(name.str().c_str(),std::ios_base::out);
    if (not controlfile_)
      dserror("could not open control file '%s' for writing", name.str().c_str());

    time_t time_value;
    time_value = time(NULL);

    char hostname[31];
    struct passwd *user_entry;
#ifndef WIN_MUENCH
    user_entry = getpwuid(getuid());
    gethostname(hostname, 30);
#else
    strcpy(hostname, "unknown host");
#endif

    controlfile_ << "# baci output control file\n"
                 << "# created by "
#if !defined(WIN_MUENCH) && !defined(HPUX_GNU)
                 << user_entry->pw_name
#else
                 << "unknown"
#endif
                 << " on " << hostname << " at " << ctime(&time_value)
                 << "# using code revision " CHANGEDREVISION " from " CHANGEDDATE " \n\n"
                 << "input_file = \"" << inputfile_ << "\"\n"
                 << "problem_type = \"" << problemtype_ << "\"\n"
                 << "ndim = " << ndim_ << "\n"
                 << "\n";

    // insert back reference
    if (restart)
    {
      controlfile_ << "restarted_run = \""
                   << restartname
                   << "\"\n\n";
    }

    controlfile_ << std::flush;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::OutputControl::NewResultFile(int numb_run)
{
  if (numb_run > 1)
  {
    unsigned pos = filename_.rfind("_");
    if (pos==string::npos)
      dserror("inconsistent file name");
    filename_ = filename_.substr(0, pos);
  }

  std::stringstream name;
  name << filename_ << "_"<< numb_run;
  filename_ = name.str();
  name << ".control";
  controlfile_.close();
  controlfile_.open(name.str().c_str(),std::ios_base::out);
  if (not controlfile_)
    dserror("could not open control file '%s' for writing", name.str().c_str());

  time_t time_value;
  time_value = time(NULL);

  char hostname[31];
  struct passwd *user_entry;
#ifndef WIN_MUENCH
  user_entry = getpwuid(getuid());
  gethostname(hostname, 30);
#else
  strcpy(hostname, "unknown host");
#endif

  controlfile_ << "# baci output control file\n"
               << "# created by "
#if !defined(WIN_MUENCH) && !defined(HPUX_GNU)
               << user_entry->pw_name
#else
               << "unknown"
#endif
               << " on " << hostname << " at " << ctime(&time_value)
               << "# using code revision " CHANGEDREVISION " from " CHANGEDDATE " \n\n"
               << "input_file = \"" << inputfile_ << "\"\n"
               << "problem_type = \"" << problemtype_ << "\"\n"
               << "ndim = " << ndim_ << "\n"
               << "\n";

  controlfile_ << std::flush;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::InputControl::InputControl(std::string filename, const bool serial)
  : filename_(filename)
{
  std::stringstream name;
  name << filename << ".control";

  if (!serial)
    parse_control_file(&table_, name.str().c_str());
  else
    parse_control_file_serial(&table_, name.str().c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::InputControl::~InputControl()
{
  destroy_map(&table_);
}

#endif
