/*----------------------------------------------------------------------*/
/*!
\file io_drt_micro.cpp

\brief output context of micro discretization

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/wiechert
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#ifndef WIN_MUENCH
#include <pwd.h>
#endif

#include "compile_settings.h"
#include "io_drt_micro.H"
#include "io_drt.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_globalproblem.H"

#define BINIO_VERSION "0.3"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
  *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
  *----------------------------------------------------------------------*/
extern struct _FILES           allfiles;

using namespace std;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::MicroDiscretizationWriter::MicroDiscretizationWriter(RefCountPtr<DRT::Discretization> dis,
                                                         int probnum,
                                                         int ele,
                                                         int gp)
  : IO::DiscretizationWriter::DiscretizationWriter(dis, probnum)
{
  // this string is needed to create the new output file name
  string outputname;
  int number = -1;
  int number_rrun = -1;
  // this string is needed to create the file name of the restarted run
  string rrun(allfiles.outputfile_kenner);
  if (genprob.restart)
  {
    // name of output file

    unsigned len = cfname_.length()-1;
    // remove trailing number
    for (; (len > 0) && isdigit(cfname_[len]); --len)
    {}
    if ((len < cfname_.length()-1) && (len > 0) && (cfname_[len] == '-'))
    {
      string tmp2 = cfname_.substr(0, len);
      outputname = tmp2;
      number = atoi(&(cfname_[len+1]));
    }
    else
      dserror("removal of trailing number for microscale (restart) failed");

    // name of restarted run

    len = rrun.length()-1;
    string r;
    // remove trailing number
    for (; (len > 0) && isdigit(rrun[len]); --len)
    {}
    if ((len < rrun.length()-1) && (len > 0) && (rrun[len] == '-'))
    {
      r = rrun.substr(0, len);
      number_rrun = atoi(&(rrun[len+1]));
    }
    else
    {
      r = rrun;
    }
    ostringstream rr;
    rr << r << "_el" << ele << "_gp" << gp;
    if (number_rrun != -1)
      rr << "-" << number_rrun;
    rrun = rr.str();
    rr << ".control";
    name_restart_ = rr.str();
  }
  else
  {
    outputname = cfname_;
  }

  ostringstream s;

  s << outputname << "_el" << ele << "_gp" << gp;
  outputname = s.str();

  if (genprob.restart)
  {
    s << "-" << number;
  }

  cfname_ = s.str();
  s << ".control";
  cf_ = fopen(s.str().c_str(), "wb");

  static CHAR* problem_names[] = PROBLEMNAMES;
  time_t time_value;
  CHAR hostname[31];
  struct passwd *user_entry;

#ifndef WIN_MUENCH
  user_entry = getpwuid(getuid());
  gethostname(hostname, 30);
#else
  strcpy(hostname, "unknown host");
#endif
  time_value = time(NULL);

  // note that currently input file name that is output is that of the
  // corresponding macroscale!!

  fprintf(cf_, "# ccarat output control file\n"
          "# created by %s on %s at %s"
          "# using code revision " CHANGEDREVISION " from " CHANGEDDATE " \n\n"
          "version = \"" BINIO_VERSION "\"\n"
          "input_file = \"%s\"\n"
          "problem_type = \"%s\"\n"
          "ndim = %d\n"
          "\n",
#if !defined(WIN_MUENCH) && !defined(HPUX_GNU)
          user_entry->pw_name,
#else
          "unknown",
#endif
          hostname,
          ctime(&time_value),
          allfiles.inputfile_name,
          problem_names[genprob.probtyp],
          genprob.ndim);

  if (genprob.restart)
  {
    fprintf(cf_,
            "restarted_run = \"%s\"\n\n",
            rrun.c_str());
  }

  fflush(cf_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::MicroDiscretizationWriter::~MicroDiscretizationWriter()
{
  fclose(cf_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::MicroDiscretizationReader::MicroDiscretizationReader(Teuchos::RCP<DRT::Discretization> dis,
                                                         int step,
                                                         string name):
 DiscretizationReader(dis, step)
{
  parse_control_file(&microfile_, name.c_str());
  FindResultGroup(step,&microfile_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::MicroDiscretizationReader::~MicroDiscretizationReader()
{
  destroy_map(&microfile_);
}


#endif
