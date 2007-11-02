/*----------------------------------------------------------------------*/
/*!
\file io_drt.cpp

\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <iostream>
#include <sstream>
#include <string>

#include "io_drt_micro.H"
#include "io_drt.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_globalproblem.H"

#define BINIO_VERSION "0.2"

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
  ostringstream s;
  s << cfname_ << "_el" << ele << "_gp" << gp;
  cfname_ = s.str();
  s << ".control";

  cf_ = fopen(s.str().c_str(), "wb");

  static CHAR* problem_names[] = PROBLEMNAMES;

  // note that currently input file name that is output is that of the
  // corresponding macroscale!!

  fprintf(cf_, "# ccarat output control file\n"
          "# using io version: $Id: io_singlefile.cpp 2881 2007-03-22 14:57:03Z kuettler $ \n\n"
          "version = \"" BINIO_VERSION "\"\n"
          "input_file = \"%s\"\n"
          "problem_type = \"%s\"\n"
          "ndim = %d\n"
          "\n",
          allfiles.inputfile_name,
          problem_names[genprob.probtyp],
          genprob.ndim);

  fflush(cf_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::MicroDiscretizationWriter::~MicroDiscretizationWriter()
{
  fclose(cf_);
}

#endif
