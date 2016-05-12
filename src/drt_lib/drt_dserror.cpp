/*!----------------------------------------------------------------------
\file drt_dserror.cpp

\brief printing error messages

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_dserror.H"

#include <stdexcept>
#include <string.h>
#include <iostream>

#include <mpi.h>
#include <execinfo.h>
#include <unistd.h>

#include <cxxabi.h>
#include <stdio.h>


static int         latest_line = -1;
static std::string latest_file = "{dserror_func call without prototype}";

/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dsassert_func(const char* file, const int line, const int test, const char* text)
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(const std::string file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dslatest(const char* file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dserror_func(const char* text, ...)
{
  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  const std::size_t BUFLEN = 8192;
  char errbuf[BUFLEN];

  va_list ap;
  va_start(ap, text);

  snprintf(errbuf,BUFLEN,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsnprintf(&errbuf[strlen(errbuf)],BUFLEN-strlen(errbuf),text,ap);

#ifdef ENABLE_STACKTR
// print stacktrace
  int nptrs;
  void *buffer[100];
  char **strings;

  nptrs   = backtrace(buffer, 100);
  strings = backtrace_symbols(buffer, nptrs);

  snprintf(&errbuf[strlen(errbuf)], BUFLEN-strlen(errbuf), "\n\n--- stacktrace ---");

  // start stack trace where we actually got the error, not this function
  int frame = 0;
  while ((frame < nptrs)
         &&
         ((std::string(strings[frame]).find ("cpp_dserror_func") != std::string::npos)
          ||
          (std::string(strings[frame]).find ("cpp_dsassert_func") != std::string::npos)))
    ++frame;

  // find name of application for getting the code location of the error
#ifdef ENABLE_ADVANCED_STACKTR
  std::string applicationname;
  for (int frame2 = nptrs-1; frame2 > 0; --frame2)
    {
      std::string entry (strings[frame2]);
      const std::size_t start = entry.find('('),
                        end   = entry.find('+');
      std::string functionname = entry.substr (start+1,end-start-1);
      if (functionname == "main")
        {
          applicationname = entry.substr(0, start);
          break;
        }
    }
#endif


  int startframe = frame;
  for ( ; frame < nptrs; ++frame)
  {
    std::string entry (strings[frame]);

#ifdef ENABLE_ADVANCED_STACKTR
    // get code location of error in call stack by running the program addr2line
    // (on linux, only when debug symbols are present)
    std::string filename;
    const std::size_t address0 = entry.find(" [");
    const std::size_t address1 = entry.find("]", address0);

    if (not applicationname.empty() &&
        address0 != std::string::npos) {
      std::string progname = "addr2line " + entry.substr(address0+2,address1-2-address0) + " -e " + applicationname;
      FILE *stream = popen( progname.c_str(), "r");
      if (stream != 0) {
        char path[4096];
        if (fgets(path, sizeof(path)-1, stream) != 0) {
          std::string pathname (path);
          if (pathname.find("?") == std::string::npos)
            filename = "  (" +
              pathname.substr(pathname.find_last_of('/')+1,
                              pathname.find('\n')-pathname.find_last_of('/')-1) + ")";
        }
      }
      pclose(stream);
    }
    else
    if (address0 != std::string::npos)
      filename = entry.substr(address0,address1-address0+1);

    // demangle class names if possible
    const std::size_t start = entry.find('('),
                      end   = entry.find('+');
    std::string functionname = entry.substr(start+1,end-start-1);
    int status;
    char *p = abi::__cxa_demangle(functionname.c_str(), 0, 0, &status);
    if (status == 0)
      snprintf(&errbuf[strlen(errbuf)], BUFLEN-strlen(errbuf), "\n[%2d]: %s%s", frame-startframe, p, filename.c_str());
    else
      snprintf(&errbuf[strlen(errbuf)], BUFLEN-strlen(errbuf), "\n[%2d]: %s%s", frame-startframe, functionname.c_str(), filename.c_str());
    free(p);

    if (functionname == "main")
      break;

#else
    snprintf(&errbuf[strlen(errbuf)], BUFLEN-strlen(errbuf), "\n[%2d]: %s", frame-startframe, entry.c_str());
#endif
  }
  snprintf(&errbuf[strlen(errbuf)], BUFLEN-strlen(errbuf), "\n------------------\n");

  free(strings);
#endif

  va_end(ap);

  if (nprocs > 1)
    throw std::runtime_error(errbuf);
  else {
    // if only one processor is running, do not throw because that destroys the call
    // history in debuggers like gdb. Instead, call abort which preserves the location
    // (but do not do it in parallel because otherwise an error leads to a call stack
    // from every processor with OpenMPI, which has a severe impact in readability).
    char line[] = "=========================================================================\n";
    std::cout << "\n\n"
              << line
              << errbuf
              << "\n"
              << line
              << "\n" << std::endl;
    std::abort();
  }
} /* end of dserror_func */

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dserror_func(const std::string text, ...)
{
  cpp_dserror_func(text.c_str());
} /* end of dserror_func */
