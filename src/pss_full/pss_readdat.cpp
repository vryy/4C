
#ifdef CCADISCRET

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

// include all C stuff from ccarat
extern "C" {

#include "../headers/standardtypes.h"

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>

*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

}


using namespace std;

/*----------------------------------------------------------------------*/
/*!
  \brief Section start positions of excluded section in input file

  Another global variable!

  \author u.kue
  \date 03/07
 */
/*----------------------------------------------------------------------*/
map<string,ifstream::pos_type> ExcludedSectionPositions;


/*----------------------------------------------------------------------*/
/*!
  \brief remove all leading and trailing whitespaces from a string.

  Note: consecutive whitespaces inside the string will be reduced to a
  single space.

  \author u.kue
  \date 03/07
*/
/*----------------------------------------------------------------------*/
static std::string trim(const std::string& line)
{
  std::istringstream t;
  std::string s;
  std::string newline;
  t.str(line);
  while (t >> s)
  {
    newline.append(s);
    newline.append(" ");
  }
  if (newline.size()>0)
    newline.resize(newline.size()-1);
  return newline;
}


/*----------------------------------------------------------------------*/
/*!
  \brief read dat file and allocate the global variables on proc 0

  Returns the number of allocated bytes. The internal buffer now has a
  variable number of lines. So some more work is needed (and some
  memory saved.)

  Some special sections are excluded. These have to be read directly
  from the file later on.

  No include file support here!

  \author u.kue
  \date 03/07
*/
/*----------------------------------------------------------------------*/
extern "C" unsigned ReadDat(char* filename)
{
  vector<string> exclude;

  exclude.push_back("--NODE COORDS");
  exclude.push_back("--STRUCTURE ELEMENTS");
  exclude.push_back("--FLUID ELEMENTS");
  exclude.push_back("--ALE ELEMENTS");
  //exclude.push_back("--DNODE-NODE TOPOLOGY");
  //exclude.push_back("--DLINE-NODE TOPOLOGY");
  //exclude.push_back("--DSURF-NODE TOPOLOGY");
  //exclude.push_back("--DVOL-NODE TOPOLOGY");

  ifstream file(filename);
  list<string> content;
  bool ignoreline = false;
  string line;
  string::size_type arraysize = 0;

  // loop all input lines
  while (getline(file, line))
  {
    // remove comments
    string::size_type loc = line.find("//");
    if (loc != string::npos)
    {
      line = line.substr(0,loc);
    }

    // remove trailing and leading whitespaces
    // compact internal whitespaces
    line = trim(line);

    // exclude all special sections
    // this includes the section header and all lines in that section
    if (ignoreline)
    {
      if (line.find("--")==0)
      {
        ignoreline = false;
      }
    }

    // Two sections to be ignored can follow each other. We need
    // independent tests.
    if (!ignoreline)
    {
      for (vector<int>::size_type i=0; i<exclude.size(); ++i)
      {
        if (line.find(exclude[i]) != string::npos)
        {
          ExcludedSectionPositions[exclude[i]] = file.tellg();
          ignoreline = true;
          break;
        }
      }
    }

    // remember line
    if (!ignoreline && line.length() > 0)
    {
      content.push_back(line);
      arraysize += line.length()+1;
    }
  }

  // setup global variables
  allfiles.numrows = static_cast<int>(content.size());

  /* allocate space for copy of file */
  allfiles.input_file_hook = static_cast<char*>(CCACALLOC(arraysize,sizeof(char)));

  // CAUTION: We allocate one more row pointer that necessary. This
  // pointer will be used to point to temporary lines when the
  // excluded section are read (on proc 0). Of course that's just
  // another EVIL HACK. Don't tell anybody.
  allfiles.input_file = static_cast<char**>(CCACALLOC(content.size()+1,sizeof(char*)));

  // fill file buffer
  char* ptr = allfiles.input_file_hook;
  int counter = 0;
  for (list<string>::iterator i=content.begin(); i!=content.end(); ++i)
  {
    strcpy(ptr,i->c_str());
    allfiles.input_file[counter] = ptr;
    ptr += i->length()+1;
    counter += 1;
  }

  if (ptr - allfiles.input_file_hook != static_cast<int>(arraysize))
    dserror("internal error in file read");

  return arraysize;
}

#endif
