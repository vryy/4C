/*----------------------------------------------------------------------*/
/*!
\file input_parameterlist.cc

\brief read an entire Teuchos::ParameterList from the input file

This requires a slight extension to the normal input file format.
The Teuchos::ParameterList is a hierarchical list and allows
whitespaces in both keys and values.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef TRILINOS_PACKAGE

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>

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
  \brief read all PL sections from the input file into the given list

  Caution: To keep things simple we use the global input variables
  directly and do not call the fortran style fr-functions. So this is
  closely tied to the underlying input data structure.

  \author u.kue
  \date 03/07
 */
/*----------------------------------------------------------------------*/
void input_ReadParameterList(Teuchos::ParameterList& list)
{
  for (int line = 0; line<allfiles.numrows; ++line)
  {
    // a section those name starts with "List" is part of the
    // parameter list. The remainder of the section name gives the
    // actual sublist.
    char* directory = strstr(allfiles.input_file[line],"--List");
    if (directory!=NULL)
    {
      directory += 6;
      std::string sublistname = directory;

      // parse sublist name, '/' separates sublists
      std::vector<std::string::size_type> pos;
      pos.push_back(sublistname.find('/'));
      while (pos.back()!=std::string::npos)
      {
        pos.push_back(sublistname.find('/',pos.back()+1));
      }

      // find sublist
      Teuchos::ParameterList* sublist = &list;
      for (std::vector<std::string::size_type>::size_type i=0;
           i<pos.size()-1;
           ++i)
      {
        std::string::size_type start = pos[i]+1;
        std::string::size_type end   = pos[i+1];

        if (end!=std::string::npos)
        {
          sublist = &sublist->sublist(sublistname.substr(start,end-start));
        }
        else
        {
          sublist = &sublist->sublist(trim(sublistname.substr(start)));
        }
      }

      // read section content and put it into the sublist
      for (line += 1; line<allfiles.numrows; ++line)
      {
        // done with this section if the line starts with "--"
        if (allfiles.input_file[line][0]=='-' &&
            allfiles.input_file[line][1]=='-')
        {
          line -= 1;
          break;
        }

        std::string strline = allfiles.input_file[line];

        // we expect a line: key = value
        // The first = in the line will be taken for the
        // separator. Thus we cannot have a = in a key.
        std::string::size_type delim = strline.find('=');
        if (delim==std::string::npos)
          dserror("no key=value pair in line %d: %s", line, strline.c_str());

        std::string key   = strline.substr(0,delim);
        std::string value = strline.substr(delim+1);

        // remove whitespaces. The simple way.
        // Note: this reduces all secutive white
        key   = trim(key);
        value = trim(value);

        // Now parse the value. Find integers and doubles if there are
        // any.
        const char* v = value.c_str();
        char* endptr = NULL;

        // Try converging to int first. If the end pointer points to
        // the trailing zero, we are done.
        long int iv = strtol(v, &endptr, 10);
        if (*endptr=='\0')
        {
          sublist->set(key,static_cast<int>(iv));
        }
        else
        {
          double dv = strtod(v, &endptr);
          if (*endptr=='\0')
          {
            sublist->set(key,dv);
          }
          else
          {
            if (value=="True" || value=="true" || value=="TRUE")
              sublist->set(key,true);
            else if (value=="False" || value=="false" || value=="FALSE")
              sublist->set(key,false);
            else
              sublist->set(key,value);
          }
        }
      }
    }
  }

  // debug
  //list.print(std::cout);
}

// Now here comes the crude C binding. A global variable! I know this
// is bad. This way we can call one simple function from C to read the
// whole thing.

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Teuchos::ParameterList> globalparameterlist;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern "C" void input_ReadGlobalParameterList()
{
  globalparameterlist = Teuchos::rcp(new Teuchos::ParameterList("Global List"));
  input_ReadParameterList(*globalparameterlist);
}

#endif
