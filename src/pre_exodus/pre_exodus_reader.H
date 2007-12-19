/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_reader.H

\brief preprocessor reader for exodusII format 

<pre>
Maintainer: Moritz & Georg
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with the exodus format and the accessible data
is handed to a c++ object mesh.
*/
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_READER_H_
#define PRE_EXODUS_READER_H_

#include <string>
#include <iostream>

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;

class mesh
{
public:
  //! constructor
  mesh(const string &exofilename);
  
  //! destructor
  virtual ~mesh();
};



#endif /*PRE_EXODUS_READER_H_*/
