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

#include "pre_exodus_reader.H"

using namespace std;

mesh::mesh(const string &exofilename)
{
  cout << "mesh" << endl;
  cout << exofilename.c_str() << endl;
  
  return;
}

mesh::~mesh()
{
  return;
}

