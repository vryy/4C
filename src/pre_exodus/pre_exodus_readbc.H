/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_readbc.cpp

\brief pre_exodus bc-file reader 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here is everything related with reading a bc file
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifndef PRE_EXODUS_READBC_H_
#define PRE_EXODUS_READBC_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"

using namespace std;
using namespace Teuchos;

namespace EXODUS
{
// differentiate between underlying mesh entity
enum mesh_entity{bceb,bcns,bcss};

// differentiate between corresponding condition_type
enum cond_type{element,dvol,dsurf,dline,dpoint,empty,invalid};

struct bc_entity{
  int id;
  mesh_entity me;
  cond_type ct;
  vector<string> specs;
};

//! Read BC-file
vector<map<int,bc_entity> > ReadBCFile(const string bcfile);

//! Read bc_entity specifications
vector<string> ReadBCEntity(const string actcond);

//! Check condition type against valid types
inline EXODUS::cond_type CheckCondType(const string buffer);

//! Conversion
inline string CondTypeToString(const EXODUS::cond_type);

//! Print bc_entity
void PrintBCEntity(ostream& os, const EXODUS::bc_entity ent);

} // of namespace EXODUS


#endif /*PRE_EXODUS_READBC_H_*/
#endif /*EXODUS*/

