/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_writedat.H

\brief pre_exodus .dat-file writer 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here is everything related with writing a dat-file
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifndef PRE_EXODUS_WRITEDAT_H_
#define PRE_EXODUS_WRITEDAT_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_readbc.H"

using namespace std;
using namespace Teuchos;

/*!
\brief EXODUS: namespace of the ExodusII-to-BACI preprocessing tool

*/
namespace EXODUS
{
//! write datfile
int WriteDatFile(const string& datfile, const EXODUS::Mesh& mymesh,
    const string& headfile, const vector<EXODUS::elem_def>& eledefs, const vector<EXODUS::cond_def>& condefs);

//! datfile Intro
string WriteDatIntro(const EXODUS::Mesh& mymesh);

//! headfile content into datfile
string WriteDatHead(const string& headfile);

//! design description
string WriteDatDesign(const vector<EXODUS::cond_def>& condefs);

//! conditions
string WriteDatConditions(const vector<EXODUS::cond_def>& condefs);

//! DesignNode - Node Topology
string WriteDatDesignTopology(const vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh);

//! return set of nodes in bc_entity
const set<int> GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m);

//! Nodes into datfile
string WriteDatNodes(const EXODUS::Mesh& mymesh);

//! Elements into datfile
string WriteDatEles(const vector<EXODUS::elem_def>& eledefs, const EXODUS::Mesh& mymesh);

//! Elements from eblock
void DatEles(const EXODUS::ElementBlock& eb, const EXODUS::elem_def& acte, int& struele, stringstream& dat);

//! as e.g. also ElementBlocks can hold 'BACI-Conditions' we find everything what will be a 'BACI-Condition"
vector<EXODUS::bc_entity> FindBACIConditions(const vector<map<int,EXODUS::bc_entity> >& bc_specs);

inline string CondGeomTypeToString(const EXODUS::cond_def& def)
{
  switch (def.gtype){
  case DRT::Condition::Volume:  return "DVOL  "; break; 
  case DRT::Condition::Surface: return "DSURF "; break; 
  case DRT::Condition::Line:    return "DLINE "; break; 
  case DRT::Condition::Point:   return "DPOINT"; break;
  case DRT::Condition::NoGeom:  return ""; break;
  default: dserror("Unknown Condition GeometryType");
  }
  return "";
}

} // of namespace EXODUS


#endif /*PRE_EXODUS_WRITEDAT_H_*/
#endif /*EXODUS*/

