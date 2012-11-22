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

#include "pre_exodus_readbc.H"


/*!
\brief EXODUS: namespace of the ExodusII-to-BACI preprocessing tool

*/
namespace EXODUS
{
//! write datfile
int WriteDatFile(const string& datfile, const EXODUS::Mesh& mymesh,
    const string& headfile, const vector<EXODUS::elem_def>& eledefs, const vector<EXODUS::cond_def>& condefs,
    const map<int,map<int,std::vector<std::vector<double> > > >& elecenterlineinfo);

//! datfile Intro
void WriteDatIntro(const string& headfile,const EXODUS::Mesh& mymesh, ostream& dat);

//! headfile content into datfile
void WriteDatHead(const string& headfile, ostream& dat);

//! remove a specific input file section
void RemoveDatSection(const string secname, string& headstring);

//! design description
void WriteDatDesign(const vector<EXODUS::cond_def>& condefs, ostream& dat);

//! conditions
void WriteDatConditions(const vector<EXODUS::cond_def>& condefs,const EXODUS::Mesh& mymesh, ostream& dat);

//! special condition SurfLocsys: calculates Normal and Tangent
vector<double> CalcNormalSurfLocsys(const int ns_id,const EXODUS::Mesh& m);

//! DesignNode - Node Topology
void WriteDatDesignTopology(const vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, ostream& dat);

//! return set of nodes in bc_entity
const std::set<int> GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m);

//! Nodes into datfile
void WriteDatNodes(const EXODUS::Mesh& mymesh, ostream& dat);

//! Elements into datfile
void WriteDatEles(const vector<EXODUS::elem_def>& eledefs, const EXODUS::Mesh& mymesh, ostream& dat,
    const map<int,map<int,std::vector<std::vector<double> > > >& elecenterlineinfo);

//! Elements from eblock
void DatEles(Teuchos::RCP< const EXODUS::ElementBlock> eb, const EXODUS::elem_def& acte, int& struele, ostream& dat,
    const map<int,map<int,std::vector<std::vector<double> > > >& elecenterlineinfo,const int eb_id);

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

