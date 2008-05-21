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
#include "pre_exodus_readbc.H"
#include "pre_exodus_reader.H"


using namespace std;
using namespace Teuchos;

void EXODUS::ReadBCFile(const string& bcfile, vector<EXODUS::elem_def>& eledefs, vector<EXODUS::cond_def>& condefs)
{
  // first we read the whole file into one stream/string
  stringstream bcstream;
  const char *bcfilechar;
  bcfilechar = bcfile.c_str();
  ifstream bcfstream(bcfilechar, ifstream::in);
  if (!bcfstream.good()){
    cout << endl << "Unable to open file: " << bcfile.c_str() << endl;
    dserror("Unable to open bc-file");
  }
  while (bcfstream.good()) bcstream << (char) bcfstream.get();
  bcfstream.close();
  // string which contains the whole file
  string allconds = bcstream.str();
  allconds.erase(allconds.end()-1);  //delete last 'whatisthis'-char
  
  // get rid of first part
  size_t found;
  found = allconds.find("BCSPECS");
  if (found==string::npos) dserror ("No specifications found in bcfile");
  allconds.erase(allconds.begin(),allconds.begin()+found);
  
  // get rid of 'validconditions' part
  found = allconds.find("VALIDCONDITIONS");
  allconds.erase(allconds.begin()+found,allconds.end());
  
  // define markers
  const string ebmarker("*eb");
  const string nsmarker("*ns");
  const string ssmarker("*ss");
  const int markerlength=3;
  const string marker("*");
  
  // necessary counters
  int E_id = 0; //the 'E num -' in the datfile
  int ndp = 0; int ndl = 0; int nds = 0; int ndv = 0;
  
  // map to avoid double assignment
  map<int,int> dp2Eid;
  map<int,int> dl2Eid;
  map<int,int> ds2Eid;
  map<int,int> dv2Eid;
  
  found = allconds.find_first_of(marker);  
  while (found != string::npos){
    int startpos=found;
    found = allconds.find(marker,found+1);  //step forward to find next match
    // get actual condition
    string actcond = allconds.substr(startpos,found-startpos);
    
    
    // find out what mesh_entity type we have
    string mesh_entity = actcond.substr(0,3);

    // get its id
    size_t found2 = actcond.find_first_of("=");
    string buffer = actcond.substr(markerlength,found2-markerlength);
    // convert string to int
    istringstream bufferstream(buffer);
    int id;
    bufferstream >> id;

    // condition type
    size_t left = actcond.find_first_of("\"",0);
    size_t right = actcond.find_first_of("\"",left+1);
    string type = actcond.substr(left+1,right-left-1);
    
    if (mesh_entity.compare(ebmarker)==0) {
      // in case of eb we differntiate between 'element' or 'condition'
      if (type.compare("ELEMENT")==0) {
        EXODUS::elem_def edef = EXODUS::ReadEdef(mesh_entity, id, actcond);
        eledefs.push_back(edef);
      }
      else if (type.compare("CONDITION")==0) {
        EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
        switch (cdef.gtype) {
        case DRT::Condition::Point: {
          if (dp2Eid.find(id) != dp2Eid.end())
            E_id = dp2Eid.find(id)->second;
          else {
            ++ndp;
            E_id = ndp;
            dp2Eid.insert(pair<int,int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Line: {
          if (dl2Eid.find(id) != dl2Eid.end())
            E_id = dl2Eid.find(id)->second;
          else {
            ++ndl;
            E_id = ndl;
            dl2Eid.insert(pair<int,int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Surface: {
          if (ds2Eid.find(id) != ds2Eid.end())
            E_id = ds2Eid.find(id)->second;
          else {
            ++nds;
            E_id = nds;
            ds2Eid.insert(pair<int,int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Volume: {
          if (dv2Eid.find(id) != dv2Eid.end())
            E_id = dv2Eid.find(id)->second;
          else {
            ++ndv;
            E_id = ndv;
            dv2Eid.insert(pair<int,int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::NoGeom:
          E_id = 0;
          break;
        default:
          dserror("geometry type unspecified");
        }
        cdef.e_id = E_id;
        condefs.push_back(cdef);
      }
      else {
        cout << "Undefined type for eb"<<id<<": "<<type<<endl;
        dserror("Undefined type!");
      }
    }
    else if (mesh_entity.compare(nsmarker)==0) {
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
      switch (cdef.gtype) {
      case DRT::Condition::Point: {
        if (dp2Eid.find(id) != dp2Eid.end())
          E_id = dp2Eid.find(id)->second;
        else {
          ++ndp;
          E_id = ndp;
          dp2Eid.insert(pair<int,int>(id, E_id));
        }
        break;
      }
      case DRT::Condition::Line: {
        if (dl2Eid.find(id) != dl2Eid.end())
          E_id = dl2Eid.find(id)->second;
        else {
          ++ndl;
          E_id = ndl;
          dl2Eid.insert(pair<int,int>(id, E_id));
        }
        break;
      }
      case DRT::Condition::Surface: {
        if (ds2Eid.find(id) != ds2Eid.end())
          E_id = ds2Eid.find(id)->second;
        else {
          ++nds;
          E_id = nds;
          ds2Eid.insert(pair<int,int>(id, E_id));
        }
        break;
      }
      case DRT::Condition::Volume: {
        if (dv2Eid.find(id) != dv2Eid.end())
          E_id = dv2Eid.find(id)->second;
        else {
          ++ndv;
          E_id = ndv;
          dv2Eid.insert(pair<int,int>(id, E_id));
        }
        break;
      }
      case DRT::Condition::NoGeom:
        E_id = 0;
        break;
      default:
        dserror("geometry type unspecified");
      }
      cdef.e_id = E_id;
      condefs.push_back(cdef);
    }
    else if (mesh_entity.compare(ssmarker)==0) {
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
      condefs.push_back(cdef);
    }
    else
      dserror("Undefined mesh_type");

  }
  
  return;
}

EXODUS::elem_def EXODUS::ReadEdef(const string& mesh_entity,const int id, const string& actcond)
{
  EXODUS::elem_def edef;
  edef.id = id;
  edef.me = EXODUS::bceb;

  // read sectionname
  size_t left = actcond.find("sectionname=\"");  // 13 chars
  size_t right = actcond.find_first_of("\"",left+13);
  edef.sec = actcond.substr(left+13,right-(left+13));

  // read description
  left = actcond.find("description=\"");  // 13 chars
  right = actcond.find_first_of("\"",left+13);
  edef.desc = actcond.substr(left+13,right-(left+13));

  // read ename
  left = actcond.find("elementname=\"");  // 13 chars
  right = actcond.find_first_of("\"",left+13);
  edef.ename = actcond.substr(left+13,right-(left+13));
  
  return edef;
}

EXODUS::cond_def EXODUS::ReadCdef(const string& mesh_entity,const int id, const string& actcond)
{
  EXODUS::cond_def cdef;
  cdef.id = id;
  if (mesh_entity.compare(1,2,"eb")==0) cdef.me = EXODUS::bceb;
  else if (mesh_entity.compare(1,2,"ns")==0) cdef.me = EXODUS::bcns;
  else if (mesh_entity.compare(1,2,"ss")==0) cdef.me = EXODUS::bcss;
  else dserror("Undefined mesh_type");
  
  // read sectionname
  size_t left = actcond.find("sectionname=\"");  // 13 chars
  size_t right = actcond.find_first_of("\"",left+13);
  string secname = actcond.substr(left+13,right-(left+13));
  cdef.sec = secname;

  // read description
  left = actcond.find("description=\"");  // 13 chars
  right = actcond.find_first_of("\"",left+13);
  string description = actcond.substr(left+13,right-(left+13));
  cdef.desc = description;

  // figure out geometry type
  cdef.gtype = DRT::Condition::NoGeom;  // default
  size_t found = secname.find("POINT");
  if (found!=string::npos){
    cdef.gtype = DRT::Condition::Point;
    return cdef;
  }
  found = secname.find("LINE");
  if (found!=string::npos){
    cdef.gtype = DRT::Condition::Line;
    return cdef;
  }
  found = secname.find("SURF");
  if (found!=string::npos){
    cdef.gtype = DRT::Condition::Surface;
    return cdef;
  }
  found = secname.find("VOL");
  if (found!=string::npos) cdef.gtype = DRT::Condition::Volume;
  
  return cdef;
}


void EXODUS::PrintBCDef(ostream& os, const EXODUS::elem_def& def)
{
  string mesh_entity;
  if (def.me==EXODUS::bceb) mesh_entity = "ElementBlock";
  else if (def.me==EXODUS::bcns) mesh_entity = "NodeSet";
  else if (def.me==EXODUS::bcss) mesh_entity = "SideSet";
  os << "The ELEMENT definition " << def.id << " refers to a " << mesh_entity << endl;
  os << "Sectionname: " << def.sec << endl;
  os << "Description: " << def.desc << endl;
  os << endl;
}
void EXODUS::PrintBCDef(ostream& os, const EXODUS::cond_def& def)
{
  string mesh_entity;
  if (def.me==EXODUS::bceb) mesh_entity = "ElementBlock";
  else if (def.me==EXODUS::bcns) mesh_entity = "NodeSet";
  else if (def.me==EXODUS::bcss) mesh_entity = "SideSet";
  os << "The CONDITION definition " << def.id << " refers to a " << mesh_entity << endl;
  os << "Sectionname: " << def.sec << endl;
  os << "Description: " << def.desc << endl;
  os << endl;
}

#endif //D_EXODUS
