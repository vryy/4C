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
    
    if (mesh_entity.compare(ebmarker)==0){
      // in case of eb we differntiate between 'element' or 'condition'
      if (type.compare("ELEMENT")==0){
        EXODUS::elem_def edef = EXODUS::ReadEdef(mesh_entity,id,actcond);
        eledefs.push_back(edef);
      } else if (type.compare("CONDITION")==0){
        EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity,id,actcond);
        condefs.push_back(cdef);
      } else {
        cout << "Undefined type for eb"<<id<<": "<<type<<endl;
        dserror("Undefined type!");
      }
    } else if (mesh_entity.compare(nsmarker)==0){
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity,id,actcond);
      condefs.push_back(cdef);
    } else if (mesh_entity.compare(ssmarker)==0){
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity,id,actcond);
      condefs.push_back(cdef);
    } else dserror("Undefined mesh_type");
      
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
  if (found!=string::npos) cdef.gtype = DRT::Condition::Point;
  found = secname.find("LINE");
  if (found!=string::npos) cdef.gtype = DRT::Condition::Line;
  found = secname.find("SURF");
  if (found!=string::npos) cdef.gtype = DRT::Condition::Surface;
  found = secname.find("VOL");
  if (found!=string::npos) cdef.gtype = DRT::Condition::Volume;
  
  // figure out number of 'E' topo-entity for datfile
  left = description.find_first_of("E");
  right = description.find_first_of("-");
  string Enum = description.substr(left+1,right-(left+1));
  // convert string to int
  istringstream Enumstream(Enum);
  Enumstream >> cdef.e_id;
  
  return cdef;
}

vector<string> EXODUS::ReadBCEntity(const string actcond)
{
  // vector of all specs of current condition 
  vector<string> specs;

  // condition number
  size_t left = 0;
  size_t right = actcond.find_first_of("\"",left+1);
  string buffer = actcond.substr(0,right);
  specs.push_back(buffer);
  
  // condition kind
  left = actcond.find("boundr_cond=\"");  // 13 chars
  right = actcond.find_first_of("\"",left+13);
  buffer = actcond.substr(left+13,right-(left+13));
  specs.push_back(buffer);
  
  // condition description
  left = actcond.find_first_of("\"",right+1);
  right = actcond.find_first_of("\"",left+1);
  //buffer = actcond.substr(right+2,actcond.find_first_of("\"",right+2)-(right+2));
  buffer = actcond.substr(left+1,right-(left+1));
  specs.push_back(buffer);
  
  // element type
  left = actcond.find("type=\"");  // 6 chars
  right = actcond.find_first_of("\"",left+6);
  buffer = actcond.substr(left+6,right-(left+6));
  specs.push_back(buffer);
  
  // element name
  left = actcond.find_first_of("\"",right+1);
  right = actcond.find_first_of("\"",left+1);
  //buffer = actcond.substr(right+2,actcond.find_first_of("\"",right+2)-(right+2));
  buffer = actcond.substr(left+1,right-(left+1));
  specs.push_back(buffer);
  
  // element shape
  left = actcond.find_first_of("\"",right+1);
  right = actcond.find_first_of("\"",left+1);
  //buffer = actcond.substr(right+2,actcond.find_first_of("\"",right+2)-(right+2));
  buffer = actcond.substr(left+1,right-(left+1));
  specs.push_back(buffer);

  // element properties
  left = actcond.find("prop=\"");  // 6 chars
  right = actcond.find_first_of("\"",left+6);
  buffer = actcond.substr(left+6,right-(left+6));
  specs.push_back(buffer);
  
  return specs;
}

inline EXODUS::cond_type EXODUS::CheckCondType(const string buffer)
{
  if      (buffer.compare(0,7,"ELEMENT")==0)   return EXODUS::element;
  else if (buffer.compare(0,4,"DVOL")==0)      return EXODUS::dvol;
  else if (buffer.compare(0,5,"DSURF")==0)     return EXODUS::dsurf;
  else if (buffer.compare(0,5,"DLINE")==0)     return EXODUS::dline;
  else if (buffer.compare(0,6,"DPOINT")==0)    return EXODUS::dpoint;
  else if (buffer.size()==0)                   return EXODUS::empty;
  else {
    cout << "Warning! No valid condition type '";
    cout << buffer << "' specified for entity " << endl;
    dserror ("ConditionType has to be ELEMENT, DVOL, DSURF, DLINE, DPOINT, or empty"); 
  }
  return EXODUS::invalid; // we should never come until here
}

inline string EXODUS::CondTypeToString(const EXODUS::cond_type cond)
{
  switch(cond){
  case EXODUS::element: return "ELEMENT";    break;
  case EXODUS::dvol:    return "DVOL   ";    break;
  case EXODUS::dsurf:   return "DSURF  ";    break;
  case EXODUS::dline:   return "DLINE  ";    break;
  case EXODUS::dpoint:  return "DPOINT ";    break;
  case EXODUS::empty:   return "";           break;
  default: dserror("Unknown CondType");
  }
  return "";
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
