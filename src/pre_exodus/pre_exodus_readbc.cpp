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

vector<map<int,EXODUS::bc_entity> > EXODUS::ReadBCFile(const string bcfile)
{
  map<int,EXODUS::bc_entity> eb_bcs; // all bc_entities for ElementBlocks
  map<int,EXODUS::bc_entity> ns_bcs; // all bc_entities for NodeSets
  map<int,EXODUS::bc_entity> ss_bcs; // all bc_entities for SideSets

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
    
    // store current entity
    EXODUS::bc_entity actentity;
    
    // find out what mesh_entity type we have
    string mesh_entity = actcond.substr(0,3);
    if (mesh_entity.compare(ebmarker)==0) actentity.me = EXODUS::bceb;
    else if (mesh_entity.compare(nsmarker)==0) actentity.me = EXODUS::bcns;
    else if (mesh_entity.compare(ssmarker)==0) actentity.me = EXODUS::bcss;
    else {
      cout << "Don't understand mesh_entity specifier: " << mesh_entity;
      cout << "Has to be " << ebmarker <<", "<<nsmarker<<" or "<<ssmarker<<endl;
      dserror("Mesh_entity not found");
    }
    
    // get its id
    size_t found2 = actcond.find_first_of("=");
    string buffer = actcond.substr(markerlength,found2-markerlength);
    // convert string to int
    istringstream bufferstream(buffer);
    int id;
    bufferstream >> id;
    actentity.id = id;
    
    // condition type
    size_t left = actcond.find_first_of("\"",0);
    size_t right = actcond.find_first_of("\"",left+1);
    buffer = actcond.substr(left+1,right-left-1);
    actentity.ct = CheckCondType(buffer);
    
    // reduce actcond by mesh_entity, id, and condition_type
    left = actcond.find_first_of("\"",right+1);
    right = actcond.find_last_of("\"");
    buffer = actcond.substr(left+1);//,right-(left+1));
    
    // read the rest of actcond into vector of strings
    actentity.specs = ReadBCEntity(buffer);
    
    //cout << "ID: " << id << endl;
    //EXODUS::PrintBCEntity(cout,actentity);
    
    // insert into corresponding map
    if (actentity.me==EXODUS::bceb) eb_bcs.insert(pair<int,EXODUS::bc_entity>(id,actentity));
    else if (actentity.me==EXODUS::bcns) ns_bcs.insert(pair<int,EXODUS::bc_entity>(id,actentity));
    else if (actentity.me==EXODUS::bcss) ss_bcs.insert(pair<int,EXODUS::bc_entity>(id,actentity));
    else dserror("Problems with bc entity");
  }
  
  vector<map<int,EXODUS::bc_entity> > allbcs(3);
  allbcs[0] = eb_bcs;
  allbcs[1] = ns_bcs;
  allbcs[2] = ss_bcs;
  
  return allbcs;
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

void EXODUS::PrintBCEntity(ostream& os, const EXODUS::bc_entity ent)
{
  string mesh_entity;
  if (ent.me==EXODUS::bceb) mesh_entity = "ElementBlock";
  else if (ent.me==EXODUS::bcns) mesh_entity = "NodeSet";
  else if (ent.me==EXODUS::bcss) mesh_entity = "SideSet";
  os << "The BC-Entity refers to a " << mesh_entity;
  os << ", is of type " << EXODUS::CondTypeToString(ent.ct) << " and has specs:" << endl;
  const vector<string> specs = ent.specs;
  vector<string>::const_iterator it;
  for (it=specs.begin(); it != specs.end(); ++it) os << *it << endl;
  os << endl;
}

#endif //D_EXODUS
