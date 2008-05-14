/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_writedat.cpp

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
#include "pre_exodus_writedat.H"
#include "pre_exodus_reader.H"
#include "../drt_lib/drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"


using namespace std;
using namespace Teuchos;

int EXODUS::WriteDatFile(const string& datfile, const EXODUS::Mesh& mymesh,
    const string& headfile, const vector<EXODUS::elem_def>& eledefs, const vector<EXODUS::cond_def>& condefs)
{
  // open datfile
  ofstream dat(datfile.c_str());
  if (!dat) dserror("failed to open file: %s", datfile.c_str());

  // write dat-file intro
  string datintro = EXODUS::WriteDatIntro(mymesh);
  dat << datintro;
  
  // write "header"
  string dathead = EXODUS::WriteDatHead(headfile);
  dat << dathead;
  
  // write "design description"
  string datdesign = EXODUS::WriteDatDesign(condefs);
  dat << datdesign;
  
  // write conditions
  string datconditions = EXODUS::WriteDatConditions(condefs);
  dat << datconditions;
  
  // write design-topology
  string datdesigntopo = EXODUS::WriteDatDesignTopology(condefs,mymesh);
  dat << datdesigntopo;
  
  // write nodal coordinates
  string datnodes = EXODUS::WriteDatNodes(mymesh);
  dat << datnodes;
  
  // write elements
  string dateles = EXODUS::WriteDatEles(eledefs,mymesh);
  dat << dateles;

  // write END
  dat << "---------------------------------------------------------------END\n"\
         "// END\n";
  
  // close datfile
  if (dat.is_open()) dat.close();

  return 0;
}


string EXODUS::WriteDatIntro(const EXODUS::Mesh& mymesh)
{
  stringstream dat;
  dat <<"==================================================================\n" \
        "        General Data File CCARAT\n" \
  "==================================================================\n" \
  "-------------------------------------------------------------TITLE\n" \
  "created by pre_exodus \n" \
  "------------------------------------------------------PROBLEM SIZE\n";
  dat << "ELEMENTS " << '\t' << mymesh.GetNumEle() << endl;
  dat << "NODES    " << '\t' << mymesh.GetNumNodes() << endl;
  dat << "DIM      " << '\t' << mymesh.GetNumDim() << endl;
  dat << "MATERIALS" << '\t' << "1" << endl;
  dat << "NUMDF    " << '\t' << "6" << endl;

  return dat.str();
}

string EXODUS::WriteDatHead(const string& headfile)
{
  stringstream head;
  const char *headfilechar;
  headfilechar = headfile.c_str();
  ifstream header(headfilechar, ifstream::in);
  while (header.good()) head << (char) header.get();
  //while (!header.eof()) head << (char) header.get();
  header.close();
  string headstring = head.str();
  headstring.erase(headstring.end()-1);
  return headstring;
}

vector<EXODUS::bc_entity> EXODUS::FindBACIConditions(const vector<map<int,EXODUS::bc_entity> >& bc_specs)
{
  vector<EXODUS::bc_entity> bacis;
  map<int,EXODUS::bc_entity>::const_iterator it;
  for (unsigned int i = 0; i < bc_specs.size(); ++i) {
    map<int,EXODUS::bc_entity> actmap = bc_specs.at(i);
    for(it=actmap.begin();it!=actmap.end();++it){
      EXODUS::bc_entity actentity = it->second;
      if(actentity.ct==EXODUS::dvol || actentity.ct==EXODUS::dsurf 
          || actentity.ct==EXODUS::dline || actentity.ct==EXODUS::dpoint){
        bacis.push_back(actentity);
      }
    }
  }
  return bacis;
}

string EXODUS::WriteDatDesign(const vector<EXODUS::cond_def>& condefs)
{
  stringstream dat;
  vector<EXODUS::cond_def>::const_iterator it;
  int ndp=0; int ndl=0; int nds=0; int ndv=0;

  for (it=condefs.begin(); it != condefs.end(); ++it){
    EXODUS::cond_def acte = *it;
    if (acte.gtype == DRT::Condition::Volume)        ndv++;
    else if (acte.gtype == DRT::Condition::Surface)  nds++;
    else if (acte.gtype == DRT::Condition::Line)     ndl++;
    else if (acte.gtype == DRT::Condition::Point)    ndp++;
    else if (acte.gtype == DRT::Condition::NoGeom) break;
    else dserror ("Cannot identify Condition GeometryType");
  }
  dat << "------------------------------------------------DESIGN DESCRIPTION" << endl;
  dat << "NDPOINT " << ndp << endl;
  dat << "NDLINE  " << ndl << endl;
  dat << "NDSURF  " << nds << endl;
  dat << "NDVOL   " << ndv << endl;
  
  dat << "-----------------------------------------------------DESIGN POINTS" << endl;
  dat << "------------------------------------------------------DESIGN LINES" << endl;
  dat << "---------------------------------------------------DESIGN SURFACES" << endl;
  dat << "----------------------------------------------------DESIGN VOLUMES" << endl;
  
  return dat.str();
}

string EXODUS::WriteDatConditions(const vector<EXODUS::cond_def>& condefs)
{
  stringstream dat;
  
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > condlist = DRT::INPUT::ValidConditions();
  vector<string> allsectionnames;
  
  // count how often we have one specific condition
  map<string,vector<int> > count_cond;
  map<string,vector<int> >::const_iterator count;
  vector<int>::const_iterator i_c;
  vector<EXODUS::cond_def>::const_iterator i_cond;
  for (i_cond = condefs.begin(); i_cond != condefs.end(); ++i_cond)
    (count_cond[(*i_cond).sec]).push_back((*i_cond).id);
  
  
  
  for (unsigned int i=0; i<(*condlist).size(); ++i)
  {
    size_t linelength = 66;
    string sectionname = (*condlist)[i]->SectionName();
    allsectionnames.push_back(sectionname);
    string dash(linelength-sectionname.size(),'-');
    dat << dash <<sectionname << endl;
    DRT::Condition::GeometryType gtype = (*condlist)[i]->GeometryType();
    string geo;
    switch (gtype)
    {
    case DRT::Condition::Point:   geo = "DPOINT "; break;
    case DRT::Condition::Line:    geo = "DLINE  ";  break;
    case DRT::Condition::Surface: geo = "DSURF  ";  break;
    case DRT::Condition::Volume:  geo = "DVOL   ";   break;
    default:
      dserror("geometry type unspecified");
    }
    count = count_cond.find(sectionname);
    if (count == count_cond.end()) dat << geo << "0" << endl;
    else  {
      dat << geo << (count->second).size() << endl;
      for (i_c=(count->second).begin();i_c!=(count->second).end();++i_c){
        EXODUS::cond_def actcon = condefs[*i_c];
        dat << "E " << actcon.e_id << " - " << actcon.desc << endl;
      }
    }
  }
  
  
  return dat.str();
}

string EXODUS::WriteDatDesignTopology(const vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh)
{
  stringstream dat;
  
  // sort baciconds w.r.t. underlying topology
  vector<EXODUS::cond_def> dpoints;
  vector<EXODUS::cond_def> dlines;
  vector<EXODUS::cond_def> dsurfs;
  vector<EXODUS::cond_def> dvols;
  
  vector<EXODUS::cond_def>::const_iterator it;
  for (it=condefs.begin();it!=condefs.end();++it){
    EXODUS::cond_def acte = *it;
    if (acte.gtype==DRT::Condition::Point) dpoints.push_back(acte);
    else if (acte.gtype==DRT::Condition::Line) dlines.push_back(acte);
    else if (acte.gtype==DRT::Condition::Surface) dsurfs.push_back(acte);
    else if (acte.gtype==DRT::Condition::Volume) dvols.push_back(acte);
    else if (acte.gtype==DRT::Condition::NoGeom) break;
    else dserror ("Cannot identify Condition GeometryType");
  }
  
  dat << "-----------------------------------------------DNODE-NODE TOPOLOGY"<<endl;
  for (it=dpoints.begin();it!=dpoints.end();++it){
    EXODUS::cond_def acte = *it;
    const set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DNODE " << acte.e_id << endl;
      //dat << EXODUS::CondGeomTypeToString(acte) << " " << acte.e_id << endl;
    }
  }
  dat << "-----------------------------------------------DLINE-NODE TOPOLOGY"<<endl;
  for (it=dlines.begin();it!=dlines.end();++it){
    EXODUS::cond_def acte = *it;
    const set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DNODE " << acte.e_id << endl;
    }
  }
  dat << "-----------------------------------------------DSURF-NODE TOPOLOGY"<<endl;
  for (it=dsurfs.begin();it!=dsurfs.end();++it){
    EXODUS::cond_def acte = *it;
    const set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DSURFACE " << acte.e_id << endl;
    }
  }
  dat << "------------------------------------------------DVOL-NODE TOPOLOGY"<<endl;
  for (it=dvols.begin();it!=dvols.end();++it){
    EXODUS::cond_def acte = *it;
    const set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DNODE " << acte.e_id << endl;
    }
  }

  return dat.str();
}

const set<int> EXODUS::GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m)
{
  if (e.me==EXODUS::bcns){
    EXODUS::NodeSet ns = m.GetNodeSet(e.id);
    return ns.GetNodeSet();
  } else if (e.me==EXODUS::bceb){
    set<int> allnodes;
    EXODUS::ElementBlock& eb = m.GetElementBlock(e.id);
    RCP<const map<int,vector<int> > > eles = eb.GetEleConn();
    map<int,vector<int> >::const_iterator i_ele;
    for(i_ele=eles->begin(); i_ele != eles->end(); ++i_ele){
      const vector<int> nodes = i_ele->second;
      vector<int>::const_iterator i;
      for(i=nodes.begin();i!=nodes.end();++i) allnodes.insert(*i);
    }
    return allnodes;
  } else if (e.me==EXODUS::bcss){
    set<int> allnodes;
    EXODUS::SideSet ss = m.GetSideSet(e.id);
    const map<int,vector<int> > eles = ss.GetSideSet();
    map<int,vector<int> >::const_iterator i_ele;
    for(i_ele=eles.begin(); i_ele != eles.end(); ++i_ele){
      const vector<int> nodes = i_ele->second;
      vector<int>::const_iterator i;
      for(i=nodes.begin();i!=nodes.end();++i) allnodes.insert(*i);
    }
    return allnodes;
  } else dserror("Cannot identify mesh_entity");
  set<int> n;
  return n;
}

string EXODUS::WriteDatNodes(const EXODUS::Mesh& mymesh)
{
  stringstream dat;
  dat << "-------------------------------------------------------NODE COORDS" << endl;
  dat.precision(8);
  map<int,vector<double> > nodes = mymesh.GetNodes();
  map<int,vector<double> >::const_iterator i_node;
  for (i_node = nodes.begin(); i_node != nodes.end(); ++i_node) {
    vector<double> coords = i_node->second;
    dat << "NODE " << i_node->first+1 << "  " << '\t' << "COORD" << '\t';
    for(unsigned int i=0; i<coords.size(); ++i) dat << scientific << coords[i] << '\t';
    dat << endl;
  }
  return dat.str();
}

string EXODUS::WriteDatEles(const vector<elem_def>& eledefs, const EXODUS::Mesh& mymesh)
{
  stringstream dat;
  dat << "----------------------------------------------------------ELEMENTS" << endl;
  
  // sort elements w.r.t. structure, fluid, ale, etc.
  vector<EXODUS::elem_def> strus;
  vector<EXODUS::elem_def> fluids;
  vector<EXODUS::elem_def> ales;
  vector<EXODUS::elem_def> levels;
  vector<EXODUS::elem_def>::const_iterator i_et;
  
  for(i_et=eledefs.begin();i_et!=eledefs.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    if (acte.sec.compare("STRUCTURE")==0) strus.push_back(acte);
    else if (acte.sec.compare("FLUID")==0) fluids.push_back(acte);
    else if (acte.sec.compare("ALE")==0) ales.push_back(acte);
    else if (acte.sec.compare("LEVELSET")==0) levels.push_back(acte);
    else{
      cout << "Unknown ELEMENT sectionname in eb" << acte.id << ": '" << acte.sec << "'!" << endl;
      dserror("Unknown ELEMENT sectionname");
    }
  }
  
  int ele = 1; // BACI-Dat eles start with 1
  
  // print structure elements
  dat << "------------------------------------------------STRUCTURE ELEMENTS" << endl;
  for(i_et=strus.begin();i_et!=strus.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    EXODUS::ElementBlock eb = mymesh.GetElementBlock(acte.id);
    string dateles = EXODUS::DatEles(eb,acte,ele);
    dat << dateles;
  }

  // print fluid elements
  dat << "----------------------------------------------------FLUID ELEMENTS" << endl;
  for(i_et=fluids.begin();i_et!=fluids.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    EXODUS::ElementBlock eb = mymesh.GetElementBlock(acte.id);
    string dateles = EXODUS::DatEles(eb,acte,ele);
    dat << dateles;
  }

  // print ale elements
  dat << "------------------------------------------------------ALE ELEMENTS" << endl;
  for(i_et=ales.begin();i_et!=ales.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    EXODUS::ElementBlock eb = mymesh.GetElementBlock(acte.id);
    string dateles = EXODUS::DatEles(eb,acte,ele);
    dat << dateles;
  }

  // print levelset elements
  dat << "-------------------------------------------------LEVELSET ELEMENTS" << endl;
  for(i_et=levels.begin();i_et!=levels.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    EXODUS::ElementBlock eb = mymesh.GetElementBlock(acte.id);
    string dateles = EXODUS::DatEles(eb,acte,ele);
    dat << dateles;
  }

  return dat.str();
}

string EXODUS::DatEles(const EXODUS::ElementBlock& eb, const EXODUS::elem_def& acte, int& struele)
{
  stringstream dat;
  RCP<const map<int,vector<int> > > eles = eb.GetEleConn();
  map<int,vector<int> >::const_iterator i_ele;
  for (i_ele=eles->begin();i_ele!=eles->end();++i_ele){
    const vector<int> nodes = i_ele->second;
    vector<int>::const_iterator i_n;
    dat << "   " << struele;
    dat << " " << acte.ename;                       // e.g. "SOLIDH8"
    dat << " " << DistypeToString(PreShapeToDrt(eb.GetShape()));
    dat << "  ";
    for(i_n=nodes.begin();i_n!=nodes.end();++i_n) dat << *i_n << " ";
    dat << "   " << acte.desc << endl;              // e.g. "MAT 1"
    struele ++;
  }
  return dat.str();
}


#endif //D_EXODUS
