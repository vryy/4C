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
#include "pre_exodus_soshextrusion.H" // to calculate normal
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"


using namespace Teuchos;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::WriteDatFile(const string& datfile, const EXODUS::Mesh& mymesh,
    const string& headfile, const vector<EXODUS::elem_def>& eledefs, const vector<EXODUS::cond_def>& condefs,
    const map<int,map<int,vector<vector<double> > > >& elecenterlineinfo)
{
  // open datfile
  std::ofstream dat(datfile.c_str());
  if (!dat) dserror("failed to open file: %s", datfile.c_str());

  // write dat-file intro
  EXODUS::WriteDatIntro(headfile,mymesh,dat);

  // write "header"
  EXODUS::WriteDatHead(headfile,dat);

  // write "design description"
  EXODUS::WriteDatDesign(condefs,dat);

  // write conditions
  EXODUS::WriteDatConditions(condefs,mymesh,dat);

  // write design-topology
  EXODUS::WriteDatDesignTopology(condefs,mymesh,dat);

  // write nodal coordinates
  EXODUS::WriteDatNodes(mymesh,dat);

  // write elements
  EXODUS::WriteDatEles(eledefs,mymesh,dat,elecenterlineinfo);

  // write END
  dat << "---------------------------------------------------------------END\n"\
         "// END\n";

  // close datfile
  if (dat.is_open()) dat.close();

  return 0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatIntro(const string& headfile, const EXODUS::Mesh& mymesh, ostream& dat)
{
  dat <<"==================================================================\n" \
        "                  General Data File BACI\n" \
  "==================================================================\n" \
  "-------------------------------------------------------------TITLE\n" \
  "created by pre_exodus \n" \
  "------------------------------------------------------PROBLEM SIZE\n";
  // print number of elements and nodes just as an comment instead of
  // a valid parameter (prevents possible misuse of these parameters in BACI)
  dat << "//ELEMENTS    " << mymesh.GetNumEle() << endl;
  dat << "//NODES       " << mymesh.GetNumNodes() << endl;
  // parameter for the number of spatial dimensions
  dat << "DIM           " << mymesh.GetBACIDim() << endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatHead(const string& headfile, ostream& dat)
{
  std::stringstream head;
  const char *headfilechar = headfile.c_str();
  std::ifstream header(headfilechar, std::ifstream::in);
  if (not header.good()){
    cout << endl << "Unable to open file: " << headfilechar << endl;
    dserror("Unable to open head-file");
  }
  while (header.good()) head << (char) header.get();
  //while (!header.eof()) head << (char) header.get();
  header.close();
  string headstring = head.str();

  // delete sections which will be written by WriteDatDesign()
  RemoveDatSection("PROBLEM SIZE", headstring);
  RemoveDatSection("DESIGN DESCRIPTION", headstring);

  // delete very first line with comment "//"
  if (headstring.find("//")== 0)
    headstring.erase(headstring.find("//"),headstring.find("\n")+1);//-headstring.find("//"));

  size_t comment = headstring.find("\n//");
  while (comment != string::npos){
    headstring.erase(comment+1,headstring.find("\n",comment+1)-comment);
    comment = headstring.find("\n//",comment);
  }

  // remove eof character
  headstring.erase(headstring.end()-1);

  // now put everything to the input file
  dat<<headstring<<endl;
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::RemoveDatSection(const string secname, string& headstring)
{
  const size_t secpos = headstring.find(secname);
  if(secpos!=string::npos)
  {
    // where does this section line actually start?
    const size_t endoflastline = headstring.substr(0,secpos).rfind("\n");
    // want to keep the newline character of line before
    const size_t sectionbegin = endoflastline+1;
    // now we remove the whole section
    const size_t sectionend = headstring.find("---", secpos);
    headstring.erase(sectionbegin,sectionend-sectionbegin);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatDesign(const vector<EXODUS::cond_def>& condefs, ostream& dat)
{
  vector<EXODUS::cond_def>::const_iterator it;
  int ndp=0; int ndl=0; int nds=0; int ndv=0;

  for (it=condefs.begin(); it != condefs.end(); ++it){
    EXODUS::cond_def acte = *it;
    if (acte.gtype == DRT::Condition::Volume)        ndv++;
    else if (acte.gtype == DRT::Condition::Surface)  nds++;
    else if (acte.gtype == DRT::Condition::Line)     ndl++;
    else if (acte.gtype == DRT::Condition::Point)    ndp++;
    else if (acte.gtype == DRT::Condition::NoGeom);
    else dserror ("Cannot identify Condition GeometryType");
  }
  dat << "------------------------------------------------DESIGN DESCRIPTION" << endl;
  dat << "NDPOINT " << ndp << endl;
  dat << "NDLINE  " << ndl << endl;
  dat << "NDSURF  " << nds << endl;
  dat << "NDVOL   " << ndv << endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatConditions(const vector<EXODUS::cond_def>& condefs,const EXODUS::Mesh& mymesh, ostream& dat)
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > condlist = DRT::INPUT::ValidConditions();

  // count how often we have one specific condition
  map<string,vector<int> > count_cond;
  map<string,vector<int> >::const_iterator count;
  vector<int>::const_iterator i_c;
  for (unsigned int i_cond = 0; i_cond < condefs.size(); ++i_cond)
    (count_cond[condefs.at(i_cond).sec]).push_back(i_cond);

  // loop all valid conditions that BACI knows
  for (unsigned int i=0; i<(*condlist).size(); ++i)
  {
    size_t linelength = 66;
    string sectionname = (*condlist)[i]->SectionName();

    // ignore conditions occurring zero times
    count = count_cond.find(sectionname);
    if (count == count_cond.end())
      continue;

    // only write conditions provided by user
    string dash(linelength-sectionname.size(),'-');
    dat << dash << sectionname << endl;
    DRT::Condition::GeometryType gtype = (*condlist)[i]->GeometryType();
    string geo;
    switch (gtype)
    {
    case DRT::Condition::Point:   geo = "DPOINT "; break;
    case DRT::Condition::Line:    geo = "DLINE  "; break;
    case DRT::Condition::Surface: geo = "DSURF  "; break;
    case DRT::Condition::Volume:  geo = "DVOL   "; break;
    default:
      dserror("geometry type unspecified");
    }

    dat << geo << (count->second).size() << endl;
    for (i_c=(count->second).begin();i_c!=(count->second).end();++i_c)
    {
      EXODUS::cond_def actcon = condefs[*i_c];
      string name;
      string pname;
      if (actcon.me==EXODUS::bcns){
        name = (mymesh.GetNodeSet(actcon.id).GetName());
        pname = (mymesh.GetNodeSet(actcon.id).GetPropName());
      } else if (actcon.me==EXODUS::bceb){
        name = (mymesh.GetElementBlock(actcon.id)->GetName());
      } else if (actcon.me==EXODUS::bcss){
        name = (mymesh.GetSideSet(actcon.id).GetName());
      } else dserror ("Unidentified Actcon");
      if((name!="")){
        dat << "// " << name;
        if (pname!="none"){ dat << " " << pname;}
        dat << endl;
      } else if (pname!="none") dat << "// " << pname << endl;

      // write the condition
      if (actcon.desc != ""){
        dat << "E " << actcon.e_id << " - " << actcon.desc << endl;
      } else if ((actcon.sec == "DESIGN SURF LOCSYS CONDITIONS") && (actcon.me==EXODUS::bcns)) {
        // special case for locsys conditions: calculate normal
        vector<double> normtang = EXODUS::CalcNormalSurfLocsys(actcon.id,mymesh);
        dat << "E " << actcon.e_id << " - ";
        for (unsigned int i = 0; i < normtang.size() ; ++i) dat <<  std::setprecision (10) << std::fixed << normtang[i] << " ";
        dat << endl;
      } else
        dat << "E " << actcon.e_id << " - " << actcon.desc << endl;
    }
    // remove sectionname from map, since writing is done
    count_cond.erase(sectionname);
  }
  if (count_cond.size() > 0) // there are conditions left that were not recognized!!
  {
    cout<<endl<<endl;
    for (count = count_cond.begin(); count != count_cond.end(); ++count)
      cout<<"Section name  "<< count->first <<"  is not valid. Typo?"<<endl;

    dserror("There are invalid condition names in your bc file (see list above)");
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
vector<double> EXODUS::CalcNormalSurfLocsys(const int ns_id,const EXODUS::Mesh& m)
{
  vector<double> normaltangent;
  EXODUS::NodeSet ns = m.GetNodeSet(ns_id);

  std::set<int> nodes_from_nodeset = ns.GetNodeSet();
  std::set<int>::iterator it;

  // compute normal
  std::set<int>::iterator surfit = nodes_from_nodeset.begin();
  int origin = *surfit;      // get first set node
  ++surfit;
  int head1 = *surfit;  // get second set node
  ++surfit;

  std::set<int>::iterator thirdnode;

  // find third node such that a proper normal can be computed
  for(it=surfit;it!=nodes_from_nodeset.end();++it){
    thirdnode = it;
    normaltangent = EXODUS::Normal(head1,origin,*thirdnode,m);
    if (normaltangent.size() != 1) break;
  }
  if (normaltangent.size()==1){
    dserror("Warning! No normal defined for SurfLocsys within nodeset '%s'!",(ns.GetName()).c_str());
  }

  // find tangent by Gram-Schmidt
  vector<double> t(3);
  t.at(0) = 1.0; t.at(1) = 0.0; t.at(2) = 0.0; // try this one
  double sp = t[0]*normaltangent[0] + t[1]*normaltangent[1] + t[2]*normaltangent[2]; // scalar product
  // subtract projection
  t.at(0) -= normaltangent[0]*sp;
  t.at(1) -= normaltangent[1]*sp;
  t.at(2) -= normaltangent[2]*sp;

  // very unlucky case
  if (t.at(0) < 1.0E-14){
    t.at(0) = 0.0; t.at(1) = 1.0; t.at(2) = 0.0; // rather use this
  }

  normaltangent.push_back(t.at(0));
  normaltangent.push_back(t.at(1));
  normaltangent.push_back(t.at(2));

  return normaltangent;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatDesignTopology(const vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, ostream& dat)
{
  // sort baciconds w.r.t. underlying topology
  map<int,EXODUS::cond_def> dpoints;
  map<int,EXODUS::cond_def> dlines;
  map<int,EXODUS::cond_def> dsurfs;
  map<int,EXODUS::cond_def> dvols;

  map<int,EXODUS::cond_def>::const_iterator it;
  vector<EXODUS::cond_def>::const_iterator i;
  for (i=condefs.begin();i!=condefs.end();++i){
    EXODUS::cond_def acte = *i;
    if (acte.gtype==DRT::Condition::Point) dpoints.insert(std::pair<int,EXODUS::cond_def>(acte.e_id,acte));
    else if (acte.gtype==DRT::Condition::Line) dlines.insert(std::pair<int,EXODUS::cond_def>(acte.e_id,acte));
    else if (acte.gtype==DRT::Condition::Surface) dsurfs.insert(std::pair<int,EXODUS::cond_def>(acte.e_id,acte));
    else if (acte.gtype==DRT::Condition::Volume) dvols.insert(std::pair<int,EXODUS::cond_def>(acte.e_id,acte));
    else if (acte.gtype==DRT::Condition::NoGeom);
    else dserror ("Cannot identify Condition GeometryType");
  }

  dat << "-----------------------------------------------DNODE-NODE TOPOLOGY"<<endl;
  for (it=dpoints.begin();it!=dpoints.end();++it){
    EXODUS::cond_def acte = it->second;
    const std::set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    std::set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DNODE " << acte.e_id << endl;
      //dat << EXODUS::CondGeomTypeToString(acte) << " " << acte.e_id << endl;
    }
  }
  dat << "-----------------------------------------------DLINE-NODE TOPOLOGY"<<endl;
  for (it=dlines.begin();it!=dlines.end();++it){
    EXODUS::cond_def acte = it->second;
    const std::set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    std::set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DLINE " << acte.e_id << endl;
    }
  }
  dat << "-----------------------------------------------DSURF-NODE TOPOLOGY"<<endl;
  for (it=dsurfs.begin();it!=dsurfs.end();++it){
    EXODUS::cond_def acte = it->second;
    const std::set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    std::set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DSURFACE " << acte.e_id << endl;
    }
  }
  dat << "------------------------------------------------DVOL-NODE TOPOLOGY"<<endl;
  for (it=dvols.begin();it!=dvols.end();++it){
    EXODUS::cond_def acte = it->second;
    const std::set<int> nodes = EXODUS::GetNsFromBCEntity(acte,mymesh);
    std::set<int>::const_iterator i;
    for(i=nodes.begin();i!=nodes.end();++i){
      dat << "NODE    " << *i << " " << "DVOL " << acte.e_id << endl;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const std::set<int> EXODUS::GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m)
{
  if (e.me==EXODUS::bcns){
    EXODUS::NodeSet ns = m.GetNodeSet(e.id);
    return ns.GetNodeSet();
  } else if (e.me==EXODUS::bceb){
    std::set<int> allnodes;
    RCP<EXODUS::ElementBlock> eb = m.GetElementBlock(e.id);
    RCP<const map<int,vector<int> > > eles = eb->GetEleConn();
    map<int,vector<int> >::const_iterator i_ele;
    for(i_ele=eles->begin(); i_ele != eles->end(); ++i_ele){
      const vector<int> nodes = i_ele->second;
      vector<int>::const_iterator i;
      for(i=nodes.begin();i!=nodes.end();++i) allnodes.insert(*i);
    }
    return allnodes;
  } else if (e.me==EXODUS::bcss){
    std::set<int> allnodes;
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
  std::set<int> n;
  return n;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatNodes(const EXODUS::Mesh& mymesh, ostream& dat)
{
  dat << "-------------------------------------------------------NODE COORDS" << endl;
  dat.precision(16);
  RCP<map<int,vector<double> > > nodes = mymesh.GetNodes();
  map<int,vector<double> >::const_iterator i_node;
  for (i_node = nodes->begin(); i_node != nodes->end(); ++i_node)
  {
    vector<double> coords = i_node->second;
    dat << "NODE " << i_node->first << "  " << '\t' << "COORD" << '\t';
    for(unsigned int i=0; i<coords.size(); ++i) dat << std::scientific << coords[i] << '\t';
    dat << endl;
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatEles(const vector<elem_def>& eledefs, const EXODUS::Mesh& mymesh, ostream& dat,
     const map<int,map<int,vector<vector<double> > > >& elecenterlineinfo)
{
  // sort elements w.r.t. structure, fluid, ale, scalar transport, thermo, etc.
  vector<EXODUS::elem_def> strus;
  vector<EXODUS::elem_def> fluids;
  vector<EXODUS::elem_def> ales;
  vector<EXODUS::elem_def> transport;
  vector<EXODUS::elem_def> thermo;
  vector<EXODUS::elem_def>::const_iterator i_et;

  for(i_et=eledefs.begin();i_et!=eledefs.end();++i_et){
    EXODUS::elem_def acte = *i_et;
    if (acte.sec.compare("STRUCTURE")==0) strus.push_back(acte);
    else if (acte.sec.compare("FLUID")==0) fluids.push_back(acte);
    else if (acte.sec.compare("ALE")==0) ales.push_back(acte);
    else if (acte.sec.compare("TRANSPORT")==0) transport.push_back(acte);
    else if (acte.sec.compare("THERMO")==0) thermo.push_back(acte);
    else if (acte.sec.compare("")==0);
    else{
      cout << "Unknown ELEMENT sectionname in eb" << acte.id << ": '" << acte.sec << "'!" << endl;
      dserror("Unknown ELEMENT sectionname");
    }
  }

  int startele = 1; // BACI-Dat eles start with 1, this int is adapted for more than one element section

  // print structure elements
  dat << "------------------------------------------------STRUCTURE ELEMENTS" << endl;
  for(i_et=strus.begin();i_et!=strus.end();++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb,acte,startele,dat,elecenterlineinfo,acte.id);
  }

  // print fluid elements
  dat << "----------------------------------------------------FLUID ELEMENTS" << endl;
  for(i_et=fluids.begin();i_et!=fluids.end();++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb,acte,startele,dat,elecenterlineinfo,acte.id);
  }

  // print ale elements
  dat << "------------------------------------------------------ALE ELEMENTS" << endl;
  for(i_et=ales.begin();i_et!=ales.end();++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb,acte,startele,dat,elecenterlineinfo,acte.id);
  }

  // print transport elements
  dat << "------------------------------------------------TRANSPORT ELEMENTS" << endl;
  for(i_et=transport.begin();i_et!=transport.end();++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb,acte,startele,dat,elecenterlineinfo,acte.id);
  }

  // print thermo elements
  dat << "---------------------------------------------------THERMO ELEMENTS" << endl;
  for(i_et=thermo.begin();i_et!=thermo.end();++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb,acte,startele,dat,elecenterlineinfo,acte.id);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::DatEles(RCP< const EXODUS::ElementBlock> eb, const EXODUS::elem_def& acte, int& startele, ostream& datfile,
    const map<int,map<int,vector<vector<double> > > >& elescli,const int eb_id)
{
  RCP<const map<int,vector<int> > > eles = eb->GetEleConn();
  map<int,vector<int> >::const_iterator i_ele;
  for (i_ele=eles->begin();i_ele!=eles->end();++i_ele)
  {
    std::stringstream dat; // first build up the string for actual element line
    const vector<int> nodes = i_ele->second;
    vector<int>::const_iterator i_n;
    dat << "   " << startele;
    dat << " " << acte.ename;                       // e.g. "SOLIDH8"
    dat << " " << DistypeToString(PreShapeToDrt(eb->GetShape()));
    dat << "  ";
    for(i_n=nodes.begin();i_n!=nodes.end();++i_n) dat << *i_n << " ";
    dat << "   " << acte.desc;              // e.g. "MAT 1"
    if(elescli.size()!=0){
    	// quick check wether elements in ele Block have a fiber direction
    	if(elescli.find(eb_id)!=elescli.end())
    	{
    		// write local cosy from centerline to each element
    		vector<vector<double> > ecli = (elescli.find(eb_id)->second).find(i_ele->first)->second;
		  dat << " RAD " << std::fixed << std::setprecision(8) << ecli[0][0] << " " << ecli[0][1] << " " << ecli[0][2];
		  dat << " AXI " << ecli[1][0] << " " << ecli[1][1] << " " << ecli[1][2];
		  dat << " CIR " << ecli[2][0] << " " << ecli[2][1] << " " << ecli[2][2];
    	}
    }
    dat << endl;  // finish this element line

    startele ++;
    datfile<<dat.str(); // only one access to the outfile (saves system time)
  }
  return;
}


#endif //D_EXODUS
