/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus .dat-file writer

\level 1


Here is everything related with writing a dat-file
 */
/*----------------------------------------------------------------------*/
#include "pre_exodus_writedat.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_soshextrusion.H"  // to calculate normal
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::WriteDatFile(const std::string& datfile, const EXODUS::Mesh& mymesh,
    const std::string& headfile, const std::vector<EXODUS::elem_def>& eledefs,
    const std::vector<EXODUS::cond_def>& condefs,
    const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo)
{
  // open datfile
  std::ofstream dat(datfile.c_str());
  if (!dat) dserror("failed to open file: %s", datfile.c_str());

  // write dat-file intro
  EXODUS::WriteDatIntro(headfile, mymesh, dat);

  // write "header"
  EXODUS::WriteDatHead(headfile, dat);

  // write "design description"
  EXODUS::WriteDatDesign(condefs, dat);

  // write conditions
  EXODUS::WriteDatConditions(condefs, mymesh, dat);

  // write design-topology
  EXODUS::WriteDatDesignTopology(condefs, mymesh, dat);

  // write nodal coordinates
  EXODUS::WriteDatNodes(mymesh, dat);

  // write elements
  EXODUS::WriteDatEles(eledefs, mymesh, dat, elecenterlineinfo);

  // write END
  dat << "---------------------------------------------------------------END\n"
         "// END\n";

  // close datfile
  if (dat.is_open()) dat.close();

  return 0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatIntro(
    const std::string& headfile, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  dat << "==================================================================\n"
         "                  General Data File BACI\n"
         "==================================================================\n"
         "-------------------------------------------------------------TITLE\n"
         "created by pre_exodus\n"
         "------------------------------------------------------PROBLEM SIZE\n";
  // print number of elements and nodes just as an comment instead of
  // a valid parameter (prevents possible misuse of these parameters in BACI)
  dat << "//ELEMENTS    " << mymesh.GetNumEle() << std::endl;
  dat << "//NODES       " << mymesh.GetNumNodes() << std::endl;
  // parameter for the number of spatial dimensions
  dat << "DIM           " << mymesh.GetBACIDim() << std::endl;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatHead(const std::string& headfile, std::ostream& dat)
{
  std::stringstream head;
  const char* headfilechar = headfile.c_str();
  std::ifstream header(headfilechar, std::ifstream::in);
  if (not header.good())
  {
    std::cout << std::endl << "Unable to open file: " << headfilechar << std::endl;
    dserror("Unable to open head-file");
  }
  while (header.good()) head << (char)header.get();
  // while (!header.eof()) head << (char) header.get();
  header.close();
  std::string headstring = head.str();

  // delete sections which will be written by WriteDatDesign()
  RemoveDatSection("PROBLEM SIZE", headstring);
  RemoveDatSection("DESIGN DESCRIPTION", headstring);

  // delete very first line with comment "//"
  if (headstring.find("//") == 0)
    headstring.erase(
        headstring.find("//"), headstring.find('\n') + 1);  //-headstd::string.find("//"));

  size_t comment = headstring.find("\n//");
  while (comment != std::string::npos)
  {
    headstring.erase(comment + 1, headstring.find('\n', comment + 1) - comment);
    comment = headstring.find("\n//", comment);
  }

  // remove eof character
  headstring.erase(headstring.end() - 1);

  // now put everything to the input file
  dat << headstring << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::RemoveDatSection(const std::string& secname, std::string& headstring)
{
  const size_t secpos = headstring.find(secname);
  if (secpos != std::string::npos)
  {
    // where does this section line actually start?
    const size_t endoflastline = headstring.substr(0, secpos).rfind('\n');
    // want to keep the newline character of line before
    const size_t sectionbegin = endoflastline + 1;
    // now we remove the whole section
    const size_t sectionend = headstring.find("---", secpos);
    headstring.erase(sectionbegin, sectionend - sectionbegin);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatDesign(const std::vector<EXODUS::cond_def>& condefs, std::ostream& dat)
{
  int ndp = 0;
  int ndl = 0;
  int nds = 0;
  int ndv = 0;

  for (const auto& condition_definition : condefs)
  {
    switch (condition_definition.gtype)
    {
      case DRT::Condition::Volume:
        ++ndv;
        break;
      case DRT::Condition::Surface:
        ++nds;
        break;
      case DRT::Condition::Line:
        ++ndl;
        break;
      case DRT::Condition::Point:
        ++ndp;
        break;
      case DRT::Condition::NoGeom:
        // do nothing
        break;
      default:
        dserror("Cannot identify Condition GeometryType");
        break;
    }
  }
  dat << "------------------------------------------------DESIGN DESCRIPTION" << std::endl;
  dat << "NDPOINT " << ndp << std::endl;
  dat << "NDLINE  " << ndl << std::endl;
  dat << "NDSURF  " << nds << std::endl;
  dat << "NDVOL   " << ndv << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatConditions(
    const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> condlist =
      DRT::INPUT::ValidConditions();

  // count how often we have one specific condition
  std::map<std::string, std::vector<int>> count_cond;
  std::map<std::string, std::vector<int>>::const_iterator count;
  std::vector<int>::const_iterator i_c;
  for (int i_cond = 0; i_cond < static_cast<int>(condefs.size()); ++i_cond)
    (count_cond[condefs.at(i_cond).sec]).push_back(i_cond);

  // loop all valid conditions that BACI knows
  for (auto& condition : *condlist)
  {
    size_t linelength = 66;
    std::string sectionname = condition->SectionName();

    // ignore conditions occurring zero times
    count = count_cond.find(sectionname);
    if (count == count_cond.end()) continue;

    // only write conditions provided by user
    std::string dash(linelength - sectionname.size(), '-');
    dat << dash << sectionname << std::endl;
    std::string geo;
    switch (condition->GeometryType())
    {
      case DRT::Condition::Point:
        geo = "DPOINT ";
        break;
      case DRT::Condition::Line:
        geo = "DLINE  ";
        break;
      case DRT::Condition::Surface:
        geo = "DSURF  ";
        break;
      case DRT::Condition::Volume:
        geo = "DVOL   ";
        break;
      default:
        dserror("geometry type unspecified");
    }

    dat << geo << (count->second).size() << std::endl;
    for (i_c = (count->second).begin(); i_c != (count->second).end(); ++i_c)
    {
      EXODUS::cond_def actcon = condefs[*i_c];
      std::string name;
      std::string pname;
      if (actcon.me == EXODUS::bcns)
      {
        name = (mymesh.GetNodeSet(actcon.id).GetName());
        pname = (mymesh.GetNodeSet(actcon.id).GetPropName());
      }
      else if (actcon.me == EXODUS::bceb)
      {
        name = (mymesh.GetElementBlock(actcon.id)->GetName());
      }
      else if (actcon.me == EXODUS::bcss)
      {
        name = (mymesh.GetSideSet(actcon.id).GetName());
      }
      else
        dserror("Unidentified Actcon");
      if ((name != ""))
      {
        dat << "// " << name;
        if (pname != "none")
        {
          dat << " " << pname;
        }
        dat << std::endl;
      }
      else if (pname != "none")
        dat << "// " << pname << std::endl;

      // write the condition
      if ((actcon.sec == "DESIGN SURF LOCSYS CONDITIONS") && (actcon.me == EXODUS::bcns))
      {
        // special case for locsys conditions: calculate normal
        std::vector<double> normal_tangent = EXODUS::CalcNormalSurfLocsys(actcon.id, mymesh);
        dat << "E " << actcon.e_id << " - ";
        for (double normtang : normal_tangent)
          dat << std::setprecision(10) << std::fixed << normtang << " ";
        dat << std::endl;
      }
      else
        dat << "E " << actcon.e_id << " - " << actcon.desc << std::endl;
    }
    // remove sectionname from map, since writing is done
    count_cond.erase(sectionname);
  }
  if (count_cond.size() > 0)  // there are conditions left that were not recognized!!
  {
    std::cout << std::endl << std::endl;
    for (count = count_cond.begin(); count != count_cond.end(); ++count)
      std::cout << "Section name  " << count->first << "  is not valid. Typo?" << std::endl;

    dserror("There are invalid condition names in your bc file (see list above)");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::CalcNormalSurfLocsys(const int ns_id, const EXODUS::Mesh& m)
{
  std::vector<double> normaltangent;
  EXODUS::NodeSet ns = m.GetNodeSet(ns_id);

  std::set<int> nodes_from_nodeset = ns.GetNodeSet();
  std::set<int>::iterator it;

  // compute normal
  auto surfit = nodes_from_nodeset.begin();
  int origin = *surfit;  // get first set node
  ++surfit;
  int head1 = *surfit;  // get second set node
  ++surfit;

  std::set<int>::iterator thirdnode;

  // find third node such that a proper normal can be computed
  for (it = surfit; it != nodes_from_nodeset.end(); ++it)
  {
    thirdnode = it;
    normaltangent = EXODUS::Normal(head1, origin, *thirdnode, m);
    if (normaltangent.size() != 1) break;
  }
  if (normaltangent.size() == 1)
  {
    dserror(
        "Warning! No normal defined for SurfLocsys within nodeset '%s'!", (ns.GetName()).c_str());
  }

  // find tangent by Gram-Schmidt
  std::vector<double> t(3);
  t.at(0) = 1.0;
  t.at(1) = 0.0;
  t.at(2) = 0.0;  // try this one
  double sp = t[0] * normaltangent[0] + t[1] * normaltangent[1] +
              t[2] * normaltangent[2];  // scalar product
  // subtract projection
  t.at(0) -= normaltangent[0] * sp;
  t.at(1) -= normaltangent[1] * sp;
  t.at(2) -= normaltangent[2] * sp;

  // very unlucky case
  if (t.at(0) < 1.0E-14)
  {
    t.at(0) = 0.0;
    t.at(1) = 1.0;
    t.at(2) = 0.0;  // rather use this
  }

  normaltangent.push_back(t.at(0));
  normaltangent.push_back(t.at(1));
  normaltangent.push_back(t.at(2));

  return normaltangent;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatDesignTopology(
    const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  // sort baciconds w.r.t. underlying topology
  std::map<int, EXODUS::cond_def> dpoints;
  std::map<int, EXODUS::cond_def> dlines;
  std::map<int, EXODUS::cond_def> dsurfs;
  std::map<int, EXODUS::cond_def> dvols;

  for (const auto& conditiondefinition : condefs)
  {
    switch (conditiondefinition.gtype)
    {
      case DRT::Condition::Point:
        dpoints.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case DRT::Condition::Line:
        dlines.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case DRT::Condition::Surface:
        dsurfs.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case DRT::Condition::Volume:
        dvols.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case DRT::Condition::NoGeom:
        // do nothing
        break;
      default:
        dserror("Cannot identify Condition GeometryType");
        break;
    }
  }

  dat << "-----------------------------------------------DNODE-NODE TOPOLOGY" << std::endl;
  for (const auto& dpoint : dpoints)
  {
    const auto nodes = EXODUS::GetNsFromBCEntity(dpoint.second, mymesh);
    for (auto node : nodes)
    {
      dat << "NODE    " << node << " "
          << "DNODE " << dpoint.second.e_id << std::endl;
    }
  }
  dat << "-----------------------------------------------DLINE-NODE TOPOLOGY" << std::endl;
  for (const auto& dline : dlines)
  {
    const auto nodes = EXODUS::GetNsFromBCEntity(dline.second, mymesh);
    for (auto node : nodes)
    {
      dat << "NODE    " << node << " "
          << "DLINE " << dline.second.e_id << std::endl;
    }
  }
  dat << "-----------------------------------------------DSURF-NODE TOPOLOGY" << std::endl;
  for (const auto& dsurf : dsurfs)
  {
    const auto nodes = EXODUS::GetNsFromBCEntity(dsurf.second, mymesh);
    for (auto node : nodes)
    {
      dat << "NODE    " << node << " "
          << "DSURFACE " << dsurf.second.e_id << std::endl;
    }
  }
  dat << "------------------------------------------------DVOL-NODE TOPOLOGY" << std::endl;
  for (const auto& dvol : dvols)
  {
    const auto nodes = EXODUS::GetNsFromBCEntity(dvol.second, mymesh);
    for (auto node : nodes)
    {
      dat << "NODE    " << node << " "
          << "DVOL " << dvol.second.e_id << std::endl;
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const std::set<int> EXODUS::GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m)
{
  if (e.me == EXODUS::bcns)
  {
    EXODUS::NodeSet ns = m.GetNodeSet(e.id);
    return ns.GetNodeSet();
  }
  else if (e.me == EXODUS::bceb)
  {
    std::set<int> allnodes;
    Teuchos::RCP<EXODUS::ElementBlock> eb = m.GetElementBlock(e.id);
    Teuchos::RCP<const std::map<int, std::vector<int>>> eles = eb->GetEleConn();
    for (const auto& ele : *eles)
    {
      const std::vector<int> nodes = ele.second;
      for (auto node : nodes) allnodes.insert(node);
    }
    return allnodes;
  }
  else if (e.me == EXODUS::bcss)
  {
    std::set<int> allnodes;
    EXODUS::SideSet ss = m.GetSideSet(e.id);
    const std::map<int, std::vector<int>> eles = ss.GetSideSet();
    for (const auto& ele : eles)
    {
      const std::vector<int> nodes = ele.second;
      for (auto node : nodes) allnodes.insert(node);
    }
    return allnodes;
  }
  else
    dserror("Cannot identify mesh_entity");
  std::set<int> n;
  return n;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatNodes(const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  dat << "-------------------------------------------------------NODE COORDS" << std::endl;
  dat.precision(16);
  Teuchos::RCP<std::map<int, std::vector<double>>> nodes = mymesh.GetNodes();

  for (const auto& node : *nodes)
  {
    std::vector<double> coords = node.second;
    dat << "NODE " << std::setw(9) << node.first << " COORD";
    for (double coord : coords) dat << " " << std::setw(23) << std::scientific << coord;
    dat << std::endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::WriteDatEles(const std::vector<elem_def>& eledefs, const EXODUS::Mesh& mymesh,
    std::ostream& dat,
    const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo)
{
  // sort elements w.r.t. structure, fluid, ale, scalar transport, thermo, etc.
  std::vector<EXODUS::elem_def> strus;
  std::vector<EXODUS::elem_def> fluids;
  std::vector<EXODUS::elem_def> ales;
  std::vector<EXODUS::elem_def> lubrication;
  std::vector<EXODUS::elem_def> transport;
  std::vector<EXODUS::elem_def> transport2;
  std::vector<EXODUS::elem_def> thermo;
  std::vector<EXODUS::elem_def> cell;
  std::vector<EXODUS::elem_def> cellscatra;
  std::vector<EXODUS::elem_def> elemag;
  std::vector<EXODUS::elem_def> artery;
  std::vector<EXODUS::elem_def>::const_iterator i_et;

  for (i_et = eledefs.begin(); i_et != eledefs.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    if (acte.sec == "STRUCTURE")
      strus.push_back(acte);
    else if (acte.sec == "FLUID")
      fluids.push_back(acte);
    else if (acte.sec == "ALE")
      ales.push_back(acte);
    else if (acte.sec == "LUBRICATION")
      lubrication.push_back(acte);
    else if (acte.sec == "TRANSPORT")
      transport.push_back(acte);
    else if (acte.sec == "TRANSPORT2")
      transport2.push_back(acte);
    else if (acte.sec == "THERMO")
      thermo.push_back(acte);
    else if (acte.sec == "CELL")
      cell.push_back(acte);
    else if (acte.sec == "CELLSCATRA")
      cellscatra.push_back(acte);
    else if (acte.sec == "ELECTROMAGNETIC")
      elemag.push_back(acte);
    else if (acte.sec == "ARTERY")
      artery.push_back(acte);
    else if (acte.sec == "")
      ;
    else
    {
      std::cout << "Unknown ELEMENT sectionname in eb" << acte.id << ": '" << acte.sec << "'!"
                << std::endl;
      dserror("Unknown ELEMENT sectionname");
    }
  }

  int startele =
      1;  // BACI-Dat eles start with 1, this int is adapted for more than one element section

  // print structure elements
  dat << "------------------------------------------------STRUCTURE ELEMENTS" << std::endl;
  for (i_et = strus.begin(); i_et != strus.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print fluid elements
  dat << "----------------------------------------------------FLUID ELEMENTS" << std::endl;
  for (i_et = fluids.begin(); i_et != fluids.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print ale elements
  dat << "------------------------------------------------------ALE ELEMENTS" << std::endl;
  for (i_et = ales.begin(); i_et != ales.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print Lubrication elements
  dat << "------------------------------------------------LUBRICATION ELEMENTS" << std::endl;
  for (i_et = lubrication.begin(); i_et != lubrication.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print transport elements
  dat << "------------------------------------------------TRANSPORT ELEMENTS" << std::endl;
  for (i_et = transport.begin(); i_et != transport.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print transport2 elements
  dat << "------------------------------------------------TRANSPORT2 ELEMENTS" << std::endl;
  for (i_et = transport2.begin(); i_et != transport2.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print thermo elements
  dat << "---------------------------------------------------THERMO ELEMENTS" << std::endl;
  for (i_et = thermo.begin(); i_et != thermo.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print cell elements
  dat << "---------------------------------------------------CELL ELEMENTS" << std::endl;
  for (i_et = cell.begin(); i_et != cell.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print cellscatra elements
  dat << "---------------------------------------------------CELLSCATRA ELEMENTS" << std::endl;
  for (i_et = cellscatra.begin(); i_et != cellscatra.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print electromagnetic elements
  dat << "---------------------------------------------------ELECTROMAGNETIC ELEMENTS" << std::endl;
  for (i_et = elemag.begin(); i_et != elemag.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }

  // print artery elements
  dat << "-------------------------------------------------------ARTERY ELEMENTS" << std::endl;
  for (i_et = artery.begin(); i_et != artery.end(); ++i_et)
  {
    EXODUS::elem_def acte = *i_et;
    Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(acte.id);
    EXODUS::DatEles(eb, acte, startele, dat, elecenterlineinfo, acte.id);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::DatEles(Teuchos::RCP<const EXODUS::ElementBlock> eb, const EXODUS::elem_def& acte,
    int& startele, std::ostream& datfile,
    const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elescli, const int eb_id)
{
  auto eles = eb->GetEleConn();
  for (const auto& ele : *eles)
  {
    std::stringstream dat;  // first build up the std::string for actual element line
    const std::vector<int> nodes = ele.second;
    std::vector<int>::const_iterator i_n;
    dat << "   " << startele;
    dat << " " << acte.ename;  // e.g. "SOLIDH8"
    dat << " " << DistypeToString(PreShapeToDrt(eb->GetShape()));
    dat << "  ";
    for (auto node : nodes) dat << node << " ";
    dat << "   " << acte.desc;  // e.g. "MAT 1"
    if (elescli.size() != 0)
    {
      // quick check wether elements in ele Block have a fiber direction
      if (elescli.find(eb_id) != elescli.end())
      {
        // write local cosy from centerline to each element
        std::vector<std::vector<double>> ecli =
            (elescli.find(eb_id)->second).find(ele.first)->second;
        dat << " RAD " << std::fixed << std::setprecision(8) << ecli[0][0] << " " << ecli[0][1]
            << " " << ecli[0][2];
        dat << " AXI " << ecli[1][0] << " " << ecli[1][1] << " " << ecli[1][2];
        dat << " CIR " << ecli[2][0] << " " << ecli[2][1] << " " << ecli[2][2];
      }
    }
    dat << std::endl;  // finish this element line

    startele++;
    datfile << dat.str();  // only one access to the outfile (saves system time)
  }
}
