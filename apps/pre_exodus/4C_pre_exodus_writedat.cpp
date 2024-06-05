/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus .dat-file writer

\level 1


Here is everything related with writing a dat-file
 */
/*----------------------------------------------------------------------*/
#include "4C_pre_exodus_writedat.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_inpar_validconditions.hpp"
#include "4C_pre_exodus_reader.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::WriteDatFile(const std::string& datfile, const EXODUS::Mesh& mymesh,
    const std::string& headfile, const std::vector<EXODUS::ElemDef>& eledefs,
    const std::vector<EXODUS::CondDef>& condefs)
{
  // open datfile
  std::ofstream dat(datfile.c_str());
  if (!dat) FOUR_C_THROW("failed to open file: %s", datfile.c_str());

  // write dat-file intro
  EXODUS::WriteDatIntro(headfile, mymesh, dat);

  // write "header"
  EXODUS::WriteDatHead(headfile, dat);

  // write conditions
  EXODUS::WriteDatConditions(condefs, mymesh, dat);

  // write design-topology
  EXODUS::WriteDatDesignTopology(condefs, mymesh, dat);

  // write nodal coordinates
  EXODUS::WriteDatNodes(mymesh, dat);

  // write elements
  EXODUS::WriteDatEles(eledefs, mymesh, dat);

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
         "                   General Data File 4C\n"
         "==================================================================\n"
         "-------------------------------------------------------------TITLE\n"
         "created by pre_exodus\n"
         "------------------------------------------------------PROBLEM SIZE\n";
  // print number of elements and nodes just as an comment instead of
  // a valid parameter (prevents possible misuse of these parameters in 4C)
  dat << "//ELEMENTS    " << mymesh.GetNumEle() << std::endl;
  dat << "//NODES       " << mymesh.GetNumNodes() << std::endl;
  // parameter for the number of spatial dimensions
  dat << "DIM           " << mymesh.GetFourCDim() << std::endl;
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
    FOUR_C_THROW("Unable to open head-file");
  }
  while (header.good()) head << (char)header.get();
  // while (!header.eof()) head << (char) header.get();
  header.close();
  std::string headstring = head.str();

  // delete sections which has been written by WriteDatIntro already
  RemoveDatSection("PROBLEM SIZE", headstring);

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
void EXODUS::WriteDatConditions(
    const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  using namespace FourC;

  Teuchos::RCP<std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>> condlist =
      INPUT::ValidConditions();

  // count how often we have one specific condition
  std::map<std::string, std::vector<int>> count_cond;
  std::map<std::string, std::vector<int>>::const_iterator count;
  std::vector<int>::const_iterator i_c;
  for (int i_cond = 0; i_cond < static_cast<int>(condefs.size()); ++i_cond)
    (count_cond[condefs.at(i_cond).sec]).push_back(i_cond);

  // loop all valid conditions that 4C knows
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
      case CORE::Conditions::geometry_type_point:
        geo = "DPOINT ";
        break;
      case CORE::Conditions::geometry_type_line:
        geo = "DLINE  ";
        break;
      case CORE::Conditions::geometry_type_surface:
        geo = "DSURF  ";
        break;
      case CORE::Conditions::geometry_type_volume:
        geo = "DVOL   ";
        break;
      default:
        FOUR_C_THROW("geometry type unspecified");
    }

    dat << geo << (count->second).size() << std::endl;
    for (i_c = (count->second).begin(); i_c != (count->second).end(); ++i_c)
    {
      EXODUS::CondDef actcon = condefs[*i_c];
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
        FOUR_C_THROW("Unidentified Actcon");
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
      if (actcon.desc == "" and actcon.sec == "DESIGN SURF LOCSYS CONDITIONS" and
          actcon.me == EXODUS::bcns)
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

    FOUR_C_THROW("There are invalid condition names in your bc file (see list above)");
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

  const auto compute_normal = [](int head1, int origin, int head2, const EXODUS::Mesh& basemesh)
  {
    std::vector<double> normal(3);
    std::vector<double> h1 = basemesh.GetNode(head1);
    std::vector<double> h2 = basemesh.GetNode(head2);
    std::vector<double> o = basemesh.GetNode(origin);

    normal[0] = ((h1[1] - o[1]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[1] - o[1]));
    normal[1] = -((h1[0] - o[0]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[0] - o[0]));
    normal[2] = ((h1[0] - o[0]) * (h2[1] - o[1]) - (h1[1] - o[1]) * (h2[0] - o[0]));

    double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    const double epsilon = 1.E-4;
    if (length > epsilon)
    {
      normal[0] = normal[0] / length;
      normal[1] = normal[1] / length;
      normal[2] = normal[2] / length;
    }
    else
    {  // normal is undefined, vectors seem collinear
      normal.resize(1);
      normal[0] = 0.0;
    }

    return normal;
  };

  // find third node such that a proper normal can be computed
  for (it = surfit; it != nodes_from_nodeset.end(); ++it)
  {
    thirdnode = it;
    normaltangent = compute_normal(head1, origin, *thirdnode, m);
    if (normaltangent.size() != 1) break;
  }
  if (normaltangent.size() == 1)
  {
    FOUR_C_THROW(
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
    const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  using namespace FourC;

  // sort 4C conditions w.r.t. underlying topology
  std::map<int, EXODUS::CondDef> dpoints;
  std::map<int, EXODUS::CondDef> dlines;
  std::map<int, EXODUS::CondDef> dsurfs;
  std::map<int, EXODUS::CondDef> dvols;

  for (const auto& conditiondefinition : condefs)
  {
    switch (conditiondefinition.gtype)
    {
      case CORE::Conditions::geometry_type_point:
        dpoints.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case CORE::Conditions::geometry_type_line:
        dlines.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case CORE::Conditions::geometry_type_surface:
        dsurfs.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case CORE::Conditions::geometry_type_volume:
        dvols.insert(std::make_pair(conditiondefinition.e_id, conditiondefinition));
        break;
      case CORE::Conditions::geometry_type_no_geom:
        // do nothing
        break;
      default:
        FOUR_C_THROW("Cannot identify Condition GeometryType");
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
std::set<int> EXODUS::GetNsFromBCEntity(const EXODUS::CondDef& e, const EXODUS::Mesh& m)
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
    FOUR_C_THROW("Cannot identify mesh_entity");
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
void EXODUS::WriteDatEles(
    const std::vector<ElemDef>& eledefs, const EXODUS::Mesh& mymesh, std::ostream& dat)
{
  // sort elements w.r.t. structure, fluid, ale, scalar transport, thermo, etc.
  std::vector<EXODUS::ElemDef> structure_elements;
  std::vector<EXODUS::ElemDef> fluid_elements;
  std::vector<EXODUS::ElemDef> ale_elements;
  std::vector<EXODUS::ElemDef> lubrication_elements;
  std::vector<EXODUS::ElemDef> transport_elements;
  std::vector<EXODUS::ElemDef> transport2_elements;
  std::vector<EXODUS::ElemDef> thermo_elements;
  std::vector<EXODUS::ElemDef> cell_elements;
  std::vector<EXODUS::ElemDef> cellscatra_elements;
  std::vector<EXODUS::ElemDef> elemag_elements;
  std::vector<EXODUS::ElemDef> artery_elements;

  for (const auto& element_definition : eledefs)
  {
    if (element_definition.sec == "STRUCTURE")
      structure_elements.push_back(element_definition);
    else if (element_definition.sec == "FLUID")
      fluid_elements.push_back(element_definition);
    else if (element_definition.sec == "ALE")
      ale_elements.push_back(element_definition);
    else if (element_definition.sec == "LUBRICATION")
      lubrication_elements.push_back(element_definition);
    else if (element_definition.sec == "TRANSPORT")
      transport_elements.push_back(element_definition);
    else if (element_definition.sec == "TRANSPORT2")
      transport2_elements.push_back(element_definition);
    else if (element_definition.sec == "THERMO")
      thermo_elements.push_back(element_definition);
    else if (element_definition.sec == "CELL")
      cell_elements.push_back(element_definition);
    else if (element_definition.sec == "CELLSCATRA")
      cellscatra_elements.push_back(element_definition);
    else if (element_definition.sec == "ELECTROMAGNETIC")
      elemag_elements.push_back(element_definition);
    else if (element_definition.sec == "ARTERY")
      artery_elements.push_back(element_definition);
    else if (element_definition.sec == "")
      ;
    else
    {
      std::cout << "Unknown ELEMENT sectionname in eb" << element_definition.id << ": '"
                << element_definition.sec << "'!" << std::endl;
      FOUR_C_THROW("Unknown ELEMENT sectionname");
    }
  }

  // element ids in 4C dat files start with 1, this int is adapted for more than one element section
  int startele = 1;

  const auto printElementSection =
      [&](const std::vector<ElemDef>& ele_vector, const std::string& section_name)
  {
    const unsigned padding_length = 66;
    // we need at least 2 dashes at the beginning to be recognizable to the dat file reader
    const unsigned min_num_preceding_dashes = 2;

    if (section_name.length() > padding_length - min_num_preceding_dashes)
      FOUR_C_THROW("The section name you chose exceeds padding length");

    std::string padded_section_name(section_name);
    padded_section_name.insert(
        padded_section_name.begin(), padding_length - padded_section_name.length(), '-');

    dat << padded_section_name << std::endl;

    for (const auto& ele : ele_vector)
    {
      Teuchos::RCP<EXODUS::ElementBlock> eb = mymesh.GetElementBlock(ele.id);
      EXODUS::DatEles(eb, ele, startele, dat, ele.id);
    }
  };

  // print structure elements
  if (!structure_elements.empty()) printElementSection(structure_elements, "STRUCTURE ELEMENTS");

  // print fluid elements
  if (!fluid_elements.empty()) printElementSection(fluid_elements, "FLUID ELEMENTS");

  // print ale elements
  if (!ale_elements.empty()) printElementSection(ale_elements, "ALE ELEMENTS");

  // print Lubrication elements
  if (!lubrication_elements.empty())
    printElementSection(lubrication_elements, "LUBRICATION ELEMENTS");

  // print transport elements
  if (!transport_elements.empty()) printElementSection(transport_elements, "TRANSPORT ELEMENTS");

  // print transport2 elements
  if (!transport2_elements.empty()) printElementSection(transport2_elements, "TRANSPORT2 ELEMENTS");

  // print thermo elements
  if (!thermo_elements.empty()) printElementSection(thermo_elements, "THERMO ELEMENTS");

  // print cell elements
  if (!cell_elements.empty()) printElementSection(cell_elements, "CELL ELEMENTS");

  // print cellscatra elements
  if (!cellscatra_elements.empty()) printElementSection(cellscatra_elements, "CELLSCATRA ELEMENTS");

  // print electromagnetic elements
  if (!elemag_elements.empty()) printElementSection(elemag_elements, "ELECTROMAGNETIC ELEMENTS");

  // print artery elements
  if (!artery_elements.empty()) printElementSection(artery_elements, "ARTERY ELEMENTS");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::DatEles(Teuchos::RCP<const EXODUS::ElementBlock> eb, const EXODUS::ElemDef& acte,
    int& startele, std::ostream& datfile, const int eb_id)
{
  auto eles = eb->GetEleConn();
  for (const auto& ele : *eles)
  {
    std::stringstream dat;  // first build up the std::string for actual element line
    const std::vector<int> nodes = ele.second;
    std::vector<int>::const_iterator i_n;
    dat << "   " << startele;
    dat << " " << acte.ename;  // e.g. "SOLIDH8"
    dat << " " << CORE::FE::CellTypeToString(PreShapeToDrt(eb->GetShape()));
    dat << "  ";
    for (auto node : nodes) dat << node << " ";
    dat << "   " << acte.desc;  // e.g. "MAT 1"
    dat << std::endl;           // finish this element line

    startele++;
    datfile << dat.str();  // only one access to the outfile (saves system time)
  }
}

FOUR_C_NAMESPACE_CLOSE
