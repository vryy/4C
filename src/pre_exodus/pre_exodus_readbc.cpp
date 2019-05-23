/*----------------------------------------------------------------------*/
/*!

\brief pre_exodus bc-file reader

\maintainer Martin Kronbichler

\level 1

Here is everything related with reading a bc file
*/
/*----------------------------------------------------------------------*/
#include "pre_exodus_readbc.H"
#include "pre_exodus_reader.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ReadBCFile(const std::string& bcfile, std::vector<EXODUS::elem_def>& eledefs,
    std::vector<EXODUS::cond_def>& condefs)
{
  // first we read the whole file into one stream/std::string
  std::stringstream bcstream;
  const char* bcfilechar = bcfile.c_str();

  std::ifstream bcfstream(bcfilechar, std::ifstream::in);
  if (!bcfstream.good())
  {
    std::cout << std::endl << "Unable to open file: " << bcfile << std::endl;
    dserror("Unable to open bc-file");
  }

  while (bcfstream.good()) bcstream << (char)bcfstream.get();

  bcfstream.close();

  // std::string which contains the whole file
  std::string allconds = bcstream.str();
  allconds.erase(allconds.end() - 1);  // delete last 'whatisthis'-char

  // get rid of first part
  size_t found;
  found = allconds.find("BCSPECS");
  if (found == std::string::npos)
    dserror("No specifications found in bcfile. BCSPECS section is missing!");
  allconds.erase(allconds.begin(), allconds.begin() + found);

  // get rid of 'validconditions' part, if it can be found
  found = allconds.find("VALIDCONDITIONS");
  if (found != std::string::npos) allconds.erase(allconds.begin() + found, allconds.end());

  // erase all comments
  size_t comment_start = allconds.find_first_of("//");
  size_t comment_end = 0;
  while (comment_start != std::string::npos)
  {
    if (allconds[comment_start - 1] != '\n')
    {
      comment_end = comment_start + 2;
      comment_start = allconds.find_first_of("//", comment_end);
      continue;
    }

    comment_end = allconds.find_first_of('\n', comment_start);
    allconds.erase(comment_start, comment_end - comment_start);
    comment_start = allconds.find_first_of("//", comment_start - 1);
  }

  // define markers
  const std::string ebmarker("*eb");
  const std::string nsmarker("*ns");
  const std::string ssmarker("*ss");
  const int markerlength = 3;
  const std::string marker("*");

  // necessary counters
  int E_id = 0;  // the 'E num -' in the datfile
  int ndp = 0;
  int ndl = 0;
  int nds = 0;
  int ndv = 0;

  // map to avoid double assignment
  std::map<int, int> eb_dp2Eid;  // ele block point E id
  std::map<int, int> eb_dl2Eid;  // ele block line E id
  std::map<int, int> eb_ds2Eid;  // ele block surf E id
  std::map<int, int> eb_dv2Eid;  // ele block vol E id
  std::map<int, int> ns_dp2Eid;  // node set point E id
  std::map<int, int> ns_dl2Eid;  // node set line E id
  std::map<int, int> ns_ds2Eid;  // node set surf E id
  std::map<int, int> ns_dv2Eid;  // node set vol E id

  if (allconds.find("**") != std::string::npos)
    dserror("String '**' detected. More * than one not allowed due to usage as a marker.");

  found = allconds.find_first_of(marker);
  while (found != std::string::npos)
  {
    int startpos = found;
    found = allconds.find(marker, found + 1);  // step forward to find next match
    // get actual condition
    std::string actcond = allconds.substr(startpos, found - startpos);

    // ensure substd::string has minimum length!
    if (actcond.size() < 3) dserror("Substd::string is too short");

    // find out what mesh_entity type we have
    std::string mesh_entity = actcond.substr(0, 3);

    // get its id
    size_t found2 = actcond.find_first_of("=");
    std::string buffer = actcond.substr(markerlength, found2 - markerlength);
    // convert std::string to int
    std::istringstream bufferstream(buffer);
    int id;
    bufferstream >> id;

    // condition type
    size_t left = actcond.find_first_of("\"", 0);
    size_t right = actcond.find_first_of("\"", left + 1);
    std::string type = actcond.substr(left + 1, right - left - 1);

    if (mesh_entity.compare(ebmarker) == 0)
    {
      // in case of eb we differntiate between 'element' or 'condition'
      if (type.compare("ELEMENT") == 0)
      {
        EXODUS::elem_def edef = EXODUS::ReadEdef(mesh_entity, id, actcond);
        eledefs.push_back(edef);
      }
      else if (type.compare("CONDITION") == 0)
      {
        // the geometry type is figured out by finding the identifier e.g. "SURF"
        // within the condition name (sectionname) which therefore should carry it!
        EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
        switch (cdef.gtype)
        {
          case DRT::Condition::Point:
          {
            if (eb_dp2Eid.find(id) != eb_dp2Eid.end())
              E_id = eb_dp2Eid.find(id)->second;
            else
            {
              ++ndp;
              E_id = ndp;
              eb_dp2Eid.insert(std::pair<int, int>(id, E_id));
            }
            break;
          }
          case DRT::Condition::Line:
          {
            if (eb_dl2Eid.find(id) != eb_dl2Eid.end())
              E_id = eb_dl2Eid.find(id)->second;
            else
            {
              ++ndl;
              E_id = ndl;
              eb_dl2Eid.insert(std::pair<int, int>(id, E_id));
            }
            break;
          }
          case DRT::Condition::Surface:
          {
            if (eb_ds2Eid.find(id) != eb_ds2Eid.end())
              E_id = eb_ds2Eid.find(id)->second;
            else
            {
              ++nds;
              E_id = nds;
              eb_ds2Eid.insert(std::pair<int, int>(id, E_id));
            }
            break;
          }
          case DRT::Condition::Volume:
          {
            if (eb_dv2Eid.find(id) != eb_dv2Eid.end())
              E_id = eb_dv2Eid.find(id)->second;
            else
            {
              ++ndv;
              E_id = ndv;
              eb_dv2Eid.insert(std::pair<int, int>(id, E_id));
            }
            break;
          }
          case DRT::Condition::NoGeom:
            E_id = 0;
            break;
          default:
            dserror("geometry type unspecified");
            break;
        }
        cdef.e_id = E_id;
        condefs.push_back(cdef);
      }
      else
      {
        std::cout << "Undefined type for eb" << id << ": " << type << std::endl;
        dserror("Undefined type!");
      }
    }
    else if (mesh_entity.compare(nsmarker) == 0)
    {
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
      switch (cdef.gtype)
      {
        case DRT::Condition::Point:
        {
          if (ns_dp2Eid.find(id) != ns_dp2Eid.end())
            E_id = ns_dp2Eid.find(id)->second;
          else
          {
            ++ndp;
            E_id = ndp;
            ns_dp2Eid.insert(std::pair<int, int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Line:
        {
          if (ns_dl2Eid.find(id) != ns_dl2Eid.end())
            E_id = ns_dl2Eid.find(id)->second;
          else
          {
            ++ndl;
            E_id = ndl;
            ns_dl2Eid.insert(std::pair<int, int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Surface:
        {
          if (ns_ds2Eid.find(id) != ns_ds2Eid.end())
            E_id = ns_ds2Eid.find(id)->second;
          else
          {
            ++nds;
            E_id = nds;
            ns_ds2Eid.insert(std::pair<int, int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::Volume:
        {
          if (ns_dv2Eid.find(id) != ns_dv2Eid.end())
            E_id = ns_dv2Eid.find(id)->second;
          else
          {
            ++ndv;
            E_id = ndv;
            ns_dv2Eid.insert(std::pair<int, int>(id, E_id));
          }
          break;
        }
        case DRT::Condition::NoGeom:
          E_id = 0;
          break;
        default:
          dserror("geometry type unspecified");
          break;
      }
      cdef.e_id = E_id;
      condefs.push_back(cdef);
    }
    else if (mesh_entity.compare(ssmarker) == 0)
    {
      EXODUS::cond_def cdef = EXODUS::ReadCdef(mesh_entity, id, actcond);
      condefs.push_back(cdef);
    }
    else
      dserror(
          "Cannot identify marker '%s'. Use *el (element block), *ns (nodeset) or *ss (sideset)",
          mesh_entity.c_str());
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::elem_def EXODUS::ReadEdef(
    const std::string& mesh_entity, const int id, const std::string& actcond)
{
  EXODUS::elem_def edef;
  edef.id = id;
  edef.me = EXODUS::bceb;

  // read sectionname
  size_t left = actcond.find("sectionname=\"");  // 13 chars
  size_t right = actcond.find_first_of("\"", left + 13);
  edef.sec = actcond.substr(left + 13, right - (left + 13));

  // read description
  left = actcond.find("description=\"");  // 13 chars
  right = actcond.find_first_of("\"", left + 13);
  edef.desc = actcond.substr(left + 13, right - (left + 13));

  // read ename
  left = actcond.find("elementname=\"");  // 13 chars
  right = actcond.find_first_of("\"", left + 13);
  edef.ename = actcond.substr(left + 13, right - (left + 13));

  return edef;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::cond_def EXODUS::ReadCdef(
    const std::string& mesh_entity, const int id, const std::string& actcond)
{
  EXODUS::cond_def cdef;
  cdef.id = id;
  if (mesh_entity.compare(1, 2, "eb") == 0)
    cdef.me = EXODUS::bceb;
  else if (mesh_entity.compare(1, 2, "ns") == 0)
    cdef.me = EXODUS::bcns;
  else if (mesh_entity.compare(1, 2, "ss") == 0)
    cdef.me = EXODUS::bcss;
  else
    dserror("Cannot identify marker. Use *el (element block), *ns (nodeset) or *ss (sideset)");

  // read sectionname
  size_t left = actcond.find("sectionname=\"");  // 13 chars
  size_t right = actcond.find_first_of("\"", left + 13);
  std::string secname = actcond.substr(left + 13, right - (left + 13));
  cdef.sec = secname;

  // read description
  left = actcond.find("description=\"");  // 13 chars
  right = actcond.find_first_of("\"", left + 13);
  std::string description = actcond.substr(left + 13, right - (left + 13));
  cdef.desc = description;

  // figure out geometry type
  cdef.gtype = DRT::Condition::NoGeom;  // default
  size_t found = secname.find("POINT");
  if (found != std::string::npos)
  {
    cdef.gtype = DRT::Condition::Point;
    return cdef;
  }
  found = secname.find("LINE");
  if (found != std::string::npos)
  {
    cdef.gtype = DRT::Condition::Line;
    return cdef;
  }
  found = secname.find("SURF");
  if (found != std::string::npos)
  {
    cdef.gtype = DRT::Condition::Surface;
    return cdef;
  }
  found = secname.find("VOL");
  if (found != std::string::npos) cdef.gtype = DRT::Condition::Volume;

  return cdef;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintBCDef(std::ostream& os, const EXODUS::elem_def& def)
{
  std::string mesh_entity;
  if (def.me == EXODUS::bceb)
    mesh_entity = "ElementBlock";
  else if (def.me == EXODUS::bcns)
    mesh_entity = "NodeSet";
  else if (def.me == EXODUS::bcss)
    mesh_entity = "SideSet";
  os << "The ELEMENT definition " << def.id << " refers to a " << mesh_entity << std::endl;
  os << "Sectionname: " << def.sec << std::endl;
  os << "Description: " << def.desc << std::endl;
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintBCDef(std::ostream& os, const EXODUS::cond_def& def)
{
  std::string mesh_entity;
  if (def.me == EXODUS::bceb)
    mesh_entity = "ElementBlock";
  else if (def.me == EXODUS::bcns)
    mesh_entity = "NodeSet";
  else if (def.me == EXODUS::bcss)
    mesh_entity = "SideSet";
  os << "The CONDITION definition " << def.id << " refers to a " << mesh_entity << std::endl;
  os << "Sectionname: " << def.sec << std::endl;
  os << "Description: " << def.desc << std::endl;
  os << std::endl;
}


/*----------------------------------------------------------------------*
 * check if periodic boundary conditions are defined        u.may 02/10 *
 *----------------------------------------------------------------------*/
bool EXODUS::PeriodicBoundaryConditionsFound(std::vector<EXODUS::cond_def> condefs)
{
  bool pbc_defined = false;
  for (unsigned int i = 0; i < condefs.size(); i++)
  {
    // check if pbc found
    if (condefs[i].sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
        (condefs[i].desc.find("ANGLE 0.0") != std::string::npos ||
            condefs[i].desc.find("ANGLE 0") != std::string::npos))
    {
      pbc_defined = true;
      break;
    }
  }
  return pbc_defined;
}



/*----------------------------------------------------------------------*
 * correct nodal coordinates for                                        *
 * periodic boundary conditions are defined                 u.may 02/10 *
 *----------------------------------------------------------------------*/
void EXODUS::CorrectNodalCoordinatesForPeriodicBoundaryConditions(
    EXODUS::Mesh& mesh, std::vector<EXODUS::cond_def> condefs)
{
  CorrectYZPlaneForPeriodicBoundaryConditions(mesh, condefs);
  CorrectXZPlaneForPeriodicBoundaryConditions(mesh, condefs);
  CorrectXYPlaneForPeriodicBoundaryConditions(mesh, condefs);
  return;
}



/*----------------------------------------------------------------------*
 * correct nodal coordinates in yz plane for                            *
 * periodic boundary conditions are defined                 u.may 02/10 *
 *----------------------------------------------------------------------*/
void EXODUS::CorrectYZPlaneForPeriodicBoundaryConditions(
    EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs)
{
  // loop over all conditions
  for (unsigned int i = 0; i < condefs.size(); i++)
  {
    // yz plane : y and z coordinates of two matching nodes have to be equal
    // yz plane:  x coordinates differ by a constant length
    if (condefs[i].sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
        condefs[i].desc.find("Master PLANE yz") != std::string::npos &&
        (condefs[i].desc.find("ANGLE 0.0") != std::string::npos ||
            condefs[i].desc.find("ANGLE 0") != std::string::npos))
    {
      // store master pbc
      const EXODUS::cond_def master_con = condefs[i];
      const int master_nodeset_id = master_con.id;

      // get condition id for master pbc
      size_t end_m = master_con.desc.find("Master");  // master
      std::string string_mconditionid = master_con.desc.substr(0, end_m - 1);
      // convert std::string to int
      std::istringstream string_mconditionidstream(string_mconditionid);
      int master_conditionID = -1;
      string_mconditionidstream >> master_conditionID;

      // find matching slave pbc node set
      EXODUS::cond_def slave_con;
      int slave_conditionID = -1;
      bool matchingIDfound = false;
      for (unsigned int i_slave = 0; i_slave < condefs.size(); i_slave++)
      {
        slave_con = condefs[i_slave];
        if (slave_con.sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
            slave_con.desc.find("Slave PLANE yz") != std::string::npos &&
            (slave_con.desc.find("ANGLE 0.0") != std::string::npos ||
                slave_con.desc.find("ANGLE 0") != std::string::npos))
        {
          size_t end_s = slave_con.desc.find("Slave");                            // slave
          std::string string_sconditionid = slave_con.desc.substr(0, end_s - 1);  // one whitespaces
          // convert std::string to int
          std::istringstream string_sconditionidstream(string_sconditionid);
          string_sconditionidstream >> slave_conditionID;

          if (slave_conditionID == master_conditionID) matchingIDfound = true;
        }
        if (matchingIDfound) break;
      }
      if (matchingIDfound == false)  // this indicates that ANGLE != 0.0 for slave plane
        continue;                    // in this case do not proceed, but leave here!!

      // read tolerance
      size_t start_tol = slave_con.desc.find("ABSTREETOL");
      size_t tol_length = std::string("ABSTREETOL").length();
      start_tol = start_tol + tol_length + 1;
      std::string string_tol =
          slave_con.desc.substr(start_tol, slave_con.desc.length() - start_tol);  // two whitespaces
      // convert std::string to int
      std::istringstream string_tolstream(string_tol);
      double abstol = -1.0;
      string_tolstream >> abstol;

      const int slave_nodeset_id = slave_con.id;
      const EXODUS::NodeSet master_nodeset = mesh.GetNodeSet(master_nodeset_id);
      const EXODUS::NodeSet slave_nodeset = mesh.GetNodeSet(slave_nodeset_id);

      if (slave_nodeset.GetNumNodes() != master_nodeset.GetNumNodes())
      {
        std::cout << "yz num master nodes = " << master_nodeset.GetNumNodes() << std::endl;
        std::cout << "yz num slave nodes = " << slave_nodeset.GetNumNodes() << std::endl;
        dserror("num master nodes != num slave nodes before adjusting coords");
      }

      const std::set<int> master_nodeset_ids = master_nodeset.GetNodeSet();
      const std::set<int> slave_nodeset_ids = slave_nodeset.GetNodeSet();

      // reference values for x on the master side and x on the slave side
      const double x_master = (mesh.GetNode(*master_nodeset_ids.begin()))[0];
      const double x_slave = (mesh.GetNode(*slave_nodeset_ids.begin()))[0];

      // loop over all master nodes
      for (std::set<int>::const_iterator m_node_id = master_nodeset_ids.begin();
           m_node_id != master_nodeset_ids.end(); m_node_id++)
      {
        // get master node coordinates
        std::vector<double> m_coord = mesh.GetNode(*m_node_id);
        m_coord[0] = x_master;
        mesh.SetNode(*m_node_id, m_coord);

        int count_slavenodes = 0;
        // loop over all slave nodes and find matching node
        for (std::set<int>::const_iterator s_node_id = slave_nodeset_ids.begin();
             s_node_id != slave_nodeset_ids.end(); s_node_id++)
        {
          count_slavenodes++;
          // get master node coordinates
          std::vector<double> s_coord = mesh.GetNode(*s_node_id);

          // compare coordinates
          bool equal = false;
          if ((fabs(m_coord[1] - s_coord[1]) < abstol) && (fabs(m_coord[2] - s_coord[2]) < abstol))
            equal = true;

          if ((count_slavenodes == (int)slave_nodeset_ids.size()) && !equal)
            dserror("no matching slave node %d found, adjust your tolerance", *m_node_id);

          if (!equal)
            continue;
          else
          {
            // if equal copy y and z
            s_coord[0] = x_slave;
            s_coord[1] = m_coord[1];
            s_coord[2] = m_coord[2];
            mesh.SetNode(*s_node_id, s_coord);
            break;
          }
        }  // loop all slave nodes
      }    // loop all master nodes
      if (mesh.GetNodeSet(master_nodeset_id).GetNumNodes() !=
          mesh.GetNodeSet(slave_nodeset_id).GetNumNodes())
      {
        std::cout << "yz num master nodes = " << mesh.GetNodeSet(master_nodeset_id).GetNumNodes()
                  << std::endl;
        std::cout << "yz num slave nodes = " << mesh.GetNodeSet(slave_nodeset_id).GetNumNodes()
                  << std::endl;
        dserror("num master nodes != num slave nodes after adjusting coords");
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 * correct nodal coordinates in xz plane for                            *
 * periodic boundary conditions are defined                 u.may 02/10 *
 *----------------------------------------------------------------------*/
void EXODUS::CorrectXZPlaneForPeriodicBoundaryConditions(
    EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs)
{
  // loop over all conditions
  for (unsigned int i = 0; i < condefs.size(); i++)
  {
    // xz plane : x and z coordinates of two matching nodes have to be equal
    // xz plane:  y coordinates differ by a constant length
    if (condefs[i].sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
        condefs[i].desc.find("Master PLANE xz") != std::string::npos &&
        (condefs[i].desc.find("ANGLE 0.0") != std::string::npos ||
            condefs[i].desc.find("ANGLE 0") != std::string::npos))
    {
      // store master pbc
      const EXODUS::cond_def master_con = condefs[i];
      const int master_nodeset_id = master_con.id;

      // get condition id for master pbc
      size_t end_m = master_con.desc.find("Master");  // master
      std::string string_mconditionid = master_con.desc.substr(0, end_m - 1);
      // convert std::string to int
      std::istringstream string_mconditionidstream(string_mconditionid);
      int master_conditionID = -1;
      string_mconditionidstream >> master_conditionID;

      // find matching slave pbc node set
      EXODUS::cond_def slave_con;
      int slave_conditionID = -1;
      bool matchingIDfound = false;

      for (unsigned int i_slave = 0; i_slave < condefs.size(); i_slave++)
      {
        slave_con = condefs[i_slave];
        if (slave_con.sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
            slave_con.desc.find("Slave PLANE xz") != std::string::npos &&
            (slave_con.desc.find("ANGLE 0.0") != std::string::npos ||
                slave_con.desc.find("ANGLE 0") != std::string::npos))
        {
          size_t end_s = slave_con.desc.find("Slave");                            // slave
          std::string string_sconditionid = slave_con.desc.substr(0, end_s - 1);  // one whitespaces
          // convert std::string to int
          std::istringstream string_sconditionidstream(string_sconditionid);
          string_sconditionidstream >> slave_conditionID;

          if (slave_conditionID == master_conditionID) matchingIDfound = true;
        }
        if (matchingIDfound) break;
      }
      if (matchingIDfound == false)  // this indicates that ANGLE != 0.0 for slave plane
        continue;                    // in this case do not proceed, but leave here!!

      // read tolerance
      size_t start_tol = slave_con.desc.find("ABSTREETOL");
      size_t tol_length = std::string("ABSTREETOL").length();
      start_tol = start_tol + tol_length + 1;
      std::string string_tol =
          slave_con.desc.substr(start_tol, slave_con.desc.length() - start_tol);  // two whitespaces
      // convert std::string to int
      std::istringstream string_tolstream(string_tol);
      double abstol = -1.0;
      string_tolstream >> abstol;

      const int slave_nodeset_id = slave_con.id;
      const EXODUS::NodeSet master_nodeset = mesh.GetNodeSet(master_nodeset_id);
      const EXODUS::NodeSet slave_nodeset = mesh.GetNodeSet(slave_nodeset_id);

      if (slave_nodeset.GetNumNodes() != master_nodeset.GetNumNodes())
      {
        std::cout << "xz num master nodes = " << master_nodeset.GetNumNodes() << std::endl;
        std::cout << "xz num slave nodes = " << slave_nodeset.GetNumNodes() << std::endl;
        dserror("xz num master nodes != num slave nodes before adjusting coords");
      }

      const std::set<int> master_nodeset_ids = master_nodeset.GetNodeSet();
      const std::set<int> slave_nodeset_ids = slave_nodeset.GetNodeSet();

      // reference values for x on the master side and x on the slave side
      const double y_master = (mesh.GetNode(*master_nodeset_ids.begin()))[1];  // get y-coord
      const double y_slave = (mesh.GetNode(*slave_nodeset_ids.begin()))[1];    // get y-coord

      // loop over all master nodes
      for (std::set<int>::const_iterator m_node_id = master_nodeset_ids.begin();
           m_node_id != master_nodeset_ids.end(); m_node_id++)
      {
        // get master node coordinates
        std::vector<double> m_coord = mesh.GetNode(*m_node_id);
        m_coord[1] = y_master;
        mesh.SetNode(*m_node_id, m_coord);

        int count_slavenodes = 0;
        // loop over all slave nodes and find matching node
        for (std::set<int>::const_iterator s_node_id = slave_nodeset_ids.begin();
             s_node_id != slave_nodeset_ids.end(); s_node_id++)
        {
          count_slavenodes++;
          // get master node coordinates
          std::vector<double> s_coord = mesh.GetNode(*s_node_id);

          // compare coordinates
          bool equal = false;
          if ((fabs(m_coord[0] - s_coord[0]) < abstol) && (fabs(m_coord[2] - s_coord[2]) < abstol))
            equal = true;

          if ((count_slavenodes == (int)slave_nodeset_ids.size()) && !equal)
            dserror("no matching slave node %d found, adjust your tolerance", *m_node_id);

          if (!equal)
            continue;
          else
          {
            // if equal copy x and z
            s_coord[0] = m_coord[0];
            s_coord[1] = y_slave;
            s_coord[2] = m_coord[2];
            mesh.SetNode(*s_node_id, s_coord);
            break;
          }
        }  // loop all slave nodes
      }    // loop all master nodes
      if (mesh.GetNodeSet(master_nodeset_id).GetNumNodes() !=
          mesh.GetNodeSet(slave_nodeset_id).GetNumNodes())
      {
        std::cout << "xz num master nodes = " << mesh.GetNodeSet(master_nodeset_id).GetNumNodes()
                  << std::endl;
        std::cout << "xz num slave nodes = " << mesh.GetNodeSet(slave_nodeset_id).GetNumNodes()
                  << std::endl;
        dserror("xz num master nodes != num slave nodes after adjusting coords");
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 * correct nodal coordinates in yz plane for                            *
 * periodic boundary conditions are defined                 u.may 02/10 *
 *----------------------------------------------------------------------*/
void EXODUS::CorrectXYPlaneForPeriodicBoundaryConditions(
    EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs)
{
  // loop over all conditions
  for (unsigned int i = 0; i < condefs.size(); i++)
  {
    // xy plane : x and y coordinates of two matching nodes have to be equal
    // xy plane:  z coordinates differ by a constant length
    if (condefs[i].sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
        condefs[i].desc.find("Master PLANE xy") != std::string::npos &&
        (condefs[i].desc.find("ANGLE 0.0") != std::string::npos ||
            condefs[i].desc.find("ANGLE 0") != std::string::npos))
    {
      // store master pbc
      const EXODUS::cond_def master_con = condefs[i];
      const int master_nodeset_id = master_con.id;

      // get condition id for master pbc
      size_t end_m = master_con.desc.find("Master");  // master
      std::string string_mconditionid = master_con.desc.substr(0, end_m - 1);
      // convert std::string to int
      std::istringstream string_mconditionidstream(string_mconditionid);
      int master_conditionID = -1;
      string_mconditionidstream >> master_conditionID;

      // find matching slave pbc node set
      EXODUS::cond_def slave_con;
      int slave_conditionID = -1;
      bool matchingIDfound = false;
      for (unsigned int i_slave = 0; i_slave < condefs.size(); i_slave++)
      {
        slave_con = condefs[i_slave];
        if (slave_con.sec.find("PERIODIC BOUNDARY CONDITIONS") != std::string::npos &&
            slave_con.desc.find("Slave PLANE xy") != std::string::npos &&
            (slave_con.desc.find("ANGLE 0.0") != std::string::npos ||
                slave_con.desc.find("ANGLE 0") != std::string::npos))
        {
          size_t end_s = slave_con.desc.find("Slave");                            // slave
          std::string string_sconditionid = slave_con.desc.substr(0, end_s - 1);  // one whitespaces
          // convert std::string to int
          std::istringstream string_sconditionidstream(string_sconditionid);
          string_sconditionidstream >> slave_conditionID;

          if (slave_conditionID == master_conditionID) matchingIDfound = true;
        }
        if (matchingIDfound) break;
      }
      if (matchingIDfound == false)  // this indicates that ANGLE != 0.0 for slave plane
        continue;                    // in this case do not proceed, but leave here!!

      // read tolerance
      size_t start_tol = slave_con.desc.find("ABSTREETOL");
      size_t tol_length = std::string("ABSTREETOL").length();
      start_tol = start_tol + tol_length + 1;
      std::string string_tol =
          slave_con.desc.substr(start_tol, slave_con.desc.length() - start_tol);  // two whitespaces
      // convert std::string to int
      std::istringstream string_tolstream(string_tol);
      double abstol = -1.0;
      string_tolstream >> abstol;

      const int slave_nodeset_id = slave_con.id;
      const EXODUS::NodeSet master_nodeset = mesh.GetNodeSet(master_nodeset_id);
      const EXODUS::NodeSet slave_nodeset = mesh.GetNodeSet(slave_nodeset_id);

      if (slave_nodeset.GetNumNodes() != master_nodeset.GetNumNodes())
      {
        std::cout << "xy num master nodes = " << master_nodeset.GetNumNodes() << std::endl;
        std::cout << "xy num slave nodes = " << slave_nodeset.GetNumNodes() << std::endl;
        dserror("xy num master nodes != num slave nodes before adjusting coords");
      }

      const std::set<int> master_nodeset_ids = master_nodeset.GetNodeSet();
      const std::set<int> slave_nodeset_ids = slave_nodeset.GetNodeSet();

      // reference values for z on the master side and z on the slave side
      const double z_master = (mesh.GetNode(*master_nodeset_ids.begin()))[2];
      const double z_slave = (mesh.GetNode(*slave_nodeset_ids.begin()))[2];

      // loop over all master nodes
      for (std::set<int>::const_iterator m_node_id = master_nodeset_ids.begin();
           m_node_id != master_nodeset_ids.end(); m_node_id++)
      {
        // get master node coordinates
        std::vector<double> m_coord = mesh.GetNode(*m_node_id);
        m_coord[2] = z_master;
        mesh.SetNode(*m_node_id, m_coord);

        int count_slavenodes = 0;
        // loop over all slave nodes and find matching node
        for (std::set<int>::const_iterator s_node_id = slave_nodeset_ids.begin();
             s_node_id != slave_nodeset_ids.end(); s_node_id++)
        {
          count_slavenodes++;
          // get master node coordinates
          std::vector<double> s_coord = mesh.GetNode(*s_node_id);

          // compare coordinates x and y
          bool equal = false;
          if ((fabs(m_coord[0] - s_coord[0]) < abstol) && (fabs(m_coord[1] - s_coord[1]) < abstol))
            equal = true;

          if ((count_slavenodes == (int)slave_nodeset_ids.size()) && !equal)
            dserror("no matching slave node %d found, adjust your tolerance", *m_node_id);

          if (!equal)
            continue;
          else
          {
            // if equal copy x and y
            s_coord[0] = m_coord[0];
            s_coord[1] = m_coord[1];
            s_coord[2] = z_slave;
            mesh.SetNode(*s_node_id, s_coord);
            break;
          }
        }  // loop all slave nodes
      }    // loop all master nodes
      if (mesh.GetNodeSet(master_nodeset_id).GetNumNodes() !=
          mesh.GetNodeSet(slave_nodeset_id).GetNumNodes())
      {
        std::cout << "xy num master nodes = " << mesh.GetNodeSet(master_nodeset_id).GetNumNodes()
                  << std::endl;
        std::cout << "xy num slave nodes = " << mesh.GetNodeSet(slave_nodeset_id).GetNumNodes()
                  << std::endl;
        dserror("xy num master nodes != num slave nodes after adjusting coords");
      }
    }
  }
  return;
}
