/*----------------------------------------------------------------------*/
/*! \file

\brief preprocessor reader for exodusII format

\maintainer Martin Kronbichler

\level 1

Here everything related with the exodus format and the accessible data
is handed to a c++ object mesh.
*/
/*----------------------------------------------------------------------*/
#include "pre_exodus_reader.H"
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "pre_exodus_soshextrusion.H"  //for gmsh plot

#include <exodusII.h>


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
EXODUS::Mesh::Mesh(const std::string exofilename)
{
  int error;
  int CPU_word_size, IO_word_size;
  float exoversion;               /* version of exodus */
  CPU_word_size = sizeof(double); /* size of a double */
  IO_word_size = 0;               /* use what is stored in file */

  const char* exofilenamechar = exofilename.c_str();

  // open EXODUS II file
  exoid_ = ex_open(exofilenamechar, EX_READ, &CPU_word_size, &IO_word_size, &exoversion);
  if (exoid_ <= 0) dserror("Error while opening EXODUS II file %s", exofilenamechar);

  // print version
  std::cout << "File " << exofilename << " was created with EXODUS II library version "
            << exoversion << std::endl;

  // read database parameters
  int num_elem_blk, num_node_sets, num_side_sets, num_nodes;
  char title[MAX_LINE_LENGTH + 1];
  error = ex_get_init(exoid_, title, &baci_dim_, &num_nodes, &num_elem_, &num_elem_blk,
      &num_node_sets, &num_side_sets);
  title_ = std::string(title);

  num_dim_ = 3;

  // get nodal coordinates
  {
    std::vector<double> x(num_nodes);
    std::vector<double> y(num_nodes);
    std::vector<double> z(num_nodes);
    error = ex_get_coord(exoid_, &x[0], &y[0], &z[0]);
    if (error != 0) dserror("exo error returned");

    // store nodes in map
    nodes_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
    for (int i = 0; i < num_nodes; ++i)
    {
      std::vector<double> coords;
      coords.push_back(x[i]);
      coords.push_back(y[i]);
      coords.push_back(z[i]);
      nodes_->insert(std::pair<int, std::vector<double>>(
          i + 1, coords));  // to store the EXO-ID starting with 1
    }
  }  // free coordinate vectors x, y ,z

  // Get all ElementBlocks
  {
    std::vector<int> epropID(num_elem_blk);
    std::vector<int> ebids(num_elem_blk);
    error = ex_get_elem_blk_ids(exoid_, &(ebids[0]));
    if (error != 0) dserror("exo error returned");
    error = ex_get_prop_array(exoid_, EX_ELEM_BLOCK, "ID", &(epropID[0]));
    if (error != 0) dserror("exo error returned");
    for (int i = 0; i < num_elem_blk; ++i)
    {
      // Read Element Blocks into Map
      char mychar[MAX_STR_LENGTH + 1];
      int num_el_in_blk, num_nod_per_elem, num_attr;
      // error = ex_get_elem_block (exoid_, epropID[i], mychar, &num_el_in_blk, &num_nod_per_elem,
      // &num_attr);
      error =
          ex_get_elem_block(exoid_, ebids[i], mychar, &num_el_in_blk, &num_nod_per_elem, &num_attr);
      if (error != 0) dserror("exo error returned");
      // prefer std::string to store element type
      std::string ele_type(mychar);

      // get ElementBlock name
      error = ex_get_name(exoid_, EX_ELEM_BLOCK, ebids[i], mychar);
      if (error != 0) dserror("exo error returned");
      // prefer std::string to store name
      std::string blockname(mychar);

      // get element connectivity
      std::vector<int> allconn(num_nod_per_elem * num_el_in_blk);
      error = ex_get_elem_conn(exoid_, ebids[i], &allconn[0]);
      if (error != 0) dserror("exo error returned");
      Teuchos::RCP<std::map<int, std::vector<int>>> eleconn =
          Teuchos::rcp(new std::map<int, std::vector<int>>);
      for (int j = 0; j < num_el_in_blk; ++j)
      {
        std::vector<int> actconn;
        actconn.reserve(num_nod_per_elem);
        for (int k = 0; k < num_nod_per_elem; ++k)
        {
          actconn.push_back(allconn[k + j * num_nod_per_elem]);
        }
        eleconn->insert(std::pair<int, std::vector<int>>(j, actconn));
      }
      Teuchos::RCP<ElementBlock> actEleBlock =
          Teuchos::rcp(new ElementBlock(StringToShape(ele_type), eleconn, blockname));

      // Add this ElementBlock into Mesh map
      elementBlocks_.insert(std::pair<int, Teuchos::RCP<ElementBlock>>(ebids[i], actEleBlock));
    }
  }  // end of element section

  // get all NodeSets
  {
    std::map<int, NodeSet> prelimNodeSets;  // prelim due to possible prop names
    std::vector<int> npropID(num_node_sets);
    error = ex_get_prop_array(exoid_, EX_NODE_SET, "ID", &(npropID[0]));
    for (int i = 0; i < num_node_sets; ++i)
    {
      // Read NodeSet params
      int num_nodes_in_set, num_df_in_set;
      error = ex_get_node_set_param(exoid_, npropID[i], &num_nodes_in_set, &num_df_in_set);

      // get NodeSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exoid_, EX_NODE_SET, npropID[i], mychar);
      // prefer std::string to store name
      std::string nodesetname(mychar);

      // get nodes in node set
      std::vector<int> node_set_node_list(num_nodes_in_set);
      error = ex_get_node_set(exoid_, npropID[i], &(node_set_node_list[0]));
      if (error > 0)
        std::cout << "'ex_get_node_set' returned warning while reading node set " << npropID[i]
                  << std::endl;
      else if (error < 0)
        dserror("error reading node set");
      std::set<int> nodes_in_set;
      for (int j = 0; j < num_nodes_in_set; ++j) nodes_in_set.insert(node_set_node_list[j]);
      NodeSet actNodeSet(nodes_in_set, nodesetname, "none");

      // Add this NodeSet into Mesh map (here prelim due to pro names)
      prelimNodeSets.insert(std::pair<int, NodeSet>(npropID[i], actNodeSet));
    }

    /* Read NodeSet property names ***********************************************
     * They are assigned by ICEM and provide recognition */
    int num_props;
    float fdum;
    char cdum;  // dummy argument
    error = ex_inquire(exoid_, EX_INQ_NS_PROP, &num_props, &fdum, &cdum);
    // allocate memory for NodeSet property names
    char** prop_names = new char*[num_props];
    for (int i = 0; i < num_props; ++i)
    {
      prop_names[i] = new char[MAX_STR_LENGTH + 1];
    }
    // get prop names of node sets
    error = ex_get_prop_names(exoid_, EX_NODE_SET, prop_names);

    // Add prop names to final Mesh NodeSet if available
    std::map<int, NodeSet>::const_iterator i_ns;
    if ((num_props - 1) == num_node_sets)
    {
      int i = 1;  // id of propname, starts with 1 because 0 is "ID"
      for (i_ns = prelimNodeSets.begin(); i_ns != prelimNodeSets.end(); ++i_ns)
      {
        std::string propname(prop_names[i]);
        const NodeSet actNodeSet = i_ns->second;
        std::string ns_name = actNodeSet.GetName();
        if (ns_name.size() == 0) ns_name = propname;
        const NodeSet newNodeSet(actNodeSet.GetNodeSet(), ns_name, propname);
        nodeSets_.insert(std::pair<int, NodeSet>(i_ns->first, newNodeSet));
        ++i;  // next propname refers to next NodeSet
      }
    }
    else
    {
      // this is the standard case without prop names
      nodeSets_ = prelimNodeSets;
    }

    // clean up node set names
    for (int i = 0; i < num_props; i++)
    {
      delete[] prop_names[i];
    }
    delete[] prop_names;
  }  // end of nodeset section
  // ***************************************************************************

  // get all SideSets
  if (num_side_sets > 0)
  {
    std::vector<int> spropID(num_side_sets);
    error = ex_get_prop_array(exoid_, EX_SIDE_SET, "ID", &(spropID[0]));
    for (int i = 0; i < num_side_sets; ++i)
    {
      // get SideSet name
      char mychar[MAX_STR_LENGTH + 1];
      error = ex_get_name(exoid_, EX_SIDE_SET, spropID[i], mychar);
      // prefer std::string to store name
      std::string sidesetname(mychar);

      // Read SideSet params
      int num_side_in_set, num_dist_fact_in_set;
      error = ex_get_side_set_param(exoid_, spropID[i], &num_side_in_set, &num_dist_fact_in_set);

      // get SideSet
      std::vector<int> side_set_elem_list(num_side_in_set);
      std::vector<int> side_set_side_list(num_side_in_set);
      error =
          ex_get_side_set(exoid_, spropID[i], &(side_set_elem_list[0]), &(side_set_side_list[0]));
      if (error != 0) dserror("error reading side set");
      std::map<int, std::vector<int>> sides_in_set;
      for (int j = 0; j < num_side_in_set; ++j)
      {
        std::vector<int> side(2);  // first entry is element, second side
        side[0] = side_set_elem_list[j];
        side[1] = side_set_side_list[j];
        sides_in_set.insert(std::pair<int, std::vector<int>>(j, side));
      }

      SideSet actSideSet(sides_in_set, sidesetname);

      // Add this SideSet into Mesh map
      sideSets_.insert(std::pair<int, SideSet>(spropID[i], actSideSet));
    }
  }  // end of sideset section

  // close ExoFile
  CloseExo();

  return;
}

EXODUS::Mesh::Mesh()
{
  nodes_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
  num_dim_ = 3;
  baci_dim_ = 3;
  num_elem_ = 0;
  exoid_ = 0;
  title_ = "emptymesh";
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
EXODUS::Mesh::~Mesh()
{
  // CloseExo();
  return;
}

/*----------------------------------------------------------------------*
 |  Extension constructor (public)                             maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::Mesh::Mesh(const EXODUS::Mesh& basemesh,
    const Teuchos::RCP<std::map<int, std::vector<double>>> extNodes,
    const std::map<int, Teuchos::RCP<ElementBlock>>& extBlocks,
    const std::map<int, NodeSet>& extNodesets, const std::map<int, SideSet>& extSidesets,
    const std::string newtitle)
    : title_(newtitle.c_str())
{
  // get all data from basemesh
  const int basedim = basemesh.GetNumDim();
  const int bacidim = basemesh.GetBACIDim();
  const int basenumele = basemesh.GetNumEle();
  Teuchos::RCP<std::map<int, std::vector<double>>> baseNodes =
      Teuchos::rcp(new std::map<int, std::vector<double>>);
  baseNodes = basemesh.GetNodes();
  // Teuchos::RCP<std::map<int,std::vector<double> > > baseNodes = basemesh.GetNodes();
  std::map<int, Teuchos::RCP<ElementBlock>> baseEblocks = basemesh.GetElementBlocks();
  std::map<int, NodeSet> baseNodesets = basemesh.GetNodeSets();
  std::map<int, SideSet> baseSidesets = basemesh.GetSideSets();

  //  // get infos from extension
  //  int extnumele = extBlocks.size();

  /********************* merge everything into new mesh ***********************/
  num_dim_ = basedim;
  baci_dim_ = bacidim;
  int total_num_elem = basenumele;
  //  num_elem_ = basenumele; // + extnumele;
  exoid_ =
      basemesh
          .GetExoId();  // basefile still used for writing minor infos, e.g. qa record or coordnames

  // merge nodes
  std::map<int, std::vector<double>>::const_iterator i_node;
  for (i_node = extNodes->begin(); i_node != extNodes->end(); ++i_node)
  {
    std::pair<std::map<int, std::vector<double>>::iterator, bool> check;
    check = baseNodes->insert(std::pair<int, std::vector<double>>(i_node->first, i_node->second));
    // happens when concatenating: if (check.second == false)  dserror("Extension node already
    // exists!");
  }
  nodes_ = baseNodes;

  // merge ElementBlocks
  std::map<int, Teuchos::RCP<ElementBlock>>::const_iterator i_block;
  for (i_block = baseEblocks.begin(); i_block != baseEblocks.end(); ++i_block)
  {
    elementBlocks_.insert(
        std::pair<int, Teuchos::RCP<ElementBlock>>(i_block->first, i_block->second));
  }
  for (i_block = extBlocks.begin(); i_block != extBlocks.end(); ++i_block)
  {
    std::pair<std::map<int, Teuchos::RCP<ElementBlock>>::iterator, bool> check;
    check = elementBlocks_.insert(
        std::pair<int, Teuchos::RCP<ElementBlock>>(i_block->first, i_block->second));
    if (check.second == false)
      dserror("Extension ElementBlock already exists!");
    else
      total_num_elem += i_block->second->GetNumEle();
  }
  num_elem_ = total_num_elem;

  // merge NodeSets
  std::map<int, NodeSet>::const_iterator i_ns;
  for (i_ns = baseNodesets.begin(); i_ns != baseNodesets.end(); ++i_ns)
  {
    nodeSets_.insert(std::pair<int, NodeSet>(i_ns->first, i_ns->second));
  }
  for (i_ns = extNodesets.begin(); i_ns != extNodesets.end(); ++i_ns)
  {
    std::pair<std::map<int, NodeSet>::iterator, bool> check;
    check = nodeSets_.insert(std::pair<int, NodeSet>(i_ns->first, i_ns->second));
    if (check.second == false) dserror("Extension NodeSet already exists!");
  }

  // merge SideSets
  std::map<int, SideSet>::const_iterator i_ss;
  for (i_ss = baseSidesets.begin(); i_ss != baseSidesets.end(); ++i_ss)
  {
    sideSets_.insert(std::pair<int, SideSet>(i_ss->first, i_ss->second));
  }
  for (i_ss = extSidesets.begin(); i_ss != extSidesets.end(); ++i_ss)
  {
    std::pair<std::map<int, SideSet>::iterator, bool> check;
    check = sideSets_.insert(std::pair<int, SideSet>(i_ss->first, i_ss->second));
    if (check.second == false) dserror("Extension SideSet already exists!");
  }
}


/*----------------------------------------------------------------------*
 |  Close corresponding Exofile(public)                        maf 12/07|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::CloseExo() const
{
  // close exodus II file
  int exoid = GetExoId();
  int error = ex_close(exoid);
  if (error < 0) dserror("error while closing exodus II file");
}


/*----------------------------------------------------------------------*
 |  Print method (public)                                      maf 12/07|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::Print(std::ostream& os, bool verbose) const
{
  os << "Mesh consists of ";
  os << GetNumNodes() << " Nodes, ";
  os << num_elem_ << " Elements, organized in " << std::endl;
  os << GetNumElementBlocks() << " ElementBlocks, ";
  os << GetNumNodeSets() << " NodeSets, ";
  os << GetNumSideSets() << " SideSets ";
  os << std::endl << std::endl;
  if (verbose)
  {
    os << "ElementBlocks" << std::endl;
    std::map<int, Teuchos::RCP<ElementBlock>>::const_iterator it;
    std::map<int, Teuchos::RCP<ElementBlock>> eleBlocks = GetElementBlocks();
    for (it = eleBlocks.begin(); it != eleBlocks.end(); it++)
    {
      os << it->first << ": ";
      it->second->Print(os);
    }
    os << std::endl << "NodeSets" << std::endl;
    std::map<int, NodeSet>::const_iterator it2;
    std::map<int, NodeSet> nodeSets = GetNodeSets();
    for (it2 = nodeSets.begin(); it2 != nodeSets.end(); it2++)
    {
      os << "NodeSet " << it2->first << ": ";
      it2->second.Print(os);
    }
    os << std::endl << "SideSets" << std::endl;
    os << "Warning: SideSets are not yet fully supported by PreExodus!" << std::endl;
    std::map<int, SideSet>::const_iterator it3;
    std::map<int, SideSet> sideSets = GetSideSets();
    for (it3 = sideSets.begin(); it3 != sideSets.end(); it3++)
    {
      os << "SideSet " << it3->first << ": ";
      it3->second.Print(os);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string EXODUS::Mesh::GetTitle() const
{
  std::string title(title_, int(MAX_LINE_LENGTH + 1));
  return title;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<EXODUS::ElementBlock> EXODUS::Mesh::GetElementBlock(const int id) const
{
  if (elementBlocks_.find(id) == elementBlocks_.end())
  {
    std::cout << "ElementBlock " << id << " not found!" << std::endl;
    dserror("ElementBlock not found");
  }
  return (elementBlocks_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet EXODUS::Mesh::GetNodeSet(const int id) const
{
  if (nodeSets_.find(id) == nodeSets_.end())
  {
    std::cout << "NodeSet " << id << " not found!" << std::endl;
    dserror("NodeSet not found");
  }
  return (nodeSets_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::SideSet EXODUS::Mesh::GetSideSet(const int id) const
{
  if (sideSets_.find(id) == sideSets_.end())
  {
    std::cout << "SideSet " << id << " not found!" << std::endl;
    dserror("SideSet not found");
  }
  return (sideSets_.find(id))->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::PrintNodes(std::ostream& os, bool storeid) const
{
  std::map<int, std::vector<double>>::const_iterator it;
  for (it = nodes_->begin(); it != nodes_->end(); it++)
  {
    if (storeid) os << "MapID: " << it->first;
    int exoid = it->first + 1;
    os << " ExoID: " << exoid << " : ";
    const std::vector<double> mycoords = it->second;
    for (int i = 0; i < signed(mycoords.size()); i++)
    {
      os << mycoords[i] << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::GetNode(const int NodeID) const
{
  std::map<int, std::vector<double>>::const_iterator it = nodes_->find(NodeID);
  return it->second;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::SetNode(const int NodeID, const std::vector<double> coord)
{
  // if entry exits already , delete it first and the insert the new value
  // other wise nothing is inserted
  if (nodes_->find(NodeID) != nodes_->end()) nodes_->erase(NodeID);

  nodes_->insert(std::make_pair(NodeID, coord));
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::Mesh::SetNsd(const int nsd)
{
  if (nsd != 2 && nsd != 3) return;

  baci_dim_ = nsd;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::vector<int>> EXODUS::Mesh::GetSideSetConn(const SideSet sideset) const
{
  std::cout << "Creating SideSet Connectivity... " << std::endl;
  fflush(stdout);

  Epetra_SerialComm Comm;
  Epetra_Time time(Comm);
  Teuchos::RCP<Teuchos::Time> timetot;
  timetot = Teuchos::TimeMonitor::getNewTimer("Side Set Connect total");
  Teuchos::RCP<Teuchos::Time> time1 = Teuchos::TimeMonitor::getNewTimer("One Side Set");
  Teuchos::RCP<Teuchos::Time> time2 = Teuchos::TimeMonitor::getNewTimer("Get one Element Block");
  Teuchos::RCP<Teuchos::Time> time3 = Teuchos::TimeMonitor::getNewTimer("Get Ele Conn");
  Teuchos::RCP<Teuchos::Time> time4 = Teuchos::TimeMonitor::getNewTimer("Get one Ele");
  Teuchos::RCP<Teuchos::Time> time5 = Teuchos::TimeMonitor::getNewTimer("Build one Side Conn");
  Teuchos::RCP<Teuchos::Time> time6 = Teuchos::TimeMonitor::getNewTimer("One Side Set");
  Teuchos::RCP<Teuchos::Time> time7 =
      Teuchos::TimeMonitor::getNewTimer("Get all Eblocks and Econns");
  Teuchos::RCP<Teuchos::TimeMonitor> tm_total = Teuchos::rcp(new Teuchos::TimeMonitor(*timetot));

  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>> mysides = sideset.GetSideSet();
  std::map<int, std::vector<int>>::iterator i_side;

  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;
  // Also we once get all EBlocks and EConns to enable quick access
  std::vector<EXODUS::ElementBlock> eblocks;
  std::vector<std::map<int, std::vector<int>>> econns;
  Teuchos::RCP<Teuchos::TimeMonitor> tm7 = Teuchos::rcp(new Teuchos::TimeMonitor(*time7));
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->GetNumEle();
    glob_eb_erange.push_back(rangebreak);
    eblocks.push_back(*i_ebs->second);
    econns.push_back(*(i_ebs->second->GetEleConn()));
  }
  tm7 = Teuchos::null;

  // fill SideSet Connectivity
  // int perc = 1;
  int tetc = 0, hexc = 0, pyrc = 0, wedgc = 0;
  for (i_side = mysides.begin(); i_side != mysides.end(); ++i_side)
  {
    Teuchos::RCP<Teuchos::TimeMonitor> tm1 = Teuchos::rcp(new Teuchos::TimeMonitor(*time1));
    int actele = i_side->second.at(0) -
                 1;  // ExoIds start from 1, but we from 0 //ToDo: remove -1 idconfusion
    int actface = i_side->second.at(1) - 1;  // ExoIds start from 1, but we from 0
    // find actual EBlock where actele lies in
    int actebid = -1;
    for (unsigned int i = 0; i < glob_eb_erange.size(); ++i)
      if (actele < glob_eb_erange[i])
      {
        actebid = i - 1;
        break;
      }
    Teuchos::RCP<Teuchos::TimeMonitor> tm2 = Teuchos::rcp(new Teuchos::TimeMonitor(*time2));
    // EXODUS::ElementBlock acteb = ebs.find(actebid)->second;
    tm2 = Teuchos::null;
    // EXODUS::ElementBlock::Shape actshape = acteb.GetShape();
    if (actebid < 0) dserror("invalid element block id");
    EXODUS::ElementBlock::Shape actshape = eblocks[actebid].GetShape();
    Teuchos::RCP<Teuchos::TimeMonitor> tm3 = Teuchos::rcp(new Teuchos::TimeMonitor(*time3));
    // map<int,std::vector<int> > acteconn = acteb.GetEleConn();
    tm3 = Teuchos::null;
    // get act parent ele from actual Side
    int parent_ele_id = actele - glob_eb_erange[actebid];
    Teuchos::RCP<Teuchos::TimeMonitor> tm4 = Teuchos::rcp(new Teuchos::TimeMonitor(*time4));
    // vector<int> parent_ele = acteconn.find(parent_ele_id)->second;
    std::vector<int> parent_ele = econns[actebid].find(parent_ele_id)->second;
    tm4 = Teuchos::null;
    // Face to ElementNode Map
    //// **** temporary hex map due to conflicts between side numbering exo<->baci
    Teuchos::RCP<Teuchos::TimeMonitor> tm5 = Teuchos::rcp(new Teuchos::TimeMonitor(*time5));
    switch (actshape)
    {
      case ElementBlock::tet4:
      {
        tetc++;
        break;
      }
      case ElementBlock::hex8:
      {
        actface = HexSideNumberExoToBaci(actface);
        hexc++;
        break;
      }
      case ElementBlock::pyramid5:
      {
        //      vector<std::vector<int> > test =
        //      DRT::UTILS::getEleNodeNumberingSurfaces(PreShapeToDrt(actshape)); for(unsigned int
        //      j=0; j<test.size(); ++j) PrintVec(std::cout,test[j]);
        actface = PyrSideNumberExoToBaci(actface);
        pyrc++;
        break;
      }
      case ElementBlock::wedge6:
      {
        //      vector<std::vector<int> > test =
        //      DRT::UTILS::getEleNodeNumberingSurfaces(PreShapeToDrt(actshape)); for(unsigned int
        //      j=0; j<test.size(); ++j) PrintVec(std::cout,test[j]);
        wedgc++;
        break;
      }
      default:
      {
        std::cout << ShapeToString(actshape) << ":" << std::endl;
        dserror("Parent Element Type not supported");
        break;
      }
    }
    std::vector<int> childmap =
        DRT::UTILS::getEleNodeNumberingSurfaces(PreShapeToDrt(actshape))[actface];
    // child gets its node ids
    std::vector<int> child;
    for (unsigned int j = 0; j < childmap.size(); ++j) child.push_back(parent_ele[childmap[j]]);
    //    PrintVec(std::cout,childmap);
    //    PrintVec(std::cout,child);
    // some checking
    if ((child.size() != 3) && (child.size() != 4))
    {
      PrintVec(std::cout, child);
      PrintVec(std::cout, childmap);
      PrintVec(std::cout, parent_ele);
      std::cout << ShapeToString(actshape) << ",Face: " << actface << ",childsize:" << child.size()
                << std::endl;
      dserror("Child Ele error");
    }
    // insert child into SideSet Connectivity
    conn.insert(std::pair<int, std::vector<int>>(i_side->first, child));
    tm5 = Teuchos::null;

    //    // progress output
    //    tm1 = Teuchos::null;
    //    if (signed(i_side->first) == perc * signed(mysides.size())/100){
    //      Teuchos::TimeMonitor::summarize();
    //      std::cout << perc << " % of " << mysides.size() << " Sides done" << std::endl;
    //      perc ++;
    //      fflush(stdout);
    //    }
  }
  tm_total = Teuchos::null;
  //  Teuchos::TimeMonitor::summarize();
  std::cout << "...done" << std::endl;
  fflush(stdout);

  return conn;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::vector<int>> EXODUS::Mesh::GetSideSetConn(
    const SideSet sideset, bool checkoutside) const
{
  std::cout << "Creating SideSet Connectivity with outside-check... " << std::endl;
  fflush(stdout);

  Epetra_SerialComm Comm;
  Epetra_Time time(Comm);
  Teuchos::RCP<Teuchos::Time> timetot;
  timetot = Teuchos::TimeMonitor::getNewTimer("Side Set Connect total");
  Teuchos::RCP<Teuchos::Time> time1 = Teuchos::TimeMonitor::getNewTimer("One Side Set");
  Teuchos::RCP<Teuchos::TimeMonitor> tm_total = Teuchos::rcp(new Teuchos::TimeMonitor(*timetot));

  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>> mysides = sideset.GetSideSet();
  std::map<int, std::vector<int>>::iterator i_side;

  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;
  // Also we once get all EBlocks and EConns to enable quick access
  std::vector<EXODUS::ElementBlock> eblocks;
  std::vector<std::map<int, std::vector<int>>> econns;
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->GetNumEle();
    glob_eb_erange.push_back(rangebreak);
    eblocks.push_back(*i_ebs->second);
    econns.push_back(*(i_ebs->second->GetEleConn()));
  }

  // fill SideSet Connectivity
  int tetc = 0, hexc = 0, pyrc = 0, wedgc = 0;
  for (i_side = mysides.begin(); i_side != mysides.end(); ++i_side)
  {
    Teuchos::RCP<Teuchos::TimeMonitor> tm1 = Teuchos::rcp(new Teuchos::TimeMonitor(*time1));
    int actele = i_side->second.at(0) -
                 1;  // ExoIds start from 1, but we from 0 //ToDo: remove -1 idconfusion
    int actface = i_side->second.at(1) - 1;  // ExoIds start from 1, but we from 0
    // find actual EBlock where actele lies in
    int actebid = -1;
    for (unsigned int i = 0; i < glob_eb_erange.size(); ++i)
      if (actele < glob_eb_erange[i])
      {
        actebid = i - 1;
        break;
      }
    if (actebid < 0) dserror("invalid element block id");
    EXODUS::ElementBlock::Shape actshape = eblocks[actebid].GetShape();

    // get act parent ele from actual Side
    int parent_ele_id = actele - glob_eb_erange[actebid];
    std::vector<int> parent_ele = econns[actebid].find(parent_ele_id)->second;

    // Face to ElementNode Map
    //// **** temporary hex map due to conflicts between side numbering exo<->baci
    switch (actshape)
    {
      case ElementBlock::tet4:
      {
        tetc++;
        break;
      }
      case ElementBlock::hex8:
      {
        actface = HexSideNumberExoToBaci(actface);
        hexc++;
        break;
      }
      case ElementBlock::pyramid5:
      {
        actface = PyrSideNumberExoToBaci(actface);
        pyrc++;
        break;
      }
      case ElementBlock::wedge6:
      {
        wedgc++;
        break;
      }
      default:
      {
        std::cout << ShapeToString(actshape) << ":" << std::endl;
        dserror("Parent Element Type not supported");
        break;
      }
    }
    std::vector<int> childmap =
        DRT::UTILS::getEleNodeNumberingSurfaces(PreShapeToDrt(actshape))[actface];

    std::vector<int> child;
    if (checkoutside)
      child = OutsideOrientedSide(parent_ele, childmap);
    else
      for (unsigned int j = 0; j < childmap.size(); ++j) child.push_back(parent_ele[childmap[j]]);

    // some checking
    if ((child.size() != 3) && (child.size() != 4))
    {
      PrintVec(std::cout, child);
      PrintVec(std::cout, childmap);
      PrintVec(std::cout, parent_ele);
      std::cout << ShapeToString(actshape) << ",Face: " << actface << ",childsize:" << child.size()
                << std::endl;
      dserror("Child Ele error");
    }

    // insert child into SideSet Connectivity
    conn.insert(std::pair<int, std::vector<int>>(i_side->first, child));
    tm1 = Teuchos::null;
  }
  tm_total = Teuchos::null;
  //  Teuchos::TimeMonitor::summarize();
  std::cout << "...done" << std::endl;
  fflush(stdout);

  return conn;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::Mesh::OutsideOrientedSide(
    const std::vector<int> parentele, const std::vector<int> sidemap) const
{
  // first guess of child
  std::vector<int> child;
  // set simplifies later inverse search
  std::set<int> childset;
  for (unsigned int j = 0; j < sidemap.size(); ++j)
  {
    child.push_back(parentele[sidemap[j]]);
    childset.insert(parentele[sidemap[j]]);
  }

  // set of parentele for later inverse search
  std::set<int> parentset;
  for (unsigned int i = 0; i < parentele.size(); ++i) parentset.insert(parentele.at(i));

  // find parentele node not within side
  int insidenode = -1;
  std::set<int>::iterator it;
  for (it = parentset.begin(); it != parentset.end(); ++it)
  {
    if (childset.find(*it) == childset.end())
    {
      insidenode = *it;
      break;
    }
  }

  // build normal at first side node
  std::vector<double> sidenormal = Normal(child.back(), child.front(), child.at(1));
  // build vector from first side node to inside element node
  std::vector<double> insidevec = NodeVec(child.front(), insidenode);

  // scalar-product
  double scp =
      sidenormal[0] * insidevec[0] + sidenormal[1] * insidevec[1] + sidenormal[2] * insidevec[2];

  std::vector<int> out_side;
  if (scp < 0)
  {
    std::vector<int> reversechild;
    std::vector<int>::reverse_iterator rit;
    for (rit = child.rbegin(); rit < child.rend(); ++rit) out_side.push_back(*rit);
  }
  else
    out_side = child;

  return out_side;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::Normal(const int head1, const int origin, const int head2) const
{
  std::vector<double> normal(3);
  std::vector<double> h1 = GetNode(head1);
  std::vector<double> h2 = GetNode(head2);
  std::vector<double> o = GetNode(origin);

  normal[0] = ((h1[1] - o[1]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[1] - o[1]));
  normal[1] = -((h1[0] - o[0]) * (h2[2] - o[2]) - (h1[2] - o[2]) * (h2[0] - o[0]));
  normal[2] = ((h1[0] - o[0]) * (h2[1] - o[1]) - (h1[1] - o[1]) * (h2[0] - o[0]));

  double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  normal[0] = normal[0] / length;
  normal[1] = normal[1] / length;
  normal[2] = normal[2] / length;

  return normal;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> EXODUS::Mesh::NodeVec(const int tail, const int head) const
{
  std::vector<double> nv(3);
  std::vector<double> t = GetNode(tail);
  std::vector<double> h = GetNode(head);
  nv[0] = h[0] - t[0];
  nv[1] = h[1] - t[1];
  nv[2] = h[2] - t[2];
  double length = sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
  nv[0] = nv[0] / length;
  nv[1] = nv[1] / length;
  nv[2] = nv[2] / length;
  return nv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<EXODUS::ElementBlock> EXODUS::Mesh::SideSetToEBlocks(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sideconn) const
{
  std::vector<ElementBlock> eblocks;
  // map<int,std::vector<int> > sideconn = sideset.GetSideSet();
  std::map<int, std::vector<int>>::const_iterator i_ele;
  Teuchos::RCP<std::map<int, std::vector<int>>> quadconn =
      Teuchos::rcp(new std::map<int, std::vector<int>>);
  int quadcounter = 0;
  Teuchos::RCP<std::map<int, std::vector<int>>> triconn =
      Teuchos::rcp(new std::map<int, std::vector<int>>);
  int tricounter = 0;
  for (i_ele = sideconn.begin(); i_ele != sideconn.end(); ++i_ele)
  {
    int numnodes = i_ele->second.size();
    if (numnodes == 4)
    {
      quadconn->insert(std::pair<int, std::vector<int>>(quadcounter, i_ele->second));
      quadcounter++;
    }
    else if (numnodes == 3)
    {
      triconn->insert(std::pair<int, std::vector<int>>(tricounter, i_ele->second));
      tricounter++;
    }
    else
      dserror("Number of basenodes for conversion from SideSet to EBlock not supported");
  }
  if (quadcounter > 0)
  {
    std::ostringstream quadblockname;
    quadblockname << sideset.GetName() << "quad";
    EXODUS::ElementBlock neweblock(ElementBlock::quad4, quadconn, quadblockname.str());
    eblocks.push_back(neweblock);
  }
  if (tricounter > 0)
  {
    std::ostringstream triblockname;
    triblockname << sideset.GetName() << "tri";
    EXODUS::ElementBlock neweblock(ElementBlock::tri3, triconn, triblockname.str());
    eblocks.push_back(neweblock);
  }

  return eblocks;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet EXODUS::Mesh::SideSetToNodeSet(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sideconn) const
{
  std::map<int, std::vector<int>>::const_iterator i_side;
  std::vector<int>::const_iterator i_node;
  std::set<int> nodes;
  for (i_side = sideconn.begin(); i_side != sideconn.end(); ++i_side)
    for (i_node = i_side->second.begin(); i_node != i_side->second.end(); ++i_node)
      nodes.insert(*i_node);  // nodes.insert(i_side->second.at(i_node));
  std::ostringstream nodesetname;
  nodesetname << "nodes";  // sideset.GetName() << "nodes";
  std::string propname = "";
  EXODUS::NodeSet nodeset(nodes, nodesetname.str(), propname);

  return nodeset;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::set<int> EXODUS::Mesh::GetSideSetNodes(
    const EXODUS::SideSet& sideset, const std::map<int, std::vector<int>>& sideconn) const
{
  std::map<int, std::vector<int>>::const_iterator i_side;
  std::vector<int>::const_iterator i_node;
  std::set<int> nodes;
  for (i_side = sideconn.begin(); i_side != sideconn.end(); ++i_side)
    for (i_node = i_side->second.begin(); i_node != i_side->second.end(); ++i_node)
      nodes.insert(*i_node);  // nodes.insert(i_side->second.at(i_node));
  return nodes;
}

/*----------------------------------------------------------------------*
 |  Write Mesh into exodus file (public)                       maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::WriteMesh(const std::string newexofilename) const
{
  std::cout << "Writing Mesh into: " << newexofilename << std::endl;
  // std::string newexofile(newexofilename);
  // newexofile += ".exo";
  const char* newexofilechar = newexofilename.c_str();

  int CPU_word_size, IO_word_size, exoid, error;
  CPU_word_size = sizeof(double); /* use double*/
  IO_word_size = 8;               /* store variables as doubles */

  /* create EXODUS II file */
  exoid = ex_create(newexofilechar, /* filename path */
      EX_CLOBBER,                   /* create mode */
      &CPU_word_size,               /* CPU double word size in bytes */
      &IO_word_size);               /* I/O double word size in bytes */

  int num_elem_blk = GetNumElementBlocks();
  int num_node_sets = GetNumNodeSets();
  int num_side_sets = GetNumSideSets();
  int num_nodes = GetNumNodes();
  /* initialize file with parameters */
  const char* title = title_.c_str();
  error = ex_put_init(
      exoid, title, num_dim_, num_nodes, num_elem_, num_elem_blk, num_node_sets, num_side_sets);
  if (error != 0) dserror("error in exfile init");

  /* Write QA record based on original exofile */
  int num_qa_rec;
  char* qa_record[MAX_STR_LENGTH][4];  // should be MAX_QA_REC][4], but this is nowhere defined!;
  char cdum;                           // dummy variable
  float fdum;
  /* read QA records */
  ex_inquire(exoid_, EX_INQ_QA, &num_qa_rec, &fdum, &cdum); /* write QA records */
  for (int i = 0; i < num_qa_rec; i++)
    for (int j = 0; j < 4; j++) qa_record[i][j] = (char*)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  // for (int j=0; j<4; j++) qa_record[i][j] = new (char)[MAX_STR_LENGTH+1];
  error = ex_get_qa(exoid_, qa_record);
  error = ex_put_qa(exoid, num_qa_rec, qa_record);
  for (int i = 0; i < num_qa_rec; i++)
    for (int j = 0; j < 4; j++)
      // delete [] qa_record[i][j];
      free(qa_record[i][j]);

  // Write coord names based on original exofile
  char* coord_names[3];
  for (int i = 0; i < num_dim_; i++)
    coord_names[i] = (char*)calloc((MAX_STR_LENGTH + 1), sizeof(char));
  error = ex_get_coord_names(exoid_, coord_names);
  error = ex_put_coord_names(exoid, coord_names);
  for (int i = 0; i < num_dim_; i++) free(coord_names[i]);

  // Write nodal coordinates
  std::vector<double> xc(num_nodes);
  std::vector<double> yc(num_nodes);
  std::vector<double> zc(num_nodes);
  std::map<int, std::vector<double>>::const_iterator it;
  Teuchos::RCP<std::map<int, std::vector<double>>> nodes = GetNodes();
  for (it = nodes->begin(); it != nodes->end(); ++it)
  {
    xc[it->first - 1] = it->second[0];  // vector starts with 0
    yc[it->first - 1] = it->second[1];  // vector starts with 0
    zc[it->first - 1] = it->second[2];  // vector starts with 0
  }
  error = ex_put_coord(exoid, &xc[0], &yc[0], &zc[0]);

  // Write NodeSets ************************************************************
  std::map<int, NodeSet>::const_iterator ins;
  const std::map<int, NodeSet> nss = GetNodeSets();
  for (ins = nss.begin(); ins != nss.end(); ++ins)
  {
    const int nsID = ins->first;
    const NodeSet ns = ins->second;
    const int num_nodes_in_set = ns.GetNumNodes();
    if (num_nodes_in_set > 0)  // do not bother if nodeset empty
    {
      const std::string name = ns.GetName();
      const char* nsname = name.c_str();
      const std::string propname = ns.GetPropName();
      const std::set<int> nodes = ns.GetNodeSet();
      error = ex_put_node_set_param(exoid,  // of write file
          nsID,                             // node set id
          num_nodes_in_set,
          0);  // yet no distribution factors
      if (error != 0) dserror("error writing node set params");
      std::vector<int> nodelist(num_nodes_in_set);
      ns.FillNodelistArray(&nodelist[0]);
      error = ex_put_node_set(exoid, nsID, &nodelist[0]);
      if (error != 0) dserror("error writing node set \"%s\" ", nsname);
      error = ex_put_name(exoid, EX_NODE_SET, nsID, nsname);
      if (error != 0) dserror("error writing node set name");
    }
  }

  // Write ElementBlocks  ******************************************************
  std::map<int, Teuchos::RCP<ElementBlock>>::const_iterator iebs;
  const std::map<int, Teuchos::RCP<ElementBlock>> ebs = GetElementBlocks();
  for (iebs = ebs.begin(); iebs != ebs.end(); iebs++)
  {
    const int blockID = iebs->first;
    const ElementBlock eb = (*iebs->second);
    const ElementBlock::Shape shape = eb.GetShape();
    const std::string shapestring = ShapeToString(shape);
    const std::vector<int> exampleconn = eb.GetEleNodes(0);  // iebs->first);
    const int num_nod_per_elem = exampleconn.size();
    const int numele = eb.GetNumEle();
    const char* elem_type = shapestring.c_str();
    error = ex_put_elem_block(exoid,  // of write file
        blockID,                      // element block id
        elem_type,                    // its name
        numele,                       // number of element in block
        num_nod_per_elem,             // num of nodes per ele
        1);                           // num of attributes, not supported yet ->1
    if (error != 0) dserror("error writing element block");
    // Write Element Connectivity
    std::vector<int> conn(num_nod_per_elem * numele);
    eb.FillEconnArray(&conn[0]);
    error = ex_put_elem_conn(exoid, blockID, &conn[0]);
    if (error != 0) dserror("error writing element block conns");
    // write block name
    const std::string bname = eb.GetName();
    const char* blockname = bname.c_str();
    error = ex_put_name(exoid, EX_ELEM_BLOCK, blockID, blockname);
    if (error != 0) dserror("error writing element block name");
  }

  // Write SideSets ************************************************************
  std::map<int, SideSet>::const_iterator iss;
  const std::map<int, SideSet> sss = GetSideSets();
  for (iss = sss.begin(); iss != sss.end(); ++iss)
  {
    const int ssID = iss->first;
    const SideSet ss = iss->second;
    const int num_side_in_set = ss.GetNumSides();
    const std::string name = ss.GetName();
    const char* ssname = name.c_str();
    error = ex_put_side_set_param(exoid,  // of write file
        ssID,                             // side set id
        num_side_in_set,
        0);  // yet no distribution factors
    if (error != 0) dserror("error writing side set params");
    std::vector<int> side_set_elem_list(num_side_in_set);
    std::vector<int> side_set_side_list(num_side_in_set);
    // in case the sideset is newly created we have to adjust element ids to global numbering
    std::map<int, std::vector<int>> globalsides;
    if (iss->second.GetFirstSideSet().size() == 3)
    {
      globalsides = GlobalifySSeleids(ssID);
      ss.FillSideLists(&side_set_elem_list[0], &side_set_side_list[0], globalsides);
    }
    else
      ss.FillSideLists(&side_set_elem_list[0], &side_set_side_list[0]);

    error = ex_put_side_set(exoid, ssID, &side_set_elem_list[0], &side_set_side_list[0]);
    if (error != 0) dserror("error writing side set");
    error = ex_put_name(exoid, EX_SIDE_SET, ssID, ssname);
    if (error != 0) dserror("error writing sideset name");
  }

  // ***************************************************************************

  // close file
  error = ex_close(exoid);
  if (error != 0) dserror("error closing exodus file");
  std::cout << ".. finished" << std::endl;
}

/*----------------------------------------------------------------------*
 |  Add Element Block to mesh(public)                          maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::AddElementBlock(const Teuchos::RCP<EXODUS::ElementBlock> eblock) const
{
  std::map<int, Teuchos::RCP<ElementBlock>> eblocks = GetElementBlocks();
  eblocks.insert(std::pair<int, Teuchos::RCP<ElementBlock>>(GetNumElementBlocks() + 1, eblock));
}

/*----------------------------------------------------------------------*
 |  Erase Element Block from mesh(public)                      maf 07/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::EraseElementBlock(const int id)
{
  int red_numele = GetElementBlock(id)->GetNumEle();
  elementBlocks_.erase(id);
  num_elem_ = num_elem_ - red_numele;
}

/*----------------------------------------------------------------------*
 |  Erase SideSet from mesh(public)                            maf 07/08|
 *----------------------------------------------------------------------*/
void EXODUS::Mesh::EraseSideSet(const int id) { sideSets_.erase(id); }

/*------------------------------------------------------------------------*
 | - calculates the midpoint of each element                               |
 | - returns map <midpoint-ID,std::pair<eblock-ID,element-ID> >         SP 06/08|
 *------------------------------------------------------------------------*/
std::map<int, std::pair<int, int>> EXODUS::Mesh::createMidpoints(
    std::map<int, std::vector<double>>& midpoints, const std::vector<int>& eb_ids) const
{
  // map that will be returned
  std::map<int, std::pair<int, int>> conn_mpID_elID;

  //  //initialising midpoints
  //  this->midpoints_ = Teuchos::rcp(new std::map<int,std::vector<double> >);

  // auxiliary variables
  int counter_elements = 0;
  int nodes_per_element = 0;

  std::vector<double> sumVector(3, 0);
  std::vector<double> midPoint(3, 0);

  std::map<int, Teuchos::RCP<ElementBlock>> EBlocks;
  std::map<int, Teuchos::RCP<ElementBlock>>::const_iterator it;
  // work only on eblocks in eb_ids
  std::vector<int>::const_iterator id;
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs;
  for (id = eb_ids.begin(); id != eb_ids.end(); ++id)
  {
    Teuchos::RCP<EXODUS::ElementBlock> acteb = this->GetElementBlock(*id);
    EBlocks.insert(std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(*id, acteb));
  }

  std::map<int, std::vector<int>> EleConn;
  std::map<int, std::vector<int>>::const_iterator it_2;

  std::vector<int>::const_iterator it_3;

  // loop over ElementBlocks
  for (it = EBlocks.begin(); it != EBlocks.end(); ++it)
  {
    EleConn = *(it->second->GetEleConn());

    // loop over element connectivity
    for (it_2 = EleConn.begin(); it_2 != EleConn.end(); ++it_2)
    {
      counter_elements++;
      nodes_per_element = 0;

      sumVector[0] = 0;
      sumVector[1] = 0;
      sumVector[2] = 0;

      // loop over each Node in one element connectivity
      for (it_3 = it_2->second.begin(); it_3 != it_2->second.end(); ++it_3)
      {
        nodes_per_element++;

        // sum of two vectors
        sumVector[0] += GetNode(*it_3)[0];
        sumVector[1] += GetNode(*it_3)[1];
        sumVector[2] += GetNode(*it_3)[2];
      }

      // midpoint of element i
      midPoint[0] = sumVector[0] / nodes_per_element;
      midPoint[1] = sumVector[1] / nodes_per_element;
      midPoint[2] = sumVector[2] / nodes_per_element;

      // insert calculated midpoint in midpoints_
      midpoints.insert(std::pair<int, std::vector<double>>(counter_elements, midPoint));
      // conn_mpID_elID = (midpoint-ID, eblock-ID, element-ID)
      std::pair<int, int> eb_e = std::make_pair(it->first, it_2->first);
      conn_mpID_elID.insert(std::pair<int, std::pair<int, int>>(counter_elements, eb_e));
    }
  }
  // EXODUS::PrintMap(std::cout,conn_mpID_elID);
  return conn_mpID_elID;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::vector<int>> EXODUS::Mesh::GlobalifySSeleids(const int ssid) const
{
  SideSet ss = GetSideSet(ssid);

  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;

  // Range Vector for global eleID identification in SideSet
  std::vector<int> glob_eb_erange(1, 0);
  int rangebreak = 0;

  std::map<int, int> ebid_rangepos;
  int rangepos = 0;
  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    rangebreak += i_ebs->second->GetNumEle();
    glob_eb_erange.push_back(rangebreak);
    ebid_rangepos.insert(std::pair<int, int>(i_ebs->first, rangepos));
    ++rangepos;
  }

  std::map<int, std::vector<int>> sideset = ss.GetSideSet();
  std::map<int, std::vector<int>>::iterator i_ss;

  for (i_ss = sideset.begin(); i_ss != sideset.end(); ++i_ss)
  {
    if (i_ss->second.size() != 3) dserror("Problem in new SideSet!");  // double check
    int ebid = i_ss->second.at(2);
    int lowerbound = glob_eb_erange.at(ebid_rangepos.find(ebid)->second);
    i_ss->second.at(0) = i_ss->second.at(0) + lowerbound + 1;
  }

  return sideset;
}


/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize mesh                             MF 07/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::PlotElementBlocksGmsh(const std::string fname, const EXODUS::Mesh& mymesh) const
{
  Teuchos::RCP<std::map<int, std::vector<double>>> nodes = mymesh.GetNodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Mesh \" {" << std::endl;

  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = mymesh.GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator eb_it;
  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>>::const_iterator it;

  for (eb_it = ebs.begin(); eb_it != ebs.end(); ++eb_it)
  {
    Teuchos::RCP<std::map<int, std::vector<int>>> actconn = eb_it->second->GetEleConn();
    for (it = actconn->begin(); it != actconn->end(); ++it)
    {
      int eleid = it->first;
      const std::vector<int> elenodes = it->second;
      int numnodes = elenodes.size();
      if (numnodes == 6)
        gmshfilecontent << "SI(";
      else if (numnodes == 8)
        gmshfilecontent << "SH(";
      for (unsigned int i = 0; i < elenodes.size(); ++i)
      {
        gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
        if (i == (elenodes.size() - 1))
          gmshfilecontent << ")";
        else
          gmshfilecontent << ",";
      }
      gmshfilecontent << "{";
      for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
      gmshfilecontent << eleid << "};" << std::endl;
    }
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize mesh                             MF 07/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::PlotElementBlocksGmsh(
    const std::string fname, const EXODUS::Mesh& mymesh, const std::vector<int>& ebids) const
{
  Teuchos::RCP<std::map<int, std::vector<double>>> nodes = mymesh.GetNodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Mesh \" {" << std::endl;

  std::vector<int>::const_iterator id;
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs;
  for (id = ebids.begin(); id != ebids.end(); ++id)
  {
    Teuchos::RCP<EXODUS::ElementBlock> acteb = mymesh.GetElementBlock(*id);
    ebs.insert(std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(*id, acteb));
  }

  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator eb_it;
  std::map<int, std::vector<int>> conn;
  std::map<int, std::vector<int>>::const_iterator it;

  for (eb_it = ebs.begin(); eb_it != ebs.end(); ++eb_it)
  {
    Teuchos::RCP<std::map<int, std::vector<int>>> actconn = eb_it->second->GetEleConn();
    for (it = actconn->begin(); it != actconn->end(); ++it)
    {
      int eleid = it->first;
      const std::vector<int> elenodes = it->second;
      int numnodes = elenodes.size();
      if (numnodes == 6)
        gmshfilecontent << "SI(";
      else if (numnodes == 8)
        gmshfilecontent << "SH(";
      for (unsigned int i = 0; i < elenodes.size(); ++i)
      {
        gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
        gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
        if (i == (elenodes.size() - 1))
          gmshfilecontent << ")";
        else
          gmshfilecontent << ",";
      }
      gmshfilecontent << "{";
      for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
      gmshfilecontent << eleid << "};" << std::endl;
    }
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize all nodes                        SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::PlotNodesGmsh() const
{
  std::ofstream f_system("mesh_all_nodes.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Nodes \" {" << std::endl;

  std::map<int, std::vector<double>>::const_iterator it;
  // loop over all nodes
  for (it = nodes_->begin(); it != nodes_->end(); it++)
  {
    const std::vector<double> mycoords = it->second;

    // writing of coordinates of each node
    gmshfilecontent << "SP(";
    gmshfilecontent << mycoords[0] << ",";
    gmshfilecontent << mycoords[1] << ",";
    gmshfilecontent << mycoords[2];
    gmshfilecontent << ")";

    // writing of node-ID
    gmshfilecontent << "{";
    gmshfilecontent << it->first << "," << it->first << "," << it->first << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize connectivity                     MF 12/08|
 *------------------------------------------------------------------------*/
void EXODUS::Mesh::PlotConnGmsh(const std::string fname, const EXODUS::Mesh& mymesh,
    const std::map<int, std::vector<int>>& conn) const
{
  Teuchos::RCP<std::map<int, std::vector<double>>> nodes = mymesh.GetNodes();
  std::ofstream f_system(fname.c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Connectivity \" {" << std::endl;

  std::map<int, std::vector<int>>::const_iterator it;

  for (it = conn.begin(); it != conn.end(); ++it)
  {
    int eleid = it->first;
    const std::vector<int> elenodes = it->second;
    int numnodes = elenodes.size();
    if (numnodes == 6)
      gmshfilecontent << "SI(";
    else if (numnodes == 8)
      gmshfilecontent << "SH(";
    else if (numnodes == 3)
      gmshfilecontent << "ST(";
    else if (numnodes == 4)
      gmshfilecontent << "SQ(";
    for (unsigned int i = 0; i < elenodes.size(); ++i)
    {
      gmshfilecontent << nodes->find(elenodes.at(i))->second[0] << ",";
      gmshfilecontent << nodes->find(elenodes.at(i))->second[1] << ",";
      gmshfilecontent << nodes->find(elenodes.at(i))->second[2];
      if (i == (elenodes.size() - 1))
        gmshfilecontent << ")";
      else
        gmshfilecontent << ",";
    }
    gmshfilecontent << "{";
    for (unsigned int i = 0; i < (elenodes.size() - 1); ++i) gmshfilecontent << eleid << ",";
    gmshfilecontent << eleid << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::ElementBlock::ElementBlock(ElementBlock::Shape Distype,
    Teuchos::RCP<std::map<int, std::vector<int>>>& eleconn, std::string name)
    : distype_(Distype), eleconn_(eleconn), name_(name.c_str())
{
  // do a sanity check
  for (std::map<int, std::vector<int>>::const_iterator elem = eleconn->begin();
       elem != eleconn->end(); ++elem)
  {
    if (DRT::UTILS::getNumberOfElementNodes(PreShapeToDrt(Distype)) != (int)elem->second.size())
    {
      dserror("number of read nodes does not fit the distype");
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::ElementBlock::~ElementBlock() { return; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::ElementBlock::GetEleNodes(int i) const
{
  std::map<int, std::vector<int>>::const_iterator it = eleconn_->find(i);
  return it->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::ElementBlock::GetEleNode(int ele, int node) const
{
  std::map<int, std::vector<int>>::const_iterator it = eleconn_->find(ele);
  if (it == eleconn_->end()) dserror("Element not found");
  std::vector<int> elenodes = GetEleNodes(ele);
  return elenodes[node];
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ElementBlock::FillEconnArray(int* connarray) const
{
  const std::map<int, std::vector<int>>::const_iterator iele;
  int numele = eleconn_->size();
  for (int i = 0; i < numele; ++i)
  {
    std::vector<int> ele = GetEleNodes(i);
    int num_nod_per_elem = ele.size();
    for (int j = 0; j < num_nod_per_elem; ++j)
    {
      connarray[i * num_nod_per_elem + j] = GetEleNode(i, j);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ElementBlock::Print(std::ostream& os, bool verbose) const
{
  os << "Element Block, named: " << name_ << std::endl
     << "of Shape: " << ShapeToString(distype_) << std::endl
     << "has " << GetNumEle() << " Elements" << std::endl;
  if (verbose)
  {
    std::map<int, std::vector<int>>::const_iterator it;
    for (it = eleconn_->begin(); it != eleconn_->end(); it++)
    {
      os << "Ele " << it->first << ": ";
      const std::vector<int> myconn = it->second;  // GetEleNodes(int(it));
      for (int i = 0; i < signed(myconn.size()); i++)
      {
        os << myconn[i] << ",";
      }
      os << std::endl;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet::NodeSet(
    const std::set<int>& nodeids, const std::string& name, const std::string& propname)
    : nodeids_(nodeids), name_(name.c_str()), propname_(propname.c_str())
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::NodeSet::~NodeSet() { return; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::NodeSet::Print(std::ostream& os, bool verbose) const
{
  os << "Node Set, named: " << name_ << std::endl
     << "Property Name: " << propname_ << std::endl
     << "has " << GetNumNodes() << " Nodes" << std::endl;
  if (verbose)
  {
    os << "Contains Nodes:" << std::endl;
    std::set<int>::iterator it;
    for (it = nodeids_.begin(); it != nodeids_.end(); it++) os << *it << ",";
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::NodeSet::FillNodelistArray(int* nodelist) const
{
  std::set<int> nlist = GetNodeSet();
  std::set<int>::iterator it;
  int i = 0;
  for (it = nlist.begin(); it != nlist.end(); ++it)
  {
    nodelist[i] = (*it);
    ++i;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::SideSet::SideSet(const std::map<int, std::vector<int>>& sides, const std::string& name)
    : sides_(sides), name_(name.c_str())
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::SideSet::FillSideLists(int* elemlist, int* sidelist) const
{
  std::map<int, std::vector<int>> sides = GetSideSet();
  std::map<int, std::vector<int>>::iterator it;
  int i = 0;
  for (it = sides.begin(); it != sides.end(); ++it)
  {
    elemlist[i] = it->second[0];
    sidelist[i] = it->second[1];
    ++i;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::SideSet::FillSideLists(
    int* elemlist, int* sidelist, const std::map<int, std::vector<int>>& sides) const
{
  std::map<int, std::vector<int>>::const_iterator it;
  int i = 0;
  for (it = sides.begin(); it != sides.end(); ++it)
  {
    elemlist[i] = it->second[0];
    sidelist[i] = it->second[1];
    ++i;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::SideSet::Print(std::ostream& os, bool verbose) const
{
  os << "SideSet, named: " << name_ << std::endl
     << "has " << GetNumSides() << " Sides" << std::endl;
  if (verbose)
  {
    std::map<int, std::vector<int>>::const_iterator it;
    for (it = sides_.begin(); it != sides_.end(); it++)
    {
      os << "Side " << it->first << ": ";
      os << "Ele: " << it->second.at(0) << ", Side: " << it->second.at(1) << std::endl;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::SideSet::~SideSet() { return; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::Mesh EXODUS::QuadtoTri(EXODUS::Mesh& basemesh)
{
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> neweblocks;  // here the new EBlocks are stored
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = basemesh.GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;

  for (i_ebs = ebs.begin(); i_ebs != ebs.end(); ++i_ebs)
  {
    Teuchos::RCP<EXODUS::ElementBlock> quadblock = i_ebs->second;
    EXODUS::ElementBlock::Shape quadshape = quadblock->GetShape();
    if ((quadshape != EXODUS::ElementBlock::quad4) && (quadshape != EXODUS::ElementBlock::shell4))
    {
      // dserror("Only quad4 or shell4 in quad->tri conversion");
      std::cout << "Warning! Only quad4 or shell4 in quad->tri conversion. Skipping EBlock"
                << std::endl;
    }
    else
    {
      Teuchos::RCP<std::map<int, std::vector<int>>> quad_conn = quadblock->GetEleConn();
      std::map<int, std::vector<int>>::const_iterator i_quad;
      Teuchos::RCP<std::map<int, std::vector<int>>> triconn =
          Teuchos::rcp(new std::map<int, std::vector<int>>);

      for (i_quad = quad_conn->begin(); i_quad != quad_conn->end(); ++i_quad)
      {
        std::vector<int> quad = i_quad->second;
        std::vector<int> tri1(3);
        std::vector<int> tri2(3);
        tri1[0] = quad[0];
        tri1[1] = quad[1];
        tri1[2] = quad[2];
        tri2[0] = quad[2];
        tri2[1] = quad[3];
        tri2[2] = quad[0];
        int tri1_id = 2 * i_quad->first;
        int tri2_id = 2 * i_quad->first + 1;
        triconn->insert(std::pair<int, std::vector<int>>(tri1_id, tri1));
        triconn->insert(std::pair<int, std::vector<int>>(tri2_id, tri2));
      }

      Teuchos::RCP<EXODUS::ElementBlock> triblock = Teuchos::rcp(
          new EXODUS::ElementBlock(EXODUS::ElementBlock::tri3, triconn, quadblock->GetName()));
      neweblocks.insert(std::pair<int, Teuchos::RCP<EXODUS::ElementBlock>>(i_ebs->first, triblock));
      basemesh.EraseElementBlock(i_ebs->first);
    }
  }

  std::string newtitle = "trimesh";
  std::map<int, EXODUS::NodeSet> emptynodeset;
  std::map<int, EXODUS::SideSet> emptysideset;
  if (basemesh.GetNumSideSets() > 0)
    std::cout << "Warning! SideSets will not be transferred by quad->tri!" << std::endl;
  Teuchos::RCP<std::map<int, std::vector<double>>> emptynodes =
      Teuchos::rcp(new std::map<int, std::vector<double>>);
  EXODUS::Mesh trimesh(basemesh, emptynodes, neweblocks, emptynodeset, emptysideset, newtitle);
  return trimesh;
}



/*----------------------------------------------------------------------*
  TINY HELPER FUNCTIONS
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, std::vector<int>> mymap)
{
  std::map<int, std::vector<int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::vector<int> actvec = iter->second;
    std::vector<int>::iterator i;
    for (i = actvec.begin(); i < actvec.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, std::set<int>> mymap)
{
  std::map<int, std::set<int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::set<int> actset = iter->second;
    std::set<int>::iterator i;
    for (i = actset.begin(); i != actset.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, std::vector<double>> mymap)
{
  std::map<int, std::vector<double>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::vector<double> actvec = iter->second;
    std::vector<double>::iterator i;
    for (i = actvec.begin(); i < actvec.end(); ++i)
    {
      os << *i << ",";
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<double, int> mymap)
{
  std::map<double, int>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": " << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, std::map<int, int>> mymap)
{
  std::map<int, std::map<int, int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::map<int, int> actmap = iter->second;
    std::map<int, int>::iterator i;
    for (i = actmap.begin(); i != actmap.end(); ++i)
    {
      os << i->first << ": " << i->second;
    }
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, std::pair<int, int>> mymap)
{
  std::map<int, std::pair<int, int>>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    std::pair<int, int> actpair = iter->second;
    os << actpair.first << " <=> " << actpair.second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, int> mymap)
{
  std::map<int, int>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    os << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintMap(std::ostream& os, const std::map<int, double> mymap)
{
  std::map<int, double>::const_iterator iter;
  for (iter = mymap.begin(); iter != mymap.end(); ++iter)
  {
    os << iter->first << ": ";
    os << iter->second;
    os << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintVec(std::ostream& os, const std::vector<int> actvec)
{
  std::vector<int>::const_iterator i;
  for (i = actvec.begin(); i < actvec.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintVec(std::ostream& os, const std::vector<double> actvec)
{
  std::vector<double>::const_iterator i;
  for (i = actvec.begin(); i < actvec.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::PrintSet(std::ostream& os, const std::set<int> actset)
{
  std::set<int>::iterator i;
  for (i = actset.begin(); i != actset.end(); ++i)
  {
    os << *i << ",";
  }
  os << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::HexSideNumberExoToBaci(const int exoface)
{
  const int map[6] = {1, 2, 3, 4, 0, 5};
  return map[exoface];
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::PyrSideNumberExoToBaci(const int exoface)
{
  const int map[5] = {1, 2, 3, 4, 0};
  return map[exoface];
}
