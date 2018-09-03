/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_centerline.cpp

\brief pre_exodus centerline definition

\level 1

\maintainer Martin Kronbichler

 */
/*----------------------------------------------------------------------*/

#include "pre_exodus_centerline.H"
#include <iostream>
#include "../drt_mat/matpar_bundle.H"  //
#include "../drt_lib/drt_colors.H"
#include "pre_exodus_soshextrusion.H"  // for Gmsh plot

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, std::map<int, std::vector<std::vector<double>>>> EXODUS::EleCenterlineInfo(
    std::string& cline, EXODUS::Mesh& mymesh, const std::vector<double> coordcorr)
{
  if (cline == "mesh")
  {  // create Centerline object from NodeSet
    int centerlineid = -1;

    std::map<int, EXODUS::NodeSet> nss = mymesh.GetNodeSets();
    std::map<int, EXODUS::NodeSet>::const_iterator i_ns;
    // check for Centerline or Centerpoint Nodeset
    bool clbool = true;
    for (i_ns = nss.begin(); i_ns != nss.end(); ++i_ns)
    {
      const std::string myname = i_ns->second.GetName();
      if (myname.find("centerline") != std::string::npos)
      {
        centerlineid = i_ns->first;
        clbool = true;
      }
      else if (myname.find("centerpoint") != std::string::npos)
      {
        clbool = false;
        centerlineid = i_ns->first;
      }
    }
    if (centerlineid == -1) dserror("Have not found centerline NodeSet");

    EXODUS::Centerline myCLine(nss.find(centerlineid)->second, mymesh.GetNodes());
    myCLine.PlotCL_Gmsh();  // generation of accordant Gmsh-file

    // get rid of helper eb where the centerline ns was based on
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = mymesh.GetElementBlocks();
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_eb;
    std::vector<int> eb_ids;
    // check for Centerline ElementBlock
    for (i_eb = ebs.begin(); i_eb != ebs.end(); ++i_eb)
    {
      const std::string myname = i_eb->second->GetName();
      if (myname.find("centerline") != std::string::npos) mymesh.EraseElementBlock(i_eb->first);
      // check for CenterPoints block
      else if (myname.find("centerpoint") != std::string::npos)
        mymesh.EraseElementBlock(i_eb->first);
      else
        eb_ids.push_back(i_eb->first);
    }

    // generation of coordinate systems
    std::map<int, std::map<int, std::vector<std::vector<double>>>> centlineinfo;
    if (clbool) centlineinfo = EXODUS::element_cosys(myCLine, mymesh, eb_ids);
    // generation of degenerated coordinate systems
    else
      centlineinfo = EXODUS::element_degcosys(myCLine, mymesh, eb_ids);
    if (clbool) EXODUS::PlotCosys(myCLine, mymesh, eb_ids);  // generation of accordant Gmsh-file
    // plot mesh to gmsh
    std::string meshname = "centerlinemesh.gmsh";
    mymesh.PlotElementBlocksGmsh(meshname, mymesh);


    return centlineinfo;
  }
  else if (cline.find(".exo") != std::string::npos)
  {
    // read centerline from another exodus file
    EXODUS::Mesh centerlinemesh(cline);

    int centerlineid = -1;

    std::map<int, EXODUS::NodeSet> nss = centerlinemesh.GetNodeSets();
    std::map<int, EXODUS::NodeSet>::const_iterator i_ns;
    // check for Centerline or Centerpoint Nodeset
    bool clbool = true;
    for (i_ns = nss.begin(); i_ns != nss.end(); ++i_ns)
    {
      const std::string myname = i_ns->second.GetName();
      if (myname.find("centerline") != std::string::npos)
      {
        centerlineid = i_ns->first;
        clbool = true;
      }
      else if (myname.find("centerpoint") != std::string::npos)
      {
        clbool = false;
        centerlineid = i_ns->first;
      }
    }
    if (centerlineid == -1) dserror("Have not found centerline NodeSet");
    EXODUS::Centerline myCLine(nss.find(centerlineid)->second, centerlinemesh.GetNodes());
    myCLine.PlotCL_Gmsh();  // generation of accordant Gmsh-file

    // get rid of helper eb where the centerline ns was based on
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = mymesh.GetElementBlocks();
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_eb;
    std::vector<int> eb_ids;
    // check for Centerline ElementBlock
    for (i_eb = ebs.begin(); i_eb != ebs.end(); ++i_eb)
    {
      const std::string myname = i_eb->second->GetName();
      if (myname.find("centerline") != std::string::npos)
        centerlinemesh.EraseElementBlock(i_eb->first);
      // check for CenterPoints block
      else if (myname.find("centerpoint") != std::string::npos)
        centerlinemesh.EraseElementBlock(i_eb->first);
      else
        eb_ids.push_back(i_eb->first);
    }

    // generation of coordinate systems
    std::map<int, std::map<int, std::vector<std::vector<double>>>> centlineinfo;
    if (clbool)
    {
      // switch with which method normals are computed
      // get nodeset with all surface nodes
      int surfnodes_id = -1;
      std::map<int, EXODUS::NodeSet> nss = mymesh.GetNodeSets();
      std::map<int, EXODUS::NodeSet>::const_iterator i_ns;
      // check for Surface Nodeset for normal calculation
      for (i_ns = nss.begin(); i_ns != nss.end(); ++i_ns)
      {
        const std::string myname = i_ns->second.GetName();
        if (myname.find("roof_surface_all") != std::string::npos)
        {
          surfnodes_id = i_ns->first;
        }
      }
      if (surfnodes_id == -1)
      {
        centlineinfo = EXODUS::element_cosys(myCLine, mymesh, eb_ids);
      }
      // get all node ids from this nodeset
      else
      {
        std::cout << RED << "Found surface nodes NodeSet! Fiber directions computed accordingly"
                  << END_COLOR << std::endl;

        std::set<int> all_surfnodes = nss.find(surfnodes_id)->second.GetNodeSet();

        // int node_id_first_node = *all_surfnodes.begin();
        // select only one element block for fiber direction calculation
        eb_ids.clear();
        for (i_eb = ebs.begin(); i_eb != ebs.end(); ++i_eb)
        {
          // search for node_id in all el_block
          // i_eb->second->GetEleConn()->find(node_id_first_node);
          // set key featuer o be the name for now
          const std::string myname = i_eb->second->GetName();
          if (myname.find("exth3") != std::string::npos)
          {
            eb_ids.push_back((i_eb)->first);
          }
        }
        // check if there is more than one element block
        if (eb_ids.size() != 1)
          dserror("ERROR: COMPUTATION OF FIBER DIRECTIONS IS ONLY SUPPORTED FOR ONE ELEMENT BLOCK");
        // check if Elementblock is of shape HEX8
        if (mymesh.GetElementBlock(eb_ids[0])->GetShape() != EXODUS::ElementBlock::hex8)
          dserror("ERROR: COMPUATION OF FIBER DIRECTION SUPPORTS ONLY HEX 8 ELEMENTS");
        centlineinfo = EXODUS::element_cosys(myCLine, mymesh, eb_ids, all_surfnodes);
      }
    }
    // generation of degenerated coordinate systems
    else
      centlineinfo = EXODUS::element_degcosys(myCLine, mymesh, eb_ids);
    if (clbool) EXODUS::PlotCosys(myCLine, mymesh, eb_ids);  // generation of accordant Gmsh-file
    // plot mesh to gmsh
    std::string meshname = "centerlinemesh.gmsh";
    centerlinemesh.PlotElementBlocksGmsh(meshname, mymesh);


    return centlineinfo;
  }
  else
  {  // creation of a Centerline object from file
    std::cout << "Reading centerline..." << std::endl;
    EXODUS::Centerline myCLine(cline, coordcorr);
    std::cout << "...done" << std::endl;

    // myCLine.PrintPoints();

    // get ids of the eblocks you want to calculate the locsys's
    std::string identifier = "ext";
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = mymesh.GetElementBlocks();
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_eb;
    std::vector<int> eb_ids;

    for (i_eb = ebs.begin(); i_eb != ebs.end(); ++i_eb)
    {
      std::string actname = i_eb->second->GetName();
      size_t found;
      found = actname.find(identifier);
      if (found != std::string::npos) eb_ids.push_back(i_eb->first);
    }

    std::cout << "Generating local cosys..." << std::endl;
    myCLine.PlotCL_Gmsh();  // generation of accordant Gmsh-file
    // generation of coordinate systems
    std::map<int, std::map<int, std::vector<std::vector<double>>>> centlineinfo =
        EXODUS::element_cosys(myCLine, mymesh, eb_ids);
    std::cout << "...done" << std::endl;

    //    std::cout << "Generating gmsh plots..." << std::endl;
    //    EXODUS::PlotCosys(myCLine,mymesh,eb_ids);       //generation of accordant Gmsh-file
    //
    //    // plot mesh to gmsh
    //    std::string meshname = "centerlinemesh.gmsh";
    //    mymesh.PlotElementBlocksGmsh(meshname,mymesh,eb_ids);
    //    std::cout << "...done" << std::endl;


    return centlineinfo;
  }

  // weirdo impossible case
  std::map<int, std::map<int, std::vector<std::vector<double>>>> mymap;
  return mymap;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int, double> EXODUS::NdCenterlineThickness(std::string cline, const std::set<int>& nodes,
    const std::map<int, std::vector<int>>& conn, const EXODUS::Mesh& mesh, const double ratio,
    const std::vector<double> coordcorr)
{
  std::map<int, double> ndthick;
  std::set<int>::const_iterator i_node;

  // create centerline object with diameter info
  std::cout << "Reading centerline and generating variable nodal thicknesses..." << std::endl;
  EXODUS::Centerline mycline(cline, coordcorr);
  std::map<int, std::vector<double>> clpoints = *(mycline.GetPoints());
  std::map<int, double> cldiams = *(mycline.GetDiams());

  mycline
      .PlotCL_Gmsh();  // plot centerline with coloured diameter (change comment in plot function!)
  mesh.PlotConnGmsh("extrusionsideset.gmsh", mesh, conn);

  // go through all nodes and search for closest cline-point to get respective thickness
  for (i_node = nodes.begin(); i_node != nodes.end(); ++i_node)
  {
    std::vector<double> nodecoords = mesh.GetNode(*i_node);

    double min_distance = EXODUS::distance3d(nodecoords, clpoints.begin()->second);
    int clID = -1;

    // loop over all points of the centerline to find closest
    for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
         cl_iter != clpoints.end(); ++cl_iter)
    {
      double dist = EXODUS::distance3d(nodecoords, cl_iter->second);
      if (min_distance > dist)
      {
        min_distance = dist;
        clID = cl_iter->first;
      }
    }

    // calculate thickness at current node
    double radius = 0.5 * cldiams.find(clID)->second;
    // ratio = thickness / radius
    double thickness = ratio * radius;

    ndthick.insert(std::pair<int, double>(*i_node, thickness));
  }
  std::cout << "...done" << std::endl;

  return ndthick;
}



/*------------------------------------------------------------------------*
 |Ctor                                                            SP 06/08|
 *------------------------------------------------------------------------*/
EXODUS::Centerline::Centerline(std::string filename, std::vector<double> coordcorr)
{
  // initialization of points_
  points_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
  diam_ = Teuchos::rcp(new std::map<int, double>);

  // routine to open file
  std::ifstream infile;
  infile.open(filename.c_str(), std::ifstream::in);

  // check
  if (!infile)
  {
    std::cout << "Could not open Centerline file: " << filename << std::endl;
    dserror("Could not open Centerline file!");
  }

  // read in the whole file into a "table"
  // for large file this might be memory intensive!
  typedef std::vector<float> Row;
  std::vector<Row> table;

  while (infile)
  {
    std::string line;
    getline(infile, line);
    std::istringstream is(line);
    Row row;
    while (is)
    {
      float data;
      is >> data;
      row.push_back(data);
    }
    table.push_back(row);
  }
  infile.close();

  // sort the table
  int clp_id = 0;
  for (unsigned i = 0; i < table.size(); i++)
  {
    Row row = table[i];
    if (row.size() > 1)
    {  // if true we have a "number-row"
      std::vector<double> clp(3);
      for (int j = 0; j < 3; ++j)
      {
        // correct possible offset of coords
        // clp[j] = table[i][j] + coordcorr[j];  // first 3 numbers are CL coords
        clp[j] = row[j] + coordcorr[j];  // first 3 numbers are CL coords
      }
      points_->insert(std::pair<int, std::vector<double>>(clp_id, clp));  // fill std::map
      // diam_->insert(std::pair<int,double>(clp_id,row[7])); // fitted diameter is 7th entry in
      // cline table
      diam_->insert(
          std::pair<int, double>(clp_id, row[7]));  // minimum diameter is 8th entry in cline table
      ++clp_id;
    }
  }

  // PrintMap(std::cout,*points_);

  /* Stefans old code to read matlab file and shift coords
  //auxiliary variables
  double d;
  int i=0,j=0;
  std::vector<double> CLPoint(3,0);

  while(infile.read((char*) &d, sizeof(d)))  //reads coordinates of points in d
  {
    //and stores all three coordinates of one point in CLPoint

    // displacement of coordinate systems of the centerline and the mesh must be considered
    // therefore we transform the centerline to fit into the mesh
    double delta = 205;
    double scale = -1.0;
    // mind that coordinates also have switched! x->y y->x
    switch(i%3)
    {
    case 0: CLPoint[1] = scale*d + delta; break;
    case 1: CLPoint[0] = d; break;
    case 2:
      {
        CLPoint[2] = d;
        //fill points_ with current point of centerline
        points_->insert(std::pair<int,std::vector<double> >(j,CLPoint));
        ++j;
        break;
      }
    default: break;
    }
    ++i;
  }
  infile.close();
  */
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EXODUS::Centerline::Centerline(
    const EXODUS::NodeSet& ns, const Teuchos::RCP<std::map<int, std::vector<double>>> nodes)
{
  // initialization of points_
  points_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
  diam_ = Teuchos::rcp(new std::map<int, double>);  // stays empty for mesh-based clines
  std::set<int> nodeset = ns.GetNodeSet();
  std::set<int>::const_iterator it;
  int id = 0;
  // ns.Print(std::cout,true);
  for (it = nodeset.begin(); it != nodeset.end(); ++it)
  {
    std::vector<double> node = nodes->find(*it)->second;
    points_->insert(std::pair<int, std::vector<double>>(id, node));
    ++id;
  }
}

/*------------------------------------------------------------------------*
 |Dtor                                                            SP 06/08|
 *------------------------------------------------------------------------*/
EXODUS::Centerline::~Centerline() {}

/*------------------------------------------------------------------------*
 |displays points_ on console                                      SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::Centerline::PrintPoints()
{
  for (std::map<int, std::vector<double>>::const_iterator it = points_->begin();
       it != points_->end(); ++it)
  {
    std::cout << it->first << ": " << it->second[0] << " " << it->second[1] << " " << it->second[2]
              << std::endl;
  }
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize points                            SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::Centerline::PlotCL_Gmsh()
{
  std::ofstream gmshFile("centerline.gmsh");
  gmshFile << "View \" Centerline \" {" << std::endl;

  for (std::map<int, std::vector<double>>::const_iterator it = points_->begin();
       it != points_->end(); ++it)
  {
    gmshFile << "SP(" << it->second[0] << "," << it->second[1] << "," << it->second[2] << "){";
    // gmshFile << it->first << "};" << std::endl;  // id as color
    gmshFile << diam_->find(it->first)->second << "};" << std::endl;  // diameter as color
  }

  gmshFile << "};";
  gmshFile.close();
}

/*------------------------------------------------------------------------*
 |calculates distance of two 3-dim vectors                         SP 06/08|
 *------------------------------------------------------------------------*/
double EXODUS::distance3d(std::vector<double> v1, std::vector<double> v2)
{
  double distance;
  distance = sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) +
                  (v1[2] - v2[2]) * (v1[2] - v2[2]));
  return distance;
}

/*------------------------------------------------------------------------*
 |calculates difference of two 3-dim vectors                       SP 06/08|
 *------------------------------------------------------------------------*/
std::vector<double> EXODUS::substract3d(std::vector<double> v1, std::vector<double> v2)
{
  std::vector<double> result(3, 0);

  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];

  return result;
}

/*------------------------------------------------------------------------*
 |calculates sum of two 3-dim vectors                              SP 06/08|
 *------------------------------------------------------------------------*/
std::vector<double> EXODUS::add3d(std::vector<double> v1, std::vector<double> v2)
{
  std::vector<double> result(3, 0);

  result[0] = v1[0] + v2[0];
  result[1] = v1[1] + v2[1];
  result[2] = v1[2] + v2[2];

  return result;
}


/*------------------------------------------------------------------------*
 |calculates cross product of two 3-dim vectors                    SP 06/08|
 *------------------------------------------------------------------------*/
std::vector<double> EXODUS::cross_product3d(std::vector<double> v1, std::vector<double> v2)
{
  std::vector<double> result(3, 0);

  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return result;
}

/*------------------------------------------------------------------------*
 |calculates scalar product of two 3-dim vectors                   JB 06/08|
 *------------------------------------------------------------------------*/
double EXODUS::scalar_product3d(std::vector<double> v1, std::vector<double> v2)
{
  double result;

  result = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  return result;
}

/*------------------------------------------------------------------------*
 |normalizes a 3-dim vector                                       SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::normalize3d(std::vector<double>& v)
{
  double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] = v[0] / d;
  v[1] = v[1] / d;
  v[2] = v[2] / d;
}

/*------------------------------------------------------------------------*
 |- creates local coordinate systems for each element             JB 12/09|
 |- based on surface normals calculated from surface nodes
 |- returns a map containing calc. directions referred to each element    |
 *------------------------------------------------------------------------*/
std::map<int, std::map<int, std::vector<std::vector<double>>>> EXODUS::element_cosys(
    EXODUS::Centerline& mycline, const EXODUS::Mesh& mymesh, const std::vector<int>& eb_ids,
    std::set<int>& mysurfnodes)
{
  // std::map which stores the surface normals for each element
  std::map<int, std::vector<double>> ele_normals;

  std::map<std::pair<int, int>, std::vector<double>> np_eb_el;
  // Get all Element Blocks
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> ebs = mymesh.GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator i_ebs;
  // vector to store the  surface nodes of the each element
  std::vector<std::vector<double>> ele_surf_nodes;
  std::vector<std::vector<double>>::iterator counter_ele_surfnodes;
  std::map<int, std::vector<int>> EleConn;
  std::map<int, std::vector<int>>::const_iterator it_2;
  std::vector<int>::const_iterator it_3;
  // store node connecting vectors
  std::vector<double> ele_vec_12, ele_vec_14, ele_vec_23, ele_vec_21, ele_vec_34, ele_vec_32,
      ele_vec_41, ele_vec_43;
  // vectors to store the nodal normal vectors
  std::vector<double> ele_vec_norm_1, ele_vec_norm_2, ele_vec_norm_3, ele_vec_norm_4;
  std::vector<double> ele_vec_norm(3, 0);
  // int counter_ele_surfnodes =0;
  int ele_counter = 0;


  EleConn = *(mymesh.GetElementBlock(eb_ids[0])->GetEleConn());
  // start loop over all elements
  // Loop over all elements in Block
  for (it_2 = EleConn.begin(); it_2 != EleConn.end(); ++it_2)
  {
    ele_counter++;
    counter_ele_surfnodes = ele_surf_nodes.begin();
    ele_surf_nodes.clear();
    // Loop over a single element
    for (it_3 = it_2->second.begin(); it_3 != it_2->second.end(); ++it_3)
    {
      // Find the nodes that are in surface nodes nodeset

      // if(all_surfnodes.find(*it_3)!=all_surfnodes.end())
      if (mysurfnodes.find(*it_3) != mysurfnodes.end())
      {
        ele_surf_nodes.insert(counter_ele_surfnodes, (mymesh.GetNode(*it_3)));
        counter_ele_surfnodes = ele_surf_nodes.begin();
        ;
      }
    }  // loop over element

    // Calculate  element normal vector from element surface nodes
    // only if more than 2 nodes lie on the surface
    if ((ele_surf_nodes.size()) > 2)
    {
      // Get two vectors in element plane node 1
      ele_vec_12 = substract3d(ele_surf_nodes.at(1), ele_surf_nodes.at(0));
      ele_vec_14 = substract3d(ele_surf_nodes.at(3), ele_surf_nodes.at(0));
      // Get two vectors in element plane node 2
      ele_vec_23 = substract3d(ele_surf_nodes.at(2), ele_surf_nodes.at(1));
      ele_vec_21 = substract3d(ele_surf_nodes.at(0), ele_surf_nodes.at(1));
      // Get two vectors in element plane node 3
      ele_vec_34 = substract3d(ele_surf_nodes.at(3), ele_surf_nodes.at(2));
      ele_vec_32 = substract3d(ele_surf_nodes.at(1), ele_surf_nodes.at(2));
      // Get two vectors in element plane node 4
      ele_vec_41 = substract3d(ele_surf_nodes.at(0), ele_surf_nodes.at(3));
      ele_vec_43 = substract3d(ele_surf_nodes.at(2), ele_surf_nodes.at(3));
      // create nodal normal vectors
      ele_vec_norm_1 = cross_product3d(ele_vec_12, ele_vec_14);
      ele_vec_norm_2 = cross_product3d(ele_vec_23, ele_vec_21);
      ele_vec_norm_3 = cross_product3d(ele_vec_34, ele_vec_32);
      ele_vec_norm_4 = cross_product3d(ele_vec_41, ele_vec_43);
      // build average element normal vector
      ele_vec_norm[0] =
          0.25 * (ele_vec_norm_1[0] + ele_vec_norm_2[0] + ele_vec_norm_3[0] + ele_vec_norm_4[0]);
      ele_vec_norm[1] =
          0.25 * (ele_vec_norm_1[1] + ele_vec_norm_2[1] + ele_vec_norm_3[1] + ele_vec_norm_4[1]);
      ele_vec_norm[2] =
          0.25 * (ele_vec_norm_1[2] + ele_vec_norm_2[2] + ele_vec_norm_3[2] + ele_vec_norm_4[2]);
      // normalize normal vector
      normalize3d(ele_vec_norm);
      ele_normals.insert(std::pair<int, std::vector<double>>(ele_counter, ele_vec_norm));
      // Write data in std::map which stores ele Block data as well
      // pair<int,int> eb_e = make_pair(i_ebs->first,it_2->first);
      std::pair<int, int> eb_e = std::make_pair(eb_ids[0], it_2->first);
      // np_eb_el contains ((eblock-ID, element-ID),ele_normal)
      np_eb_el.insert(std::pair<std::pair<int, int>, std::vector<double>>(eb_e, ele_vec_norm));
    }
  }

  std::map<int, std::vector<double>> midpoints;  // here midpoints are stored
  // mp_eb_el contains (midpoint-ID, (eblock-ID, element-ID))
  std::map<int, std::pair<int, int>> mp_eb_el = mymesh.createMidpoints(midpoints, eb_ids);
  // conn_mp_cp will contain (midpoint-ID, centerpoint-ID_1, centerpoint-ID_2)
  std::map<int, std::vector<int>> conn_mp_cp;
  // auxiliary variables
  int clID = -1;
  int clID_2 = -1;
  std::vector<int> ids(2, 0);
  std::vector<double> mean_cl_dir(3, 0);
  std::list<std::pair<int, double>> distances;
  std::list<std::pair<double, std::vector<double>>> cl_direction;
  std::map<int, std::vector<double>> clpoints = *(mycline.GetPoints());

  // calculate mean direction of centerline to determine fiber direction
  // loop over all points of the centerline to find the two nearest
  for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
       cl_iter != clpoints.end(); ++cl_iter)
  {
    cl_direction.push_back(
        std::make_pair(EXODUS::distance3d((*clpoints.begin()).second, cl_iter->second),
            EXODUS::substract3d((*clpoints.begin()).second, cl_iter->second)));
  }
  cl_direction.sort();
  mean_cl_dir = cl_direction.back().second;
  normalize3d(mean_cl_dir);

  // loop over all midpoints of all elements
  for (std::map<int, std::vector<double>>::const_iterator el_iter = midpoints.begin();
       el_iter != midpoints.end(); ++el_iter)
  {
    distances.clear();

    // loop over all points of the centerline to find the two nearest
    for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
         cl_iter != clpoints.end(); ++cl_iter)
    {
      distances.push_back(
          std::make_pair(cl_iter->first, EXODUS::distance3d(el_iter->second, cl_iter->second)));
    }
    // sort distance map not by key but by value of second entry (here distance)
    distances.sort(MyDataSortPredicate);

    // get id of closest cl point
    clID = distances.front().first;
    // remove first element form list
    distances.pop_front();
    // get id of second closest cl point
    clID_2 = distances.front().first;

    ids[0] = clID;
    ids[1] = clID_2;
    conn_mp_cp[el_iter->first] = ids;
  }

  // in this section the three directions of all local coordinate systems are calculated
  // with the aid of conn_mp_cp
  // map that will be returned containing (eblock-ID, element-ID, directions of local coordinate
  // systems)
  std::map<int, std::map<int, std::vector<std::vector<double>>>> ebID_elID_local_cosy;

  std::vector<double> r_0(3, 0);
  std::vector<double> r_1(3, 0);
  std::vector<double> r_2(3, 0);
  std::vector<double> r_3(3, 0);
  std::vector<double> r_4(3, 0);
  std::vector<std::vector<double>> directions;

  // loop over conn_mp_cp
  for (std::map<int, std::vector<int>>::const_iterator it = conn_mp_cp.begin();
       it != conn_mp_cp.end(); ++it)
  {
    directions.clear();
    // position vector from centerline point 1 to midpoint of element
    r_0 = EXODUS::substract3d(
        midpoints.find(it->first)->second, mycline.GetPoints()->find(it->second[0])->second);
    // introduce vector to check normal direction
    r_4 = r_0;
    normalize3d(r_0);
    // take the normal vector calculated from surface nodes if existent
    std::pair<int, int> eb_el = mp_eb_el.find(it->first)->second;
    std::map<std::pair<int, int>, std::vector<double>>::iterator normal;
    normal = np_eb_el.find(eb_el);

    if (np_eb_el.find(eb_el) != np_eb_el.end())
    {
      r_0 = np_eb_el.find(eb_el)->second;
      // std::cout << RED << r_0[0] << r_0[1] << r_0[2] << END_COLOR << std::endl;
    }

    // position vector from centerline point 1 to centerline point 2 (axial direction)
    r_1 = EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,
        mycline.GetPoints()->find(it->second[0])->second);
    normalize3d(r_1);

    // check direction of fiber with mean_cl_dir by scalar product
    if (EXODUS::scalar_product3d(mean_cl_dir, r_1) < 0.0)
    {
      r_1[0] = -r_1[0];
      r_1[1] = -r_1[1];
      r_1[2] = -r_1[2];
    }

    // r_2 = r_0 x r_1 (calculate circumferential direction)
    r_2 = EXODUS::cross_product3d(r_0, r_1);
    normalize3d(r_2);

    // r_3 = r_1 x r_2 (radial direction)
    r_3 = EXODUS::cross_product3d(r_1, r_2);
    normalize3d(r_3);
    // take calculated normal vector instead
    if (ele_normals.find(it->first) != ele_normals.end())
    {
      r_3[0] = ele_normals.find(it->first)->second[0];
      r_3[1] = ele_normals.find(it->first)->second[1];
      r_3[2] = ele_normals.find(it->first)->second[2];
    }

    directions.push_back(r_3);
    directions.push_back(r_1);
    directions.push_back(r_2);

    // ebID_elID_local_cosy(ebID,elID,directions)
    // std::pair<int,int> eb_el = mp_eb_el.find(it->first)->second;
    ebID_elID_local_cosy[eb_el.first][eb_el.second] = directions;
  }
  return ebID_elID_local_cosy;
}

// function for sorting pairs by value of -> second
bool MyDataSortPredicate(std::pair<int, double> lhs, std::pair<int, double> rhs)
{
  return lhs.second < rhs.second;
}


/*------------------------------------------------------------------------*
 |- creates local coordinate systems for each element             SP 06/08|
 |- returns a map containing calc. directions referred to each element    |
 *------------------------------------------------------------------------*/
std::map<int, std::map<int, std::vector<std::vector<double>>>> EXODUS::element_cosys(
    EXODUS::Centerline& mycline, const EXODUS::Mesh& mymesh, const std::vector<int>& eb_ids)
{
  std::map<int, std::vector<double>> midpoints;  // here midpoints are stored
  // mp_eb_el contains (midpoint-ID, (eblock-ID, element-ID))
  std::map<int, std::pair<int, int>> mp_eb_el = mymesh.createMidpoints(midpoints, eb_ids);
  // conn_mp_cp will contain (midpoint-ID, centerpoint-ID_1, centerpoint-ID_2)
  std::map<int, std::vector<int>> conn_mp_cp;
  // auxiliary variables
  int clID = -1;
  int clID_2 = -1;
  std::vector<int> ids(2, 0);
  double min_distance, temp;

  std::map<int, std::vector<double>> clpoints = *(mycline.GetPoints());

  // this search should later be replaced by a nice search-tree!
  // in this section for each element the nearest point on the centerline is searched
  // and the ids of each element midpoint and the accordant centerline points are stored
  //
  // loop over all midpoints of all elements
  for (std::map<int, std::vector<double>>::const_iterator el_iter = midpoints.begin();
       el_iter != midpoints.end(); ++el_iter)
  {
    min_distance = -1;

    // loop over all points of the centerline to find nearest
    for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
         cl_iter != clpoints.end(); ++cl_iter)
    {
      temp = EXODUS::distance3d(el_iter->second, cl_iter->second);

      if (min_distance == -1)  // just for the first step
      {
        min_distance = temp;
        clID = cl_iter->first;
      }
      else
      {
        if (min_distance > temp)
        {
          min_distance = temp;
          clID = cl_iter->first;
        }
      }
    }
    // storage of IDs in conn_mp_cp
    if (clID == ((int)mycline.GetPoints()->size()) - 1)
      clID_2 = clID - 1;
    else
      clID_2 = clID + 1;
    ids[0] = clID;
    ids[1] = clID_2;
    conn_mp_cp[el_iter->first] = ids;
  }

  // in this section the three directions of all local coordinate systems are calculated
  // with the aid of conn_mp_cp
  //
  // map that will be returned containing (eblock-ID, element-ID, directions of local coordinate
  // systems)
  std::map<int, std::map<int, std::vector<std::vector<double>>>> ebID_elID_local_cosy;

  std::vector<double> r_0(3, 0);
  std::vector<double> r_1(3, 0);
  std::vector<double> r_2(3, 0);
  std::vector<double> r_3(3, 0);
  std::vector<std::vector<double>> directions;

  // loop over conn_mp_cp
  for (std::map<int, std::vector<int>>::const_iterator it = conn_mp_cp.begin();
       it != conn_mp_cp.end(); ++it)
  {
    directions.clear();

    // position vector from centerline point 1 to midpoint of element
    r_0 = EXODUS::substract3d(
        midpoints.find(it->first)->second, mycline.GetPoints()->find(it->second[0])->second);
    normalize3d(r_0);
    // position vector from centerline point 1 to centerline point 2 (axial direction)
    r_1 = EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,
        mycline.GetPoints()->find(it->second[0])->second);
    normalize3d(r_1);

    // if last CLPoint has been reached
    if (it->second[0] == ((int)mycline.GetPoints()->size()) - 1)
    {
      r_1[0] = -r_1[0];
      r_1[1] = -r_1[1];
      r_1[2] = -r_1[2];
    }

    // r_2 = r_0 x r_1 (circumferential direction)
    r_2 = EXODUS::cross_product3d(r_0, r_1);
    normalize3d(r_2);

    // r_3 = r_1 x r_2 (radial direction)
    r_3 = EXODUS::cross_product3d(r_1, r_2);
    normalize3d(r_3);

    // directions = {r_3,r_1,r_2}
    directions.push_back(r_3);
    directions.push_back(r_1);
    directions.push_back(r_2);

    // ebID_elID_local_cosy(ebID,elID,directions)
    std::pair<int, int> eb_el = mp_eb_el.find(it->first)->second;
    ebID_elID_local_cosy[eb_el.first][eb_el.second] = directions;
  }
  return ebID_elID_local_cosy;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
std::map<int, std::map<int, std::vector<std::vector<double>>>> EXODUS::element_degcosys(
    EXODUS::Centerline& mycline, const EXODUS::Mesh& mymesh, const std::vector<int>& eb_ids)
{
  std::map<int, std::vector<double>> midpoints;  // here midpoints are stored
  // mp_eb_el contains (midpoint-ID, (eblock-ID, element-ID))
  std::map<int, std::pair<int, int>> mp_eb_el = mymesh.createMidpoints(midpoints, eb_ids);
  // conn_mp_cp will contain (midpoint-ID, centerpoint-ID)
  std::map<int, int> conn_mp_cp;
  // auxiliary variables
  int clID = -1;
  double min_distance, temp;

  std::map<int, std::vector<double>> clpoints = *(mycline.GetPoints());

  // this search should later be replaced by a nice search-tree!
  // in this section for each element the nearest point on the centerline is searched
  // and the ids of each element midpoint and the accordant centerline points are stored
  //
  // loop over all midpoints of all elements
  for (std::map<int, std::vector<double>>::const_iterator el_iter = midpoints.begin();
       el_iter != midpoints.end(); ++el_iter)
  {
    min_distance = -1;

    // loop over all points of the centerline to find nearest
    for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
         cl_iter != clpoints.end(); ++cl_iter)
    {
      temp = EXODUS::distance3d(el_iter->second, cl_iter->second);

      if (min_distance == -1)  // just for the first step
      {
        min_distance = temp;
        clID = cl_iter->first;
      }
      else
      {
        if (min_distance > temp)
        {
          min_distance = temp;
          clID = cl_iter->first;
        }
      }
    }
    conn_mp_cp[el_iter->first] = clID;
  }

  // in this section the three directions of all local coordinate systems are calculated
  // with the aid of conn_mp_cp
  //
  // map that will be returned containing (eblock-ID, (element-ID, directions of local coordinate
  // systems))
  std::map<int, std::map<int, std::vector<std::vector<double>>>> ebID_elID_local_cosy;

  std::vector<double> r_0(3, 0);
  //  std::vector<double> r_1(3,0);
  //  std::vector<double> r_2(3,0);
  //  std::vector<double> r_3(3,0);
  std::vector<std::vector<double>> directions;

  // loop over conn_mp_cp
  for (std::map<int, int>::const_iterator it = conn_mp_cp.begin(); it != conn_mp_cp.end(); ++it)
  {
    directions.clear();

    // position vector from centerline point 1 to midpoint of element
    r_0 = EXODUS::substract3d(
        midpoints.find(it->first)->second, mycline.GetPoints()->find(it->second)->second);
    normalize3d(r_0);
    //    //position vector from centerline point 1 to centerline point 2 (axial direction)
    //    r_1 =
    //    EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,mycline.GetPoints()->find(it->second[0])->second);
    //    normalize3d(r_1);
    //
    //    //if last CLPoint has been reached
    //    if (it->second[0] == ((int) mycline.GetPoints()->size())-1 )
    //    {
    //      r_1[0]=-r_1[0];
    //      r_1[1]=-r_1[1];
    //      r_1[2]=-r_1[2];
    //    }
    //
    //    //r_2 = r_0 x r_1 (circumferential direction)
    //    r_2 = EXODUS::cross_product3d(r_0,r_1);
    //    normalize3d(r_2);
    //
    //    //r_3 = r_1 x r_2 (radial direction)
    //    r_3 = EXODUS::cross_product3d(r_1,r_2);
    //    normalize3d(r_3);

    // directions = {r_0,r_0,r_0}
    directions.push_back(r_0);
    directions.push_back(r_0);
    directions.push_back(r_0);

    // ebID_elID_local_cosy(ebID,elID,directions)
    std::pair<int, int> eb_el = mp_eb_el.find(it->first)->second;
    ebID_elID_local_cosy[eb_el.first][eb_el.second] = directions;
  }
  return ebID_elID_local_cosy;
}

/*------------------------------------------------------------------------*
 |- creates local coordinate systems for each element like element_cosys   |
 |- generates gmsh-file to visualize coordinate systems            SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::PlotCosys(
    EXODUS::Centerline& mycline, const EXODUS::Mesh& mymesh, const std::vector<int>& eb_ids)
{
  std::map<int, std::vector<double>> midpoints;  // here midpoints are stored
  // mp_eb_el contains (midpoint-ID, eblock-ID, element-ID)
  std::map<int, std::pair<int, int>> mp_eb_el = mymesh.createMidpoints(midpoints, eb_ids);
  // conn_mp_cp will contain (midpoint-ID, centerpoint-ID_1, centerpoint-ID_2)
  std::map<int, std::vector<int>> conn_mp_cp;
  // auxiliary variables
  int clID = -1;
  int clID_2 = -1;
  std::vector<int> ids(2, 0);
  double min_distance, temp;

  std::map<int, std::vector<double>> clpoints = *(mycline.GetPoints());

  // in this section for each element the nearest point on the centerline is searched
  // and the ids of each element midpoint and the accordant centerline points are stored
  //
  // loop over all midpoints of all elements
  for (std::map<int, std::vector<double>>::const_iterator el_iter = midpoints.begin();
       el_iter != midpoints.end(); ++el_iter)
  {
    min_distance = -1;

    // loop over all points of the centerline
    for (std::map<int, std::vector<double>>::const_iterator cl_iter = clpoints.begin();
         cl_iter != clpoints.end(); ++cl_iter)
    {
      temp = EXODUS::distance3d(el_iter->second, cl_iter->second);

      if (min_distance == -1)  // just for the first step
      {
        min_distance = temp;
        clID = cl_iter->first;
      }
      else
      {
        if (min_distance > temp)
        {
          min_distance = temp;
          clID = cl_iter->first;
        }
      }
    }
    // storage of IDs in conn_mp_cp
    if (clID == ((int)mycline.GetPoints()->size()) - 1)
      clID_2 = clID - 1;
    else
      clID_2 = clID + 1;
    ids[0] = clID;
    ids[1] = clID_2;
    conn_mp_cp[el_iter->first] = ids;
  }

  // in this section the three directions of all local coordinate systems are calculated
  // with the aid of conn_mp_cp
  // auxiliary variables
  std::vector<double> r_0(3, 0);
  std::vector<double> r_1(3, 0);
  std::vector<double> r_2(3, 0);
  std::vector<double> r_3(3, 0);
  std::vector<std::vector<double>> directions;
  std::map<int, std::vector<std::vector<double>>> mpID_directions;

  // loop over conn_mp_cp
  for (std::map<int, std::vector<int>>::const_iterator it = conn_mp_cp.begin();
       it != conn_mp_cp.end(); ++it)
  {
    directions.clear();

    // position vector from centerline point 1 to midpoint of element
    r_0 = EXODUS::substract3d(
        midpoints.find(it->first)->second, mycline.GetPoints()->find(it->second[0])->second);
    normalize3d(r_0);
    // position vector from centerline point 1 to centerline point 2 (axial direction)
    r_1 = EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,
        mycline.GetPoints()->find(it->second[0])->second);
    normalize3d(r_1);

    // if last CLPoint has been reached
    if (it->second[0] == ((int)mycline.GetPoints()->size()) - 1)
    {
      r_1[0] = -r_1[0];
      r_1[1] = -r_1[1];
      r_1[2] = -r_1[2];
    }

    // r_2 = r_0 x r_1 (circumferential direction)
    r_2 = EXODUS::cross_product3d(r_0, r_1);
    normalize3d(r_2);

    // r_3 = r_1 x r_2 (radial direction)
    r_3 = EXODUS::cross_product3d(r_1, r_2);
    normalize3d(r_3);

    // directions = {r_3,r_1,r_2}
    directions.push_back(r_3);
    directions.push_back(r_1);
    directions.push_back(r_2);

    mpID_directions[it->first] = directions;
  }

  //
  // generation of gmsh-file to visualize local coordinate systems
  //
  std::ofstream gmshFile("local_coordinate_systems.gmsh");

  gmshFile << "View \" local coordinate systems \" {" << std::endl;

  for (std::map<int, std::vector<std::vector<double>>>::iterator iti = mpID_directions.begin();
       iti != mpID_directions.end(); ++iti)
  {
    std::vector<double> mp = midpoints.find(iti->first)->second;
    std::vector<double> r1 = iti->second[0];
    std::vector<double> r2 = iti->second[1];
    std::vector<double> r3 = iti->second[2];
    gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << r1[0] << "," << r1[1]
             << "," << r1[2] << "};" << std::endl;
    gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << 2. * r2[0] << ","
             << 2. * r2[1] << "," << 2. * r2[2] << "};" << std::endl;
    gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << 3. * r3[0] << ","
             << 3. * r3[1] << "," << 3. * r3[2] << "};" << std::endl;
    //    //VL(mp,mp,mp,mp+r_3,mp+r_3,mp+r_3){1,1,1,1,1,1};
    //    gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," <<
    //    midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] <<
    //    "," << midpoints.find(iti->first)->second[0] + iti->second[0][0] << "," <<
    //    midpoints.find(iti->first)->second[1] + iti->second[0][1] << "," <<
    //    midpoints.find(iti->first)->second[2] + iti->second[0][2] << "){1,1,1,1,1,1};" <<
    //    std::endl;
    //
    //    //VL(mp,mp,mp,mp+r_1,mp+r_1,mp+r_1){2,2,2,2,2,2};
    //    gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," <<
    //    midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] <<
    //    "," << midpoints.find(iti->first)->second[0] + iti->second[1][0] << "," <<
    //    midpoints.find(iti->first)->second[1] + iti->second[1][1] << "," <<
    //    midpoints.find(iti->first)->second[2] + iti->second[1][2] << "){2,2,2,2,2,2};" <<
    //    std::endl;
    //
    //    //VL(mp,mp,mp,mp+r_2,mp+r_2,mp+r_2){3,3,3,3,3,3};
    //    gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," <<
    //    midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] <<
    //    "," << midpoints.find(iti->first)->second[0] + iti->second[2][0] << "," <<
    //    midpoints.find(iti->first)->second[1] + iti->second[2][1] << "," <<
    //    midpoints.find(iti->first)->second[2] + iti->second[2][2] << "){3,3,3,3,3,3};" <<
    //    std::endl;
  }
  gmshFile << "};";
  gmshFile.close();

  return;
}
