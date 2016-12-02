/*!-----------------------------------------------------------------------------------------------*
\file tetrahedradecomposition.cpp

  \brief tetrahedralization procedure based on level set function

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

<pre>
\level 3
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "element_coordtrafo.H"
#include "integrationcell.H"
#include "tetrahedradecomposition.H"
#include "intersection_service_templates.H"
#include "../drt_combust/combust_refinementcell.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/standardtypes_cpp.H"


/*----------------------------------------------------------------------------------*
 | Decompose into tetrahedra and intersect with interface            sons 01/11     |
 | (entry point for this class)                                                     |
 *----------------------------------------------------------------------------------*/
GEO::TetrahedraDecomposition::TetrahedraDecomposition(
    const COMBUST::RefinementCell* cell,
    GEO::BoundaryIntCells & listBoundaryIntCellsperEle,
    GEO::DomainIntCells & listDomainIntCellsperEle
) :
cell_(cell),
num_degeneratedTets_(0)
{
  // vector for created tetrahedra
  std::vector<GEO::DomainIntCell> tets;
  // gvalues for tets
  LINALG::SerialDenseMatrix tets_gvalues(6, 4);
  // decompose refinementcells into 6 tetrahedra and add to tets
  DecomposeIntoTetrahedra(tets, tets_gvalues);
  // intersect 6 tetrahedra with interface
  IntersectTetrahdera(tets, listBoundaryIntCellsperEle, tets_gvalues);
  // add saved tetrahedra to listDomainIntCellsperEle
  for(std::size_t i = 0;i < tets.size();i++){
    listDomainIntCellsperEle.push_back(tets[i]);
  }
  // delete temporary tetrahedra
  tets.clear();

  if(num_degeneratedTets_ != 0)
    std::cout << "Warning: " << num_degeneratedTets_ << " degenerated tet" << (num_degeneratedTets_>1 ? "s" : "") << "!" << std::endl;
}


/*----------------------------------------------------------------------------------*
 |Decompose into tetrahedra                                          sons 09/10     |
 *----------------------------------------------------------------------------------*/
void GEO::TetrahedraDecomposition::DecomposeIntoTetrahedra(std::vector<GEO::DomainIntCell>& tets, LINALG::SerialDenseMatrix& tets_gvalues)
{
  unsigned int i = 0;
  unsigned int j = 0;

  // Get element vertices
  const std::vector<std::vector<double> >& vertices = cell_->GetVertexCoord();

  // get G-function values at vertices from refinement cell
  const std::vector<double>& gfuncvalues = cell_->GetGfuncValues();
  std::vector<double> gfuncvalues_tet;


  // Get element and nodes
  const DRT::Element* ele = cell_->Ele();
  const DRT::Node** nodes = (const DRT::Node**)ele->Nodes();


  int num_nodes = ele->NumNode();

  LINALG::SerialDenseMatrix tet_coords(3, num_nodes);
  LINALG::SerialDenseMatrix tet_pcoords(3, num_nodes);

  // There are multiple ways to do that, this is one example
  const int tetindex[6][4] = {  {0, 5, 7, 4},
      {0, 1, 7, 5},
      {1, 6, 7, 5},
      {0, 7, 2, 3},
      {0, 7, 1, 2},
      {1, 7, 6, 2} };
  //  const int tetindex[6][4] = {  {0, 2, 3, 6},
  //                                {0, 3, 7, 6},
  //                                {0, 7, 4, 6},
  //                                {0, 5, 6, 4},
  //                                {1, 5, 6, 0},
  //                                {1, 6, 2, 0} };
  //  const int tetindex[5][4] = {  {0, 1, 2, 5},
  //                                {0, 2, 7, 5},
  //                                {0, 2, 3, 7},
  //                                {0, 5, 7, 4},
  //                                {2, 7, 5, 6}};

  // Divide hex8 element into 6 tet4 elements
  for(i = 0; i < 6; i++)
  {
    gfuncvalues_tet.clear();

    for(j = 0; j < 4; j++)
    {
      for(int k = 0; k < 3; k++)
      {
        tet_coords(k,j) = vertices[tetindex[i][j]][k];
        tet_pcoords(k,j) = nodes[tetindex[i][j]]->X()[k];
        tets_gvalues(i,j) = gfuncvalues[tetindex[i][j]];
        gfuncvalues_tet.push_back(gfuncvalues[tetindex[i][j]]);
      }


    }

    bool inGplus = GEO::TetrahedraDecomposition::GetIntCellDomainInElement(tet_coords,
        gfuncvalues_tet, DRT::Element::hex8, DRT::Element::tet4);

    if(!GEO::checkDegenerateTet(0, tet_coords, tet_pcoords))
      tets.push_back(GEO::DomainIntCell(DRT::Element::tet4, tet_coords, tet_pcoords, inGplus));

  }
}

/*----------------------------------------------------------------------------------*
 | Create tetrahedra (domain integration cells)                          sons 01/11 |
 *----------------------------------------------------------------------------------*/
void GEO::TetrahedraDecomposition::createTetrahedra(
    unsigned int num_tetrahedra,
    LINALG::SerialDenseMatrix & tet_coords,
    LINALG::SerialDenseMatrix & vertices,
    const int* tetindex,
    LINALG::SerialDenseMatrix & tet_pcoords,
    LINALG::SerialDenseMatrix & nodes, bool & intersected,
    std::vector<GEO::DomainIntCell> & tets,
    int & sign_node_four,
    bool tetgplus[2]
)
{
  for(unsigned int j = 0;j < num_tetrahedra;j++)
  {
    for(unsigned int k = 0;k < 4;k++)
    {
      for(unsigned int d = 0;d < 3;d++)
      {
        tet_coords(d, k) = vertices(d, tetindex[j*4+k]);
        tet_pcoords(d, k) = nodes(d, tetindex[j*4+k]);
      }
    }

    // set intersected flag, so we know we must delete the original tetrahedra
    intersected = true;

    // Intersection with line 0,2,3 will result in wrong tet numbering, checkDegenerateTet will repair it
    if(!GEO::checkDegenerateTet(0, tet_coords, tet_pcoords)){
      tets.push_back(GEO::DomainIntCell(DRT::Element::tet4, tet_coords, tet_pcoords, (sign_node_four == -1) ? !tetgplus[j] : tetgplus[j]));
    }else{
      //std::cout << "WARNING: degenerated tet4!" << std::endl;
      intersected = false;
      num_degeneratedTets_++;
    }
  }

}
/*----------------------------------------------------------------------------------*
 | Create triangles (boundary integration cells)                     sons 01/11     |
 *----------------------------------------------------------------------------------*/
void GEO::TetrahedraDecomposition::createTriangles(unsigned int num_triangles, LINALG::SerialDenseMatrix trianglecoord, LINALG::SerialDenseMatrix vertices, const int* triindex, LINALG::SerialDenseMatrix phystrianglecoord, LINALG::SerialDenseMatrix nodes, LINALG::Matrix<3,1> temp, int comparison_vertex, int sign_node_four, GEO::BoundaryIntCells & tris, bool intersected)
{
  if(!intersected)
  {
    intersected = true;
    return;
  }


  for(unsigned int k = 0;k < num_triangles;k++)
  {
    // create and add boundary cells (2 tris)
    for(unsigned int j = 0;j < 3;j++)
    {
      for(unsigned int dim = 0;dim < 3;dim++)
      {
        trianglecoord(dim, j) = vertices(dim, triindex[k*3+j]);
        phystrianglecoord(dim, j) = nodes(dim, triindex[k*3+j]);
      }
    }

    for(unsigned int dim = 0;dim < 3;dim++)
      temp(dim) = nodes(dim, comparison_vertex);

    checkTriangleOrientation(trianglecoord, phystrianglecoord, temp, sign_node_four);
    tris.push_back(GEO::BoundaryIntCell(DRT::Element::tri3, -1, trianglecoord, Teuchos::null, phystrianglecoord, true));
  }
}

/*----------------------------------------------------------------------------------*
 |Intersect tetrahedra with interface                                sons 09/10     |
 *----------------------------------------------------------------------------------*/
void GEO::TetrahedraDecomposition::IntersectTetrahdera(std::vector<GEO::DomainIntCell>& tets, GEO::BoundaryIntCells& tris, LINALG::SerialDenseMatrix& tets_gvalues)
{
  std::map<int, std::vector<double> > ip;
  std::map<int, std::vector<double> >::iterator it;
  std::vector<double> tet_gvalues;

  LINALG::SerialDenseMatrix vertices(3, 8);
  LINALG::SerialDenseMatrix nodes(3, 8);

  LINALG::SerialDenseMatrix trianglecoord(3,3);
  LINALG::SerialDenseMatrix phystrianglecoord(3,3);

  std::vector<int> intersected_tets;

  std::vector<double> tempVertex;

  LINALG::Matrix<3,1> temp;
  LINALG::Matrix<3,1> tetcoord;

  int sign_node_four = 0;
  int num_tets = 0;

  bool intersected = false;

  std::vector<double> ip_coords;

  // 4 nodes, three dimensions
  LINALG::SerialDenseMatrix tet_coords(3, 4);
  LINALG::SerialDenseMatrix tet_pcoords(3, 4);

  // get number of tets (decomposition may lead to 5/6 tets)
  num_tets = tets.size();

  // loop over all 6 tets
  for(int i = 0; i < num_tets; i++)
  {
    tet_gvalues.clear();
    ip_coords.clear();
    tet_coords.Zero();
    tet_pcoords.Zero();

    intersected = false;

    for(unsigned int j = 0; j < 4; j++)
      tet_gvalues.push_back(tets_gvalues(i, j));

    // find intersection points
    int num_ips = FindIntersectionPointsTet4(tets[i], tet_gvalues, vertices, nodes, sign_node_four);

    // tetrahedra is not intersected
    if(num_ips == 0) {

      // number of Zero-G-points (saved in dummy vertex (0,0)) must be 3 (touched)
      if(vertices(0,0) != 3.0)
        continue;

      int triindex[3] = {1, 2, 3};

      // create 1 boundary triangle
      createTriangles(1, trianglecoord, vertices, triindex, phystrianglecoord, nodes, temp, 7, sign_node_four, tris, true);

      continue;
    }

    // 1 intersection point
    else if(num_ips == 1) {

      const int tetindex[8] =  {0, 7, 5, 6,
          0, 5, 7, 4};

      bool tetgplus[2] = {false, true};

      // create 2 tets
      createTetrahedra(2, tet_coords, vertices, tetindex, tet_pcoords, nodes, intersected, tets, sign_node_four, tetgplus);

      const int triindex[3] = {0, 5, 7};

      // create 1 boundary triangle
      createTriangles(1, trianglecoord, vertices, triindex, phystrianglecoord, nodes, temp, 4, sign_node_four, tris, intersected);
    }

    // 2 intersection points, 1 tet + 1 pyramide = 3 tets
    else if(num_ips == 2) {
      // There are multiple ways to do that, this is one example
      const int tetindex[12] =  {0, 1, 5, 4,
          0, 1, 5, 7,
          0, 5, 7, 6};

      // g-plus/g-minus
      bool tetgplus[3] = {true, false, false};

      // create 3 tets
      createTetrahedra(3, tet_coords, vertices, tetindex, tet_pcoords, nodes, intersected, tets, sign_node_four, tetgplus);

      const int triindex[3] = {0, 1, 5};
      // create 1 boundary triangle
      createTriangles(1, trianglecoord, vertices, triindex, phystrianglecoord, nodes, temp, 4, sign_node_four, tris, intersected);
    }

    // 3 intersection points -> 1 tet + 1 wedge = 4 tets
    else if(num_ips == 3)
    {
      // There are multiple ways to do that, this is one example
      const int tetindex[16] =   {0, 1, 2, 4,
          0, 1, 2, 7,
          0, 1, 7, 5,
          0, 5, 7, 6};

      // g-plus/g-minus
      bool tetgplus[4] = {true, false, false, false};

      // create 4 tets
      createTetrahedra(4, tet_coords, vertices, tetindex, tet_pcoords, nodes, intersected, tets, sign_node_four, tetgplus);

      int triindex[3] = {0, 1, 2};

      // create 1 boundary triangle
      createTriangles(1, trianglecoord, vertices, triindex, phystrianglecoord, nodes, temp, 4, sign_node_four, tris, intersected);

    }

    // we have 4 intersection points -> wedge + wedge -> 3 tets + 3 tets = 6 tets
    else if(num_ips == 4)
    {
      // There are multiple ways to do that, this is one example
      const int tetindex[24] =  {0, 2, 4, 5,
          0, 2, 5, 3,
          0, 3, 5, 1,
          0, 1, 6, 7,
          0, 1, 7, 3,
          0, 3, 7, 2};

      // g-plus/g-minus
      bool tetgplus[6] = {true, true, true, false, false, false};

      // create and add 6 tets
      createTetrahedra(6, tet_coords, vertices, tetindex, tet_pcoords, nodes, intersected, tets, sign_node_four, tetgplus);

      const int triindex[6] = {0, 1, 2,
          1, 3, 2};
      // create 2 boundary triangles
      createTriangles(2, trianglecoord, vertices, triindex, phystrianglecoord, nodes, temp, 4, sign_node_four, tris, intersected);
    }

    else
    {
      std::cout << "Invalid number of intersection points: " << num_ips << "(Touched?)" << std::endl;
      continue;
    }

    // mark tet as intersected
    if(intersected)
      intersected_tets.push_back(i);

  }

  int counter = 0;

  //   Delete tets that have been intersected (we just created new tets for those)
  for(unsigned int i = 0; i < intersected_tets.size(); i++)
  {
    tets.erase(tets.begin()+intersected_tets[i]+counter);
    counter--;
  }

}

/*----------------------------------------------------------------------------------*
 |Find intersection points (tet4)                                    sons 09/10     |
 *----------------------------------------------------------------------------------*/
int GEO::TetrahedraDecomposition::FindIntersectionPointsTet4(GEO::DomainIntCell& tet, std::vector<double>& tet_gvalues, LINALG::SerialDenseMatrix& vertices, LINALG::SerialDenseMatrix& nodes, int& sign_node_four)
{
  //-------------------------------------------------
  // get vector of lines with corresponding vertices
  //-------------------------------------------------
  // lines: edge numbers and corresponding vertices (hex8)
  std::vector<std::vector<int> > lines;

  if(tet.Shape() != DRT::Element::tet4)
    return -1;

  int num_ips = 0;

  int ip_line[6] = {0};
  int saved_nodes[4] = {0};

  LINALG::Matrix<3,1> coord;

  LINALG::SerialDenseMatrix xyze(3,8);
  for(int inode=0;inode<8;inode++)
  {
    xyze(0,inode) = cell_->Ele()->Nodes()[inode]->X()[0];
    xyze(1,inode) = cell_->Ele()->Nodes()[inode]->X()[1];
    xyze(2,inode) = cell_->Ele()->Nodes()[inode]->X()[2];
  }

  const LINALG::SerialDenseMatrix vertexcoord = tet.CellNodalPosXiDomain();

  lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::tet4);

  //-------------------------------
  // determine intersection points
  //-------------------------------
  // loop edges of refinement cell
  for(std::size_t iline=0; iline<lines.size(); iline++)
  {
    // get G-function value of the two vertices defining an edge
    double gfuncval1 = tet_gvalues[lines[iline][0]];
    double gfuncval2 = tet_gvalues[lines[iline][1]];

    std::vector<double> coordinates(3);

    // check for change of sign along edge
    if (gfuncval1*gfuncval2 < 0.0)
    {
      ip_line[iline] = 1;

      for (int dim = 0; dim < 3; dim++)
      {
        // vertices have one coordinate in common
        if(vertexcoord(dim,lines[iline][0]) == vertexcoord(dim,lines[iline][1]))
        {
          // intersection point has the same coordinate for that direction
          coordinates[dim] = vertexcoord(dim,lines[iline][0]);
        }
        else // compute intersection point
        {
          // linear interpolation (for hex8)
          // x = x1 + (phi(=0) - phi1)/(phi2 - phi1)*(x2 - x1)
          // store intersection point coordinate (local element coordinates) for every component dim
          coordinates[dim] = vertexcoord(dim,lines[iline][0]) - gfuncval1 / (gfuncval2 - gfuncval1)
              * (vertexcoord(dim,lines[iline][1]) - vertexcoord(dim,lines[iline][0]));

        }

        // save intersection points as 0,1,2,3 (if available)
        vertices(dim, num_ips) = coordinates[dim];
        coord(dim) = coordinates[dim];

      }

      // save nodes 4 and 6
      // nodes on first intersected line
      if(num_ips == 0)
      {
        for(int dim = 0; dim < 3; dim++)
        {
          vertices(dim, 4) = vertexcoord(dim,lines[iline][0]);
          vertices(dim, 6) = vertexcoord(dim,lines[iline][1]);
        }

        // save sign of the g-value of node 4 (so we know which tet is g-plus/g-minus)
        sign_node_four = SIGN(gfuncval1);

        saved_nodes[lines[iline][0]] = 1;
        saved_nodes[lines[iline][1]] = 1;
      }
      num_ips++;
    }
    else
      ip_line[iline] = 0; // mark line as not intersected
  }

  // 0 intersection points, save number of "touched" points in vertices(0,0)
  if(num_ips == 0)
  {
    int counter = 0;

    // save remaining nodes
    for(int i = 0; i < 4; i++) {
      if(tet_gvalues[i] == 0)
      {
        for(int dim = 0; dim < 3; dim++)
          vertices(dim, counter+1) = vertexcoord(dim, i);
        counter++;
      }
      else if(tet_gvalues[i] > 0) {
        sign_node_four = 1;
        for(int dim = 0; dim < 3; dim++)
          vertices(dim, 7) = vertexcoord(dim, i);
      }
      else {
        sign_node_four = -1;
        for(int dim = 0; dim < 3; dim++)
          vertices(dim, 7) = vertexcoord(dim, i);
      }
    }
    vertices(0, 0) = counter;
  }

  // 1 intersection point
  else if(num_ips == 1)
  {
    bool first = true;
    bool debug = false;
    for(int i = 0; i < 4; i++)
    {
      // "touched" points will be 5 and 7
      if(tet_gvalues[i] == 0)
      {
        if(first)
        {
          for(int dim = 0; dim < 3; dim++)
            vertices(dim, 5) = vertexcoord(dim, i);

          first = false;
        }
        else
        {
          for(int dim = 0; dim < 3; dim++)
            vertices(dim, 7) = vertexcoord(dim, i);
          debug = true;
        }
      }
    }
    if(!debug)
      std::cout << "WARNING: Unknown case!" << std::endl;
  }

  // 2 intersection points
  else if(num_ips == 2)
  {

    bool got_node_five = false; // just for debugging

    // save remaining 2 nodes
    for(int i = 0; i < 4; i++)
    {
      if(saved_nodes[i] == 1)
        continue;

      // 6 is "single" point
      if((tet_gvalues[i] != 0) && (SIGN(tet_gvalues[i]) == sign_node_four))
      {
        // need to change sign of node 4 as we swap 4 and 6
        sign_node_four = -sign_node_four;

        for(int dim = 0; dim < 3; dim++) {
          // swap 4 and 6
          coord(dim) = vertices(dim, 4);
          vertices(dim, 4) = vertices(dim, 6);
          vertices(dim, 6) = coord(dim);
        }
      }

      // g-value == 0: node 5
      if(tet_gvalues[i] == 0) {
        for(int dim = 0; dim < 3; dim++)
          vertices(dim, 5) = vertexcoord(dim, i);

        got_node_five = true;
      }
      // != 0: node 7
      else {
        for(int dim = 0; dim < 3; dim++)
          vertices(dim, 7) = vertexcoord(dim, i);
      }
    }
    if(!got_node_five)
      std::cout << "ERROR: Didn't assign node 5" << std::endl;
  }

  // 3 intersection points
  if(num_ips == 3)
  {
    bool first = true;
    // save remaining 2 nodes
    for(int i = 0; i < 4; i++)
    {
      if(saved_nodes[i] == 1)
        continue;

      // 4 is single point
      if(SIGN(tet_gvalues[i]) != sign_node_four)
      {
        if(first)
          for(int dim = 0; dim < 3; dim++)
            vertices(dim, 5) = vertexcoord(dim, i);
        else
          for(int dim = 0; dim < 3; dim++)
            vertices(dim, 7) = vertexcoord(dim, i);
        first = false;
      }
      // 6 is single point
      else
      {
        if(first)
        {
          for(int dim = 0; dim < 3; dim++)
          {
            // need to change sign of node 4 as we swap 4 and 6
            sign_node_four = -sign_node_four;

            // swap 4 and 6
            coord(dim) = vertices(dim, 4);
            vertices(dim, 4) = vertices(dim, 6);
            vertices(dim, 6) = coord(dim);

            // save node 5
            vertices(dim, 5) = vertexcoord(dim, i);
          }
        }
        else
          for(int dim = 0; dim < 3; dim++)
            vertices(dim, 7) = vertexcoord(dim, i);
        first = false;
      }

    }
  }

  // 4 intersection points
  else if(num_ips == 4)
  {
    // save remaining 2 nodes
    for(int i = 0; i < 4; i++)
    {
      // already saved this node (4 or 6)
      if(saved_nodes[i] == 1)
        continue;

      for(int dim = 0; dim < 3; dim++)
      {
        // node 5
        if(SIGN(tet_gvalues[i]) == sign_node_four)
          vertices(dim, 5) = vertexcoord(dim, i);
        //node 7
        else
          vertices(dim, 7) = vertexcoord(dim, i);
      }
    }
    // swap points 2 and 3 when L1,2,3,4 is intersected (to get correct node orientation)
    if((ip_line[1] == 1) && (ip_line[2] == 1) && (ip_line[3] == 1) && (ip_line[4] == 1))
    {
      for(int dim = 0; dim < 3; dim++)
      {
        coord(dim) = vertices(dim, 2);
        vertices(dim, 2) = vertices(dim, 3);
        vertices(dim, 3) = coord(dim);
      }
    }
    // swap 4, 6 and 5,7
    if((ip_line[0] == 1) && (ip_line[2] == 1))
    {
      for(int dim = 0; dim < 3; dim++)
      {
        sign_node_four = -sign_node_four;
        // swap 4 and 6
        coord(dim) = vertices(dim, 4);
        vertices(dim, 4) = vertices(dim, 6);
        vertices(dim, 6) = coord(dim);

        // swap 5 and 7
        coord(dim) = vertices(dim, 5);
        vertices(dim, 5) = vertices(dim, 7);
        vertices(dim, 7) = coord(dim);

      }

    }
  }

  // calculate physical coordinates and save in nodes
  for(int i = 0; i < 8; i++)
  {
    for(int dim = 0; dim < 3; dim++)
      coord(dim) = vertices(dim, i);

    GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, coord);

    for(int dim = 0; dim < 3; dim++)
      nodes(dim, i) = coord(dim);
  }

  return num_ips;

}

/*----------------------------------------------------------------------------------*
 | Check orientation of triangles and reorientate if necessary       sons 10/10     |
 *----------------------------------------------------------------------------------*/
void GEO::TetrahedraDecomposition::checkTriangleOrientation(LINALG::SerialDenseMatrix& coords, LINALG::SerialDenseMatrix& pcoords, LINALG::Matrix<3,1>& node_four, int sign_node_four)
{
  // create planes consisting of 3 nodes each
  LINALG::Matrix<3,1> v01;
  LINALG::Matrix<3,1> v02;

  LINALG::Matrix<3,1> coord;

  for (int dim = 0; dim < 3; ++dim)  v01(dim,0) = pcoords(dim,1)-pcoords(dim,0);
  for (int dim = 0; dim < 3; ++dim)  v02(dim,0) = pcoords(dim,2)-pcoords(dim,0);

  LINALG::Matrix<3,1> nplane = GEO::computeCrossProduct(v01,v02);

  double value = 0.0;
  double d = 0.0;

  d = nplane(0)*pcoords(0,0) + nplane(1)*pcoords(1,0) + nplane(2)*pcoords(2,0);

  // check relative position of node 4
  value = nplane(0)*node_four(0) + nplane(1)*node_four(1) + nplane(2)*node_four(2) - d;

  // compare with the sign of node 4's g-value
  if(SIGN(value) == sign_node_four)
  {
    // swap nodes 0 and 1 (arbitrary)
    for(unsigned int dim = 0; dim < 3; dim++)
    {
      coord(dim) = pcoords(dim, 0);
      pcoords(dim, 0) = pcoords(dim, 1);
      pcoords(dim, 1) = coord(dim);

      coord(dim) = coords(dim, 0);
      coords(dim, 0) = coords(dim, 1);
      coords(dim, 1) = coord(dim);
    }
  }
}

/*------------------------------------------------------------------------------------------------*
 | compute average GfuncValue for integration cell                                                |
 *------------------------------------------------------------------------------------------------*/
bool GEO::TetrahedraDecomposition::GetIntCellDomainInElement(
    const LINALG::SerialDenseMatrix&       IntCellCoord,
    const std::vector<double>&             gfuncvalues_ele,
    const DRT::Element::DiscretizationType xfem_distype,
    const DRT::Element::DiscretizationType cell_distype
)
{
  /*------------------------------------------------------------------------------
   * - compute phi at each node of domain integration cell
   * - compute average G-value for each domain integration cell
   * -> >=0 in plus == true
   * _> <0  in minus == false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;

  if (cell_distype != DRT::Element::tet4) dserror("Tetrahedra expected!");
  const int numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;

  int numelenodes = 0;
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  }
  else {
    dserror("Discretization Type (Ele) not supported yet!");
  }

  std::vector<double> gvalcellnodes (numcellnodes);
  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    gvalcellnodes[icellnode] = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
  {
    //calculate shape function at IntCell node
    Epetra_SerialDenseVector  funct(numelenodes);
    DRT::UTILS::shape_function_3D(funct,IntCellCoord(0,icellnode),IntCellCoord(1,icellnode),IntCellCoord(2,icellnode),DRT::Element::hex8);
    /*
     *         ___
     * gval(x)=\   N(x)*gvali
     *         /
     *         ___
     */
    for (int ielenode=0; ielenode<numelenodes; ielenode++)
    {
      gvalcellnodes[icellnode] += funct(ielenode) * gfuncvalues_ele[ielenode];
    }
  }

  //calculate average Gfunc value
  double averageGvalue = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    averageGvalue += gvalcellnodes[icellnode];

  if (numcellnodes == 0.0)
    dserror("division by zero: number of vertices per cell is 0!");
  averageGvalue /= numcellnodes;

  //determine DomainIntCell position
  // remark: '>=' because of the approximation of the interface,
  // also averageGvalue=0 is possible; as always approximation errors are made
  if(averageGvalue>=0.0)
    inGplus = true;

  return inGplus;
}

