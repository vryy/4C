/*!----------------------------------------------------------------------
\file combust3_surface_evaluate.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element (tri or quad)

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "combust3_sysmat.H"
#include "combust_refinementcell.H"
#include "../drt_cut/cut_position.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/element_normals.H"
#include "../drt_geometry/position_array.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3Surface::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
    DRT::ELEMENTS::Combust3Surface::ActionType act = Combust3Surface::none;
    std::string action = params.get<std::string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action == "integrate_Shapefunction")
        act = Combust3Surface::integrate_Shapefunction;
    else if (action == "calc_flux")
        act = Combust3Surface::calc_flux;
    else if (action == "calc_Neumann_inflow")
        act = Combust3Surface::calc_Neumann_inflow;
    else dserror("Unknown type of action for Combust3_Surface");

    switch(act)
    {
      case integrate_Shapefunction:
      {
        Teuchos::RCP<const Epetra_Vector> dispnp;
        std::vector<double> mydispnp(lm.size(),0.0);

        dispnp = discretization.GetState("dispnp");
        if (dispnp!=Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
        IntegrateShapeFunction(params,discretization,lm,elevec1,mydispnp);
        break;
      }
      case calc_flux:
      {
        const Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
        IntegrateSurfaceFlow(params,discretization,lm,elevec1,myvelnp);
        break;
      }
      case calc_Neumann_inflow:
      {
        //std::cout << "ich werte jetzt den Neumann inflow Term aus!" << std::endl;
        //this->Print(std::cout);
        Epetra_Vector* phinp = ParentElement()->Phinp();

        if (this->ParentElement()->Shape() != DRT::Element::hex8)
          dserror("Neumann inflow term evaluation only implemented for hex8 elements");

        // remark: surface (2D) elements have been build on entering this function
        //         according to Neumann inflow condition on nodes

        //---------------------------------------------------------------
        // generate boundary integration cells at Neumann inflow boundary
        //---------------------------------------------------------------
        // list of boundary integration cells
        GEO::BoundaryIntCells surfaceintcelllist;
        {
          // vector holding G-function values for this surface element
          std::vector<double> gfuncvalues_surfcell;
          // extract local (element level) G-function values from global vector
          DRT::UTILS::ExtractMyNodeBasedValues(this, gfuncvalues_surfcell, *phinp);

          // vector holding G-function values for the parent (hex8) element
          std::vector<double> gfuncvalues_parent;
          // extract local (element level) G-function values from global vector
          DRT::UTILS::ExtractMyNodeBasedValues(this->ParentElement(), gfuncvalues_parent, *phinp);

          // get node coordinates of the parent (3D) element
          const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
          LINALG::Matrix<3,numnode> xyze;
          GEO::fillInitialPositionArray<DRT::Element::hex8>(this->ParentElement(), xyze);

          // get node coordinates of the surface (2D) element
          std::vector<std::vector<double> > xicoord(4, std::vector<double>(3, 0.0));

          size_t numnodesurf = this->NumNode();
          if (numnodesurf != 4)
            dserror("surfaces of hex8 elements expected (4 nodes)");

          // loop nodes of surface element
          const DRT::Node*const* nodessurf = this->Nodes();
          for (size_t inode=0; inode<numnodesurf; inode++)
          {
            // compute local element coordinates
            // remark: since we deal with nodes, we should get clean numbers (-1 or 1) for 'xicoord'
            const double* x = nodessurf[inode]->X();
            {
              LINALG::Matrix<3,1> xyznode;
              xyznode(0) = x[0];
              xyznode(1) = x[1];
              xyznode(2) = x[2];

              LINALG::Matrix<3,1> xsi;
              GEO::CUT::Position<DRT::Element::hex8> pos(xyze,xyznode);
              pos.Compute();
              xsi = pos.LocalCoordinates();
              // fill array
              xicoord[inode][0] = xsi(0);
              xicoord[inode][1] = xsi(1);
              xicoord[inode][2] = xsi(2);
            }
          }
          // create refinement cell from a fluid element -> cell will have same geometry as element!
          const Teuchos::RCP<COMBUST::RefinementCell> surfcell = Teuchos::rcp(new COMBUST::RefinementCell(ParentElement(), DRT::Element::quad4, xicoord));

          // reset G-function values close to 0 with in tolerance 1.0e-10
          // remark: facilitates handling of special cases
          for (size_t i=0; i<gfuncvalues_surfcell.size(); i++ )
          {
            if (fabs(gfuncvalues_surfcell[i])<1.0E-10)
            {
              gfuncvalues_surfcell[i] = 0.0;
              //std::cout << " G-Function value  reset to 0 " << std::endl;
            }
          }

          surfcell->SetGfuncValues(gfuncvalues_surfcell);

          // get vertex coordinates (local fluid element coordinates) from refinement cell
          const std::vector<std::vector<double> >& vertexcoord = surfcell->GetVertexCoord();

          if(surfcell->Bisected() == true)
          {
            // determine number of intersection points
            {
              // get G-function values at vertices from refinement cell
              const std::vector<double>& gfuncvalues = surfcell->GetGfuncValues();
              // get vertex coordinates (local fluid element coordinates) from refinement cell
              const std::vector<std::vector<double> >& vertexcoord = surfcell->GetVertexCoord();
              // temporary variable to store intersection points
              std::multimap<int,std::vector<double> > intersectionpoints;

              //-------------------------------------------------
              // get vector of lines with corresponding vertices
              //-------------------------------------------------
              // lines: edge numbers and corresponding vertices (hex8)
              std::vector<std::vector<int> > lines;
              //-------------------------------
              // determine intersection points
              //-------------------------------
              lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::quad4);
              // loop edges of refinement cell
              for(std::size_t iline=0; iline<lines.size(); iline++)
              {
                // get G-function value of the two vertices defining an edge
                double gfuncval1 = gfuncvalues[lines[iline][0]];
                double gfuncval2 = gfuncvalues[lines[iline][1]];

                std::vector<double> coordinates(3);

                // check for change of sign along edge
                if (gfuncval1*gfuncval2 < 0.0)
                {
                  for (int dim = 0; dim < 3; dim++)
                  {
                    // vertices have one coordinate in common
                    if(vertexcoord[lines[iline][0]][dim] == vertexcoord[lines[iline][1]][dim])
                    {
                      // intersection point has the same coordinate for that direction
                      coordinates[dim] = vertexcoord[lines[iline][0]][dim];
                    }
                    else // compute intersection point
                    {
                      // linear interpolation (for hex8)
                      // x = x1 + (phi(=0) - phi1)/(phi2 - phi1)*(x2 - x1)
                      // store intersection point coordinate (local element coordinates) for every component dim
                      coordinates[dim] = vertexcoord[lines[iline][0]][dim] - gfuncval1 / (gfuncval2 - gfuncval1)
                                * (vertexcoord[lines[iline][1]][dim] - vertexcoord[lines[iline][0]][dim]);
                    }
                  }

                  // store coordinates of intersection point for each line
                  intersectionpoints.insert(std::pair<int, std::vector<double> >( iline, coordinates));
                  //intersectionpoints[iline] = coordinates;
                }
                else
                {
                  //do nothing and go to the next line
                }
              }

              //TEST Ausgabe
              //for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
              //{
              //  std::cout<< iter->first << std::endl;
              //  std::vector<double> coord = iter->second;
              //  for (std::size_t isd=0; isd<3; isd++)
              //  {
              //    std::cout<< coord[isd] << std::endl;
              //  }
              //}

              // store intersection points in refinement cell
              surfcell->intersectionpoints_ = intersectionpoints;
            }
            // generate flame front (interface) geometry
            if (surfcell->intersectionpoints_.size()==2)
            {
              //-----------------------------------------------------
              // prepare lists, get coordinates and G-function values
              //-----------------------------------------------------
              std::vector<std::vector<double> >    pointlist;
              std::vector<std::vector<int> >       trianglelist;

              //---------------------------------------------
              // fill list of points and intersection points
              //---------------------------------------------
              int numofpoints = 0;
              int numvertex = DRT::UTILS::getNumberOfElementCornerNodes(DRT::Element::quad4);
              // get intersection points from refinement cell

              // corner points
              for (int ivertex=0; ivertex<numvertex; ivertex++)
              {
                // remark: all hex elements/cells have 8 corners;
                //         additional nodes (hex20,hex27) are inner nodes and not corners
                //         here, only the corners are written into the point list
                pointlist.push_back(vertexcoord[ivertex]);
                numofpoints++;
              }

              const std::multimap<int,std::vector<double> >& intersectionpoints = surfcell->intersectionpoints_;
              // map corresponds to intersection points, but contains the points' IDs instead of their coordinates
              // links 'pointlist' required by TetGen to 'intersectionpointlist'
              std::map<int,int> intersectionpointsids;

              // intersection points
              // std::map<ID of cut edge in element, coordinates of intersection point>
              for (std::multimap<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
              {
                pointlist.push_back(iter->second);
                intersectionpointsids.insert(std::pair<int, int>(iter->first,numofpoints));
                //intersectionpointsids[iter->first] = numofpoints;
                numofpoints++;
              }

              //TEST
              //std::cout << "number of intersection points " << intersectionpoints.size() << std::endl;
              //for (std::size_t iter=0; iter<pointlist.size(); ++iter)
              //{
              //  std::cout<< iter << std::endl;
              //  std::vector<double> coord = pointlist[iter];
              //  for (std::size_t isd=0; isd<3; isd++)
              //  {
              //    std::cout<< coord[isd] << std::endl;
              //  }
              //}

              //--------------------------------------
              // triangulate flame front within a cell
              //--------------------------------------

              // array containing the point IDs for each line
              int line[4][2] = {
                  {0,  1},
                  {1,  2},
                  {2,  3},
                  {3,  0}};

              //--------------------
              // build first polygon
              //--------------------
              std::map<int,std::vector<int> > polygonpoints;
              polygonpoints.clear();

              for (int ipoly=0; ipoly<2; ipoly++)
              {
                //contains the vertices of the polygon
                std::vector<int> polypoints;

                int iline = -1;
                int ipoint = -1;
                if (ipoly == 0) // first polygon
                {
                  // define first intersection point as starting point
                  iline = intersectionpointsids.begin()->first;
                  ipoint = intersectionpointsids.begin()->second;
                  polypoints.push_back(ipoint);
                }
                if (ipoly == 1) // second polygon
                {
                  // define second (=last) intersection point as starting point
                  iline = intersectionpointsids.rbegin()->first;
                  ipoint = intersectionpointsids.rbegin()->second;
                  polypoints.push_back(ipoint);
                }

                // add second point of intersected line
                polypoints.push_back(line[iline][1]);

                bool finished = false;
                while(!finished)
                {
                  // next line
                  iline = iline +1;
                  if (iline > 3)
                    iline = 0;
                  std::multimap<int,int>::const_iterator iter = intersectionpointsids.find(iline);
                  if (iter != intersectionpointsids.end()) // line is intersected
                  {
                    // add intersection point to list
                    polypoints.push_back(iter->second);
                    finished = true;
                  }
                  else // line is not intersected
                  {
                    // add next node to list
                    polypoints.push_back(line[iline][1]);
                  }
                }
                if (ipoly == 0) // first polygon
                  polypoints.push_back( intersectionpointsids.begin()->second);
                if (ipoly == 1) // second polygon
                  polypoints.push_back( intersectionpointsids.rbegin()->second);

                //TEST Ausgabe
                //std::cout << "polypoints" << std::endl;
                //for (int i=0; i<polypoints.size(); i++)
                //  std::cout << polypoints[i] << std::endl;

                polygonpoints[ipoly] = polypoints;
              }
              //----------------
              // build triangles
              //----------------
              for (std::size_t ipolygons=0; ipolygons<polygonpoints.size(); ipolygons++)
              {
                std::vector<int> polypoints = polygonpoints[ipolygons];
                if (polypoints.size()<4)
                  dserror("TriangulateFlameFront needs at least 3 intersectionpoints");

                //-----------------------------
                // three points form a triangle
                //-----------------------------
                if (polypoints.size()==4) //there are only three different points and three points form a triangle
                {
                  std::vector<int> trianglepoints (3);
                  for (int i=0; i<3; i++)
                  {
                    trianglepoints[i] = polypoints[i];
                  }
                  trianglelist.push_back(trianglepoints);
                }
                //----------------------------------------------------
                // add a midpoint to the pointlist to create triangles
                //----------------------------------------------------
                else
                {
                  //-----------------------------
                  //calculate midpoint of polygon
                  //-----------------------------
                  std::size_t numpoints = polypoints.size() - 1;
                  std::vector<double> midpoint (3);
                  std::vector<double> point1 (3);
                  std::vector<double> point2 (3);
                  if (numpoints%2==0) // even
                  {
                    point1 = pointlist[polypoints[0]];
                    point2 = pointlist[polypoints[numpoints/2]];
                  }
                  else // odd
                  {
                    point1 = pointlist[polypoints[0]];
                    point2 = pointlist[polypoints[(numpoints+1)/2]];
                  }

                  for (int dim=0; dim<3; dim++)
                  {
                    // compute middle of this coordinate
                    midpoint[dim] = (point2[dim] + point1[dim]) * 0.5;
                  }

                  // add midpoint to list of points defining piecewise linear complex (interface surface)
                  std::size_t midpoint_id = pointlist.size();//ids start at 0
                  //std::cout<< "pointlistsize " << pointlist.size() << "midpointid " << midpoint_id<<std::endl;
                  pointlist.push_back(midpoint);

                  //build triangles
                  for (std::size_t j=0; j<polypoints.size()-1; j++)
                  {
                    std::vector<int> trianglepoints (3);
                    trianglepoints[0] = polypoints[j];
                    trianglepoints[1] = polypoints[j+1];
                    trianglepoints[2] = midpoint_id;
                    trianglelist.push_back(trianglepoints);
                  }
                }
              }

              //-------------------------------------------
              // store triangular surface integration cells
              //-------------------------------------------
              for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
              {
                LINALG::SerialDenseMatrix trianglecoord(3,3); //3 directions, 3 nodes
                LINALG::SerialDenseMatrix phystrianglecoord(3,3);

                // average value of G-function (indicates side of interface)
                double averageGvalue = 0.0;
                bool inGplus = false;

                for (int inode=0; inode<3; inode++)
                {
                  static LINALG::Matrix<3,1> tcoord;
                  for (int dim=0; dim<3; dim++)
                  {
                    trianglecoord(dim,inode) = pointlist[trianglelist[itriangle][inode]][dim];
                    tcoord(dim) = trianglecoord(dim,inode);
                  }
                  //-------------------------------------------------------------
                  // determine which side the surface integration cell belongs to
                  //-------------------------------------------------------------
                  {
                    // number of nodes of element
                    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
                    if (gfuncvalues_parent.size() != numnode)
                      dserror("mismatch between number of phi values and number of nodes");
                    // get shape functions at vertex point
                    static LINALG::Matrix<numnode,1> funct;
                    funct.Clear();
                    DRT::UTILS::shape_function_3D(funct,tcoord(0),tcoord(1),tcoord(2),ParentElement()->Shape());

                    double phi = 0;
                    // interpolate phi value
                    for(size_t i = 0; i< numnode; i++)
                      phi += gfuncvalues_parent[i] *funct(i);

                    // calculate average Gfunc value
                    averageGvalue += phi;
                  }

                  GEO::elementToCurrentCoordinatesInPlace(ParentElement()->Shape(), xyze, tcoord);
                  for(int  dim=0; dim<3; dim++)
                    phystrianglecoord(dim,inode) = tcoord(dim);
                }
                averageGvalue /= 3;

                // determine surface cell position
                // ">=" because of the approximation of the interface,
                // also averageGvalue=0 is possible; as always approximation errors are made
                if(averageGvalue>=0.0)
                  inGplus = true;

                // store boundary integration cells in boundaryintcelllist
                // remark: boundary integration cells in bisected cells (these have always triangular boundary
                //         integration cells) are defined to belong to the "plus" domain (G>0)
                surfaceintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::tri3, -1, trianglecoord, Teuchos::null, phystrianglecoord, inGplus));
                //std::cout << "created tri3 boundary cell" << std::endl;
              }
            }
            else
            {
              std::cout << "--------------------------------" << std::endl;
              std::cout << "special case Neumann inflow term" << std::endl;
              std::cout << "--------------------------------" << std::endl;
              // print G-function values on screen
              const std::vector<double>& gfuncvalues = surfcell->GetGfuncValues();
              for (size_t i=0; i<gfuncvalues.size(); i++ )
                std::cout << "G-function value " << i << ": " << gfuncvalues[i] << std::endl;

              // remark: - case: 0 intersection points
              //                 the interface cuts through two nodes diagonally
              //                 we should build two triangular cells
              //         - case: 1 intersection point
              //                 the interface cuts through one nodes and intersects a line
              //                 we should build a triangle and 4 triangular cells
              //         - case: the interface cuts through one node
              //                 should not occur here, since this surface cell is 'touched_'
              //         - case: anything else is very strange
              std::cout << "number of intersection points: " << surfcell->intersectionpoints_.size() << std::endl;
              //dserror("exactly 2 intersection points expected");

              // compute Neumann inflow term as if this surface cell was not bisected
              // remark: this code is duplicated from the non-bisected case below
              {
                // get node coordinates of the surface (2D) element
                LINALG::SerialDenseMatrix xyzesurf(3,4);
                LINALG::SerialDenseMatrix xicoord(3,4);

                size_t numnodesurf = this->NumNode();
                if (numnodesurf != 4)
                  dserror("surfaces of hex8 elements expected (4 nodes)");

                // loop nodes of surface element
                const DRT::Node*const* nodessurf = this->Nodes();
                for (size_t inode=0; inode<numnodesurf; inode++)
                {
                  // compute local element coordinates
                  // remark: since we deal with nodes, we should get clean numbers (-1 or 1) for 'xicoord'
                  const double* x = nodessurf[inode]->X();
                  {
                    LINALG::Matrix<3,1> xyznode;
                    xyznode(0) = x[0];
                    xyznode(1) = x[1];
                    xyznode(2) = x[2];

                    LINALG::Matrix<3,1> xsi;
                    GEO::CUT::Position<DRT::Element::hex8> pos(xyze,xyznode);
                    pos.Compute();
                    xsi = pos.LocalCoordinates();
                    // fill array
                    xicoord(0,inode) = xsi(0);
                    xicoord(1,inode) = xsi(1);
                    xicoord(2,inode) = xsi(2);
                  }
                  // fill array
                  xyzesurf(0,inode) = x[0];
                  xyzesurf(1,inode) = x[1];
                  xyzesurf(2,inode) = x[2];
                }
                bool inGplus = false;
                {
                  /*------------------------------------------------------------------------------
                   * - compute phi at each node of domain integration cell
                   * - compute average G-value for each surface integration cell
                   * -> >=0 in plus == true
                   * ->  <0 in minus == false
                   * ----------------------------------------------------------------------------*/

                  int numcellnodes = 0;
                  if (surfcell->Shape() == DRT::Element::quad4)
                    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
                  else
                    dserror("quad4 surface cell expected");

                  //calculate average Gfunc value
                  double averageGvalue = 0.0;

                  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
                    averageGvalue += gfuncvalues_surfcell[icellnode];

                  if (numcellnodes == 0.0)
                    dserror("division by zero: number of vertices per cell is 0!");
                  averageGvalue /= numcellnodes;

                  // determine DomainIntCell position
                  // ">=" because of the approximation of the interface,
                  // also averageGvalue=0 is possible; as always approximation errors are made
                  if(averageGvalue>=0.0)
                    inGplus = true;
                }
                //std::cout << "created quad4 boundary cell" << std::endl;
                surfaceintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, xicoord, Teuchos::null, xyzesurf, inGplus));
              }
            }
          }
          else // non-bisected cell (touched or uncut)
          {
            // get node coordinates of the surface (2D) element
            LINALG::SerialDenseMatrix xyzesurf(3,4);
            LINALG::SerialDenseMatrix xicoord(3,4);

            size_t numnodesurf = this->NumNode();
            if (numnodesurf != 4)
              dserror("surfaces of hex8 elements expected (4 nodes)");

            // loop nodes of surface element
            const DRT::Node*const* nodessurf = this->Nodes();
            for (size_t inode=0; inode<numnodesurf; inode++)
            {
              // compute local element coordinates
              // remark: since we deal with nodes, we should get clean numbers (-1 or 1) for 'xicoord'
              const double* x = nodessurf[inode]->X();
              {
                LINALG::Matrix<3,1> xyznode;
                xyznode(0) = x[0];
                xyznode(1) = x[1];
                xyznode(2) = x[2];

                LINALG::Matrix<3,1> xsi;
                GEO::CUT::Position<DRT::Element::hex8> pos(xyze,xyznode);
                pos.Compute();
                xsi = pos.LocalCoordinates();
                // fill array
                xicoord(0,inode) = xsi(0);
                xicoord(1,inode) = xsi(1);
                xicoord(2,inode) = xsi(2);
              }
              // fill array
              xyzesurf(0,inode) = x[0];
              xyzesurf(1,inode) = x[1];
              xyzesurf(2,inode) = x[2];
            }
            bool inGplus = false;
            {
              /*------------------------------------------------------------------------------
               * - compute phi at each node of domain integration cell
               * - compute average G-value for each surface integration cell
               * -> >=0 in plus == true
               * ->  <0 in minus == false
               * ----------------------------------------------------------------------------*/

              int numcellnodes = 0;
              if (surfcell->Shape() == DRT::Element::quad4)
                numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
              else
                dserror("quad4 surface cell expected");

              //calculate average Gfunc value
              double averageGvalue = 0.0;

              for (int icellnode=0; icellnode<numcellnodes; icellnode++)
                averageGvalue += gfuncvalues_surfcell[icellnode];

              if (numcellnodes == 0.0)
                dserror("division by zero: number of vertices per cell is 0!");
              averageGvalue /= numcellnodes;

              // determine DomainIntCell position
              // ">=" because of the approximation of the interface,
              // also averageGvalue=0 is possible; as always approximation errors are made
              if(averageGvalue>=0.0)
                inGplus = true;
            }
            //std::cout << "created quad4 boundary cell" << std::endl;
            surfaceintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, xicoord, Teuchos::null, xyzesurf, inGplus));
          }
          {
            //const bool screen_out = false;
            //
            //const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("NeumannSurface", ParentElement()->Id(), 200, screen_out, 0);
            //std::ofstream gmshfilecontent(filename.c_str());
            //{
            //  gmshfilecontent << "View \" " << "SurfaceCell \" {\n";
            //  const int Id = this->Id();
            //  const int numnode = this->NumNode();
            //  const DRT::Element::DiscretizationType distype = this->Shape();
            //  const std::vector<std::vector<double> > vertexcoord = surfcell->GetVertexCoord();
            //  LINALG::SerialDenseMatrix Pos(3,numnode);
            //  for (int i=0; i<numnode; i++)
            //  {
            //    for (size_t k=0; k<3; k++)
            //    {
            //      Pos(k,i) = vertexcoord[i][k];
            //    }
            //  }
            //  IO::GMSH::cellWithScalarToStream(distype, Id, Pos, gmshfilecontent);
            //  gmshfilecontent << "};\n";
            //}
            //{
            //  gmshfilecontent << "View \" " << "SurfaceIntCells \" {\n";
            //  for(std::size_t icell = 0; icell < surfaceintcelllist.size(); icell++)
            //  {
            //    GEO::BoundaryIntCell cell = surfaceintcelllist[icell];
            //    const LINALG::SerialDenseMatrix& cellpos = cell.CellNodalPosXiDomain();//cell.CellNodalPosXYZ();
            //    const double color = 5;
            //    gmshfilecontent << IO::GMSH::cellWithScalarToString(cell.Shape(), color, cellpos) << std::endl;
            //  }
            //  gmshfilecontent << "};\n";
            //}
            //gmshfilecontent.close();
          }
        }

        //--------------------------------------------------
        // build map of surface nodes for element dofmanager
        //--------------------------------------------------
        Teuchos::RCP<XFEM::ElementDofManager> eleDofManager = ParentElement()->GetEleDofManager();

        //----------------------
        // unpack all parameters
        //----------------------
        // time integration parameters
        const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");
        const double time   = params.get<double>("time");
        const double dt     = params.get<double>("dt");
        const double theta  = params.get<double>("theta");
        const double ga_gamma  = params.get<double>("gamma");
        const double ga_alphaF = params.get<double>("alphaF");
        const double ga_alphaM = params.get<double>("alphaM");
        // generalized alpha time integration scheme
        bool genalpha = false;
        // parameter for type of linearization
        const bool newton = params.get<bool>("include reactive terms for linearisation",false);
        // instationary formulation
        bool instationary = true;
        if (timealgo == INPAR::FLUID::timeint_stationary) instationary = false;

        // general parameters for two-phase flow and premixed combustion problems
        const INPAR::COMBUST::CombustionType combusttype   = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");

        // reshape element matrices and vectors and
        const int eledim = (int)(lm).size();
        elemat1.Shape(eledim,eledim);
        elevec1.Size (eledim);
        // initialize to zero
        elemat1.Scale(0.0);
        elevec1.Scale(0.0);

        //-----------------------------
        // assemble Neumann inflow term
        //-----------------------------
        // get the list of materials
        const Teuchos::RCP<MAT::Material> material = ParentElement()->Material();

        // extract local (element level) vectors from global state vectors
        DRT::ELEMENTS::Combust3::MyStateSurface mystate(
            discretization, (lm), true, false, false, ParentElement(), phinp);

        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
            *eleDofManager, ParentElement()->NumNode(), ParentElement()->NodeIds());

        COMBUST::callSysmatNeumannInflow(
            assembly_type,
            this->ParentElement(),
            this,
            *eleDofManager,
            mystate,
            elemat1,
            elevec1,
            surfaceintcelllist,
            material,
            timealgo,
            time,
            dt,
            theta,
            ga_gamma,
            ga_alphaF,
            ga_alphaM,
            newton,
            instationary,
            genalpha,
            combusttype);

        break;
      }
      default:
        dserror("Unknown type of action for Combust3_Surface");
    } // end of switch(act)

    return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3Surface::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  //std::cout << "/!\\ warning === intersected Neumann boundary conditions in XFEM problems are not treated correctly" << std::endl;
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  const double thsl = params.get("thsl",0.0);

  const DiscretizationType distype = this->Shape();

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  // set number of nodes
  const int iel   = this->NumNode();

  DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
  case tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct       (iel);
  Epetra_SerialDenseMatrix  deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct, e0, e1, distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)
    const double fac = intpoints.qwgt[gpid] * drs * curvefac * thsl;

    // factor given by spatial function
    double functfac = 1.0;
    // determine coordinates of current Gauss point
    double coordgp[3];
    coordgp[0]=0.0;
    coordgp[1]=0.0;
    coordgp[2]=0.0;
    for (int i = 0; i< iel; i++)
    {
      coordgp[0] += xyze(0,i) * funct[i];
      coordgp[1] += xyze(1,i) * funct[i];
      coordgp[2] += xyze(2,i) * funct[i];
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        if (func) functnum = (*func)[dim];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;
        }
        elevec1[node*numdf+dim]+= funct[node]*(*onoff)[dim]*(*val)[dim]*fac*functfac;
      }
    }

  } /* end of loop over integration points gpid */

  return 0;
}


/* compute kovariant metric tensor G for fluid element        gammi 04/07

                        +-       -+
                        | g11 g12 |
                    G = |         |
                        | g12 g22 |
                        +-       -+

 where (o denotes the inner product, xyz a vector)


                            dxyz   dxyz
                    g11 =   ---- o ----
                             dr     dr

                            dxyz   dxyz
                    g12 =   ---- o ----
                             dr     ds

                            dxyz   dxyz
                    g22 =   ---- o ----
                             ds     ds


 and the square root of the first fundamental form


                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/

 they are needed for the integration over the surface element

*/
void  DRT::ELEMENTS::Combust3Surface::ComputeMetricTensorForSurface(
    const int                       numnode,
    const Epetra_SerialDenseMatrix& xyze,
    const Epetra_SerialDenseMatrix& deriv,
    LINALG::Matrix<2,2>&            metrictensor,
    double&                         detmetric
    ) const
{
  // get jacobian matrix d x / d \xi  (3x2)
  LINALG::Matrix<3,2> dxyzdrs;
  // dxyzdrs(i,j) = xyze_boundary(i,k)*deriv_boundary(j,k);
  xyze.GEMM('N','T',3,2,numnode,1.0,xyze.A(),xyze.LDA(),deriv.A(),deriv.LDA(),0.0,dxyzdrs.A(),dxyzdrs.M());

  // compute covariant metric tensor G for surface element (2x2)
  // metric = dxyzdrs(k,i)*dxyzdrs(k,j);
  metrictensor.MultiplyTN(dxyzdrs,dxyzdrs);

  detmetric = sqrt(metrictensor.Determinant());

  return;

  // this is old documentation
  // nomenclature has changed, but this might still ne helpfull to understand   henke

  /*
  |                                              0 1 2
  |                                             +-+-+-+
  |       0 1 2              0...iel-1          | | | | 0
  |      +-+-+-+             +-+-+-+-+          +-+-+-+
  |      | | | | 1           | | | | | 0        | | | | .
  |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
  |      | | | | 2           | | | | | 1        | | | | .
  |      +-+-+-+             +-+-+-+-+          +-+-+-+
  |                                             | | | | iel-1
  |                                             +-+-+-+
  |
  |       dxyzdrs             deriv              xyze^T
  |
  |
  |                                     +-            -+
  |                         | dx   dy   dz |
  |                         | --   --   -- |
  |                      | dr   dr   dr |
  |  yields               dxyzdrs =  |              |
  |                         | dx   dy   dz |
  |                         | --   --   -- |
  |                      | ds   ds   ds |
  |                                     +-            -+
  |
  */

  /*
  |
  |      +-           -+    +-            -+   +-            -+ T
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
  |      |             |    | dr   dr   dr |   | dr   dr   dr |
  |      |             |  = |              | * |              |
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
  |      |             |    | ds   ds   ds |   | ds   ds   ds |
  |      +-           -+    +-            -+   +-            -+
  |
  | the calculation of g21 is redundant since g21=g12
  */

/*
                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/
*/
}

/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)          g.bau 07/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Surface::IntegrateShapeFunction(
    Teuchos::ParameterList&          params,
    const DRT::Discretization&       discretization,
    const std::vector<int>&          lm,
    Epetra_SerialDenseVector&        elevec1,
    const std::vector<double>&       edispnp)
{
  // there are 3 velocities and 1 pressure
  return;
  const int numdf = 4;

//  const double thsl = params.get("thsl",1.0);

  const DiscretizationType distype = this->Shape();

/*  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;
*/

  // set number of nodes
  const int iel   = this->NumNode();

  DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
  case tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }

    // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector      funct       (iel);
  Epetra_SerialDenseMatrix      deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  LINALG::Matrix<2,2>           metrictensor;
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  if (edispnp.size()!=0)
  {
    for (int i=0;i<iel;i++)
    {
      xyze(0,i) += edispnp[4*i];
      xyze(1,i) += edispnp[4*i+1];
      xyze(2,i) += edispnp[4*i+2];
    }
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);

  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct, e0, e1, distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration

    ComputeMetricTensorForSurface(iel,xyze,deriv,metrictensor,drs);

    // values are multiplied by the product from inf. area element,
    // the gauss weight and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)

    //const double fac = intpoints.qwgt[gpid] * drs * thsl;
    const double fac = intpoints.qwgt[gpid] * drs;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+=
          funct[node] * fac;
      }
    }

  } /* end of loop over integration points gpid */


return;
} // DRT::ELEMENTS::Combust3Surface::IntegrateShapeFunction

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Surface::IntegrateSurfaceFlow(
    Teuchos::ParameterList&          params,
    const DRT::Discretization&       discretization,
    const std::vector<int>&          lm,
    Epetra_SerialDenseVector&        elevec1,
    const std::vector<double>&       myvelnp)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

  DRT::UTILS::GaussRule2D  gaussrule = DRT::UTILS::intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = DRT::UTILS::intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = DRT::UTILS::intrule_tri_3point;
      break;
  case tri6:
      gaussrule = DRT::UTILS::intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector      funct       (iel);
  Epetra_SerialDenseMatrix      deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);
  // node velocities
  Epetra_SerialDenseMatrix      evelnp      (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  LINALG::Matrix<2,2>           metrictensor;
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // get element velocities
  for(int i=0;i<iel;i++)
  {
    evelnp(0,i)=myvelnp[i*iel+0];
    evelnp(1,i)=myvelnp[i*iel+1];
    evelnp(2,i)=myvelnp[i*iel+2];
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);

  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    LINALG::Matrix<2,1> xi_gp;
    xi_gp(0) = intpoints.qxg[gpid][0];
    xi_gp(1) = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct, xi_gp(0), xi_gp(1), distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv, xi_gp(0), xi_gp(1), distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    ComputeMetricTensorForSurface(iel,xyze,deriv,metrictensor,drs);

    // values are multiplied by the product from inf. area element and gauss weight
    const double fac = drs * intpoints.qwgt[gpid];

    // velocity at gausspoint
    const LINALG::Matrix<3,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, funct, iel);

    // get normal vector (in x coordinates) to surface element at integration point
    LINALG::Matrix<3,1> n(true);
    GEO::computeNormalToSurfaceElement(this->Shape(), xyze, xi_gp, n);

    // flowrate = u_i * n_i
    const double flowrate = gpvelnp(0)*n(0) + gpvelnp(1)*n(1) + gpvelnp(2)*n(2);

    // store flowrate at first dof of each node
    // use negatve value so that inflow is positiv
    for (int node=0;node<iel;++node)
    {
      elevec1[node*numdf] -= funct[node] * fac * flowrate;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
  |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Surface::LocationVector(
    const Discretization&   dis,
    LocationArray&          la,
    bool                    doDirichlet,
    const std::string&      condstring,
    Teuchos::ParameterList& params
    ) const
{
  DRT::ELEMENTS::Combust3Surface::ActionType act = Combust3Surface::none;
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_Neumann_inflow")
      act = Combust3Surface::calc_Neumann_inflow;

  switch(act)
  {
  case calc_Neumann_inflow:
    // special cases: the boundary element assembles also into
    // the inner dofs of its parent element
    // note: using these actions, the element will get the parent location vector
    //       as input in the respective evaluate routines
    ParentElement()->LocationVector(dis,la,doDirichlet);
    break;
  default:
    DRT::Element::LocationVector(dis,la,doDirichlet);
    break;
  }
  return;
}
