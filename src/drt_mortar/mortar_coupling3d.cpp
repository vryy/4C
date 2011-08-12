/*!----------------------------------------------------------------------
\file mortar_coupling3d.cpp
\brief Classes for mortar coupling in 3D.

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "mortar_coupling3d.H"
#include "mortar_node.H"
#include "mortar_projector.H"
#include "mortar_integrator.H"
#include "mortar_defines.H"
#include "mortar_utils.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3d::Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                               bool auxplane, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele) :
shapefcn_(INPAR::MORTAR::shape_undefined),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3d::Coupling3d(const INPAR::MORTAR::ShapeFcn shapefcn,
                               DRT::Discretization& idiscret, int dim, bool quad,
                               bool auxplane, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele) :
shapefcn_(shapefcn),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  get communicator  (public)                                popp 06/09|
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::Coupling3d::Comm() const
{
  return idiscret_.Comm();
}

/*----------------------------------------------------------------------*
 |  Evaluate coupling (3D)                                    popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::EvaluateCoupling()
{
  // rough check whether element centers are "near"
  // whether or not quadratic 3d coupling is performed, we only
  // check the distance of the parent slave and master elements
  bool near = RoughCheckCenters();
  if (!near) return false;

  // tolerance for polygon clipping
  double tol= 0.0;

  // map to store projection parameter alpha for each master node
  // (currently only necessary for non-AUXPLANE case)
  map<int,double> projpar;

  // *******************************************************************
  // ************ Coupling with or without auxiliary plane *************
  // *******************************************************************
  if (CouplingInAuxPlane())
  {
    // compute auxiliary plane for 3D coupling
    AuxiliaryPlane();

    // rough check of orientation of element centers
    // if slave and master element center normals form an
    // angle > 90° the pair will not be considered further
    bool orient = RoughCheckOrient();
    if (!orient) return false;
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave();

    // project master element nodes onto auxiliary plane
    ProjectMaster();

    // tolerance for polygon clipping
    double sminedge = SlaveIntElement().MinEdgeSize();
    double mminedge = MasterIntElement().MinEdgeSize();
    tol = MORTARCLIPTOL * min(sminedge,mminedge);
  }

  // *******************************************************************
  else //(!CouplingInAuxPlane())
  {
    // get some data
    int nsnodes = SlaveIntElement().NumNode();
    int nmnodes = MasterIntElement().NumNode();

    // get slave vertices in slave element parameter space (direct)
    // additionally get slave vertex Ids for later linearization
    vector<vector<double> > svertices(nsnodes,vector<double>(3));
    vector<int> snodeids(1);

    for (int i=0;i<nsnodes;++i)
    {
      double xi[2] = {0.0, 0.0};
      SlaveIntElement().LocalCoordinatesOfNode(i,xi);
      svertices[i][0] = xi[0];
      svertices[i][1] = xi[1];
      svertices[i][2] = 0.0;

      // relevant ids (here only slave node id itself)
      snodeids[0] = SlaveIntElement().NodeIds()[i];

      // store into vertex data structure
      SlaveVertices().push_back(Vertex(svertices[i],Vertex::slave,snodeids,NULL,NULL,false,false,NULL,-1.0));
    }

    // get master vertices in slave element parameter space (project)
    // additionally get master vertex Ids for later linearization
    vector<vector<double> > mvertices(nmnodes,vector<double>(3));
    vector<int> mnodeids(1);
    for (int i=0;i<nmnodes;++i)
    {
      int gid = MasterIntElement().NodeIds()[i];
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MortarNode* mnode = static_cast<MortarNode*>(node);

      // do the projection
      // the third component of sxi will be the proj. parameter alpha!
      double sxi[2] = {0.0, 0.0};
      double alpha = 0.0;
      MORTAR::MortarProjector projector(3);
      //cout << "Projecting master node ID: " << mnode->Id() << endl;
      projector.ProjectElementNormal3D(*mnode,SlaveIntElement(),sxi,alpha);

      mvertices[i][0] = sxi[0];
      mvertices[i][1] = sxi[1];
      mvertices[i][2] = 0.0;

      // relevant ids (here only master node id itself)
      mnodeids[0] = gid;

      // store proj. parameter for later linearization
      projpar[gid] = alpha;

      // store into vertex data structure
      MasterVertices().push_back(Vertex(mvertices[i],Vertex::projmaster,mnodeids,NULL,NULL,false,false,NULL,-1.0));
    }

    // normal is (0,0,1) in slave element parameter space
    Auxn()[0] = 0.0; Auxn()[1] = 0.0; Auxn()[2] = 1.0;
    Lauxn() = 1.0;

    // tolerance for polygon clipping
    // minimum edge size in parameter space is 1
    tol = MORTARCLIPTOL;
  }
  // *******************************************************************

  // do polygon clipping (NEW version)
  bool clip = PolygonClippingConvexHull(SlaveVertices(),MasterVertices(),Clip(),tol);
  int clipsize = (int)(Clip().size());
  
  // within polygon clipping we may have performed a second rough check
  // if the two elements are really "near" (NOTE: this has only been done
  // if problems occured within polygon clipping, i.e. the projected master 
  // element being non-convex. In the standard case, bool clip is simply true)
  // --> this way, robustness of coupling is further increased!
  if (!clip) return false;

  // proceed only if clipping polygon is at least a triangle
  bool overlap = false;
  if (clipsize>=3) overlap = true;
  if (overlap)
  {
    // check / set  projection status of slave nodes
    HasProjStatus();

    // do linearization + triangulation of clip polygon
    Triangulation(projpar);

    // do integration of integration cells
    IntegrateCells();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with centers)           popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::RoughCheckCenters()
{
  double sme = SlaveElement().MaxEdgeSize();
  double mme = MasterElement().MaxEdgeSize();
  double near = 2.0 * max(sme,mme);

  double loccs[2] = {0.0, 0.0};
  DRT::Element::DiscretizationType dts = SlaveElement().Shape();
  if (dts==MortarElement::tri3 || dts==MortarElement::tri6)
  {
    loccs[0] = 1.0/3;
    loccs[1] = 1.0/3;
  }
  double loccm[2] = {0.0, 0.0};
  DRT::Element::DiscretizationType dtm = MasterElement().Shape();
  if (dtm==MortarElement::tri3 || dtm==MortarElement::tri6)
  {
    loccm[0] = 1.0/3;
    loccm[1] = 1.0/3;
  }

  double sc[3] = {0.0, 0.0, 0.0};
  double mc[3] = {0.0, 0.0, 0.0};
  SlaveElement().LocalToGlobal(loccs,sc,0);
  MasterElement().LocalToGlobal(loccm,mc,0);

  double cdist = sqrt((mc[0]-sc[0])*(mc[0]-sc[0])+(mc[1]-sc[1])*(mc[1]-sc[1])+(mc[2]-sc[2])*(mc[2]-sc[2]));
  if (cdist>=near) return false;
  else return true;
}

/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with master nodes)      popp 11/09|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::RoughCheckNodes()
{
  // project master nodes onto auxiliary plane
  int nnodes = MasterIntElement().NumNode();
  DRT::Node** mynodes = MasterIntElement().Nodes();
  if (!mynodes) dserror("ERROR: RoughCheckNodes: Null pointer!");
  
  // prepare check
  bool near = false;
  double sme = SlaveElement().MaxEdgeSize();
  double mme = MasterElement().MaxEdgeSize();
  double limit = 0.3 * max(sme,mme);
  
  for (int i=0;i<nnodes;++i)
  {
    MortarNode* mycnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: RoughCheckNodes: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0]-Auxc()[0])*Auxn()[0]
                + (mycnode->xspatial()[1]-Auxc()[1])*Auxn()[1]
                + (mycnode->xspatial()[2]-Auxc()[2])*Auxn()[2];

    // only near if at least one node is "close" to auxiliary plane
    // (signed distance in order to capture intrusion correctly)
    if (dist < limit) near = true;
  }

  return near;
}

/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with normals)           popp 11/10|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::RoughCheckOrient()
{
  // we first need the master element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = {0.0, 0.0};

  DRT::Element::DiscretizationType dt = MasterIntElement().Shape();
  if (dt==MortarElement::tri3 || dt==MortarElement::tri6)
  {
    loccenter[0] = 1.0/3;
    loccenter[1] = 1.0/3;
  }
  else if (dt==MortarElement::quad4 || dt==MortarElement::quad8 || dt==MortarElement::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else dserror("ERROR: RoughCheckOrient called for unknown element type");

  // compute the unit normal vector at the master element center
  double nmc[3] = {0.0, 0.0, 0.0};
  MasterIntElement().ComputeUnitNormalAtXi(loccenter,nmc);

  // check orientation with respect to slave element
  double dot = nmc[0]*Auxn()[0]+nmc[1]*Auxn()[1]+nmc[2]*Auxn()[2];
  if (dot < -1.0e-12) return true;
  else                return false;
}

/*----------------------------------------------------------------------*
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::AuxiliaryPlane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = {0.0, 0.0};

  DRT::Element::DiscretizationType dt = SlaveIntElement().Shape();
  if (dt==MortarElement::tri3 || dt==MortarElement::tri6)
  {
    loccenter[0] = 1.0/3;
    loccenter[1] = 1.0/3;
  }
  else if (dt==MortarElement::quad4 || dt==MortarElement::quad8 || dt==MortarElement::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else dserror("ERROR: AuxiliaryPlane called for unknown element type");

  // compute element center via shape fct. interpolation
  SlaveIntElement().LocalToGlobal(loccenter,Auxc(),0);

  // we then compute the unit normal vector at the element center
  Lauxn() = SlaveIntElement().ComputeUnitNormalAtXi(loccenter,Auxn());

  //cout << "Slave Element: " << SlaveIntElement().Id() << endl;
  //cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << endl;
  //cout << "->Normal: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] << endl;

  return true;
}

/*----------------------------------------------------------------------*
 |  Project slave element onto auxiliary plane (public)       popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::ProjectSlave()
{
  // project slave nodes onto auxiliary plane
  int nnodes = SlaveIntElement().NumNode();
  DRT::Node** mynodes = SlaveIntElement().Nodes();
  if (!mynodes) dserror("ERROR: ProjectSlave: Null pointer!");

  // initialize storage for slave coords + their ids
  vector<double> vertices(3);
  vector<int> snodeids(1);

  for (int i=0;i<nnodes;++i)
  {
    MortarNode* mycnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: ProjectSlave: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0]-Auxc()[0])*Auxn()[0]
                + (mycnode->xspatial()[1]-Auxc()[1])*Auxn()[1]
                + (mycnode->xspatial()[2]-Auxc()[2])*Auxn()[2];

    // compute projection
    for (int k=0;k<3;++k) vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    snodeids[0] = mycnode->Id();

    // store into vertex data structure
    SlaveVertices().push_back(Vertex(vertices,Vertex::slave,snodeids,NULL,NULL,false,false,NULL,-1.0));

    //cout << "->RealNode(S) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << endl;
    //cout << "->ProjNode(S) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " << vertices[2] << endl;
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Project master element onto auxiliary plane (public)      popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::ProjectMaster()
{
  // project master nodes onto auxiliary plane
  int nnodes = MasterIntElement().NumNode();
  DRT::Node** mynodes = MasterIntElement().Nodes();
  if (!mynodes) dserror("ERROR: ProjectMaster: Null pointer!");

  // initialize storage for master coords + their ids
  vector<double> vertices(3);
  vector<int> mnodeids(1);

  for (int i=0;i<nnodes;++i)
  {
    MortarNode* mycnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: ProjectMaster: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0]-Auxc()[0])*Auxn()[0]
                + (mycnode->xspatial()[1]-Auxc()[1])*Auxn()[1]
                + (mycnode->xspatial()[2]-Auxc()[2])*Auxn()[2];

    // compute projection
    for (int k=0;k<3;++k) vertices[k] = mycnode->xspatial()[k] - dist * Auxn()[k];

    // get node id, too
    mnodeids[0] = mycnode->Id();

    // store into vertex data structure
    MasterVertices().push_back(Vertex(vertices,Vertex::projmaster,mnodeids,NULL,NULL,false,false,NULL,-1.0));

    //cout << "->RealNode(M) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << endl;
    //cout << "->ProjNode(M) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " << vertices[2] << endl;
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Clipping of two polygons                                  popp 11/08|
 *----------------------------------------------------------------------*/
void MORTAR::Coupling3d::PolygonClipping(vector<Vertex>& poly1,
                                          vector<Vertex>& poly2,
                                          vector<Vertex>& respoly,
                                          double& tol)
{
  //**********************************************************************
  dserror("ERROR: PolygonClipping outdated, use PolygonClippingConvexHull instead!");
  //**********************************************************************
  
  // choose output
  bool out = false;
  
  //**********************************************************************
  // STEP1: Input check
  // - input polygons must consist of min. 3 vertices each
  // - rotation of poly1 must be c-clockwise w.r.t. (0,0,1) or Auxn()
  // - rotation of poly 2 changed to c-clockwise w.r.t. (0,0,1) or Auxn()
  // - both input polygons must be convex
  //**********************************************************************

  // check input variables
  if ((int)poly1.size()<3 || (int)poly2.size()<3)
    dserror("ERROR: Input Polygons must consist of min. 3 vertices each");

  // check for rotation of polygon1 (slave) and polgon 2 (master)
  // note that we implicitly already rely on convexity here!
  // first get geometric centers of polygon1 and polygon2
  double center1[3] = {0.0, 0.0, 0.0};
  double center2[3] = {0.0, 0.0, 0.0};

  for (int i=0;i<(int)poly1.size();++i)
    for (int k=0;k<3;++k)
      center1[k] += poly1[i].Coord()[k]/((int)poly1.size());

  for (int i=0;i<(int)poly2.size();++i)
    for (int k=0;k<3;++k)
      center2[k] += poly2[i].Coord()[k]/((int)poly2.size());

  if (out)
  {
    cout << "Center 1: " << center1[0] << " " << center1[1] << " " << center1[2] << endl;
    cout << "Center 2: " << center2[0] << " " << center2[1] << " " << center2[2] << endl;
  }
  
  // then we compute the counter-clockwise plane normal
  double diff1[3] = {0.0, 0.0, 0.0};
  double edge1[3] = {0.0, 0.0, 0.0};
  double diff2[3] = {0.0, 0.0, 0.0};
  double edge2[3] = {0.0, 0.0, 0.0};

  for (int k=0;k<3;++k)
  {
    diff1[k] = poly1[0].Coord()[k]-center1[k];
    edge1[k] = poly1[1].Coord()[k]-poly1[0].Coord()[k];
    diff2[k] = poly2[0].Coord()[k]-center2[k];
    edge2[k] = poly2[1].Coord()[k]-poly2[0].Coord()[k];
  }

  double cross1[3] = {0.0, 0.0, 0.0};
  double cross2[3] = {0.0, 0.0, 0.0};

  cross1[0] = diff1[1]*edge1[2]-diff1[2]*edge1[1];
  cross1[1] = diff1[2]*edge1[0]-diff1[0]*edge1[2];
  cross1[2] = diff1[0]*edge1[1]-diff1[1]*edge1[0];

  cross2[0] = diff2[1]*edge2[2]-diff2[2]*edge2[1];
  cross2[1] = diff2[2]*edge2[0]-diff2[0]*edge2[2];
  cross2[2] = diff2[0]*edge2[1]-diff2[1]*edge2[0];

  // check against auxiliary plane normal
  double check1 = cross1[0]*Auxn()[0]+cross1[1]*Auxn()[1]+cross1[2]*Auxn()[2];
  double check2 = cross2[0]*Auxn()[0]+cross2[1]*Auxn()[1]+cross2[2]*Auxn()[2];

  // check polygon 1 and throw dserror if not c-clockwise
  if (check1<=0) dserror("ERROR: Polygon 1 (slave) not ordered counter-clockwise!");

  // check polygon 2 and reorder in c-clockwise direction
  if (check2<0)
  {
    if (out) cout << "Polygon 2 (master) not ordered counter-clockwise -> reordered!" << endl;
    std::reverse(poly2.begin(), poly2.end());
  }

  // check if the two input polygons are convex
  // a polygon is convex if the scalar product of an edge normal and the
  // next edge direction is negative for all edges
  for (int i=0;i<(int)poly1.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly1.size()-1) edge[k] = poly1[i+1].Coord()[k] - poly1[i].Coord()[k];
      else edge[k] = poly1[0].Coord()[k] - poly1[i].Coord()[k];
    }

    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly1.size()-2) nextedge[k] = poly1[i+2].Coord()[k] - poly1[i+1].Coord()[k];
      else if (i==(int)poly1.size()-2) nextedge[k] = poly1[0].Coord()[k] - poly1[i+1].Coord()[k];
      else nextedge[k] = poly1[1].Coord()[k] - poly1[0].Coord()[k];
    }

    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) dserror("ERROR: Input polygon 1 not convex");
  }

  for (int i=0;i<(int)poly2.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly2.size()-1) edge[k] = poly2[i+1].Coord()[k] - poly2[i].Coord()[k];
      else edge[k] = poly2[0].Coord()[k] - poly2[i].Coord()[k];
    }

    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly2.size()-2) nextedge[k] = poly2[i+2].Coord()[k] - poly2[i+1].Coord()[k];
      else if (i==(int)poly2.size()-2) nextedge[k] = poly2[0].Coord()[k] - poly2[i+1].Coord()[k];
      else nextedge[k] = poly2[1].Coord()[k] - poly2[0].Coord()[k];
    }

    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) dserror("ERROR: Input polygon 2 not convex");
  }

  // print final input polygons to screen
  if (out)
  {
    cout << "\nInput Poylgon 1:";
    for (int i=0;i<(int)poly1.size();++i)
      cout << "\nVertex " << i << ":\t" << scientific << poly1[i].Coord()[0]
           << "\t" << poly1[i].Coord()[1] << "\t" << poly1[i].Coord()[2];
  
    cout << "\nInput Poylgon 2:";
    for (int i=0;i<(int)poly2.size();++i)
      cout << "\nVertex " << i << ":\t" << scientific << poly2[i].Coord()[0]
           << "\t" << poly2[i].Coord()[1] << "\t" << poly2[i].Coord()[2];
  
    cout << endl << endl;
  }

  //**********************************************************************
  // STEP2: Extend Vertex data structures
  // - note that poly1 is the slave element and poly2 the master element
  // - assign Next() and Prev() pointers to initialize linked structure
  //**********************************************************************

  // set previous and next Vertex pointer for all elements in lists
  for (int i=0;i<(int)poly1.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly1.size()-1)
    {
      poly1[i].AssignNext(&poly1[i+1]);
      poly1[i].AssignPrev(&poly1[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly1[i].AssignNext(&poly1[i+1]);
      poly1[i].AssignPrev(&poly1[(int)poly1.size()-1]);
    }
    // last element in list
    else
    {
      poly1[i].AssignNext(&poly1[0]);
      poly1[i].AssignPrev(&poly1[i-1]);
    }
  }
  for (int i=0;i<(int)poly2.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly2.size()-1)
    {
      poly2[i].AssignNext(&poly2[i+1]);
      poly2[i].AssignPrev(&poly2[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly2[i].AssignNext(&poly2[i+1]);
      poly2[i].AssignPrev(&poly2[(int)poly2.size()-1]);
    }
    // last element in list
    else
    {
      poly2[i].AssignNext(&poly2[0]);
      poly2[i].AssignPrev(&poly2[i-1]);
    }
  }

  //**********************************************************************
  // STEP3: Avoid degenerate cases
  // - if a point of poly1 is close (<tol) to a edge of poly2 or vice
  //   versa we move this point away from the edge by tol
  //**********************************************************************
  for (int i=0;i<(int)poly1.size();++i)
  {
    for (int j=0;j<(int)poly2.size();++j)
    {
      // we need diff vector and edge2 first
      double diff1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff1[k] = poly1[i].Coord()[k] - poly2[j].Coord()[k];
        edge2[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
      }

      // check if point of poly1 lies within [0,1] for edge2
      double checkalpha = diff1[0]*edge2[0]+diff1[1]*edge2[1]+diff1[2]*edge2[2];
      checkalpha /= (edge2[0]*edge2[0]+edge2[1]*edge2[1]+edge2[2]*edge2[2]);

      // proceed only if inside [0,1] with tolerance tol
      if (checkalpha<-tol || checkalpha>1+tol) continue;

      // compute distance from point on poly1 to edge2
      double n2[3] = {0.0, 0.0, 0.0};
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];
      double ln = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
      for (int k=0;k<3;++k) n2[k] /= ln;

      double dist = diff1[0]*n2[0]+diff1[1]*n2[1]+diff1[2]*n2[2];

      // move point away if very close to edge 2
      if (dist > -tol && dist < 0)
      {
        //cout << "Vertex " << i << " on poly1 is very close to edge " << j << " of poly2 -> moved inside!" << endl;
        poly1[i].Coord()[0] -= tol*n2[0];
        poly1[i].Coord()[1] -= tol*n2[1];
        poly1[i].Coord()[2] -= tol*n2[2];
      }
      else if (dist < tol && dist >= 0)
      {
        //cout << "Vertex " << i << " on poly1 is very close to edge " << j << " of poly2 -> moved outside!" << endl;
        poly1[i].Coord()[0] += tol*n2[0];
        poly1[i].Coord()[1] += tol*n2[1];
        poly1[i].Coord()[2] += tol*n2[2];
      }
      else
      {
        // do nothing, point is not very close
      }
    }
  }

  for (int i=0;i<(int)poly2.size();++i)
  {
    for (int j=0;j<(int)poly1.size();++j)
    {
      // we need diff vector and edge1 first
      double diff2[3] = {0.0, 0.0, 0.0};
      double edge1[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff2[k] = poly2[i].Coord()[k] - poly1[j].Coord()[k];
        edge1[k] = (poly1[j].Next())->Coord()[k] - poly1[j].Coord()[k];
      }

      // check if point of poly2 lies within [0,1] for edge1
      double checkalpha = diff2[0]*edge1[0]+diff2[1]*edge1[1]+diff2[2]*edge1[2];
      checkalpha /= (edge1[0]*edge1[0]+edge1[1]*edge1[1]+edge1[2]*edge1[2]);

      // proceed only if inside [0,1] with tolerance tol
      if (checkalpha<-tol || checkalpha>1+tol) continue;

      // compute distance from point on poly2 to edge1
      double n1[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      double ln = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
      for (int k=0;k<3;++k) n1[k] /= ln;

      double dist = diff2[0]*n1[0]+diff2[1]*n1[1]+diff2[2]*n1[2];

      // move point away if very close to edge 2
      if (dist > -tol && dist < 0)
      {
        //cout << "Vertex " << i << " on poly2 is very close to edge " << j << " of poly1 -> moved inside!" << endl;
        poly2[i].Coord()[0] -= tol*n1[0];
        poly2[i].Coord()[1] -= tol*n1[1];
        poly2[i].Coord()[2] -= tol*n1[2];
      }
      else if (dist < tol && dist >= 0)
      {
        //cout << "Vertex " << i << " on poly2 is very close to edge " << j << " of poly1 -> moved outside!" << endl;
        poly2[i].Coord()[0] += tol*n1[0];
        poly2[i].Coord()[1] += tol*n1[1];
        poly2[i].Coord()[2] += tol*n1[2];
      }
      else
      {
        // do nothing, point is not very close
      }
    }
  }

  //**********************************************************************
  // STEP4: Perform line intersection of all edge pairs
  // - this yields two new vectors of intersection vertices
  // - by default the respective edge end vertices are assumed to be
  //   the next/prev vertices and connectivity is set up accordingly
  //**********************************************************************
  vector<Vertex> intersec1;
  vector<Vertex> intersec2;

  for (int i=0;i<(int)poly1.size();++i)
  {
    for (int j=0;j<(int)poly2.size();++j)
    {
      // we need two diff vectors and edges first
      double diffp1[3] = {0.0, 0.0, 0.0};
      double diffp2[3] = {0.0, 0.0, 0.0};
      double diffq1[3] = {0.0, 0.0, 0.0};
      double diffq2[3] = {0.0, 0.0, 0.0};
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diffp1[k] = poly1[i].Coord()[k] - poly2[j].Coord()[k];
        diffp2[k] = (poly1[i].Next())->Coord()[k] - poly2[j].Coord()[k];
        diffq1[k] = poly2[j].Coord()[k] - poly1[i].Coord()[k];
        diffq2[k] = (poly2[j].Next())->Coord()[k] - poly1[i].Coord()[k];
        edge1[k] = (poly1[i].Next())->Coord()[k] - poly1[i].Coord()[k];
        edge2[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
      }

      // outward edge normals of polygon 1 and 2 edges
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];

      // check for parallelity of edges
      double parallel = edge1[0]*n2[0]+edge1[1]*n2[1]+edge1[2]*n2[2];
      if(abs(parallel)<1.0e-12)
      {
        if (out) cout << "WARNING: Detected two parallel edges! (" << i << "," << j << ")" << endl;
        continue;
      }

      if (out) cout << "Searching intersection (" << i << "," << j << ")" << endl;

      // check for intersection of non-parallel edges
      double wec_p1 = 0.0;
      double wec_p2 = 0.0;
      for (int k=0;k<3;++k)
      {
        wec_p1 += (poly1[i].Coord()[k] - poly2[j].Coord()[k]) * n2[k];
        wec_p2 += ((poly1[i].Next())->Coord()[k] - poly2[j].Coord()[k]) * n2[k];
      }

      if (wec_p1*wec_p2<=0)
      {
        double wec_q1 = 0.0;
        double wec_q2 = 0.0;
        for (int k=0;k<3;++k)
        {
          wec_q1 += (poly2[j].Coord()[k] - poly1[i].Coord()[k]) * n1[k];
          wec_q2 += ((poly2[j].Next())->Coord()[k] - poly1[i].Coord()[k]) * n1[k];
        }

        if (wec_q1*wec_q2<=0)
        {
          double alphap = wec_p1/(wec_p1-wec_p2);
          double alphaq = wec_q1/(wec_q1-wec_q2);
          vector<double> ip(3);
          vector<double> iq(3);
          for (int k=0;k<3;++k)
          {
            ip[k] = (1-alphap) * poly1[i].Coord()[k] + alphap * (poly1[i].Next())->Coord()[k];
            iq[k] = (1-alphaq) * poly2[j].Coord()[k] + alphaq * (poly2[j].Next())->Coord()[k];
            if (abs(ip[k])<1.0e-12) ip[k] = 0.0;
            if (abs(iq[k])<1.0e-12) iq[k] = 0.0;
          }

          if (out)
          {
            cout << "Found intersection! (" << i << "," << j << ") " << alphap << " " << alphaq << endl;
            cout << "On Polygon 1: " << ip[0] << " " << ip[1] << " " << ip[2] << endl;
            cout << "On Polygon 2: " << iq[0] << " " << iq[1] << " " << iq[2] << endl;
          }
          
          // generate vectors of underlying node ids for lineclip (2x slave, 2x master)
          vector<int> lcids(4);
          lcids[0] = (int)(poly1[i].Nodeids()[0]);
          lcids[1] = (int)((poly1[i].Next())->Nodeids()[0]);
          lcids[2] = (int)(poly2[j].Nodeids()[0]);
          lcids[3] = (int)((poly2[j].Next())->Nodeids()[0]);

          // store intersection points
          intersec1.push_back(Vertex(ip,Vertex::lineclip,lcids,poly1[i].Next(),&poly1[i],true,false,NULL,alphap));
          intersec2.push_back(Vertex(iq,Vertex::lineclip,lcids,poly2[j].Next(),&poly2[j],true,false,NULL,alphaq));
        }
      }
    }
  }

  /*
  // check slave points
  cout << "\nTesting slave element points" << endl;
  for (int i=0;i<(int)poly1.size();++i)
  {
    Vertex& testv = poly1[i];
    cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
    cout << "Type: " << testv.VType() << endl;
    cout << "Alpha: " << testv.Alpha() << endl;
    cout << "Node id: " << testv.Nodeids()[0] << endl << endl;
  }
  // check master points
  cout << "\nTesting master element points" << endl;
  for (int i=0;i<(int)poly2.size();++i)
  {
    Vertex& testv = poly2[i];
    cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
    cout << "Type: " << testv.VType() << endl;
    cout << "Alpha: " << testv.Alpha() << endl;
    cout << "Node id: " << testv.Nodeids()[0] << endl << endl;
  }

  // check intersection points
  cout << "\nTesting slave intersection points" << endl;
  for (int i=0;i<(int)intersec1.size();++i)
  {
    Vertex& testv = intersec1[i];
    cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
    cout << "Type: " << testv.VType() << endl;
    cout << "Alpha: " << testv.Alpha() << endl;
    cout << "Lineclip ids: " << testv.Nodeids()[0] << " " << testv.Nodeids()[1]
         << " " << testv.Nodeids()[2] << " " << testv.Nodeids()[3] << endl << endl;
  }
  // check intersection points
  cout << "\nTesting master intersection points" << endl;
  for (int i=0;i<(int)intersec2.size();++i)
  {
    Vertex& testv = intersec2[i];
    cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
    cout << "Type: " << testv.VType() << endl;
    cout << "Alpha: " << testv.Alpha() << endl;
    cout << "Lineclip ids: " << testv.Nodeids()[0] << " " << testv.Nodeids()[1]
         << " " << testv.Nodeids()[2] << " " << testv.Nodeids()[3] << endl << endl;
  }
  */

  //**********************
  // do clipping
  //**********************
  if ((int)respoly.size()!=0) dserror("ERROR: PolygonClipping: Respoly!=0 at beginning...");

  //**********************************************************************
  // STEP5: Find result polygon for no intersection case
  // - if there are no intersections polygon 1 could lie within polygon 2,
  //   polygon 2 could lie within polygon 1 or they are fully adjacent
  // - by default the respective edge end vertices are assumed to be
  //   the next/prev vertices and connectivity is set up accordingly
  //**********************************************************************
  if ((int)intersec1.size()==0 && (int)intersec2.size()==0)
  {
    // if they are nested the clip polygon is the inner input polygon
    bool poly1inner = true;
    bool poly2inner = true;

    // (A) check if poly1[0] inside poly2
    for (int i=0;i<(int)poly2.size();++i)
    {
      double edge[3] = {0.0, 0.0, 0.0};
      double diff[3] = {0.0, 0.0, 0.0};

      for (int k=0;k<3;++k)
      {
        edge[k] = (poly2[i].Next())->Coord()[k] - poly2[i].Coord()[k];
        diff[k] = poly1[0].Coord()[k] - poly2[i].Coord()[k];
      }

      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

      double check = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];

      // if check>0 then poly1[0] NOT inside poly2
      if (check>0)
      {
        poly1inner = false;
        break;
      }
    }

    if (poly1inner==true)
    {
      if (out) cout << "Polygon S lies fully inside polygon M!" << endl;
      respoly = poly1;
    }

    else
    {
      // (A) check if poly2[0] inside poly1
      for (int i=0;i<(int)poly1.size();++i)
      {
        double edge[3] = {0.0, 0.0, 0.0};
        double diff[3] = {0.0, 0.0, 0.0};

        for (int k=0;k<3;++k)
        {
          edge[k] = (poly1[i].Next())->Coord()[k] - poly1[i].Coord()[k];
          diff[k] = poly2[0].Coord()[k] - poly1[i].Coord()[k];
        }

        double n[3] = {0.0, 0.0, 0.0};
        n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
        n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
        n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

        double check = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];

        // if check>0 then poly2[0] NOT inside poly1
        if (check>0)
        {
          poly2inner = false;
          break;
        }
      }

      if (poly2inner==true)
      {
        if (out) cout << "Polygon M lies fully inside polygon S!" << endl;
        respoly = poly2;
      }

      // fully adjacent case
      else
      {
        if (out) cout << "Polygons S and M are fully adjacent!" << endl;
        vector<Vertex> empty;
        respoly = empty;
      }
    }
  }

  // invalid case
  else if ((int)intersec1.size()*(int)intersec2.size()==0)
  {
    dserror("ERROR: Found intersection points only on one input polygon...?");
  }

  //**********************************************************************
  // STEP6: Find result polygon for intersection case
  // - assign neighbor connectivity for intersection points
  // - check for edges where 2 intersection points have been found
  // - establish new connectivity (next/prev) accordingly: first for
  //   the intersection points, then for the adjacent edge nodes
  // - perform entry exit classification of intersection points
  // - build result polygon (path finding through linked data structures)
  // - check for sanity of result polygon (last==first)
  // - reorder result polygon in c-clockwise direction if necessary
  // - check if result polygon is convex
  //**********************************************************************
  else
  {
    // assign neighbor intersection nodes on the respective other polygon
    for (int i=0;i<(int)intersec1.size();++i)
      intersec1[i].AssignNeighbor(&intersec2[i]);
    for (int i=0;i<(int)intersec2.size();++i)
      intersec2[i].AssignNeighbor(&intersec1[i]);

    // check all edges for double intersections
    for (int i=0;i<(int)poly1.size();++i)
    {
      vector<Vertex*> dis;
      for (int z=0;z<(int)intersec1.size();++z)
      {
       if (intersec1[z].Next()==poly1[i].Next() && intersec1[z].Prev()==&poly1[i])
       {
         dis.push_back(&intersec1[z]);
       }
      }

      if ((int)dis.size()<2) continue;
      if ((int)dis.size()>2) dserror("ERROR: More than 2 intersections on 1 edge impossible!");

      double alpha1 = dis[0]->Alpha();
      double alpha2 = dis[1]->Alpha();

      if (alpha1<alpha2)
      {
        // ordering is poly1[i] -> dis[0] -> dis[1] -> poly1[i].Next()
        dis[0]->AssignNext(dis[1]);
        dis[1]->AssignPrev(dis[0]);
      }
      else if (alpha1==alpha2)
      {
        dserror("ERROR: Two identical intersection points on 1 edge!");
      }
      else
      {
        // ordering is poly1[i] -> dis[1] -> dis[0] -> poly1[i].Next()
        dis[1]->AssignNext(dis[0]);
        dis[0]->AssignPrev(dis[1]);
      }
    }

    for (int i=0;i<(int)poly2.size();++i)
    {
      vector<Vertex*> dis;
      for (int z=0;z<(int)intersec2.size();++z)
      {
       if (intersec2[z].Next()==poly2[i].Next() && intersec2[z].Prev()==&poly2[i])
       {
         dis.push_back(&intersec2[z]);
       }
      }

      if ((int)dis.size()<2) continue;
      if ((int)dis.size()>2) dserror("ERROR: More than 2 intersections on 1 edge impossible!");

      double alpha1 = dis[0]->Alpha();
      double alpha2 = dis[1]->Alpha();

      if (alpha1<alpha2)
      {
        // ordering is poly2[i] -> dis[0] -> dis[1] -> poly2[i].Next()
        dis[0]->AssignNext(dis[1]);
        dis[1]->AssignPrev(dis[0]);
      }
      else if (alpha1==alpha2)
      {
        dserror("ERROR: Two identical intersection points on 1 edge!");
      }
      else
      {
        // ordering is poly2[i] -> dis[1] -> dis[0] -> poly2[i].Next()
        dis[1]->AssignNext(dis[0]);
        dis[0]->AssignPrev(dis[1]);
      }
    }

    // assign new next / previous nodes for vertices near intersections
    for (int i=0;i<(int)intersec1.size();++i)
    {
      (intersec1[i].Prev())->AssignNext(&intersec1[i]);
      (intersec1[i].Next())->AssignPrev(&intersec1[i]);
    }
    for (int i=0;i<(int)intersec2.size();++i)
    {
      (intersec2[i].Prev())->AssignNext(&intersec2[i]);
      (intersec2[i].Next())->AssignPrev(&intersec2[i]);
    }

    // perform entry / exit classification of intersections
    // we move along both polygons and determine whether each intersection
    // point is an entry or exit point with respect to the other polygon.
    // this status is then stored into the vertex data structure.

    for (int i=0;i<(int)intersec1.size();++i)
    {
      // check if previous vertex is inside for first intersection
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        edge1[k] = (intersec1[i].Next())->Coord()[k] - (intersec1[i].Prev())->Coord()[k];
        edge2[k] = ((intersec1[i].Neighbor())->Next())->Coord()[k] - ((intersec1[i].Neighbor())->Prev())->Coord()[k];
      }
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];

      double check = edge1[0]*n2[0] + edge1[1]*n2[1] + edge1[2]*n2[2];
      if (check<0) intersec1[i].EntryExit()=true;
    }

    for (int i=0;i<(int)intersec1.size();++i)
    {
      // check if previous vertex is inside for first intersection
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        edge1[k] = (intersec2[i].Next())->Coord()[k] - (intersec2[i].Prev())->Coord()[k];
        edge2[k] = ((intersec2[i].Neighbor())->Next())->Coord()[k] - ((intersec2[i].Neighbor())->Prev())->Coord()[k];
      }
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];

      double check = edge1[0]*n2[0] + edge1[1]*n2[1] + edge1[2]*n2[2];
      if (check<0) intersec2[i].EntryExit()=true;
    }

    // print intersection points and their status
    if (out)
    {
      for (int i=0;i<(int)intersec1.size();++i)
      {
        cout << "Intersec1: " << i << " " << intersec1[i].Coord()[0] << " " << intersec1[i].Coord()[1] << " " << intersec1[i].Coord()[2];
        cout << " EntryExit: " << intersec1[i].EntryExit() << endl;
      }
  
      //cout << endl;
      for (int i=0;i<(int)intersec2.size();++i)
      {
        cout << "Intersec2: " << i << " " <<  intersec2[i].Coord()[0] << " " << intersec2[i].Coord()[1] << " " << intersec2[i].Coord()[2];
        cout << " EntryExit: " << intersec2[i].EntryExit() << endl;
      }
    }

    // create clipped polygon by filtering
    // We simply have to find our way through the linked data structures of
    // poly1, poly2 and intersection vertices. For this we start at an
    // intersection point and move on according to the entry / exit status.
    // When we reach the next intersection point we jump to the other polygon
    // according to the neighboring vertex pointer.
    // The result will be an ordered list of vertices of the clipped polygon!
    Vertex* current = &intersec1[0];

    // push_back start Vertex coords into result polygon
    if (out) cout << "\nStart loop on Slave at " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
    respoly.push_back(Vertex(current->Coord(),Vertex::lineclip,current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));

    do {
      // find next Vertex / Vertices (path)
      if (current->EntryExit()==true)
      {
        if (out) cout << "Intersection was Entry, so move to Next() on same polygon!" << endl;
        do {
          current = current->Next();
          if (out) cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
        } while (current->Intersect()==false);
        if (out) cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }
      else
      {
        if (out) cout << "Intersection was Exit, so move to Prev() on same polygon!" << endl;
        do {
          current = current->Prev();
          if (out) cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
        } while (current->Intersect()==false);
        if (out) cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }

      // jump to the other input polygon
      current = current->Neighbor();
      if (out) cout << "Jumping to other polygon at intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      if (out) cout << "Length of result list so far: " << (int)respoly.size() << endl;

      // check length of result polygon list
      if ((int)respoly.size()>8) dserror("ERROR: Length of result polygon > 8! Something went wrong!");
      
    } while (current!=&intersec1[0] && current!=&intersec2[0]);

    // check if last entry is identical to first entry
    // check on both intersection lists
    double fldiff[3] = {0.0, 0.0, 0.0};
    double fldiff2[3] = {0.0, 0.0, 0.0};
    bool identical = true;
    double fldist = 0.0;
    double fldist2 = 0.0;
    for (int k=0;k<3;++k)
    {
      fldiff[k] = respoly[(int)respoly.size()-1].Coord()[k] - intersec1[0].Coord()[k];
      fldist += fldiff[k]*fldiff[k];
      fldiff2[k] = respoly[(int)respoly.size()-1].Coord()[k] - intersec2[0].Coord()[k];
      fldist2 += fldiff2[k]*fldiff2[k];
    }
    fldist = sqrt(fldist);
    if (fldist>1.0e-8 && fldist2>1.0e-8) identical = false;

    // remove last entry if so, throw dserror if not so
    if (identical) respoly.pop_back();
    else
    {
      cout << "\nDifference Dists: " << fldist << " " << fldist2 << endl;
      dserror("ERROR: We did not arrive at the starting point again...?");
    }

    // collapse respoly points that are very close
    vector<Vertex> collapsedrespoly;
    for (int i=0;i<(int)respoly.size();++i)
    {
      // find distance between two consecutive points
      // first point of respoly
      if (i==0)
        collapsedrespoly.push_back(respoly[i]);

      // last point of respoly
      else if (i==(int)respoly.size()-1)
      {
        double diff[3] = {0.0, 0.0, 0.0};
        double diff2[3] = {0.0, 0.0, 0.0};

        for (int k=0;k<3;++k) diff[k] = respoly[i].Coord()[k] - respoly[i-1].Coord()[k];
        for (int k=0;k<3;++k) diff2[k] = respoly[0].Coord()[k] - respoly[i].Coord()[k];

        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        double dist2 = sqrt(diff2[0]*diff2[0]+diff2[1]*diff2[1]+diff2[2]*diff2[2]);
        double tolcollapse = 1.0e6*tol;

        if (abs(dist) >= tolcollapse && abs(dist2) >= tolcollapse)
          collapsedrespoly.push_back(respoly[i]);
        else
          if (out) cout << "Collapsed two points in result polygon!" << endl;
      }

      // standard case
      else
      {
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k) diff[k] = respoly[i].Coord()[k] - respoly[i-1].Coord()[k];

        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        double tolcollapse = 1.0e6*tol;

        if (abs(dist) >= tolcollapse)
          collapsedrespoly.push_back(respoly[i]);
        else
          if (out) cout << "Collapsed two points in result polygon!" << endl;
      }
    }

    // replace respoly by collapsed respoly
    respoly = collapsedrespoly;
    if (out) cout << "Final length of result list: " << (int)respoly.size() << endl;

    // check if respoly collapsed to nothing
    if ((int)collapsedrespoly.size()<3)
    {
      if (out) cout << "Collapsing of result polygon led to < 3 vertices -> no respoly!" << endl;
      vector<Vertex> empty;
      respoly = empty;
    }

    // check for rotation of result polygon (must be clockwise!!!)
    // first get geometric center
    double center[3] = {0.0, 0.0, 0.0};

    for (int i=0;i<(int)respoly.size();++i)
      for (int k=0;k<3;++k)
        center[k] += respoly[i].Coord()[k]/((int)respoly.size());

    if (out) cout << "\nCenter ResPoly: " << center[0] << " " << center[1] << " " << center[2] << endl;

    // then we compute the clockwise plane normal
    double diff[3] = {0.0, 0.0, 0.0};
    double edge[3] = {0.0, 0.0, 0.0};

    for (int k=0;k<3;++k)
    {
      diff[k] = respoly[0].Coord()[k]-center[k];
      edge[k] = respoly[1].Coord()[k]-respoly[0].Coord()[k];
    }

    double cross[3] = {0.0, 0.0, 0.0};

    cross[0] = diff[1]*edge[2]-diff[2]*edge[1];
    cross[1] = diff[2]*edge[0]-diff[0]*edge[2];
    cross[2] = diff[0]*edge[1]-diff[1]*edge[0];

    // check against auxiliary plane normal
    double check = cross[0]*Auxn()[0]+cross[1]*Auxn()[1]+cross[2]*Auxn()[2];

    if (check<0)
    {
      // reorder result polygon in clockwise direction
      // cout << "Result polygon not ordered counter-clockwise -> reordered!" << endl;
      std::reverse(respoly.begin(),respoly.end());
    }

    if (out)
    {
      // print final input polygons to screen
      cout << "\nResult Poylgon:";
      for (int i=0;i<(int)respoly.size();++i)
        cout << "\nVertex " << i << ":\t" << respoly[i].Coord()[0] << "\t" << respoly[i].Coord()[1] << "\t" << respoly[i].Coord()[2];
      cout << endl;
    }
    
    // check if result polygon is convex
    // a polygon is convex if the scalar product of an edge normal and the
    // next edge direction is negative for all edges
    for (int i=0;i<(int)respoly.size();++i)
    {
      // we need the edge vector first
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        if (i!=(int)respoly.size()-1) edge[k] = respoly[i+1].Coord()[k] - respoly[i].Coord()[k];
        else edge[k] = respoly[0].Coord()[k] - respoly[i].Coord()[k];
      }
      // edge normal is result of cross product
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

      // we need the next edge vector now
      double nextedge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        if (i<(int)respoly.size()-2) nextedge[k] = respoly[i+2].Coord()[k] - respoly[i+1].Coord()[k];
        else if (i==(int)respoly.size()-2) nextedge[k] = respoly[0].Coord()[k] - respoly[i+1].Coord()[k];
        else nextedge[k] = respoly[1].Coord()[k] - respoly[0].Coord()[k];
      }
      // check scalar product
      double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
      if (check>0) dserror("ERROR: Result polygon not convex!");
    }
  }
  
  if (out)
  {
    // **********************************************************************
    // STEP6: Result visualization with GMSH
    // - plot the two input polygons and their vertex numbering
    // - plot the result polygon and its vertex numbering
    // **********************************************************************
    std::ostringstream filename;
    static int gmshcount=0;
    filename << "o/gmsh_output/" << "clipping_";
    if (gmshcount<10)
      filename << 0 << 0 << 0 << 0;
    else if (gmshcount<100)
      filename << 0 << 0 << 0;
    else if (gmshcount<1000)
      filename << 0 << 0;
    else if (gmshcount<10000)
      dserror("Gmsh output implemented for a maximum of 9.999 clip polygons");
    filename << gmshcount << ".pos";
    gmshcount++;
  
    // do output to file in c-style
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream gmshfilecontent;
    gmshfilecontent << "View \" Clipping \" {" << endl;
  
    for (int i=0;i<(int)poly1.size();++i)
    {
      if (i!=(int)poly1.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << ","
                                 << poly1[i].Coord()[2] << "," << poly1[i+1].Coord()[0] << "," << poly1[i+1].Coord()[1] << ","
                                 << poly1[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << ","
                                 << poly1[i].Coord()[2] << "," << poly1[0].Coord()[0] << "," << poly1[0].Coord()[1] << ","
                                 << poly1[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << "," << poly1[i].Coord()[2] << "," << 17 << ")";
      gmshfilecontent << "{" << "S" << i << "};" << endl;
    }
  
    for (int i=0;i<(int)poly2.size();++i)
    {
      if (i!=(int)poly2.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << ","
                                 << poly2[i].Coord()[2] << "," << poly2[i+1].Coord()[0] << "," << poly2[i+1].Coord()[1] << ","
                                 << poly2[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << ","
                                 << poly2[i].Coord()[2] << "," << poly2[0].Coord()[0] << "," << poly2[0].Coord()[1] << ","
                                 << poly2[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << "," << poly2[i].Coord()[2] << "," << 17 << ")";
      gmshfilecontent << "{" << "M" << i << "};" << endl;
    }
  
    for (int i=0;i<(int)respoly.size();++i)
    {
      if (i!=(int)respoly.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << ","
                                 << respoly[i].Coord()[2] << "," << respoly[i+1].Coord()[0] << "," << respoly[i+1].Coord()[1] << ","
                                 << respoly[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << ","
                                 << respoly[i].Coord()[2] << "," << respoly[0].Coord()[0] << "," << respoly[0].Coord()[1] << ","
                                 << respoly[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << "," << respoly[i].Coord()[2] << "," << 27 << ")";
      gmshfilecontent << "{" << "R" << i << "};" << endl;
    }
  
  //  for (int i=0;i<(int)intersec1.size();++i)
  //  {
  //    gmshfilecontent << "T3(" << scientific << intersec1[i].Coord()[0] << "," << intersec1[i].Coord()[1] << "," << intersec1[i].Coord()[2] << "," << 17 << ")";
  //    if (intersec1[i].EntryExit()==true && intersec2[i].EntryExit()==true) gmshfilecontent << "{" << "SEME" << "};" << endl;
  //    else if (intersec1[i].EntryExit()==false && intersec2[i].EntryExit()==true) gmshfilecontent << "{" << "SXME" << "};" << endl;
  //    else if (intersec1[i].EntryExit()==true && intersec2[i].EntryExit()==false) gmshfilecontent << "{" << "SEMX" << "};" << endl;
  //    else gmshfilecontent << "{" << "SXMX" << "};" << endl;
  //  }
  
    gmshfilecontent << "};" << endl;
  
    // move everything to gmsh post-processing file and close it
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Clipping of two polygons (NEW version)                    popp 11/09|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::PolygonClippingConvexHull(vector<Vertex>& poly1,
                                                    vector<Vertex>& poly2,
                                                    vector<Vertex>& respoly,
                                                    double& tol)
{
  // choose output
  bool out = false;
  
  //**********************************************************************
  // STEP1: Input check
  // - input polygons must consist of min. 3 vertices each
  // - rotation of poly1 must be c-clockwise w.r.t. (0,0,1) or Auxn()
  // - rotation of poly 2 changed to c-clockwise w.r.t. (0,0,1) or Auxn()
  // - both input polygons must be convex
  //**********************************************************************

  // check input variables
  if ((int)poly1.size()<3 || (int)poly2.size()<3)
    dserror("ERROR: Input Polygons must consist of min. 3 vertices each");

  // check for rotation of polygon1 (slave) and polgon 2 (master)
  // note that we implicitly already rely on convexity here!
  // first get geometric centers of polygon1 and polygon2
  double center1[3] = {0.0, 0.0, 0.0};
  double center2[3] = {0.0, 0.0, 0.0};

  for (int i=0;i<(int)poly1.size();++i)
    for (int k=0;k<3;++k)
      center1[k] += poly1[i].Coord()[k]/((int)poly1.size());

  for (int i=0;i<(int)poly2.size();++i)
    for (int k=0;k<3;++k)
      center2[k] += poly2[i].Coord()[k]/((int)poly2.size());

  if (out)
  {
    cout << "Center 1: " << center1[0] << " " << center1[1] << " " << center1[2] << endl;
    cout << "Center 2: " << center2[0] << " " << center2[1] << " " << center2[2] << endl;
  }
  
  // then we compute the counter-clockwise plane normal
  double diff1[3] = {0.0, 0.0, 0.0};
  double edge1[3] = {0.0, 0.0, 0.0};
  double diff2[3] = {0.0, 0.0, 0.0};
  double edge2[3] = {0.0, 0.0, 0.0};

  for (int k=0;k<3;++k)
  {
    diff1[k] = poly1[0].Coord()[k]-center1[k];
    edge1[k] = poly1[1].Coord()[k]-poly1[0].Coord()[k];
    diff2[k] = poly2[0].Coord()[k]-center2[k];
    edge2[k] = poly2[1].Coord()[k]-poly2[0].Coord()[k];
  }

  double cross1[3] = {0.0, 0.0, 0.0};
  double cross2[3] = {0.0, 0.0, 0.0};

  cross1[0] = diff1[1]*edge1[2]-diff1[2]*edge1[1];
  cross1[1] = diff1[2]*edge1[0]-diff1[0]*edge1[2];
  cross1[2] = diff1[0]*edge1[1]-diff1[1]*edge1[0];

  cross2[0] = diff2[1]*edge2[2]-diff2[2]*edge2[1];
  cross2[1] = diff2[2]*edge2[0]-diff2[0]*edge2[2];
  cross2[2] = diff2[0]*edge2[1]-diff2[1]*edge2[0];

  // check against auxiliary plane normal
  double check1 = cross1[0]*Auxn()[0]+cross1[1]*Auxn()[1]+cross1[2]*Auxn()[2];
  double check2 = cross2[0]*Auxn()[0]+cross2[1]*Auxn()[1]+cross2[2]*Auxn()[2];

  // check polygon 1 and throw dserror if not c-clockwise
  if (check1<=0) dserror("ERROR: Polygon 1 (slave) not ordered counter-clockwise!");

  // check polygon 2 and reorder in c-clockwise direction
  if (check2<0)
  {
    if (out) cout << "Polygon 2 (master) not ordered counter-clockwise -> reordered!" << endl;
    std::reverse(poly2.begin(), poly2.end());
  }

  // check if the two input polygons are convex
  // a polygon is convex if the scalar product of an edge normal and the
  // next edge direction is negative for all edges
  for (int i=0;i<(int)poly1.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly1.size()-1) edge[k] = poly1[i+1].Coord()[k] - poly1[i].Coord()[k];
      else edge[k] = poly1[0].Coord()[k] - poly1[i].Coord()[k];
    }

    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly1.size()-2) nextedge[k] = poly1[i+2].Coord()[k] - poly1[i+1].Coord()[k];
      else if (i==(int)poly1.size()-2) nextedge[k] = poly1[0].Coord()[k] - poly1[i+1].Coord()[k];
      else nextedge[k] = poly1[1].Coord()[k] - poly1[0].Coord()[k];
    }

    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) dserror("ERROR: Input polygon 1 not convex");
  }

  for (int i=0;i<(int)poly2.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly2.size()-1) edge[k] = poly2[i+1].Coord()[k] - poly2[i].Coord()[k];
      else edge[k] = poly2[0].Coord()[k] - poly2[i].Coord()[k];
    }

    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];

    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly2.size()-2) nextedge[k] = poly2[i+2].Coord()[k] - poly2[i+1].Coord()[k];
      else if (i==(int)poly2.size()-2) nextedge[k] = poly2[0].Coord()[k] - poly2[i+1].Coord()[k];
      else nextedge[k] = poly2[1].Coord()[k] - poly2[0].Coord()[k];
    }

    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) 
    {
      // this may happen, so do NOT throw an error immediately
      // but instead check if the two elements to be clipped are
      // close to each other at all. If so, then really throw the
      // dserror, if not, simply continue with the next pair!
      int sid = SlaveElement().Id();
      int mid = MasterElement().Id();
      cout << "***WARNING*** Input polygon 2 not convex! (S/M-pair: " << sid << "/" << mid << ")" << endl;
      bool nearcheck = RoughCheckNodes();
      if (nearcheck && out)
      {
        //cout << "***WARNING*** Input polygon 2 not convex, but close pair!" << endl;
        std::ostringstream filename;
        static int problemcount=0;
        filename << "o/gmsh_output/" << "problem_";
        if (problemcount<10)
          filename << 0 << 0 << 0 << 0;
        else if (problemcount<100)
          filename << 0 << 0 << 0;
        else if (problemcount<1000)
          filename << 0 << 0;
        else if (problemcount<10000)
          filename << 0;
        else
          dserror("Gmsh output implemented for a maximum of 9.999 problem polygons");
        filename << problemcount << "_" << sid << "_" << mid << ".pos";
        problemcount++;
      
        // do output to file in c-style
        FILE* fp = NULL;
        fp = fopen(filename.str().c_str(), "w");
        std::stringstream gmshfilecontent;
        gmshfilecontent << "View \" ProblemClipping S" <<  sid << "/M" << mid << " \" {" << endl;
      
        // get slave and master nodes
        int nsnodes = SlaveIntElement().NumNode();
        DRT::Node** mysnodes = SlaveIntElement().Nodes();
        if (!mysnodes) dserror("ERROR: Null pointer!");
        vector<MortarNode*> mycsnodes(nsnodes);
        for (int i=0;i<nsnodes;++i)
          mycsnodes[i] = static_cast<MortarNode*> (mysnodes[i]);
        int nmnodes = MasterIntElement().NumNode();
        DRT::Node** mymnodes = MasterIntElement().Nodes();
        if (!mymnodes) dserror("ERROR: Null pointer!");
        vector<MortarNode*> mycmnodes(nmnodes);
        for (int i=0;i<nmnodes;++i)
          mycmnodes[i] = static_cast<MortarNode*> (mymnodes[i]);
          
        for (int i=0;i<nsnodes;++i)
        {
          if (i!=nsnodes-1)
          {
            gmshfilecontent << "SL(" << scientific << mycsnodes[i]->xspatial()[0] << "," << mycsnodes[i]->xspatial()[1] << ","
                                     << mycsnodes[i]->xspatial()[2] << "," << mycsnodes[i+1]->xspatial()[0] << "," << mycsnodes[i+1]->xspatial()[1] << ","
                                     << mycsnodes[i+1]->xspatial()[2] << ")";
            gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
          }
          else
          {
            gmshfilecontent << "SL(" << scientific << mycsnodes[i]->xspatial()[0] << "," << mycsnodes[i]->xspatial()[1] << ","
                                     << mycsnodes[i]->xspatial()[2] << "," << mycsnodes[0]->xspatial()[0] << "," << mycsnodes[0]->xspatial()[1] << ","
                                     <<mycsnodes[0]->xspatial()[2] << ")";
            gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
      
          }
          gmshfilecontent << "T3(" << scientific << mycsnodes[i]->xspatial()[0] << "," << mycsnodes[i]->xspatial()[1] << "," << mycsnodes[i]->xspatial()[2] << "," << 17 << ")";
          gmshfilecontent << "{" << "S" << i << "};" << endl;
        }
      
        for (int i=0;i<nmnodes;++i)
        {
          if (i!=nmnodes-1)
          {
            gmshfilecontent << "SL(" << scientific << mycmnodes[i]->xspatial()[0] << "," << mycmnodes[i]->xspatial()[1] << ","
                                     << mycmnodes[i]->xspatial()[2] << "," << mycmnodes[i+1]->xspatial()[0] << "," << mycmnodes[i+1]->xspatial()[1] << ","
                                     << mycmnodes[i+1]->xspatial()[2] << ")";
            gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
          }
          else
          {
            gmshfilecontent << "SL(" << scientific << mycmnodes[i]->xspatial()[0] << "," << mycmnodes[i]->xspatial()[1] << ","
                                     << mycmnodes[i]->xspatial()[2] << "," << mycmnodes[0]->xspatial()[0] << "," << mycmnodes[0]->xspatial()[1] << ","
                                     <<mycmnodes[0]->xspatial()[2] << ")";
            gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
      
          }
          gmshfilecontent << "T3(" << scientific << mycmnodes[i]->xspatial()[0] << "," << mycmnodes[i]->xspatial()[1] << "," << mycmnodes[i]->xspatial()[2] << "," << 17 << ")";
          gmshfilecontent << "{" << "M" << i << "};" << endl;
        }
        
        gmshfilecontent << "};" << endl;
      
        // move everything to gmsh post-processing file and close it
        fprintf(fp,gmshfilecontent.str().c_str());
        fclose(fp);
      }
      
      return false;
      
      //if (nearcheck) dserror("ERROR: Input polygon 2 not convex, but close pair!");
      //else return false;
    }
  }

  // print final input polygons to screen
  if (out)
  {
    cout << "\nInput Poylgon 1:";
    for (int i=0;i<(int)poly1.size();++i)
      cout << "\nVertex " << i << ":\t" << scientific << poly1[i].Coord()[0]
           << "\t" << poly1[i].Coord()[1] << "\t" << poly1[i].Coord()[2];
  
    cout << "\nInput Poylgon 2:";
    for (int i=0;i<(int)poly2.size();++i)
      cout << "\nVertex " << i << ":\t" << scientific << poly2[i].Coord()[0]
           << "\t" << poly2[i].Coord()[1] << "\t" << poly2[i].Coord()[2];
  
    cout << endl << endl;
  }

  //**********************************************************************
  // STEP2: Extend Vertex data structures
  // - note that poly1 is the slave element and poly2 the master element
  // - assign Next() and Prev() pointers to initialize linked structure
  //**********************************************************************

  // set previous and next Vertex pointer for all elements in lists
  for (int i=0;i<(int)poly1.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly1.size()-1)
    {
      poly1[i].AssignNext(&poly1[i+1]);
      poly1[i].AssignPrev(&poly1[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly1[i].AssignNext(&poly1[i+1]);
      poly1[i].AssignPrev(&poly1[(int)poly1.size()-1]);
    }
    // last element in list
    else
    {
      poly1[i].AssignNext(&poly1[0]);
      poly1[i].AssignPrev(&poly1[i-1]);
    }
  }
  for (int i=0;i<(int)poly2.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly2.size()-1)
    {
      poly2[i].AssignNext(&poly2[i+1]);
      poly2[i].AssignPrev(&poly2[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly2[i].AssignNext(&poly2[i+1]);
      poly2[i].AssignPrev(&poly2[(int)poly2.size()-1]);
    }
    // last element in list
    else
    {
      poly2[i].AssignNext(&poly2[0]);
      poly2[i].AssignPrev(&poly2[i-1]);
    }
  }

  //**********************************************************************
  // STEP3: Perform line intersection of all edge pairs
  // - this yields a new vector of intersection vertices
  // - by default the respective edge end vertices are assumed to be
  //   the next/prev vertices and connectivity is set up accordingly
  //**********************************************************************
  vector<Vertex> intersec;

  for (int i=0;i<(int)poly1.size();++i)
  {
    for (int j=0;j<(int)poly2.size();++j)
    {
      // we need two diff vectors and edges first
      double diffp1[3] = {0.0, 0.0, 0.0};
      double diffp2[3] = {0.0, 0.0, 0.0};
      double diffq1[3] = {0.0, 0.0, 0.0};
      double diffq2[3] = {0.0, 0.0, 0.0};
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diffp1[k] = poly1[i].Coord()[k] - poly2[j].Coord()[k];
        diffp2[k] = (poly1[i].Next())->Coord()[k] - poly2[j].Coord()[k];
        diffq1[k] = poly2[j].Coord()[k] - poly1[i].Coord()[k];
        diffq2[k] = (poly2[j].Next())->Coord()[k] - poly1[i].Coord()[k];
        edge1[k] = (poly1[i].Next())->Coord()[k] - poly1[i].Coord()[k];
        edge2[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
      }

      // outward edge normals of polygon 1 and 2 edges
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];

      // check for parallelity of edges
      double parallel = edge1[0]*n2[0]+edge1[1]*n2[1]+edge1[2]*n2[2];
      if (abs(parallel)<tol)
      {
        if (out) cout << "WARNING: Detected two parallel edges! (" << i << "," << j << ")" << endl;
        continue;
      }

      if (out) cout << "Searching intersection (" << i << "," << j << ")" << endl;

      // check for intersection of non-parallel edges
      double wec_p1 = 0.0;
      double wec_p2 = 0.0;
      for (int k=0;k<3;++k)
      {
        wec_p1 += (poly1[i].Coord()[k] - poly2[j].Coord()[k]) * n2[k];
        wec_p2 += ((poly1[i].Next())->Coord()[k] - poly2[j].Coord()[k]) * n2[k];
      }

      if (wec_p1*wec_p2<=0)
      {
        double wec_q1 = 0.0;
        double wec_q2 = 0.0;
        for (int k=0;k<3;++k)
        {
          wec_q1 += (poly2[j].Coord()[k] - poly1[i].Coord()[k]) * n1[k];
          wec_q2 += ((poly2[j].Next())->Coord()[k] - poly1[i].Coord()[k]) * n1[k];
        }

        if (wec_q1*wec_q2<=0)
        {
          double alphap = wec_p1/(wec_p1-wec_p2);
          double alphaq = wec_q1/(wec_q1-wec_q2);
          vector<double> ip(3);
          vector<double> iq(3);
          for (int k=0;k<3;++k)
          {
            ip[k] = (1-alphap) * poly1[i].Coord()[k] + alphap * (poly1[i].Next())->Coord()[k];
            iq[k] = (1-alphaq) * poly2[j].Coord()[k] + alphaq * (poly2[j].Next())->Coord()[k];
            if (abs(ip[k]) < tol) ip[k] = 0.0;
            if (abs(iq[k]) < tol) iq[k] = 0.0;
          }

          if (out)
          {
            cout << "Found intersection! (" << i << "," << j << ") " << alphap << " " << alphaq << endl;
            cout << "On Polygon 1: " << ip[0] << " " << ip[1] << " " << ip[2] << endl;
            cout << "On Polygon 2: " << iq[0] << " " << iq[1] << " " << iq[2] << endl;
          }
          
          // generate vectors of underlying node ids for lineclip (2x slave, 2x master)
          vector<int> lcids(4);
          lcids[0] = (int)(poly1[i].Nodeids()[0]);
          lcids[1] = (int)((poly1[i].Next())->Nodeids()[0]);
          lcids[2] = (int)(poly2[j].Nodeids()[0]);
          lcids[3] = (int)((poly2[j].Next())->Nodeids()[0]);

          // store intersection points
          intersec.push_back(Vertex(ip,Vertex::lineclip,lcids,poly1[i].Next(),&poly1[i],true,false,NULL,alphap));
        }
      }
    }
  }
 
  if (out)
  {
    // check intersection points
    cout << "\nTesting intersection points" << endl;
    for (int i=0;i<(int)intersec.size();++i)
    {
      Vertex& testv = intersec[i];
      cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
      cout << "Type: " << testv.VType() << endl;
      cout << "Alpha: " << testv.Alpha() << endl;
      cout << "Lineclip ids: " << testv.Nodeids()[0] << " " << testv.Nodeids()[1]
           << " " << testv.Nodeids()[2] << " " << testv.Nodeids()[3] << endl << endl;
    }
  }
    
  //**********************************************************************
  // STEP4: Collapse line intersections
  // - this yields a collapsed vector of intersection vertices
  // - those intersection points close to poly1/poly2 vertices are deleted
  //**********************************************************************
  vector<Vertex> collintersec;
  
  for (int i=0;i<(int)intersec.size();++i)
  {
    // keep track of comparisons
    bool close = false;
    
    // check against all poly1 (slave) points
    for (int j=0;j<(int)poly1.size();++j)
    {
      // distance vector
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k) diff[k] = intersec[i].Coord()[k] - poly1[j].Coord()[k];
      double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
      
      // only keep intersection point if not close
      if (dist <= tol)
      {
        close = true;
        break;
      }
    }
    
    // do only if no close poly1 (slave) point found
    if (!close)
    {
      // check against all poly2 (master) points
      for (int j=0;j<(int)poly2.size();++j)
      {
        // distance vector
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k) diff[k] = intersec[i].Coord()[k] - poly2[j].Coord()[k];
        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        
        // only keep intersection point if not close
        if (dist <= tol)
        {
          close = true;  
          break;
        }
      }
    }
    
    // keep intersection point only if not close to any poly1/poly2 point
    if (!close) collintersec.push_back(intersec[i]);
  }
  
  if (out)
  {
    // check collapsed intersection points
    cout << "\nTesting collapsed intersection points" << endl;
    for (int i=0;i<(int)collintersec.size();++i)
    {
      Vertex& testv = collintersec[i];
      cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
      cout << "Type: " << testv.VType() << endl;
      cout << "Alpha: " << testv.Alpha() << endl;
      cout << "Lineclip ids: " << testv.Nodeids()[0] << " " << testv.Nodeids()[1]
           << " " << testv.Nodeids()[2] << " " << testv.Nodeids()[3] << endl << endl;
    }
  }
  
  //**********************************************************************
  // STEP5: Create points of convex hull
  // - check all poly1 points against all poly1/poly2 edges
  // - check all poly2 points against all poly1/poly2 edges
  // - check all collintersec points against all poly1/poly2 edges
  // - points outside any poly1/poly2 edge are NOT in the convex hull
  // - as a result we obtain all points forming the convex hull
  //**********************************************************************
  vector<Vertex> convexhull;
  
  //----------------------------------------------------check poly1 points
  for (int i=0;i<(int)poly1.size();++i)
  {
    // keep track of inside / outside status
    bool outside = false;
 
    // check against all poly1 (slave) edges
    for (int j=0;j<(int)poly1.size();++j)
    {
      // we need diff vector and edge2 first
      double diff[3] = {0.0, 0.0, 0.0};
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff[k] = poly1[i].Coord()[k] - poly1[j].Coord()[k];
        edge[k] = (poly1[j].Next())->Coord()[k] - poly1[j].Coord()[k];
      }

      // compute distance from point on poly1 to edge
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
      double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      for (int k=0;k<3;++k) n[k] /= ln;

      double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
      
      // only keep point if not in outside halfspace
      if (dist > tol)
      {
        outside = true;
        break;
      }
    }
    
    // do only if not already outside w.r.t. to a poly1 (slave) edge
    if (!outside)
    {
      // check against all poly2 (master) edges
      for (int j=0;j<(int)poly2.size();++j)
      {
        // we need diff vector and edge2 first
        double diff[3] = {0.0, 0.0, 0.0};
        double edge[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k)
        {
          diff[k] = poly1[i].Coord()[k] - poly2[j].Coord()[k];
          edge[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
        }

        // compute distance from point on poly1 to edge
        double n[3] = {0.0, 0.0, 0.0};
        n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
        n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
        n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
        double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        for (int k=0;k<3;++k) n[k] /= ln;

        double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
        
        // only keep point if not in outside halfspace
        if (dist > tol)
        {
          outside = true;
          break;
        }
      }
    }
    
    // only keep point if never in outside halfspace
    if (!outside) convexhull.push_back(poly1[i]);
  }
  
  //----------------------------------------------------check poly2 points
  for (int i=0;i<(int)poly2.size();++i)
  {
    // keep track of inside / outside status
    bool outside = false;
 
    // check against all poly1 (slave) edges
    for (int j=0;j<(int)poly1.size();++j)
    {
      // we need diff vector and edge2 first
      double diff[3] = {0.0, 0.0, 0.0};
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff[k] = poly2[i].Coord()[k] - poly1[j].Coord()[k];
        edge[k] = (poly1[j].Next())->Coord()[k] - poly1[j].Coord()[k];
      }

      // compute distance from point on poly1 to edge
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
      double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      for (int k=0;k<3;++k) n[k] /= ln;

      double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
      
      // only keep point if not in outside halfspace
      if (dist > tol)
      {
        outside = true;
        break;
      }
    }
    
    // do only if not already outside w.r.t. to a poly1 (slave) edge
    if (!outside)
    {
      // check against all poly2 (master) edges
      for (int j=0;j<(int)poly2.size();++j)
      {
        // we need diff vector and edge2 first
        double diff[3] = {0.0, 0.0, 0.0};
        double edge[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k)
        {
          diff[k] = poly2[i].Coord()[k] - poly2[j].Coord()[k];
          edge[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
        }

        // compute distance from point on poly1 to edge
        double n[3] = {0.0, 0.0, 0.0};
        n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
        n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
        n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
        double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        for (int k=0;k<3;++k) n[k] /= ln;

        double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
        
        // only keep point if not in outside halfspace
        if (dist > tol)
        {
          outside = true;
          break;
        }
      }
    }
    
    // only keep point if never in outside halfspace
    if (!outside) convexhull.push_back(poly2[i]);
  }
  
  //---------------------------------------------check collintersec points
  for (int i=0;i<(int)collintersec.size();++i)
  {
    // keep track of inside / outside status
    bool outside = false;
 
    // check against all poly1 (slave) edges
    for (int j=0;j<(int)poly1.size();++j)
    {
      // we need diff vector and edge2 first
      double diff[3] = {0.0, 0.0, 0.0};
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff[k] = collintersec[i].Coord()[k] - poly1[j].Coord()[k];
        edge[k] = (poly1[j].Next())->Coord()[k] - poly1[j].Coord()[k];
      }

      // compute distance from point on poly1 to edge
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
      double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      for (int k=0;k<3;++k) n[k] /= ln;

      double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
      
      // only keep point if not in outside halfspace
      if (dist > tol)
      {
        outside = true;
        break;
      }
    }
    
    // do only if not already outside w.r.t. to a poly1 (slave) edge
    if (!outside)
    {
      // check against all poly2 (master) edges
      for (int j=0;j<(int)poly2.size();++j)
      {
        // we need diff vector and edge2 first
        double diff[3] = {0.0, 0.0, 0.0};
        double edge[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k)
        {
          diff[k] = collintersec[i].Coord()[k] - poly2[j].Coord()[k];
          edge[k] = (poly2[j].Next())->Coord()[k] - poly2[j].Coord()[k];
        }

        // compute distance from point on poly1 to edge
        double n[3] = {0.0, 0.0, 0.0};
        n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
        n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
        n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
        double ln = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        for (int k=0;k<3;++k) n[k] /= ln;

        double dist = diff[0]*n[0]+diff[1]*n[1]+diff[2]*n[2];
        
        // only keep point if not in outside halfspace
        if (dist > tol)
        {
          outside = true;
          break;
        }
      }
    }
    
    // only keep point if never in outside halfspace
    if (!outside) convexhull.push_back(collintersec[i]);
  }
  
  if (out)
  {
    // check convex hull points
    cout << "\nTesting convex hull points" << endl;
    for (int i=0;i<(int)convexhull.size();++i)
    {
      Vertex& testv = convexhull[i];
      cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
      cout << "Type: " << testv.VType() << endl;
    }
  }
  
  //**********************************************************************
  // STEP6: Collapse convex hull points
  // - this yields a collapsed vector of convex hull vertices
  // - up to now only duplicate intersection points have been eliminated
  // - this operation now removes ALL kinds of duplicate points
  // - intersection points close to poly2/poly1 points are deleted
  // - poly2 points close to poly1 vertices are deleted
  //**********************************************************************
  vector<Vertex> collconvexhull;
  
  for (int i=0;i<(int)convexhull.size();++i)
  {
    // keep track of comparisons
    bool close = false;
    
    // do not collapse poly1 (slave) points
    if (convexhull[i].VType()==MORTAR::Vertex::slave)
    {
      collconvexhull.push_back(convexhull[i]);
      continue;
    }
    
    // check remaining poly2 (master) and intersec points against poly1 (slave) points
    for (int j=0;j<(int)convexhull.size();++j)
    {
      // only collapse with poly1 (slave) points
      if (convexhull[j].VType()!=MORTAR::Vertex::slave) continue;
          
      // distance vector
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k) diff[k] = convexhull[i].Coord()[k] - convexhull[j].Coord()[k];
      double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
      
      // only keep point if not close
      if (dist <= tol)
      {
        close = true;
        break;
      }
    }
    
    // do not check poly2 (master) points
    if (convexhull[i].VType()==MORTAR::Vertex::projmaster)
    {
      if (!close) collconvexhull.push_back(convexhull[i]);
      continue;
    }
        
    // check intersec points against poly2 (master) points
    if (!close && convexhull[i].VType()==MORTAR::Vertex::lineclip)
    {
      for (int j=0;j<(int)convexhull.size();++j)
      {
        // only collapse with poly2 (master) points
        if (convexhull[j].VType()!=MORTAR::Vertex::projmaster) continue;
              
        // distance vector
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k) diff[k] = convexhull[i].Coord()[k] - convexhull[j].Coord()[k];
        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        
        // only keep intersection point if not close
        if (dist <= tol)
        {
          close = true;  
          break;
        }
      }
    }
    
    // keep intersection point only if not collapsed
    if (!close) collconvexhull.push_back(convexhull[i]);
  }
  
  if (out)
  {
    // check collapsed convex hull points
    cout << "\nTesting collapsed convex hull points" << endl;
    for (int i=0;i<(int)collconvexhull.size();++i)
    {
      Vertex& testv = collconvexhull[i];
      cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
      cout << "Type: " << testv.VType() << endl;
    }
  }
  
  //**********************************************************************
  // STEP7: Transform convex hull points to auxiliary plane
  // - x* = A * (x - p1) where p1 = translation, A = rotation
  // - this yields transformed points with coords (x_bar,y_bar,0)
  //**********************************************************************
  // number of points removed from convex hull (neede later)
  int removed = 0;
  
  // only continue if more than two points remaining
  if ((int)collconvexhull.size() < 3)
  {
    // no clip polygon if less than 3 points
    vector<Vertex> empty;
    respoly = empty;
  }
  else if ((int)collconvexhull.size()==3)
  {
    // 3 points already ARE the convex hull
    respoly = collconvexhull;
  }
  else
  {
    // do transformation to auxiliary plane
    double newzero[3] = {collconvexhull[0].Coord()[0],collconvexhull[0].Coord()[1],collconvexhull[0].Coord()[2]};
    double newxaxis[3] = {collconvexhull[1].Coord()[0]-collconvexhull[0].Coord()[0],
                          collconvexhull[1].Coord()[1]-collconvexhull[0].Coord()[1],
                          collconvexhull[1].Coord()[2]-collconvexhull[0].Coord()[2]};
    double newyaxis[3] = {Auxn()[1]*newxaxis[2]-Auxn()[2]*newxaxis[1],
                          Auxn()[2]*newxaxis[0]-Auxn()[0]*newxaxis[2],
                          Auxn()[0]*newxaxis[1]-Auxn()[1]*newxaxis[0]};
    double lnewxaxis = sqrt(newxaxis[0]*newxaxis[0]+newxaxis[1]*newxaxis[1]+newxaxis[2]*newxaxis[2]);
    double lnewyaxis = sqrt(newyaxis[0]*newyaxis[0]+newyaxis[1]*newyaxis[1]+newyaxis[2]*newyaxis[2]);
    
    // normalize
    for (int k=0;k<3;++k)
    {
      newxaxis[k] /= lnewxaxis;
      newyaxis[k] /= lnewyaxis;
    }
    
    // trafo matrix
    LINALG::Matrix<3,3> trafo;
    for (int k=0;k<3;++k)
    {
      trafo(0,k) = newxaxis[k];
      trafo(1,k) = newyaxis[k];
      trafo(2,k) = Auxn()[k]; 
    }
    
    // temporary storage for transformed points
    int np = (int)collconvexhull.size();
    Epetra_SerialDenseMatrix transformed(2,np);
        
    // transform each convex hull point
    for (int i=0;i<np;++i)
    {
      double newpoint[3] = {0.0, 0.0, 0.0};
      
      for (int j=0;j<3;++j)
        for (int k=0;k<3;++k) newpoint[j] += trafo(j,k) * (collconvexhull[i].Coord()[k] - newzero[k]);
 
      if (abs(newpoint[2]) > tol) dserror("ERROR: Transformation to aux. plane failed: z!=0 !");
      transformed(0,i) = newpoint[0];
      transformed(1,i) = newpoint[1];
    }
 
    //**********************************************************************
    // STEP8: Sort convex hull points to obtain final clip polygon
    // - this yields the final clip polygon
    // - sanity of the generated output is checked
    //**********************************************************************
    
    // (1) Find point with smallest x-value
    // (if more than 1 point with identical x-value exists, choose the one with the smallest y-value)
    
    // initialize starting point
    int startindex = 0;
    double startpoint[2] = {transformed(0,0), transformed(1,0)};
    
    for (int i=1;i<np;++i)
    {
      if (transformed(0,i) < startpoint[0])
      {
        startpoint[0] = transformed(0,i);
        startpoint[1] = transformed(1,i);
        startindex = i;
      }
      else if (transformed(0,i) == startpoint[0])
      {
        if (transformed(1,i) < startpoint[1])
        {
          startpoint[1] = transformed(1,i);
          startindex = i;
        }
      }
      else
      {
        // do nothing: starting point did not change
      }
    }
    
    if (out) cout << "Start of convex hull: Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << endl;
     
    // (2) Sort remaining points ascending w.r.t their angle with the y-axis
    // (if more than 1 point with identical angle exists, sort ascending w.r.t. their y-value)
    vector<double> cotangle(0);
    vector<double> yvalues(0);
    vector<int> sorted(0);
    
    for (int i=0;i<np;++i)
    {
      // do nothing for starting point
      if (i==startindex) continue;
      
      // compute angle and store
      double xdiff = transformed(0,i) - startpoint[0];
      double ydiff = transformed(1,i) - startpoint[1];
      
      if (xdiff < 0) dserror("ERROR: Found point with x < x_start for convex hull!");
      if (xdiff < tol) xdiff=tol;
      
      cotangle.push_back(ydiff/xdiff);
      sorted.push_back(i);
    }
    
    if (out)
    {
      cout << "Unsorted convex hull:\n";
      cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << endl;
      for (int i=0;i<np-1;++i)
        cout << "Index " << sorted[i] << "\t" << transformed(0,sorted[i]) << "\t" << transformed(1,sorted[i]) << "\t" << cotangle[sorted[i]] << endl;
    }
    
    // check if sizes are correct
    if ((int)cotangle.size() != np-1) dserror("ERROR: Size went wrong for cot angle!");   
    
    // now sort descending w.r.t cotangle = ascending w.r.t angle
    MORTAR::Sort(&cotangle[0],np-1,&sorted[0]);
    std::reverse(cotangle.begin(), cotangle.end());
    std::reverse(sorted.begin(), sorted.end());
    
    // get associated y-values
    for (int i=0;i<np-1;++i)
      yvalues.push_back(transformed(1,sorted[i]));
    
    // sort list again where two angles are identical
    for (int i=0;i<np-2;++i)
      if (cotangle[i]==cotangle[i+1])
        if (yvalues[i]>yvalues[i+1])
        {
          MORTAR::Swap(cotangle[i],cotangle[i+1]);
          MORTAR::Swap(yvalues[i],yvalues[i+1]);
          MORTAR::Swap(sorted[i],sorted[i+1]);
        }
    
    // check whether 3 points with identical angle exist (without starting point)
    int count = 0;
    
    for (int i=0;i<np-1;++i)
    {
      for (int j=0;j<np-1;++j)
      {
        // check if cotangle is identical
        if (cotangle[j]==cotangle[i]) count += 1;
      }
      
      // stop if 3 identical values found
      if (count >= 3) dserror("ERROR: 3 points on one line in convex hull!");
      else            count = 0;
    }
    
    if (out)
    {
      cout << "Sorted convex hull:\n";
      cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << endl;
      for (int i=0;i<np-1;++i)
        cout << "Index " << sorted[i] << "\t" << transformed(0,sorted[i]) << "\t" << transformed(1,sorted[i]) << "\t" << cotangle[sorted[i]] << endl;
    }
    
    // (3) Go through sorted list of points
    // (keep adding points as long as the last 3 points rotate clockwise)
    // (if 3 points rotate counter-clockwise, do NOT add current point and continue)
    
    // always push pack starting point
    Vertex* current = &collconvexhull[startindex];
    respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
    
    // boolean indicating how many points are removed
    removed = 0;
    
    // go through sorted list and check for clockwise rotation
    for (int i=0;i<np-1;++i)
    {
      double edge1[2] = {0.0, 0.0};
      double edge2[2] = {0.0, 0.0};
      
      // first triple
      if (i==0)
      {
        edge1[0] = transformed(0,sorted[0]) - startpoint[0];
        edge1[1] = transformed(1,sorted[0]) - startpoint[1];
        edge2[0] = transformed(0,sorted[1]) - transformed(0,sorted[0]);
        edge2[1] = transformed(1,sorted[1]) - transformed(1,sorted[0]);
      }
      
      // standard case
      else if (i<np-2)
      {
        edge1[0] = transformed(0,sorted[i]) - transformed(0,sorted[i-1]);
        edge1[1] = transformed(1,sorted[i]) - transformed(1,sorted[i-1]);
        edge2[0] = transformed(0,sorted[i+1]) - transformed(0,sorted[i]);
        edge2[1] = transformed(1,sorted[i+1]) - transformed(1,sorted[i]);
      }
      
      // last triple
      else /* if i = np-1 */
      {
        edge1[0] = transformed(0,sorted[i]) - transformed(0,sorted[i-1]);
        edge1[1] = transformed(1,sorted[i]) - transformed(1,sorted[i-1]);
        edge2[0] = startpoint[0] - transformed(0,sorted[i]);
        edge2[1] = startpoint[1] - transformed(1,sorted[i]);
      }
      
      // check for clockwise rotation
      double cw = edge1[0]*edge2[1]-edge1[1]*edge2[0];
      
      if (out) cout << "Check triple around point " << sorted[i] << endl;
      if (out) cout << "cw: " << cw << endl;
      
      // add point to convex hull if clockwise triple
      // increas counter "removed" if counter-clockwise triple
      if (cw <= 0)
      {
        Vertex* current = &collconvexhull[sorted[i]];
        respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
      }
      else
        removed++;
    }
    
    // (4) Keep in mind that our collconvexhull should already BE the convex hull
    // (thus only sorting should be necessary, but no point should have to be removed in (3))
    if (removed>0)
      cout << "***WARNING*** In total, " << removed << " points removed from convex hull!" << endl;
  }
  
  // **********************************************************************
  // STEP9: Result visualization with GMSH
  // - plot the two input polygons and their vertex numbering
  // - plot the result polygon and its vertex numbering
  // **********************************************************************
  if (out)
  {
    std::ostringstream filename;
    static int gmshcount=0;
    filename << "o/gmsh_output/" << "clipping_";
    if (gmshcount<10)
      filename << 0 << 0 << 0 << 0;
    else if (gmshcount<100)
      filename << 0 << 0 << 0;
    else if (gmshcount<1000)
      filename << 0 << 0;
    else if (gmshcount<10000)
      dserror("Gmsh output implemented for a maximum of 9.999 clip polygons");
    filename << gmshcount << ".pos";
    gmshcount++;
  
    // do output to file in c-style
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream gmshfilecontent;
    gmshfilecontent << "View \" Clipping \" {" << endl;
  
    for (int i=0;i<(int)poly1.size();++i)
    {
      if (i!=(int)poly1.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << ","
                                 << poly1[i].Coord()[2] << "," << poly1[i+1].Coord()[0] << "," << poly1[i+1].Coord()[1] << ","
                                 << poly1[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << ","
                                 << poly1[i].Coord()[2] << "," << poly1[0].Coord()[0] << "," << poly1[0].Coord()[1] << ","
                                 << poly1[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << poly1[i].Coord()[0] << "," << poly1[i].Coord()[1] << "," << poly1[i].Coord()[2] << "," << 17 << ")";
      gmshfilecontent << "{" << "S" << i << "};" << endl;
    }
  
    for (int i=0;i<(int)poly2.size();++i)
    {
      if (i!=(int)poly2.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << ","
                                 << poly2[i].Coord()[2] << "," << poly2[i+1].Coord()[0] << "," << poly2[i+1].Coord()[1] << ","
                                 << poly2[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << ","
                                 << poly2[i].Coord()[2] << "," << poly2[0].Coord()[0] << "," << poly2[0].Coord()[1] << ","
                                 << poly2[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << poly2[i].Coord()[0] << "," << poly2[i].Coord()[1] << "," << poly2[i].Coord()[2] << "," << 17 << ")";
      gmshfilecontent << "{" << "M" << i << "};" << endl;
    }
    
    for (int i=0;i<(int)respoly.size();++i)
    {
      if (i!=(int)respoly.size()-1)
      {
        gmshfilecontent << "SL(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << ","
                                 << respoly[i].Coord()[2] << "," << respoly[i+1].Coord()[0] << "," << respoly[i+1].Coord()[1] << ","
                                 << respoly[i+1].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      }
      else
      {
        gmshfilecontent << "SL(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << ","
                                 << respoly[i].Coord()[2] << "," << respoly[0].Coord()[0] << "," << respoly[0].Coord()[1] << ","
                                 << respoly[0].Coord()[2] << ")";
        gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
  
      }
      gmshfilecontent << "T3(" << scientific << respoly[i].Coord()[0] << "," << respoly[i].Coord()[1] << "," << respoly[i].Coord()[2] << "," << 27 << ")";
      gmshfilecontent << "{" << "R" << i << "};" << endl;
    }
    
    gmshfilecontent << "};" << endl;
  
    // move everything to gmsh post-processing file and close it
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
  }

  return true;  
}

/*----------------------------------------------------------------------*
 |  Check /set projection status of slave nodes (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::HasProjStatus()
{
  // check all nodes
  int nnodes = SlaveIntElement().NumNode();
  DRT::Node** mynodes = SlaveIntElement().Nodes();
  if (!mynodes) dserror("ERROR: HasProjStatus: Null pointer!");

  // loop over all slave nodes
  for (int i=0;i<nnodes;++i)
  {
    MortarNode* mycnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: HasProjStatus: Null pointer!");

    // loop over all vertices of clip polygon
    for (int j=0;j<(int)(Clip().size());++j)
    {
      bool identical = false;

      // check if this clip vertex is slave-type and has the
      // current slave node id
      if ((int)(Clip()[j].VType())==Vertex::slave)
        if (mycnode->Id()==Clip()[j].Nodeids()[0])
          identical = true;

      // set hasproj to true, if so
      if (identical) mycnode->HasProj()=true;
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::Triangulation(map<int,double>& projpar)
{
  // preparations
  int clipsize = (int)(Clip().size());
  vector<vector<map<int,double> > > linvertex(clipsize,vector<map<int,double> >(3));
  vector<map<int,double> > lincenter(3);

  //**********************************************************************
  // (1) Linearization of clip vertex coordinates (only non-empty for contact)
  //**********************************************************************
  VertexLinearization(linvertex,projpar);

  //**********************************************************************
  // (2) Simple cliiping polygon (triangle, quadrilateral) -> IntCells
  //**********************************************************************
  // clip polygon = triangle
  // no triangulation necessary -> 1 IntCell
  if (clipsize==3)
  {
    // IntCell vertices = clip polygon vertices
    Epetra_SerialDenseMatrix coords(3,clipsize);
    for (int i=0;i<clipsize;++i)
      for (int k=0;k<3;++k)
        coords(k,i) = Clip()[i].Coord()[k];

    // create IntCell object and push back
    Cells().push_back(rcp(new IntCell(0,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),linvertex[0],linvertex[1],linvertex[2],GetDerivAuxn())));

    // get out of here
    return true;
  }

  // clip polygon = quadrilateral
  // no triangulation necessary -> 2 IntCells
  else if (clipsize==4)
  {
    // IntCell 1 vertices = clip polygon vertices 0,1,2
    Epetra_SerialDenseMatrix coords(3,3);
    for (int k=0;k<3;++k)
    {
      coords(k,0) = Clip()[0].Coord()[k];
      coords(k,1) = Clip()[1].Coord()[k];
      coords(k,2) = Clip()[2].Coord()[k];
    }

    // create 1st IntCell object and push back
    Cells().push_back(rcp(new IntCell(0,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),linvertex[0],linvertex[1],linvertex[2],GetDerivAuxn())));

    // IntCell vertices = clip polygon vertices 2,3,0
    for (int k=0;k<3;++k)
    {
      coords(k,0) = Clip()[2].Coord()[k];
      coords(k,1) = Clip()[3].Coord()[k];
      coords(k,2) = Clip()[0].Coord()[k];
    }

    // create 2nd IntCell object and push back
    Cells().push_back(rcp(new IntCell(1,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),linvertex[2],linvertex[3],linvertex[0],GetDerivAuxn())));

    // get out of here
    return true;
  }

  //**********************************************************************
  // (3) Find center of clipping polygon (centroid formula)
  //**********************************************************************
  vector<double> clipcenter(3);
  for (int k=0;k<3;++k) clipcenter[k] = 0.0;
  double fac = 0.0;

  // first we need node averaged center
  double nac[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      nac[k] += (Clip()[i].Coord()[k] / clipsize);

  // loop over all triangles of polygon
  for (int i=0; i<clipsize; ++i)
  {
    double xi_i[3] = {0.0, 0.0, 0.0};
    double xi_ip1[3] = {0.0, 0.0, 0.0};

    // standard case
    if (i<clipsize-1)
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[i].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[i+1].Coord()[k];
    }
    // last vertex of clip polygon
    else
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[clipsize-1].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[0].Coord()[k];
    }

    // triangle area
    double diff1[3] = {0.0, 0.0, 0.0};
    double diff2[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k) diff1[k] = xi_ip1[k] - xi_i[k];
    for (int k=0;k<3;++k) diff2[k] = xi_i[k] - nac[k];

    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = diff1[1]*diff2[2] - diff1[2]*diff2[1];
    cross[1] = diff1[2]*diff2[0] - diff1[0]*diff2[2];
    cross[2] = diff1[0]*diff2[1] - diff1[1]*diff2[0];

    double Atri = 0.5 * sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

    // add contributions to clipcenter and fac
    fac += Atri;
    for (int k=0;k<3;++k) clipcenter[k] += 1.0/3.0 * (xi_i[k] + xi_ip1[k] + nac[k]) * Atri;
  }

  //final clipcenter
  for (int k=0;k<3;++k) clipcenter[k] /= fac;
  //cout << "Clipcenter: " << clipcenter[0] << " " << clipcenter[1] << " " << clipcenter[2] << endl;

  //**********************************************************************
  // (4) Linearization of clip center coordinates (only non-empty for contact)
  //**********************************************************************
  CenterLinearization(linvertex,lincenter);

  //**********************************************************************
  // (5) General clipping polygon: Triangulation -> IntCells
  //**********************************************************************
  // center-based triangulation if clip polygon > quadrilateral
  // No. of IntCells is equal to no. of clip polygon vertices
  for (int num=0;num<clipsize;++num)
  {
    // the first vertex is always the clip center
    // the second vertex is always the current clip vertex
    Epetra_SerialDenseMatrix coords(3,3);
    for (int k=0;k<3;++k)
    {
      coords(k,0) = clipcenter[k];
      coords(k,1) = Clip()[num].Coord()[k];
    }

    // the third vertex is the next vertex on clip polygon
    int numplus1 = num+1;
    if (num==clipsize-1)
    {
      for (int k=0;k<3;++k) coords(k,2) = Clip()[0].Coord()[k];
      numplus1 = 0;
    }
    else
      for (int k=0;k<3;++k) coords(k,2) = Clip()[num+1].Coord()[k];

    // create IntCell object and push back
    Cells().push_back(rcp(new IntCell(num,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),lincenter,linvertex[num],linvertex[numplus1],GetDerivAuxn())));
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3d::IntegrateCells()
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrices M (and possibly D) on the current    */
  /* integration cell of the slave / master (integration) element pair  */
  /**********************************************************************/

  // create an integrator instance with correct NumGP and Dim
  // it is sufficient to do this once as all IntCells are triangles
  MORTAR::MortarIntegrator integrator(shapefcn_,Cells()[0]->Shape());

  // loop over all integration cells
  for (int i=0;i<(int)(Cells().size());++i)
  {
    // compare intcell area with slave integration element area
    double intcellarea = Cells()[i]->Area();
    double selearea = 0.0;
    if (!CouplingInAuxPlane())
      selearea = SlaveIntElement().MoData().Area();
    else
    {
      DRT::Element::DiscretizationType dt = SlaveIntElement().Shape();
      if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        selearea = 4.0;
      else if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        selearea = 0.5;
      else dserror("ERROR: IntegrateCells: Invalid 3D slave element type");
    }

    // integrate cell only if not neglectable
    if (intcellarea < MORTARINTLIM*selearea) continue;

    // debug output of integration cells in GMSH
#ifdef MORTARGMSHCELLS
    GmshOutputCells(i);
#endif // #ifdef MORTARGMSHCELLS
    
    // set segmentation status of all slave nodes
    // (hassegment_ of a slave node is true if ANY segment/cell
    // is integrated that contributes to this slave node)
    int nnodes = SlaveIntElement().NumNode();
    DRT::Node** mynodes = SlaveIntElement().Nodes();
    if (!mynodes) dserror("ERROR: Null pointer!");
    for (int k=0;k<nnodes;++k)
    {
      MORTAR::MortarNode* mycnode = static_cast<MORTAR::MortarNode*> (mynodes[k]);
      if (!mycnode) dserror("ERROR: Null pointer!");
      mycnode->HasSegment()=true;
    }

    // *******************************************************************
    // different options for mortar integration
    // *******************************************************************
    // (1) no quadratic element(s) involved -> linear LM interpolation
    // (2) quadratic element(s) involved -> quadratic LM interpolation
    // (3) quadratic element(s) involved -> linear LM interpolation
    // (4) quadratic element(s) involved -> piecew. linear LM interpolation
    // *******************************************************************
    INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();
    
    // *******************************************************************
    // case (1)
    // *******************************************************************
    if (!Quad())
    {
      // prepare integration of M (and possibly D) on intcells
      int nrow = SlaveElement().NumNode();
      int ncol = MasterElement().NumNode();
      int ndof = static_cast<MORTAR::MortarNode*>(SlaveElement().Nodes()[0])->NumDof();
      
      // mortar matrix dimensions depend on the actual number of slave / master DOFs
      // per node, which is NOT always identical to the problem dimension
      // (e.g. fluid meshtying -> 4 DOFs per node in a 3D problem)
      // -> thus the following check has been commented out (popp 04/2011) 
      // if (ndof != Dim()) dserror("ERROR: Problem dimension and dofs per node not identical");
      
      RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,nrow*ndof));
      RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
      RCP<Epetra_SerialDenseVector> gseg = Teuchos::null;
      
      // check whether auxiliary plane coupling or not and call integrator
      if (CouplingInAuxPlane())
        integrator.IntegrateDerivCell3DAuxPlane(SlaveElement(),MasterElement(),Cells()[i],Auxn(),dseg,mseg,gseg);
      else /*(!CouplingInAuxPlane()*/
        integrator.IntegrateDerivCell3D(SlaveElement(),MasterElement(),Cells()[i],dseg,mseg,gseg);
      
      // assembly of intcell contributions to D and M
      integrator.AssembleD(Comm(),SlaveElement(),*dseg);
      integrator.AssembleM(Comm(),SlaveElement(),MasterElement(),*mseg);
    }
    
    // *******************************************************************
    // cases (2) and (3)
    // *******************************************************************
    else if (Quad() && (lmtype==INPAR::MORTAR::lagmult_quad_quad || lmtype==INPAR::MORTAR::lagmult_lin_lin))
    {
      // prepare integration of M (and possibly D) on intcells
      int nrow = SlaveElement().NumNode();
      int ncol = MasterElement().NumNode();
      int ndof = static_cast<MORTAR::MortarNode*>(SlaveElement().Nodes()[0])->NumDof();
      if (ndof != Dim()) dserror("ERROR: Problem dimension and dofs per node not identical");
      RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,nrow*ndof));
      RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
      RCP<Epetra_SerialDenseVector> gseg = Teuchos::null;
      
      // static_cast to make sure to pass in IntElement&
      MORTAR::IntElement& sintref = static_cast<MORTAR::IntElement&>(SlaveIntElement());
      MORTAR::IntElement& mintref = static_cast<MORTAR::IntElement&>(MasterIntElement());
      
      // check whether auxiliary plane coupling or not and call integrator
      if (CouplingInAuxPlane())
        integrator.IntegrateDerivCell3DAuxPlaneQuad(SlaveElement(),MasterElement(),sintref,mintref,
            Cells()[i],Auxn(),lmtype,dseg,mseg,gseg);
      else /*(!CouplingInAuxPlane()*/
        integrator.IntegrateDerivCell3DQuad(SlaveElement(),MasterElement(),sintref,mintref,
            Cells()[i],lmtype,dseg,mseg,gseg);
      
      // assembly of intcell contributions to D and M
      integrator.AssembleD(Comm(),SlaveElement(),*dseg);
      integrator.AssembleM(Comm(),SlaveElement(),MasterElement(),*mseg);
    }
    
    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (Quad() && lmtype==INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      // check for dual shape functions
      if (shapefcn_ == INPAR::MORTAR::shape_dual)
        dserror("ERROR: Piecewise linear LM interpolation not yet implemented for DUAL 3D quadratic mortar");

      // prepare integration of M (and possibly D) on intcells
      int nrow = SlaveElement().NumNode();
      int ncol = MasterElement().NumNode();
      int nintrow = SlaveIntElement().NumNode();
      int ndof = static_cast<MORTAR::MortarNode*>(SlaveElement().Nodes()[0])->NumDof();
      if (ndof != Dim()) dserror("ERROR: Problem dimension and dofs per node not identical");
      RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nintrow*ndof,nrow*ndof));
      RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nintrow*ndof,ncol*ndof));
      RCP<Epetra_SerialDenseVector> gseg = Teuchos::null;
      
      // static_cast to make sure to pass in IntElement&
      MORTAR::IntElement& sintref = static_cast<MORTAR::IntElement&>(SlaveIntElement());
      MORTAR::IntElement& mintref = static_cast<MORTAR::IntElement&>(MasterIntElement());
      
      // check whether auxiliary plane coupling or not and call integrator
      if (CouplingInAuxPlane())
        integrator.IntegrateDerivCell3DAuxPlaneQuad(SlaveElement(),MasterElement(),sintref,mintref,
            Cells()[i],Auxn(),lmtype,dseg,mseg,gseg);
      else /*(!CouplingInAuxPlane()*/
        integrator.IntegrateDerivCell3DQuad(SlaveElement(),MasterElement(),sintref,mintref,
            Cells()[i],lmtype,dseg,mseg,gseg);
      
      // assembly of intcell contributions to D and M
      // (NOTE THAT THESE ARE SPECIAL VERSIONS HERE FOR PIECEWISE LINEAR INTERPOLATION)
      integrator.AssembleD(Comm(),SlaveElement(),sintref,*dseg);
      integrator.AssembleM(Comm(),sintref,MasterElement(),*mseg);
    }
    
    // *******************************************************************
    // other cases
    // *******************************************************************
    else
    {
      dserror("ERROR: IntegrateCells: Invalid case for 3D mortar coupling LM interpolation");
    }
    // *******************************************************************
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 07/11|
 *----------------------------------------------------------------------*/
void MORTAR::Coupling3d::GmshOutputCells(int lid)
{
  // every processor writes its own cell file
  int proc = idiscret_.Comm().MyPID();
  int nproc = idiscret_.Comm().NumProc();
  
  // write each integration cell only once
  // (no overlap, only owner of slave element writes output)
  if (proc != SlaveElement().Owner()) return;
    
  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << "o/gmsh_output/" << filebase << "_cells_" << proc << ".pos";
  
  // do output to file in c-style
  FILE* fp = NULL;
  
  // static variable
  static int count = 0;
  
  // open file
  if (count==0) fp = fopen(filename.str().c_str(), "w");
  else          fp = fopen(filename.str().c_str(), "a");
  
  // plot current integration cell
  const Epetra_SerialDenseMatrix& coord = Cells()[lid]->Coords();
  
  // write output to temporary stringstream
  std::stringstream gmshfilecontent;
  
  // header and dummy elements
  if (count==0)
  {
    // header
    gmshfilecontent << "View \"Integration Cells Proc " << proc << "\" {" << endl;
    
    // dummy element 1
    gmshfilecontent << "ST(" << scientific << 0.0 << "," << 0.0 << ","
                    << 0.0 << "," << 0.0 << "," << 0.0 << ","
                    << 0.0 << "," << 0.0 << "," << 0.0 << ","
                    << 0.0 << ")";
    gmshfilecontent << "{" << scientific << 0 << "," << 0 << "," << 0 << "};" << endl;
    
    // dummy element 1
    gmshfilecontent << "ST(" << scientific << 0.0 << "," << 0.0 << ","
                    << 0.0 << "," << 0.0 << "," << 0.0 << ","
                    << 0.0 << "," << 0.0 << "," << 0.0 << ","
                    << 0.0 << ")";
    gmshfilecontent << "{" << scientific << nproc << "," << nproc << "," << nproc << "};" << endl;
  }
  
  // plot cell itself
  gmshfilecontent << "ST(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                  << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                  << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                  << coord(2,2) << ")";
  gmshfilecontent << "{" << scientific << proc << "," << proc << "," << proc << "};" << endl;
  
  // move everything to gmsh post-processing files and close them
  fprintf(fp,gmshfilecontent.str().c_str());
  fclose(fp);
  
  // increase static variable
  count += 1;
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3dQuad::Coupling3dQuad(DRT::Discretization& idiscret,
                                int dim, bool quad, bool auxplane,
                                MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                                MORTAR::IntElement& sintele,
                                MORTAR::IntElement& mintele,
                                INPAR::MORTAR::LagMultQuad& lmtype) :
MORTAR::Coupling3d(INPAR::MORTAR::shape_undefined,idiscret,dim,quad,auxplane,sele,mele),
sintele_(sintele),
mintele_(mintele),
lmtype_(lmtype)
{
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Quad())
    dserror("ERROR: Coupling3dQuad called for non-quadratic ansatz!");

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3dQuad::Coupling3dQuad(const INPAR::MORTAR::ShapeFcn shapefcn,
                                DRT::Discretization& idiscret,
                                int dim, bool quad, bool auxplane,
                                MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                                MORTAR::IntElement& sintele,
                                MORTAR::IntElement& mintele,
                                INPAR::MORTAR::LagMultQuad& lmtype) :
MORTAR::Coupling3d(shapefcn,idiscret,dim,quad,auxplane,sele,mele),
sintele_(sintele),
mintele_(mintele),
lmtype_(lmtype)
{
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Quad())
    dserror("ERROR: Coupling3dQuad called for non-quadratic ansatz!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3dManager::Coupling3dManager(DRT::Discretization& idiscret, int dim, bool quad,
                                             bool auxplane, MORTAR::MortarElement* sele,
                                             vector<MORTAR::MortarElement*> mele) :
shapefcn_(INPAR::MORTAR::shape_undefined),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // evaluate coupling
  EvaluateCoupling();

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
MORTAR::Coupling3dManager::Coupling3dManager(const INPAR::MORTAR::ShapeFcn shapefcn,
                                             DRT::Discretization& idiscret, int dim, bool quad,
                                             bool auxplane, MORTAR::MortarElement* sele,
                                             vector<MORTAR::MortarElement*> mele) :
shapefcn_(shapefcn),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // evaluate coupling
  EvaluateCoupling();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                   popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling3dManager::EvaluateCoupling()
{

  // loop over all master elements associated with this slave element
  for (int m=0;m<(int)MasterElements().size();++m)
  {
    // create Coupling3d object and push back
    Coupling().push_back(rcp(new Coupling3d(shapefcn_,idiscret_,dim_,false,auxplane_,
                                            SlaveElement(),MasterElement(m))));

    // do coupling
    Coupling()[m]->EvaluateCoupling();
  }

  return true;
}

#endif //#ifdef CCADISCRET
