/*!----------------------------------------------------------------------
\file drt_contact_coupling3d.cpp
\brief A class for mortar coupling of ONE slave element and ONE master
       element of a contact interface in 3D.

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

#include "drt_contact_coupling3d.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/09|
 *----------------------------------------------------------------------*/
CONTACT::IntElement::IntElement(int lid, int id, ElementType etype, int owner, 
                                const DRT::Element::DiscretizationType& shape, 
                                const int numnode,
                                const int* nodeids,
                                vector<DRT::Node*> nodes,
                                const bool isslave) :
CONTACT::CElement(id,etype,owner,shape,numnode,nodeids,isslave),
lid_(lid)
{
  // check for consistency of nodeids and nodes
  for (int i=0;i<NumNode();++i)
    if (nodes[i]->Id()!=nodeids[i])
      dserror("ERROR: IntElement: Inconsistency Nodes and NodeIds!");
  
  // store given nodes into class variable
  quadnode_.resize((int)nodes.size());
  for (int i=0;i<(int)nodes.size();++i)
    quadnode_[i] = nodes[i];
  
  // as discretization is already evaluated, compute area
  Area()=ComputeArea();
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Intcell::Intcell(int id, int nvertices, Epetra_SerialDenseMatrix& coords,
                          double* auxn, const DRT::Element::DiscretizationType& shape,
                          bool auxplane,
                          vector<map<int,double> >& linv1,
                          vector<map<int,double> >& linv2,
                          vector<map<int,double> >& linv3,
                          vector<map<int,double> >& linauxn) :
id_(id),
nvertices_(nvertices),
coords_(coords),
shape_(shape),
auxplane_(auxplane)
{
   // check nvertices_ and shape_
  if (nvertices_!=3) dserror("ERROR: Integration cell must have 3 vertices");
  if (shape_!=DRT::Element::tri3) dserror("ERROR: Integration cell must be tri3");
  
  // check dimensions of coords_
  if (coords_.M() != 3) dserror("ERROR: Inconsistent coord matrix");
  if (coords_.N() != nvertices_) dserror("ERROR: Inconsistent coord matrix");
  
  // store auxiliary plane normal
  for (int k=0;k<3;++k) Auxn()[k] = auxn[k];
  
  // compute area of Intcell
  double t1[3] = {0.0, 0.0, 0.0};
  double t2[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k)
  {
    t1[k]=Coords()(k,1)-Coords()(k,0);
    t2[k]=Coords()(k,2)-Coords()(k,0);
  }
  
  double t1xt2[3] = {0.0, 0.0, 0.0};
  t1xt2[0] = t1[1]*t2[2]-t1[2]*t2[1];
  t1xt2[1] = t1[2]*t2[0]-t1[0]*t2[2];
  t1xt2[2] = t1[0]*t2[1]-t1[1]*t2[0]; 
  area_ = 0.5*sqrt(t1xt2[0]*t1xt2[0]+t1xt2[1]*t1xt2[1]+t1xt2[2]*t1xt2[2]);
  
  // store vertex linearizations and auxn linearization
  linvertex_.resize(3);
  linvertex_[0] = linv1;
  linvertex_[1] = linv2;
  linvertex_[2] = linv3;
  linauxn_ = linauxn;
  
  return;
}

/*----------------------------------------------------------------------*
 |  cctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Intcell::Intcell(const Intcell& old) :
id_(old.id_),
nvertices_(old.nvertices_),
area_(old.area_),
coords_(old.coords_),
shape_(old.shape_)
{
  // empty copy constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (Intcell)        popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Intcell::LocalToGlobal(const double* xi,
                                     double* globcoord,
                                     int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");
  
  // collect fundamental data
  int nnodes = NumVertices();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  
  // Evaluate shape, get nodal coords and interpolate global coords
  EvaluateShape(xi, val, deriv);
  for (int i=0;i<3;++i) globcoord[i]=0.0;
  
  for (int i=0;i<nnodes;++i)
  {
    if (inttype==0)
    {
      // use shape function values for interpolation
      globcoord[0]+=val[i]*Coords()(0,i);
      globcoord[1]+=val[i]*Coords()(1,i);
      globcoord[2]+=val[i]*Coords()(2,i);
    }
    else if (inttype==1)
    {
      // use shape function derivatives xi for interpolation
      globcoord[0]+=deriv(i,0)*Coords()(0,i);
      globcoord[1]+=deriv(i,0)*Coords()(1,i);
      globcoord[2]+=deriv(i,0)*Coords()(2,i);
    }
    else if (inttype==2)
    {
      // use shape function derivatives eta for interpolation
      globcoord[0]+=deriv(i,1)*Coords()(0,i);
      globcoord[1]+=deriv(i,1)*Coords()(1,i);
      globcoord[2]+=deriv(i,1)*Coords()(2,i);
    }
    else
      dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions (Intcell)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Intcell::EvaluateShape(const double* xi,
    LINALG::SerialDenseVector& val, LINALG::SerialDenseMatrix& deriv)
{
  if (!xi)
    dserror("ERROR: EvaluateShape (Intcell) called with xi=NULL");
  
  // 3noded triangular element
  if(Shape()==DRT::Element::tri3)
  {
    val[0] = 1-xi[0]-xi[1]; 
    val[1] = xi[0];
    val[2] = xi[1];
    deriv(0,0) = -1.0; deriv(0,1) = -1.0;
    deriv(1,0) =  1.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  1.0;
  }
  
  // unknown case
  else dserror("ERROR: EvaluateShape (Intcell) called for type != tri3");
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant (Intcell)                   popp 11/08|
 *----------------------------------------------------------------------*/
double CONTACT::Intcell::Jacobian(double* xi)
{
  double jac = 0.0;
  vector<double> gxi(3);
  vector<double> geta(3);
    
  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
    jac = Area()*2;
  
  // unknown case
  else dserror("ERROR: Jacobian (Intcell) called for unknown ele type!");
  
  return jac;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det.              popp 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::Intcell::DerivJacobian(double* xi, vector<double>& derivjac)
{
  // initialize parameters
  int nnodes = NumVertices();
  vector<double> gxi(3);
  vector<double> geta(3);
  
  // evaluate shape functions
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  EvaluateShape(xi, val, deriv);
  
  // metrics routine gives local basis vectors
  for (int k=0;k<3;++k)
  {
    gxi[k]=Coords()(k,1)-Coords()(k,0);
    geta[k]=Coords()(k,2)-Coords()(k,0);
  }
  
  // cross product of gxi and geta
  double cross[3] = {0.0, 0.0, 0.0};
  cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
  cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
  cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
  
  double jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    
  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
  {
    // *********************************************************************
    // compute Jacobian derivative
    // *********************************************************************
    if (CouplingInAuxPlane())
      dserror("ERROR: DerivJacobian (SlaveParamSpace) called for AuxPlane case!");
    else //(!CouplingInAuxPlane())
    {
      // in this case, intcells live in slave element parameter space
      // cross[0] and cross[1] are zero!
      for (int i=0;i<nnodes;++i)
      {    
        derivjac[2*i]   += 1/jac * cross[2] * geta[1] * deriv(i,0);
        derivjac[2*i]   -= 1/jac * cross[2] * gxi[1]  * deriv(i,1);
        derivjac[2*i+1] -= 1/jac * cross[2] * geta[0] * deriv(i,0);
        derivjac[2*i+1] += 1/jac * cross[2] * gxi[0]  * deriv(i,1);
      }
    }
  }
  
  // unknown case
  else dserror("ERROR: DerivJacobian (Intcell) called for unknown ele type!");
  
  /*
  // finite difference check
  typedef map<int,double>::const_iterator CI;
  cout << "Analytical Intcell jac derivative:" << endl;
  for (CI p = derivjac.begin(); p != derivjac.end(); ++p)
  {
    cout << "dof: " << p->first << " " << p->second << endl;
  }

  double delta = 1.0e-8;
  double jacfd = 0.0;
  cout << "FD Intcell jac derivative:" << endl;
  
  for (int i=0;i<nnodes;++i)
  {
    for (int j=0;j<2;++j)
    {
      Coords()(j,i) += delta;
      
      // metrics routine gives local basis vectors
      for (int k=0;k<3;++k)
      {
        gxi[k]=Coords()(k,1)-Coords()(k,0);
        geta[k]=Coords()(k,2)-Coords()(k,0);
      }
      
      // cross product of gxi and geta
      cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
      cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
      cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
      
      jacfd = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
      cout << "dof: " << 2*i+j << " " << (jacfd-jac)/delta << endl;
      Coords()(j,i) -= delta;
    }
  }
  cout << endl;
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det. AuxPlane     popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::Intcell::DerivJacobian(double* xi, map<int,double>& derivjac)
{
  // metrics routine gives local basis vectors
  vector<double> gxi(3);
  vector<double> geta(3);
  
  for (int k=0;k<3;++k)
  {
    gxi[k]=Coords()(k,1)-Coords()(k,0);
    geta[k]=Coords()(k,2)-Coords()(k,0);
  }
  
  // cross product of gxi and geta
  double cross[3] = {0.0, 0.0, 0.0};
  cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
  cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
  cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
  
  double jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  typedef map<int,double>::const_iterator CI;
  
  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
  {
    // *********************************************************************
    // compute Jacobian derivative
    // *********************************************************************
    if (CouplingInAuxPlane())
    {
      // first vertex (Coords(k,0)) is part of gxi and geta
      for (CI p=GetDerivVertex(0)[0].begin();p!=GetDerivVertex(0)[0].end();++p)
      {
        derivjac[p->first] -= 1/jac*cross[1]*gxi[2]*(p->second);
        derivjac[p->first] += 1/jac*cross[1]*geta[2]*(p->second);
        derivjac[p->first] += 1/jac*cross[2]*gxi[1]*(p->second);
        derivjac[p->first] -= 1/jac*cross[2]*geta[1]*(p->second);
      }
      for (CI p=GetDerivVertex(0)[1].begin();p!=GetDerivVertex(0)[1].end();++p)
      {
        derivjac[p->first] += 1/jac*cross[0]*gxi[2]*(p->second);
        derivjac[p->first] -= 1/jac*cross[0]*geta[2]*(p->second);
        derivjac[p->first] -= 1/jac*cross[2]*gxi[0]*(p->second);
        derivjac[p->first] += 1/jac*cross[2]*geta[0]*(p->second);
      }
      for (CI p=GetDerivVertex(0)[2].begin();p!=GetDerivVertex(0)[2].end();++p)
      {
        derivjac[p->first] -= 1/jac*cross[0]*gxi[1]*(p->second);
        derivjac[p->first] += 1/jac*cross[0]*geta[1]*(p->second);
        derivjac[p->first] += 1/jac*cross[1]*gxi[0]*(p->second);
        derivjac[p->first] -= 1/jac*cross[1]*geta[0]*(p->second);
      }
      
      // second vertex (Coords(k,1)) is part of gxi
      for (CI p=GetDerivVertex(1)[0].begin();p!=GetDerivVertex(1)[0].end();++p)
      {
        derivjac[p->first] -= 1/jac*cross[1]*geta[2]*(p->second);
        derivjac[p->first] += 1/jac*cross[2]*geta[1]*(p->second);
      }
      for (CI p=GetDerivVertex(1)[1].begin();p!=GetDerivVertex(1)[1].end();++p)
      {
        derivjac[p->first] += 1/jac*cross[0]*geta[2]*(p->second);
        derivjac[p->first] -= 1/jac*cross[2]*geta[0]*(p->second);
      }
      for (CI p=GetDerivVertex(1)[2].begin();p!=GetDerivVertex(1)[2].end();++p)
      {
        derivjac[p->first] -= 1/jac*cross[0]*geta[1]*(p->second);
        derivjac[p->first] += 1/jac*cross[1]*geta[0]*(p->second);
      }
      
      // third vertex (Coords(k,2)) is part of geta
      for (CI p=GetDerivVertex(2)[0].begin();p!=GetDerivVertex(2)[0].end();++p)
      {
        derivjac[p->first] += 1/jac*cross[1]*gxi[2]*(p->second);
        derivjac[p->first] -= 1/jac*cross[2]*gxi[1]*(p->second);
      }
      for (CI p=GetDerivVertex(2)[1].begin();p!=GetDerivVertex(2)[1].end();++p)
      {
        derivjac[p->first] -= 1/jac*cross[0]*gxi[2]*(p->second);
        derivjac[p->first] += 1/jac*cross[2]*gxi[0]*(p->second);
      }
      for (CI p=GetDerivVertex(2)[2].begin();p!=GetDerivVertex(2)[2].end();++p)
      {
        derivjac[p->first] += 1/jac*cross[0]*gxi[1]*(p->second);
        derivjac[p->first] -= 1/jac*cross[1]*gxi[0]*(p->second);
      }
    }
    else //(!CouplingInAuxPlane())
      dserror("ERROR: DerivJacobian (AuxPlane) called for SlaveParamSpace case!"); 
  }
  
  // unknown case
  else dserror("ERROR: DerivJacobian (Intcell) called for unknown ele type!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Vertex::Vertex(vector<double> coord, Vertex::vType type, vector<int> nodeids,
                        Vertex* next, Vertex* prev, bool intersect,
                        bool entryexit, Vertex* neighbor, double alpha) :
coord_(coord),
type_(type),
nodeids_(nodeids),
next_(next),
prev_(prev),
intersect_(intersect),
entryexit_(entryexit),
neighbor_(neighbor),
alpha_(alpha)
{
   // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  cctor (public)                                            popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Vertex::Vertex(const Vertex& old) :
coord_(old.coord_),
type_(old.type_),
nodeids_(old.nodeids_),
next_(old.next_),
prev_(old.prev_),
intersect_(old.intersect_),
entryexit_(old.entryexit_),
neighbor_(old.neighbor_),
alpha_(old.alpha_)
{
   // empty copy constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3d::Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                  bool auxplane, CONTACT::CElement& sele, CONTACT::CElement& mele) :
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
 |  Evaluate coupling (3D)                                    popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::EvaluateCoupling()
{
  // rough check whether elements are "near"
  // whether or not quadratic 3d coupling is performed, we only
  // check the distance of the parent slave and master elements
  bool near = RoughCheck();
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
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave();
    
    // project master element nodes onto auxiliary plane
    ProjectMaster();
    
    // tolerance for polygon clipping
    double sminedge = SlaveIntElement().MinEdgeSize();
    double mminedge = MasterIntElement().MinEdgeSize(); 
    tol = CONTACTCLIPTOL * min(sminedge,mminedge);
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
      CNode* mnode = static_cast<CNode*>(node);
      
      // do the projection
      // the third component of sxi will be the proj. parameter alpha!
      double sxi[2] = {0.0, 0.0};
      double alpha = 0.0;
      CONTACT::Projector projector(3);
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
    tol = CONTACTCLIPTOL;
  }
  // *******************************************************************
  
  // do polygon clipping
  PolygonClipping(SlaveVertices(),MasterVertices(),Clip(),tol);
  int clipsize = (int)(Clip().size());
  
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
 |  Rough check if elements are near (3D)                     popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::RoughCheck()
{
  double sme = SlaveElement().MaxEdgeSize();
  double mme = MasterElement().MaxEdgeSize(); 
  double near = 2.0 * max(sme,mme);
  
  double loccs[2] = {0.0, 0.0};  
  DRT::Element::DiscretizationType dts = SlaveElement().Shape();
  if (dts==CElement::tri3 || dts==CElement::tri6)
  {
    loccs[0] = 1.0/3;
    loccs[1] = 1.0/3;
  }
  double loccm[2] = {0.0, 0.0};  
  DRT::Element::DiscretizationType dtm = MasterElement().Shape();
  if (dtm==CElement::tri3 || dtm==CElement::tri6)
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
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::AuxiliaryPlane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2];
    
  DRT::Element::DiscretizationType dt = SlaveIntElement().Shape();
  if (dt==CElement::tri3 || dt==CElement::tri6)
  {
    loccenter[0] = 1.0/3;
    loccenter[1] = 1.0/3;
  }
  else if (dt==CElement::quad4 || dt==CElement::quad8 || dt==CElement::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else dserror("ERROR: AuxiliaryPlane called for unknown element type");
  
  // compute element center via shape fct. interpolation
  SlaveIntElement().LocalToGlobal(loccenter,Auxc(),0);
  
  // we then compute the unit normal vector at the element center
  Lauxn() = SlaveIntElement().ComputeUnitNormalAtXi(loccenter,Auxn());
  
  // also compute linearization of the unit normal vector
  SlaveIntElement().DerivUnitNormalAtXi(loccenter,GetDerivAuxn());
  
  //cout << "Slave Element: " << SlaveIntElement().Id() << endl;
  //cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << endl;
  //cout << "->Normal: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] << endl;
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Project slave element onto auxiliary plane (public)       popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::ProjectSlave()
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
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
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
bool CONTACT::Coupling3d::ProjectMaster()
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
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
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
void CONTACT::Coupling3d::PolygonClipping(vector<Vertex>& poly1,
                                          vector<Vertex>& poly2,
                                          vector<Vertex>& respoly,
                                          double& tol)
{
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
  
  //cout << "Center 1: " << center1[0] << " " << center1[1] << " " << center1[2] << endl;
  //cout << "Center 2: " << center2[0] << " " << center2[1] << " " << center2[2] << endl;
  
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
    //cout << "Polygon 2 (master) not ordered counter-clockwise -> reordered!" << endl;
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
  //cout << "\nInput Poylgon 1:";
  //for (int i=0;i<(int)poly1.size();++i)
  //  cout << "\nVertex " << i << ":\t" << scientific << poly1[i].Coord()[0]
  //       << "\t" << poly1[i].Coord()[1] << "\t" << poly1[i].Coord()[2];
 
  //cout << "\nInput Poylgon 2:";
  //for (int i=0;i<(int)poly2.size();++i)
  //  cout << "\nVertex " << i << ":\t" << scientific << poly2[i].Coord()[0]
  //       << "\t" << poly2[i].Coord()[1] << "\t" << poly2[i].Coord()[2];
  
  //cout << endl << endl;
  
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
        //cout << "WARNING: Detected two parallel edges! (" << i << "," << j << ")" << endl;
        continue;
      }
      
      //cout << "Searching intersection (" << i << "," << j << ")" << endl;
      
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
          
          //cout << "Found intersection! (" << i << "," << j << ") " << alphap << " " << alphaq << endl;
          //cout << "On Polygon 1: " << ip[0] << " " << ip[1] << " " << ip[2] << endl;
          //cout << "On Polygon 2: " << iq[0] << " " << iq[1] << " " << iq[2] << endl;
          
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
      //cout << "Polygon S lies fully inside polygon M!" << endl; 
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
        //cout << "Polygon M lies fully inside polygon S!" << endl; 
        respoly = poly2;
      }
      
      // fully adjacent case
      else
      {
        //cout << "Polygons S and M are fully adjacent!" << endl; 
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
    //cout << endl;
    for (int i=0;i<(int)intersec1.size();++i)
    {
      //cout << "Intersec1: " << i << " " << intersec1[i].Coord()[0] << " " << intersec1[i].Coord()[1] << " " << intersec1[i].Coord()[2];
      //cout << " EntryExit: " << intersec1[i].EntryExit() << endl;
    }
    
    //cout << endl;
    for (int i=0;i<(int)intersec2.size();++i)
    {
     // cout << "Intersec2: " << i << " " <<  intersec2[i].Coord()[0] << " " << intersec2[i].Coord()[1] << " " << intersec2[i].Coord()[2];
      //cout << " EntryExit: " << intersec2[i].EntryExit() << endl;
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
    //cout << "\nStart loop on Slave at " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
    respoly.push_back(Vertex(current->Coord(),Vertex::lineclip,current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
    
    do {
      // find next Vertex / Vertices (path)
      if (current->EntryExit()==true)
      {
        //cout << "Intersection was Entry, so move to Next() on same polygon!" << endl;
        do {
          current = current->Next();
          //cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
        } while (current->Intersect()==false);
        //cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }
      else
      {
        //cout << "Intersection was Exit, so move to Prev() on same polygon!" << endl;
        do {
          current = current->Prev();
          //cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
        } while (current->Intersect()==false);
        //cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }
      
      // jump to the other input polygon
      current = current->Neighbor();
      //cout << "Jumping to other polygon at intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      //cout << "Length of result list so far: " << (int)respoly.size() << endl;
      
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
        else {}
         // cout << "Collapsed two points in result polygon!" << endl;
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
        else {}
          //cout << "Collapsed two points in result polygon!" << endl;
      }
    }
    
    // replace respoly by collapsed respoly
    respoly = collapsedrespoly;
    //cout << "Final length of result list: " << (int)respoly.size() << endl;
        
    // check if respoly collapsed to nothing
    if ((int)collapsedrespoly.size()<3)
    {
     // cout << "Collapsing of result polygon led to < 3 vertices -> no respoly!" << endl;
      vector<Vertex> empty;
      respoly = empty;
    }
 
    // check for rotation of result polygon (must be clockwise!!!)
    // first get geometric center
    double center[3] = {0.0, 0.0, 0.0};
    
    for (int i=0;i<(int)respoly.size();++i)
      for (int k=0;k<3;++k)
        center[k] += respoly[i].Coord()[k]/((int)respoly.size());
    
    //cout << "\nCenter ResPoly: " << center[0] << " " << center[1] << " " << center[2] << endl;
    
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
    
    // print final input polygons to screen
      //cout << "\nResult Poylgon:";
      //for (int i=0;i<(int)respoly.size();++i)
      //  cout << "\nVertex " << i << ":\t" << respoly[i][0] << "\t" << respoly[i][1] << "\t" << respoly[i][2];
      
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
  /*
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
  */
  
  /*
  // check result polygon points
  cout << "\nTesting result polygon points" << endl;
  for (int i=0;i<(int)respoly.size();++i)
  {
    Vertex& testv = respoly[i];
    cout << "Coords: " << testv.Coord()[0] << " " << testv.Coord()[1] << " " << testv.Coord()[2] << endl;
    cout << "Type: " << testv.VType() << endl;
    cout << "Alpha: " << testv.Alpha() << endl;
    if (testv.VType()==Vertex::lineclip)
      cout << "Lineclip ids: " << testv.Nodeids()[0] << " " << testv.Nodeids()[1]
           << " " << testv.Nodeids()[2] << " " << testv.Nodeids()[3] << endl << endl;
    else
      cout << "Node id: " << testv.Nodeids()[0] << endl << endl;
  }
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Check /set projection status of slave nodes (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::HasProjStatus()
{
  // check all nodes
  int nnodes = SlaveIntElement().NumNode();
  DRT::Node** mynodes = SlaveIntElement().Nodes();
  if (!mynodes) dserror("ERROR: HasProjStatus: Null pointer!");
    
  // loop over all slave nodes
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
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
bool CONTACT::Coupling3d::Triangulation(map<int,double>& projpar)
{
  // preparations
  int clipsize = (int)(Clip().size());
  vector<vector<map<int,double> > > linvertex(clipsize,vector<map<int,double> >(3));
  vector<map<int,double> > lincenter(3); 
  
  //**********************************************************************
  // (1) Linearization of clip vertex coordinates
  //**********************************************************************
  VertexLinearization(linvertex,projpar);
  
  //**********************************************************************
  // (2) Find center of clipping polygon (centroid formula)
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
  // (3) Linearization of clip center coordinates
  //**********************************************************************
  CenterLinearization(linvertex,lincenter);
  
  //**********************************************************************
  // (4)Triangulation -> Intcells
  //**********************************************************************
  vector<Intcell> cells;
  
  // easy if clip polygon = triangle: 1 Intcell
  if (clipsize==3)
  {
    // Intcell vertices = clip polygon vertices
    Epetra_SerialDenseMatrix coords(3,clipsize);
    for (int i=0;i<clipsize;++i)
      for (int k=0;k<3;++k)
        coords(k,i) = Clip()[i].Coord()[k];
    
    // create Intcell object and push back
    Cells().push_back(rcp(new Intcell(0,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),linvertex[0],linvertex[1],linvertex[2],GetDerivAuxn())));
    
  }
  
  // triangulation if clip polygon > triangle
  else
  {
    // No. of Intcells is equal to no. of clip polygon vertices
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
      
      // create Intcell object and push back
      Cells().push_back(rcp(new Intcell(num,3,coords,Auxn(),DRT::Element::tri3,
        CouplingInAuxPlane(),lincenter,linvertex[num],linvertex[numplus1],GetDerivAuxn())));
    }
  }
  
  return true;
}


/*----------------------------------------------------------------------*
 |  Linearization of clip polygon vertices (3D)               popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::VertexLinearization(vector<vector<map<int,double> > >& linvertex,
                                              map<int,double>& projpar, bool printderiv)
{
  typedef map<int,double>::const_iterator CI;
  
  // linearize all aux.plane slave and master nodes only ONCE
  // and use these linearizations later during lineclip linearization
  // (this speeds up the vertex linearizations in most cases, as we
  // never linearize the SAME slave or master vertex more than once)
  
  // number of nodes
  int nsrows = SlaveIntElement().NumNode();
  int nmrows = MasterIntElement().NumNode();
      
  // prepare storage for slave and master linearizations
  vector<vector<map<int,double> > > linsnodes(nsrows,vector<map<int,double> >(3));
  vector<vector<map<int,double> > > linmnodes(nmrows,vector<map<int,double> >(3));
  
  if (CouplingInAuxPlane())
  {
    // compute slave linearizations (nsrows)
    for (int i=0;i<nsrows;++i)
    {
      int sid = SlaveIntElement().NodeIds()[i];
      SlaveVertexLinearization(linsnodes[i],sid);
    }
    
    // compute master linearizations (nmrows)
    for (int i=0;i<nmrows;++i)
    {
      int mid = MasterIntElement().NodeIds()[i];
      MasterVertexLinearization(linmnodes[i],mid);
    }
  }
  
  //**********************************************************************
  // Clip polygon vertex linearization
  //**********************************************************************
  // loop over all clip polygon vertices
  for (int i=0;i<(int)Clip().size();++i)
  {
    // references to current vertex and its linearization
    Vertex& currv = Clip()[i];
    vector<map<int,double> >& currlin = linvertex[i];
    
    // decision on vertex type (slave, projmaster, linclip)
    if (currv.VType()==Vertex::slave)
    {
      if (CouplingInAuxPlane())
      {
        // get corresponding slave id
        int sid = currv.Nodeids()[0];
        
        // find corresponding slave node linearization
        int k=0;
        while (k<nsrows){
          if (SlaveIntElement().NodeIds()[k]==sid) break;
          ++k;
        }
        
        // dserror if not found
        if (k==nsrows) dserror("ERROR: Slave Id not found!");
        
        // get the correct slave node linearization
        currlin = linsnodes[k];
      }
      else //(!CouplingInAuxPlane())
      {
        // Vertex = slave node -> Linearization = 0
        // this is the easy case with nothing to do
      }
    }
    else if (currv.VType()==Vertex::projmaster)
    {
      if (CouplingInAuxPlane())
      {
        // get corresponding master id
        int mid = currv.Nodeids()[0];
              
        // find corresponding master node linearization
        int k=0;
        while (k<nmrows){
          if (MasterIntElement().NodeIds()[k]==mid) break;
          ++k;
        }
        
        // dserror if not found
        if (k==nmrows) dserror("ERROR: Master Id not found!");
        
        // get the correct master node linearization
        currlin = linmnodes[k];
      }
      else //(!CouplingInAuxPlane())
      {
        // get corresponding master id and projection alpha
        int mid = currv.Nodeids()[0];
        double alpha = projpar[mid];
        
        //cout << "Coords: " << currv.Coord()[0] << " " << currv.Coord()[1] << endl;
        
        // do master vertex linearization
        MasterVertexLinearization(currv,currlin,mid,alpha);  
      }
    }
    else if (currv.VType()==Vertex::lineclip)
    {
      // get references to the two slave vertices
      int sindex1 = -1;
      int sindex2 = -1;
      for (int j=0;j<(int)SlaveVertices().size();++j)
      {
        if (SlaveVertices()[j].Nodeids()[0]==currv.Nodeids()[0])
          sindex1 = j;
        if (SlaveVertices()[j].Nodeids()[0]==currv.Nodeids()[1])
          sindex2 = j;
      }
      if (sindex1 < 0 || sindex2 < 0 || sindex1==sindex2)
        dserror("ERROR: Lineclip linearization: (S) Something went wrong!");
      
      Vertex* sv1 = &SlaveVertices()[sindex1];
      Vertex* sv2 = &SlaveVertices()[sindex2];
      
      // get references to the two master vertices
      int mindex1 = -1;
      int mindex2 = -1;
      for (int j=0;j<(int)MasterVertices().size();++j)
      {
        if (MasterVertices()[j].Nodeids()[0]==currv.Nodeids()[2])
          mindex1 = j;
        if (MasterVertices()[j].Nodeids()[0]==currv.Nodeids()[3])
          mindex2 = j;
      }
      if (mindex1 < 0 || mindex2 < 0 || mindex1==mindex2)
        dserror("ERROR: Lineclip linearization: (M) Something went wrong!");
      
      Vertex* mv1 = &MasterVertices()[mindex1];
      Vertex* mv2 = &MasterVertices()[mindex2];
      
      // do lineclip vertex linearization
      if (CouplingInAuxPlane())
        LineclipVertexLinearization(currv,currlin,sv1,sv2,mv1,mv2,linsnodes,linmnodes);
      else
        LineclipVertexLinearization(currv,currlin,sv1,sv2,mv1,mv2,projpar);  
    }
    
    else dserror("ERROR: VertexLinearization: Invalid Vertex Type!");
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of slave vertex (3D) AuxPlane               popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::SlaveVertexLinearization(vector<map<int,double> >& currlin,
                                                   int sid)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];
    
  DRT::Element::DiscretizationType dt = SlaveIntElement().Shape();
  if (dt==CElement::tri3 || dt==CElement::tri6)
  {
    scxi[0] = 1.0/3;
    scxi[1] = 1.0/3;
  }
  else if (dt==CElement::quad4 || dt==CElement::quad8 || dt==CElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else dserror("ERROR: SlaveVertexLinearization called for unknown element type");
  
  // evlauate shape functions + derivatives at scxi
  int nrow = SlaveIntElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  SlaveIntElement().EvaluateShape(scxi,sval,sderiv,nrow);
  
  // we need all participating slave nodes
  DRT::Node** snodes = SlaveIntElement().Nodes();
  vector<CONTACT::CNode*> scnodes(nrow);
  
  for (int i=0;i<nrow;++i)
  {
    scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
    if (!scnodes[i]) dserror("ERROR: SlaveVertexLinearization: Null pointer!");
  }
  
  // we also need the corresponding slave node
  DRT::Node* snode = Discret().gNode(sid);
  if (!snode) dserror("ERROR: Cannot find node with gid %",sid);
  CNode* csnode = static_cast<CNode*>(snode);
  
  // map iterator
  typedef map<int,double>::const_iterator CI;
    
  // linearization of element center Auxc()
  vector<map<int,double> > linauxc(3);

  for (int i=0;i<nrow;++i)
  {
    linauxc[0][scnodes[i]->Dofs()[0]] += sval[i];
    linauxc[1][scnodes[i]->Dofs()[1]] += sval[i];
    linauxc[2][scnodes[i]->Dofs()[2]] += sval[i];
  }
  
  // linearization of element normal Auxn()
  vector<map<int,double> >& linauxn = GetDerivAuxn();
  
  // put everything together for slave vertex linearization
  
  // (1) slave node coordinates part
  currlin[0][csnode->Dofs()[0]] += 1.0 - Auxn()[0] * Auxn()[0];
  currlin[0][csnode->Dofs()[1]] -=       Auxn()[1] * Auxn()[0];
  currlin[0][csnode->Dofs()[2]] -=       Auxn()[2] * Auxn()[0];
  currlin[1][csnode->Dofs()[0]] -=       Auxn()[0] * Auxn()[1];
  currlin[1][csnode->Dofs()[1]] += 1.0 - Auxn()[1] * Auxn()[1];
  currlin[1][csnode->Dofs()[2]] -=       Auxn()[2] * Auxn()[1];
  currlin[2][csnode->Dofs()[0]] -=       Auxn()[0] * Auxn()[2];
  currlin[2][csnode->Dofs()[1]] -=       Auxn()[1] * Auxn()[2];
  currlin[2][csnode->Dofs()[2]] += 1.0 - Auxn()[2] * Auxn()[2];
  
  // (2) slave element center coordinates (Auxc()) part
  for (CI p=linauxc[0].begin();p!=linauxc[0].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);
  
  for (CI p=linauxc[1].begin();p!=linauxc[1].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);
  
  for (CI p=linauxc[2].begin();p!=linauxc[2].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);
  
  // (3) slave element normal (Auxn()) part
  double xdotn = (csnode->xspatial()[0]-Auxc()[0]) * Auxn()[0]
               + (csnode->xspatial()[1]-Auxc()[1]) * Auxn()[1]
               + (csnode->xspatial()[2]-Auxc()[2]) * Auxn()[2];
                    
  for (CI p=linauxn[0].begin();p!=linauxn[0].end();++p)
  {
    currlin[0][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (csnode->xspatial()[0]-Auxc()[0]) * Auxn()[k] * (p->second);
  }
  
  for (CI p=linauxn[1].begin();p!=linauxn[1].end();++p)
  {
    currlin[1][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (csnode->xspatial()[1]-Auxc()[1]) * Auxn()[k] * (p->second);
  }
  
  for (CI p=linauxn[2].begin();p!=linauxn[2].end();++p)
  {
    currlin[2][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (csnode->xspatial()[2]-Auxc()[2]) * Auxn()[k] * (p->second);
  }
 
  return true;  
}

/*----------------------------------------------------------------------*
 |  Linearization of projmaster vertex (3D) AuxPlane          popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::MasterVertexLinearization(vector<map<int,double> >& currlin,
                                                    int mid)
{
  // we first need the slave element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double scxi[2];
    
  DRT::Element::DiscretizationType dt = SlaveIntElement().Shape();
  if (dt==CElement::tri3 || dt==CElement::tri6)
  {
    scxi[0] = 1.0/3;
    scxi[1] = 1.0/3;
  }
  else if (dt==CElement::quad4 || dt==CElement::quad8 || dt==CElement::quad9)
  {
    scxi[0] = 0.0;
    scxi[1] = 0.0;
  }
  else dserror("ERROR: MasterVertexLinearization called for unknown element type");
  
  // evlauate shape functions + derivatives at scxi
  int nrow = SlaveIntElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  SlaveIntElement().EvaluateShape(scxi,sval,sderiv,nrow);
  
  // we need all participating slave nodes
  DRT::Node** snodes = SlaveIntElement().Nodes();
  vector<CONTACT::CNode*> scnodes(nrow);
  
  for (int i=0;i<nrow;++i)
  {
    scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
    if (!scnodes[i]) dserror("ERROR: MasterVertexLinearization: Null pointer!");
  }
  
  // we also need the corresponding master node
  DRT::Node* mnode = Discret().gNode(mid);
  if (!mnode) dserror("ERROR: Cannot find node with gid %",mid);
  CNode* cmnode = static_cast<CNode*>(mnode);
  
  // map iterator
  typedef map<int,double>::const_iterator CI;
    
  // linearization of element center Auxc()
  vector<map<int,double> > linauxc(3);

  for (int i=0;i<nrow;++i)
  {
    linauxc[0][scnodes[i]->Dofs()[0]] += sval[i];
    linauxc[1][scnodes[i]->Dofs()[1]] += sval[i];
    linauxc[2][scnodes[i]->Dofs()[2]] += sval[i];
  }
  
  // linearization of element normal Auxn()
  vector<map<int,double> >& linauxn = GetDerivAuxn();
    
  // put everything together for master vertex linearization
  
  // (1) master node coordinates part
  currlin[0][cmnode->Dofs()[0]] += 1.0 - Auxn()[0] * Auxn()[0];
  currlin[0][cmnode->Dofs()[1]] -=       Auxn()[1] * Auxn()[0];
  currlin[0][cmnode->Dofs()[2]] -=       Auxn()[2] * Auxn()[0];
  currlin[1][cmnode->Dofs()[0]] -=       Auxn()[0] * Auxn()[1];
  currlin[1][cmnode->Dofs()[1]] += 1.0 - Auxn()[1] * Auxn()[1];
  currlin[1][cmnode->Dofs()[2]] -=       Auxn()[2] * Auxn()[1];
  currlin[2][cmnode->Dofs()[0]] -=       Auxn()[0] * Auxn()[2];
  currlin[2][cmnode->Dofs()[1]] -=       Auxn()[1] * Auxn()[2];
  currlin[2][cmnode->Dofs()[2]] += 1.0 - Auxn()[2] * Auxn()[2];
  
  // (2) slave element center coordinates (Auxc()) part
  for (CI p=linauxc[0].begin();p!=linauxc[0].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[0] * Auxn()[k] * (p->second);
  
  for (CI p=linauxc[1].begin();p!=linauxc[1].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[1] * Auxn()[k] * (p->second);
  
  for (CI p=linauxc[2].begin();p!=linauxc[2].end();++p)
    for (int k=0;k<3;++k)
      currlin[k][p->first] += Auxn()[2] * Auxn()[k] * (p->second);
  
  // (3) slave element normal (Auxn()) part
  double xdotn = (cmnode->xspatial()[0]-Auxc()[0]) * Auxn()[0]
               + (cmnode->xspatial()[1]-Auxc()[1]) * Auxn()[1]
               + (cmnode->xspatial()[2]-Auxc()[2]) * Auxn()[2];
                    
  for (CI p=linauxn[0].begin();p!=linauxn[0].end();++p)
  {
    currlin[0][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (cmnode->xspatial()[0]-Auxc()[0]) * Auxn()[k] * (p->second);
  }
  
  for (CI p=linauxn[1].begin();p!=linauxn[1].end();++p)
  {
    currlin[1][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (cmnode->xspatial()[1]-Auxc()[1]) * Auxn()[k] * (p->second);
  }
  
  for (CI p=linauxn[2].begin();p!=linauxn[2].end();++p)
  {
    currlin[2][p->first] -= xdotn * (p->second);
    for (int k=0;k<3;++k)
      currlin[k][p->first] -= (cmnode->xspatial()[2]-Auxc()[2]) * Auxn()[k] * (p->second);
  }
    
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of projmaster vertex (3D)                   popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::MasterVertexLinearization(Vertex& currv,
                                                    vector<map<int,double> >& currlin,
                                                    int mid, double alpha)
{
  // get current vertex coordinates (in slave param. space)
  double sxi[2] = {0.0, 0.0};
  sxi[0] = currv.Coord()[0];
  sxi[1] = currv.Coord()[1];
  
  // evlauate shape functions + derivatives at sxi
  int nrow = SlaveIntElement().NumNode();
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  SlaveIntElement().EvaluateShape(sxi,sval,sderiv,nrow);
  
  // build 3x3 factor matrix L
  LINALG::Matrix<3,3> lmatrix(true);
  
  for (int z=0;z<nrow;++z)
  {
    int gid = SlaveIntElement().NodeIds()[z];
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    lmatrix(0,0) += sderiv(z,0) * snode->xspatial()[0];
    lmatrix(1,0) += sderiv(z,0) * snode->xspatial()[1];
    lmatrix(2,0) += sderiv(z,0) * snode->xspatial()[2];
    
    lmatrix(0,0) += alpha * sderiv(z,0) * snode->n()[0];
    lmatrix(1,0) += alpha * sderiv(z,0) * snode->n()[1];
    lmatrix(2,0) += alpha * sderiv(z,0) * snode->n()[2];
    
    lmatrix(0,1) += sderiv(z,1) * snode->xspatial()[0];
    lmatrix(1,1) += sderiv(z,1) * snode->xspatial()[1];
    lmatrix(2,1) += sderiv(z,1) * snode->xspatial()[2];
    
    lmatrix(0,1) += alpha * sderiv(z,1) * snode->n()[0];
    lmatrix(1,1) += alpha * sderiv(z,1) * snode->n()[1];
    lmatrix(2,1) += alpha * sderiv(z,1) * snode->n()[2];
    
    lmatrix(0,2) += sval[z] * snode->n()[0];
    lmatrix(1,2) += sval[z] * snode->n()[1];
    lmatrix(2,2) += sval[z] * snode->n()[2];        
  }
  
  // get inverse of the 3x3 matrix L (in place)
  lmatrix.Invert();
  
  // start to fill linearization maps for current vertex
  typedef map<int,double>::const_iterator CI;
  
  // (1) master node coordinates part
  DRT::Node* mnode = Discret().gNode(mid);
  if (!mnode) dserror("ERROR: Cannot find node with gid %",mid);
  CNode* cmnode = static_cast<CNode*>(mnode);
          
  currlin[0][cmnode->Dofs()[0]] += lmatrix(0,0);
  currlin[0][cmnode->Dofs()[1]] += lmatrix(0,1);
  currlin[0][cmnode->Dofs()[2]] += lmatrix(0,2);
  currlin[1][cmnode->Dofs()[0]] += lmatrix(1,0);
  currlin[1][cmnode->Dofs()[1]] += lmatrix(1,1);
  currlin[1][cmnode->Dofs()[2]] += lmatrix(1,2);
  
  // (2) all slave nodes coordinates part
  // (3) all slave nodes normals part
  for (int z=0;z<nrow;++z)
  {
    int gid = SlaveIntElement().NodeIds()[z];
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    currlin[0][snode->Dofs()[0]] -= sval[z] * lmatrix(0,0);
    currlin[0][snode->Dofs()[1]] -= sval[z] * lmatrix(0,1);
    currlin[0][snode->Dofs()[2]] -= sval[z] * lmatrix(0,2);
    currlin[1][snode->Dofs()[0]] -= sval[z] * lmatrix(1,0);
    currlin[1][snode->Dofs()[1]] -= sval[z] * lmatrix(1,1);
    currlin[1][snode->Dofs()[2]] -= sval[z] * lmatrix(1,2);
    
    // get nodal normal derivative maps (x,y and z components)
    vector<map<int,double> >& derivn = snode->GetDerivN();
    
    for (CI p=derivn[0].begin();p!=derivn[0].end();++p)
    {
      currlin[0][p->first] -= alpha * sval[z] * lmatrix(0,0) * (p->second);
      currlin[1][p->first] -= alpha * sval[z] * lmatrix(1,0) * (p->second);
    }
    for (CI p=derivn[1].begin();p!=derivn[1].end();++p)
    {
      currlin[0][p->first] -= alpha * sval[z] * lmatrix(0,1) * (p->second);
      currlin[1][p->first] -= alpha * sval[z] * lmatrix(1,1) * (p->second);
    }
    for (CI p=derivn[2].begin();p!=derivn[2].end();++p)
    {
      currlin[0][p->first] -= alpha * sval[z] * lmatrix(0,2) * (p->second);
      currlin[1][p->first] -= alpha * sval[z] * lmatrix(1,2) * (p->second);
    }
  }   
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of lineclip vertex (3D) AuxPlane            popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::LineclipVertexLinearization(Vertex& currv,
                                           vector<map<int,double> >& currlin,
                                           Vertex* sv1, Vertex* sv2, Vertex* mv1, Vertex* mv2,
                                           vector<vector<map<int,double> > >& linsnodes,
                                           vector<vector<map<int,double> > >& linmnodes)
{
  // number of nodes
  int nsrows = SlaveIntElement().NumNode();
  int nmrows = MasterIntElement().NumNode();
    
  // iterator
  typedef map<int,double>::const_iterator CI;
  
  // compute factor Z
  double crossZ[3] = {0.0, 0.0, 0.0};
  crossZ[0] = (sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])
           - (sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossZ[1] = (sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])
           - (sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]); 
  crossZ[2] = (sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])
           - (sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  double Zfac = crossZ[0]*Auxn()[0]+crossZ[1]*Auxn()[1]+crossZ[2]*Auxn()[2];
  
  // compute factor N
  double crossN[3] = {0.0, 0.0, 0.0};
  crossN[0] = (sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])
           - (sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossN[1] = (sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])
           - (sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]); 
  crossN[2] = (sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])
           - (sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  double Nfac = crossN[0]*Auxn()[0]+crossN[1]*Auxn()[1]+crossN[2]*Auxn()[2];
  
  // slave edge vector
  double sedge[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k) sedge[k] = sv2->Coord()[k] - sv1->Coord()[k];
    
  // prepare linearization derivZ
  double crossdZ1[3] = {0.0, 0.0, 0.0};
  double crossdZ2[3] = {0.0, 0.0, 0.0};
  double crossdZ3[3] = {0.0, 0.0, 0.0};
  crossdZ1[0] = (mv2->Coord()[1]-mv1->Coord()[1])*Auxn()[2]-(mv2->Coord()[2]-mv1->Coord()[2])*Auxn()[1];
  crossdZ1[1] = (mv2->Coord()[2]-mv1->Coord()[2])*Auxn()[0]-(mv2->Coord()[0]-mv1->Coord()[0])*Auxn()[2];
  crossdZ1[2] = (mv2->Coord()[0]-mv1->Coord()[0])*Auxn()[1]-(mv2->Coord()[1]-mv1->Coord()[1])*Auxn()[0];
  crossdZ2[0] = Auxn()[1]*(sv1->Coord()[2]-mv1->Coord()[2])-Auxn()[2]*(sv1->Coord()[1]-mv1->Coord()[1]);
  crossdZ2[1] = Auxn()[2]*(sv1->Coord()[0]-mv1->Coord()[0])-Auxn()[0]*(sv1->Coord()[2]-mv1->Coord()[2]);
  crossdZ2[2] = Auxn()[0]*(sv1->Coord()[1]-mv1->Coord()[1])-Auxn()[1]*(sv1->Coord()[0]-mv1->Coord()[0]);
  crossdZ3[0] = (sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])-(sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossdZ3[1] = (sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])-(sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]);
  crossdZ3[2] = (sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])-(sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  
  // prepare linearization derivN
  double crossdN1[3] = {0.0, 0.0, 0.0};
  double crossdN2[3] = {0.0, 0.0, 0.0};
  double crossdN3[3] = {0.0, 0.0, 0.0};
  crossdN1[0] = (mv2->Coord()[1]-mv1->Coord()[1])*Auxn()[2]-(mv2->Coord()[2]-mv1->Coord()[2])*Auxn()[1];
  crossdN1[1] = (mv2->Coord()[2]-mv1->Coord()[2])*Auxn()[0]-(mv2->Coord()[0]-mv1->Coord()[0])*Auxn()[2];
  crossdN1[2] = (mv2->Coord()[0]-mv1->Coord()[0])*Auxn()[1]-(mv2->Coord()[1]-mv1->Coord()[1])*Auxn()[0];
  crossdN2[0] = Auxn()[1]*(sv2->Coord()[2]-sv1->Coord()[2])-Auxn()[2]*(sv2->Coord()[1]-sv1->Coord()[1]);
  crossdN2[1] = Auxn()[2]*(sv2->Coord()[0]-sv1->Coord()[0])-Auxn()[0]*(sv2->Coord()[2]-sv1->Coord()[2]);
  crossdN2[2] = Auxn()[0]*(sv2->Coord()[1]-sv1->Coord()[1])-Auxn()[1]*(sv2->Coord()[0]-sv1->Coord()[0]);
  crossdN3[0] = (sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])-(sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossdN3[1] = (sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])-(sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]);
  crossdN3[2] = (sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])-(sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  
  // slave vertex linearization (2x)
  int sid1 = currv.Nodeids()[0];
  int sid2 = currv.Nodeids()[1];
  
  // find corresponding slave node linearizations
  int k=0;
  while (k<nsrows){
    if (SlaveIntElement().NodeIds()[k]==sid1) break;
    ++k;
  }
  
  // dserror if not found
  if (k==nsrows) dserror("ERROR: Slave Id1 not found!");
  
  // get the correct slave node linearization
  vector<map<int,double> >& slavelin0 = linsnodes[k];
  
  k=0;
  while (k<nsrows){
    if (SlaveIntElement().NodeIds()[k]==sid2) break;
    ++k;
  }
  
  // dserror if not found
  if (k==nsrows) dserror("ERROR: Slave Id2 not found!");
  
  // get the correct slave node linearization
  vector<map<int,double> >& slavelin1 = linsnodes[k];
  
  // master vertex linearization (2x)
  int mid1 = currv.Nodeids()[2];
  int mid2 = currv.Nodeids()[3];
  
  // find corresponding master node linearizations
  k=0;
  while (k<nmrows){
    if (MasterIntElement().NodeIds()[k]==mid1) break;
    ++k;
  }
  
  // dserror if not found
  if (k==nmrows) dserror("ERROR: Master Id1 not found!");
  
  // get the correct master node linearization
  vector<map<int,double> >& masterlin0 = linmnodes[k];
  
  k=0;
  while (k<nmrows){
    if (MasterIntElement().NodeIds()[k]==mid2) break;
    ++k;
  }
  
  // dserror if not found
  if (k==nmrows) dserror("ERROR: Master Id2 not found!");
  
  // get the correct master node linearization
  vector<map<int,double> >& masterlin1 = linmnodes[k];
    
  // linearization of element normal Auxn()
  vector<map<int,double> >& linauxn = GetDerivAuxn();
    
  // bring everything together -> lineclip vertex linearization
  for (int k=0;k<3;++k)
  {
    for (CI p=slavelin0[k].begin();p!=slavelin0[k].end();++p)
    {
      currlin[k][p->first] += (p->second);
      currlin[k][p->first] += Zfac/Nfac * (p->second);
      for (int dim=0;dim<3;++dim)
      {
        currlin[dim][p->first] -= sedge[dim] * 1/Nfac * crossdZ1[k] * (p->second);
        currlin[dim][p->first] -= sedge[dim] * Zfac/(Nfac*Nfac) * crossdN1[k] * (p->second);
     
      }
    }
    for (CI p=slavelin1[k].begin();p!=slavelin1[k].end();++p)
    {
      currlin[k][p->first] -= Zfac/Nfac * (p->second);
      for (int dim=0;dim<3;++dim)
      {
        currlin[dim][p->first] += sedge[dim] * Zfac/(Nfac*Nfac) * crossdN1[k] * (p->second);
      }
    }
    for (CI p=masterlin0[k].begin();p!=masterlin0[k].end();++p)
    {
      for (int dim=0;dim<3;++dim)
      {
      currlin[dim][p->first] += sedge[dim] * 1/Nfac * crossdZ1[k] * (p->second);
      currlin[dim][p->first] += sedge[dim] * 1/Nfac * crossdZ2[k] * (p->second);
      currlin[dim][p->first] -= sedge[dim] * Zfac/(Nfac*Nfac) * crossdN2[k] * (p->second);
      }
    }
    for (CI p=masterlin1[k].begin();p!=masterlin1[k].end();++p)
    {
      for (int dim=0;dim<3;++dim)
      {
      currlin[dim][p->first] -= sedge[dim] * 1/Nfac * crossdZ2[k] * (p->second);
      currlin[dim][p->first] += sedge[dim] * Zfac/(Nfac*Nfac) * crossdN2[k] * (p->second);
      }
    }
    for (CI p=linauxn[k].begin();p!=linauxn[k].end();++p)
    {
      for (int dim=0;dim<3;++dim)
      {
      currlin[dim][p->first] -= sedge[dim] * 1/Nfac * crossdZ3[k] * (p->second);
      currlin[dim][p->first] += sedge[dim] * Zfac/(Nfac*Nfac) * crossdN3[k] * (p->second);
      }
    }
  }
    
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of lineclip vertex (3D)                     popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::LineclipVertexLinearization(Vertex& currv,
                                           vector<map<int,double> >& currlin,
                                           Vertex* sv1, Vertex* sv2, Vertex* mv1, Vertex* mv2,
                                           map<int,double>& projpar)
{
  // compute factor Z
  double crossZ[3] = {0.0, 0.0, 0.0};
  crossZ[0] = (sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])
           - (sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossZ[1] = (sv1->Coord()[2]-mv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])
           - (sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]); 
  crossZ[2] = (sv1->Coord()[0]-mv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])
           - (sv1->Coord()[1]-mv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  double Zfac = crossZ[0]*Auxn()[0]+crossZ[1]*Auxn()[1]+crossZ[2]*Auxn()[2];
  
  // compute factor N
  double crossN[3] = {0.0, 0.0, 0.0};
  crossN[0] = (sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[2]-mv1->Coord()[2])
           - (sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossN[1] = (sv2->Coord()[2]-sv1->Coord()[2])*(mv2->Coord()[0]-mv1->Coord()[0])
           - (sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[2]-mv1->Coord()[2]); 
  crossN[2] = (sv2->Coord()[0]-sv1->Coord()[0])*(mv2->Coord()[1]-mv1->Coord()[1])
           - (sv2->Coord()[1]-sv1->Coord()[1])*(mv2->Coord()[0]-mv1->Coord()[0]);
  double Nfac = crossN[0]*Auxn()[0]+crossN[1]*Auxn()[1]+crossN[2]*Auxn()[2];
  
  // slave edge vector
  double sedge[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k) sedge[k] = sv1->Coord()[k] - sv2->Coord()[k];
  
  // prepare linearization derivZ
  double crossdZ1[3] = {0.0, 0.0, 0.0};
  double crossdZ2[3] = {0.0, 0.0, 0.0};
  crossdZ1[0] = Auxn()[1]*(mv2->Coord()[2]-mv1->Coord()[2])-Auxn()[2]*(mv2->Coord()[1]-mv1->Coord()[1]);
  crossdZ1[1] = Auxn()[2]*(mv2->Coord()[0]-mv1->Coord()[0])-Auxn()[0]*(mv2->Coord()[2]-mv1->Coord()[2]);
  crossdZ1[2] = Auxn()[0]*(mv2->Coord()[1]-mv1->Coord()[1])-Auxn()[1]*(mv2->Coord()[0]-mv1->Coord()[0]);
  crossdZ2[0] = Auxn()[1]*(sv1->Coord()[2]-mv1->Coord()[2])-Auxn()[2]*(sv1->Coord()[1]-mv1->Coord()[1]);
  crossdZ2[1] = Auxn()[2]*(sv1->Coord()[0]-mv1->Coord()[0])-Auxn()[0]*(sv1->Coord()[2]-mv1->Coord()[2]);
  crossdZ2[2] = Auxn()[0]*(sv1->Coord()[1]-mv1->Coord()[1])-Auxn()[1]*(sv1->Coord()[0]-mv1->Coord()[0]);
  
  // prepare linearization derivN
  double crossdN1[3] = {0.0, 0.0, 0.0};
  crossdN1[0] = Auxn()[1]*(sv2->Coord()[2]-sv1->Coord()[2])-Auxn()[2]*(sv2->Coord()[1]-sv1->Coord()[1]);
  crossdN1[1] = Auxn()[2]*(sv2->Coord()[0]-sv1->Coord()[0])-Auxn()[0]*(sv2->Coord()[2]-sv1->Coord()[2]);
  crossdN1[2] = Auxn()[0]*(sv2->Coord()[1]-sv1->Coord()[1])-Auxn()[1]*(sv2->Coord()[0]-sv1->Coord()[0]);
  
  // master vertex linearization (2x)
  vector<vector<map<int,double> > > masterlin(2,vector<map<int,double> >(3));
  
  int mid1 = currv.Nodeids()[2];
  double alpha1 = projpar[mid1];
  
  bool found1 = false;
  Vertex* masterv1 = &MasterVertices()[0];
  for (int j=0;j<(int)MasterVertices().size();++j)
  {
    if (MasterVertices()[j].Nodeids()[0]==mid1)
    {
      found1=true;
      masterv1 = &MasterVertices()[j];
      break;
    }
  }
  if (!found1) dserror("ERROR: Lineclip linearization, Master vertex 1 not found!");
  
  MasterVertexLinearization(*masterv1,masterlin[0],mid1,alpha1);
  
  int mid2 = currv.Nodeids()[3];
  double alpha2 = projpar[mid2];
  
  bool found2 = false;
  Vertex* masterv2 = &MasterVertices()[0];
  for (int j=0;j<(int)MasterVertices().size();++j)
  {
    if (MasterVertices()[j].Nodeids()[0]==mid2)
    {
      found2=true;
      masterv2 = &MasterVertices()[j];
      break;
    }
  }
  if (!found2) dserror("ERROR: Lineclip linearization, Master vertex 2 not found!");
  
  MasterVertexLinearization(*masterv2,masterlin[1],mid2,alpha2);
  
  // bring everything together -> lineclip vertex linearization
  typedef map<int,double>::const_iterator CI;
  for (int k=0;k<3;++k)
  {
    for (CI p=masterlin[0][k].begin();p!=masterlin[0][k].end();++p)
    {
      for (int dim=0;dim<2;++dim)
      {
      currlin[dim][p->first] += sedge[dim] * 1/Nfac * crossdZ1[k] * (p->second);
      currlin[dim][p->first] -= sedge[dim] * 1/Nfac * crossdZ2[k] * (p->second);
      currlin[dim][p->first] += sedge[dim] * Zfac/(Nfac*Nfac) * crossdN1[k] * (p->second);
      }
    }
    for (CI p=masterlin[1][k].begin();p!=masterlin[1][k].end();++p)
    {
      for (int dim=0;dim<2;++dim)
      {
      currlin[dim][p->first] += sedge[dim] * 1/Nfac * crossdZ2[k] * (p->second);
      currlin[dim][p->first] -= sedge[dim] * Zfac/(Nfac*Nfac) * crossdN1[k] * (p->second);
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Linearization of clip polygon center (3D)                 popp 02/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::CenterLinearization(const vector<vector<map<int,double> > >& linvertex,
                                              vector<map<int,double> >& lincenter)
{
  // preparations
  int clipsize = (int)(Clip().size());
  typedef map<int,double>::const_iterator CI;
  
  vector<double> clipcenter(3);
  for (int k=0;k<3;++k) clipcenter[k] = 0.0;
  double fac = 0.0;
  
  // first we need node averaged center
  double nac[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      nac[k] += (Clip()[i].Coord()[k] / clipsize);
  
  // loop over all triangles of polygon (1st round: preparations)
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
  
  // build factors for linearization
  double z[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k) z[k] = clipcenter[k];
  double n = fac;
  
  // first we need linearization of node averaged center
  vector<map<int,double> > linnac(3);
  
  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      for (CI p=linvertex[i][k].begin();p!=linvertex[i][k].end();++p)
        linnac[k][p->first] += 1.0/clipsize * (p->second);
    
  // loop over all triangles of polygon (2nd round: linearization)
  for (int i=0; i<clipsize; ++i)
  {
    double xi_i[3] = {0.0, 0.0, 0.0};
    double xi_ip1[3] = {0.0, 0.0, 0.0};
    int iplus1 = 0;
    
    // standard case    
    if (i<clipsize-1)
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[i].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[i+1].Coord()[k];
      iplus1 = i+1;
    }
    // last vertex of clip polygon
    else
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[clipsize-1].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[0].Coord()[k];
      iplus1 = 0;
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
    
    // linearization of cross
    vector<map<int,double> > lincross(3);
    
    for (CI p=linvertex[i][0].begin();p!=linvertex[i][0].end();++p)
    {
      lincross[1][p->first] += diff1[2] * (p->second);
      lincross[1][p->first] += diff2[2] * (p->second);
      lincross[2][p->first] -= diff1[1] * (p->second);
      lincross[2][p->first] -= diff2[1] * (p->second);
    }
    for (CI p=linvertex[i][1].begin();p!=linvertex[i][1].end();++p)
    {
      lincross[0][p->first] -= diff1[2] * (p->second);
      lincross[0][p->first] -= diff2[2] * (p->second);
      lincross[2][p->first] += diff1[0] * (p->second);
      lincross[2][p->first] += diff2[0] * (p->second);   
    }
    for (CI p=linvertex[i][2].begin();p!=linvertex[i][2].end();++p)
    {
      lincross[0][p->first] += diff1[1] * (p->second);
      lincross[0][p->first] += diff2[1] * (p->second);
      lincross[1][p->first] -= diff1[0] * (p->second);
      lincross[1][p->first] -= diff2[0] * (p->second);       
    }
    
    for (CI p=linvertex[iplus1][0].begin();p!=linvertex[iplus1][0].end();++p)
    {
      lincross[1][p->first] -= diff2[2] * (p->second);
      lincross[2][p->first] += diff2[1] * (p->second);       
    }
    for (CI p=linvertex[iplus1][1].begin();p!=linvertex[iplus1][1].end();++p)
    {
      lincross[0][p->first] += diff2[2] * (p->second);
      lincross[2][p->first] -= diff2[0] * (p->second);          
    }
    for (CI p=linvertex[iplus1][2].begin();p!=linvertex[iplus1][2].end();++p)
    {
      lincross[0][p->first] -= diff2[1] * (p->second);
      lincross[1][p->first] += diff2[0] * (p->second);           
    }
    
    for (CI p=linnac[0].begin();p!=linnac[0].end();++p)
    {
      lincross[1][p->first] -= diff1[2] * (p->second);
      lincross[2][p->first] += diff1[1] * (p->second);    
    }
    for (CI p=linnac[1].begin();p!=linnac[1].end();++p)
    {
      lincross[0][p->first] += diff1[2] * (p->second);
      lincross[2][p->first] -= diff1[0] * (p->second);      
    }
    for (CI p=linnac[2].begin();p!=linnac[2].end();++p)
    {
      lincross[0][p->first] -= diff1[1] * (p->second);
      lincross[1][p->first] += diff1[0] * (p->second);     
    }
    
    // linearization of triangle area
    map<int,double> linarea; 
    for (int k=0;k<3;++k)
      for (CI p=lincross[k].begin();p!=lincross[k].end();++p)
        linarea[p->first] += 0.25 / Atri * cross[k] * (p->second);

    // put everything together
    for (int k=0;k<3;++k)
    {
      for (CI p=linvertex[i][k].begin();p!=linvertex[i][k].end();++p)
        lincenter[k][p->first] += 1.0/(3.0*n) * Atri * (p->second);
      
      for (CI p=linvertex[iplus1][k].begin();p!=linvertex[iplus1][k].end();++p)
        lincenter[k][p->first] += 1.0/(3.0*n) * Atri * (p->second);
      
      for (CI p=linnac[k].begin();p!=linnac[k].end();++p)
        lincenter[k][p->first] += 1.0/(3.0*n) * Atri * (p->second);
      
      for (CI p=linarea.begin();p!=linarea.end();++p)
      {
        lincenter[k][p->first] += 1.0/n * 1.0/3.0 * (xi_i[k] + xi_ip1[k] + nac[k]) * (p->second);
        lincenter[k][p->first] -= z[k]/(n*n) * (p->second);
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::IntegrateCells()
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrix M and the weighted gap function g~ on  */
  /* the current integration cell of the slave / master CElement pair   */
  /**********************************************************************/
  
  // create an integrator instance with correct NumGP and Dim
  // it is sufficient to do this once as all Intcells are triangles
  CONTACT::Integrator integrator(Cells()[0]->Shape());
    
  // loop over all integration cells
  for (int i=0;i<(int)(Cells().size());++i)
  {
    // compare intcell area with slave integration element area
    double intcellarea = Cells()[i]->Area();
    double selearea = 0.0;
    if (!CouplingInAuxPlane())
      selearea = SlaveIntElement().Area();
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
    if (intcellarea<CONTACTINTLIM*selearea) continue;
    
    // do the cell integration (integrate and linearize both M and gap)
    // *******************************************************************
    // ************ Coupling with or without auxiliary plane *************
    // *******************************************************************
    int nrow = SlaveElement().NumNode();
    int ncol = MasterElement().NumNode();
    RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),ncol*Dim()));
    RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
    
    if (CouplingInAuxPlane())
    {
      if (Quad())
      {
        // static_cast to make sure to pass in IntElement&
        CONTACT::IntElement& sintref = static_cast<CONTACT::IntElement&>(SlaveIntElement());
        CONTACT::IntElement& mintref = static_cast<CONTACT::IntElement&>(MasterIntElement());
        integrator.IntegrateDerivCell3DAuxPlaneQuad(SlaveElement(),MasterElement(),
                 sintref,mintref,Cells()[i],Auxn(),mseg,gseg);
      }
      else
        integrator.IntegrateDerivCell3DAuxPlane(SlaveElement(),MasterElement(),
                                                  Cells()[i],Auxn(),mseg,gseg);
    }
    else /*(!CouplingInAuxPlane()*/
      integrator.IntegrateDerivCell3D(SlaveElement(),MasterElement(),Cells()[i],mseg,gseg);
    // *******************************************************************
    
    // do the two assemblies into the slave nodes
    // if CONTACTONEMORTARLOOP defined, then AssembleM does M AND D matrices !!!
    integrator.AssembleM(Comm(),SlaveElement(),MasterElement(),*mseg);
    integrator.AssembleG(Comm(),SlaveElement(),*gseg);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling3dQuad::Coupling3dQuad(DRT::Discretization& idiscret, int dim, bool quad,
                                bool auxplane,
                                CONTACT::CElement& sele, CONTACT::CElement& mele,
                                CONTACT::IntElement& sintele,
                                CONTACT::IntElement& mintele) :
CONTACT::Coupling3d(idiscret,dim,quad,auxplane,sele,mele),
sintele_(sintele),
mintele_(mintele)
{
  // 3D quadratic coupling only for aux. plane case
  if (!CouplingInAuxPlane())
    dserror("ERROR: Coupling3dQuad only for auxiliary plane case!");
  
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Quad())
    dserror("ERROR: Coupling3dQuad called for non-quadratic andatz!");
  
  return;
}

#endif //#ifdef CCADISCRET
