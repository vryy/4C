/*!----------------------------------------------------------------------
\file drt_celement.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_celement.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement::CElement(int id, ElementType etype, int owner, 
                            const DRT::Element::DiscretizationType& shape, 
                            const int numnode,
                            const int* nodeids,
                            const bool isslave) :
DRT::Element(id,etype,owner),
shape_(shape),
isslave_(isslave)
{
  SetNodeIds(numnode,nodeids);
  RefArea()=0.0;
  Area()=RefArea();
  searchelements_.resize(0);     //FIXME: Is this necessary???
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement::CElement(const CONTACT::CElement& old) :
DRT::Element(old),
shape_(old.shape_),
isslave_(old.isslave_),
refarea_(old.refarea_),
area_(old.area_),
searchelements_(old.searchelements_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement* CONTACT::CElement::Clone() const
{
  CONTACT::CElement* newele = new CONTACT::CElement(*this);
  return newele;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CElement& element)
{
  element.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Print(ostream& os) const
{
  os << "Contact Element ";
  DRT::Element::Print(os);
  if (isslave_) os << " Slave  Side ";
  else          os << " Master Side ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Element
  vector<char> basedata(0);
  DRT::Element::Pack(basedata);
  AddtoPack(data,basedata);
  // add shape_
  AddtoPack(data,shape_);
  // add isslave_
  AddtoPack(data,isslave_);
  // add refarea_
  AddtoPack(data,refarea_);
  // add area_
  AddtoPack(data,area_);
  // add searchelements_
  AddtoPack(data,&searchelements_,(int)searchelements_.size());

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Element::Unpack(basedata);
  // shape_
  ExtractfromPack(position,data,shape_);
  // isslave_
  ExtractfromPack(position,data,isslave_);
  // refarea_
  ExtractfromPack(position,data,refarea_);
  // area_
  ExtractfromPack(position,data,area_);
  // searchelements_
  ExtractfromPack(position,data,&searchelements_,(int)searchelements_.size());

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate element (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
int CONTACT::CElement::Evaluate(ParameterList&            params,
                                DRT::Discretization&      discretization,
                                vector<int>&              lm,
                                Epetra_SerialDenseMatrix& elemat1,
                                Epetra_SerialDenseMatrix& elemat2,
                                Epetra_SerialDenseVector& elevec1,
                                Epetra_SerialDenseVector& elevec2,
                                Epetra_SerialDenseVector& elevec3)
{
  dserror("CONTACT::CElement::Evaluate not yet impl.");
  return -1;
}

/*----------------------------------------------------------------------*
 |  Get local coordinates for local node id                   popp 12/07|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::LocalCoordinatesOfNode(int lid, double* xi)
{
  // 2D linear and quadratic case
  if (Shape()==line2 || Shape()==line3)
  {
    if (lid==0)
      xi[0]=-1.0;
    else if (lid==1)
      xi[0]= 1.0;
    else if (lid==2)
      xi[0]= 0.0;
    else
      dserror("ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range",lid,Id());
    
    // we are in the 2D case here!
    xi[1]=0.0;
  }
  
  // unknown case (3D not implemented yet)
  else
    dserror("ERROR: LocalCoordinatesOfNode called for unknown element type");
    
  return true;
}

/*----------------------------------------------------------------------*
 |  Get local numbering for global node id                    popp 12/07|
 *----------------------------------------------------------------------*/
int CONTACT::CElement::GetLocalNodeId(int nid)
{
  int lid = -1;  
  
  // look for global ID nid in element's modes
  for (int i=0;i<NumNode();++i)
  {
    if ( *(NodeIds()+i) == nid )
    {
      lid=i;
      break;
    }
  }
  
  if (lid<0)
    dserror("ERROR: Cannot find node % in segment %", nid, Id());
  
  return lid;
}

/*----------------------------------------------------------------------*
 |  Build element normal at node                              popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::BuildNormalAtNode(int nid, vector<double>& n, double& length)
{
  if (n.size()!=3) n.resize(3);
  
  // find this node in my list of nodes and get local numbering
  int lid = GetLocalNodeId(nid);
  
  // get local coordinates for this node
  double xi[2];
  LocalCoordinatesOfNode(lid,xi);
  
  // build an outward unit normal at xi and return it
  ComputeNormalAtXi(xi,n,length);

  return;
}

/*----------------------------------------------------------------------*
 |  Build element normal derivative at node                   popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DerivNormalAtNode(int nid, Epetra_SerialDenseMatrix& elens,
                                          vector<map<int,double> >& derivn)
{
  // find this node in my list of nodes and get local numbering
  int lid = GetLocalNodeId(nid);
  
  // get local coordinates for this node
  double xi[2];
  LocalCoordinatesOfNode(lid,xi);
  
  // build normal derivative at xi and return it
  DerivNormalAtXi(xi,elens,derivn);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal at loc. coord. xi                  popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::ComputeNormalAtXi(double* xi, vector<double>& n, double& length)
{
  int nnodes = NumNode();
  vector<double> val(nnodes);
  vector<double> deriv(nnodes);
  vector<double> g(3);
  
  // get shape function values and derivatives at xi
  EvaluateShape1D(xi, val, deriv, nnodes);

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord(3,nnodes);
  coord = GetNodalCoords();
  
  // build basis vector g
  for (int i=0;i<nnodes;++i)
  {
    g[0]+=deriv[i]*coord(0,i);
    g[1]+=deriv[i]*coord(1,i);
    g[2]+=deriv[i]*coord(2,i);
  }
  
  // only 2D case implemented, normal easy to compute from g
  n[0] = g[1];
  n[1] =-g[0];
  n[2] = 0.0;
  
  length = sqrt(n[0]*n[0]+n[1]*n[1]);
  if (length==0.0)
    dserror("ERROR: ComputeNormalAtXi computes normal of length zero!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal derivative at loc. coord. xi       popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DerivNormalAtXi(double* xi, Epetra_SerialDenseMatrix& elens,
                                        vector<map<int,double> >& derivn)
{
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: ComputeNormalAtXi: Null pointer!");
  
  vector<double> val(nnodes);
  vector<double> deriv(nnodes);
  vector<double> g(3);
  
  // get shape function values and derivatives at xi
  EvaluateShape1D(xi, val, deriv, nnodes);

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord(3,nnodes);
  
  // which column of elens contains this element's normal?
  int col = 2;
  for (int i=0;i<elens.N();++i)
  {
    if (elens(3,i)==Id())
    {
      col=i;
      break;
    }
  }
  if (col>1) dserror("ERROR: Something wrong with columns of elens");
  
  // ... and which does not?
  int ncol = 0;
  if (col==0) ncol=1;
  else ncol=0;
        
  // create directional derivative
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: ComputeNormalAtXi: Null pointer!");
    
    map<int,double>& derivnx = derivn[0];
    map<int,double>& derivny = derivn[1];
    
    // case1: weighted nodal normal
#ifdef CONTACTWNORMAL
    derivnx[mycnode->Dofs()[0]] +=  0;
    derivnx[mycnode->Dofs()[1]] +=  deriv[i];
    derivny[mycnode->Dofs()[0]] += -deriv[i];
    derivny[mycnode->Dofs()[1]] +=  0;
#else
    // case2: unweighted nodal normal (1 adjacent element)
    if (elens.N()==1)
    {
      derivnx[mycnode->Dofs()[0]] +=  0;
      derivnx[mycnode->Dofs()[1]] +=  deriv[i];
      derivny[mycnode->Dofs()[0]] += -deriv[i];
      derivny[mycnode->Dofs()[1]] +=  0;
    }
    // case3: unweighted nodal normal (2 adjacent elements)
    else if (elens.N()==2)
    {
      double lcol = 0.0;
      double lncol = 0.0;
     
      for (int dim=0;dim<3;++dim)
      {
        lcol  += elens(dim,col)*elens(dim,col);
        lncol += elens(dim,ncol)*elens(dim,ncol);
      }
      
      lcol  = sqrt(lcol);
      lncol = sqrt(lncol);
      
      derivnx[mycnode->Dofs()[0]] += -(elens(0,ncol)*elens(1,col)/lcol)*deriv[i];
      derivnx[mycnode->Dofs()[1]] +=  (elens(0,col)*elens(0,ncol)/lcol+lncol)*deriv[i];
      derivny[mycnode->Dofs()[0]] += -(elens(1,col)*elens(1,ncol)/lcol+lncol)*deriv[i];
      derivny[mycnode->Dofs()[1]] +=  (elens(0,col)*elens(1,ncol)/lcol)*deriv[i];
    }
    else
      dserror("ERROR: ComputeNormalAtXi: A 2D CNode can only have 1 or 2 adjacent CElements");
    
#endif // #ifdef CONTACTWNORMAL
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Get nodal coordinates of the element                      popp 01/08|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix CONTACT::CElement::GetNodalCoords()
{
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  LINALG::SerialDenseMatrix coord(3,nnodes);
  if (!mynodes)
    dserror("ERROR: GetNodalCoords: Null pointer!");
  
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode)
      dserror("ERROR: GetNodalCoords: Null pointer!");
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
  }
  
  return coord;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant - LINEAR / QUAD 1D          popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::Jacobian1D(const vector<double>& val,
                                     const vector<double>& deriv,
                                     const LINALG::SerialDenseMatrix& coord)
{
  double jac = 0.0;
  
  // 2D linear case (2noded line element)
  if (Shape()==line2)
    jac = Area()/2;
  
  // 2D quadratic case (3noded line element)
  else if (Shape()==line3)
  {
    double g[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<static_cast<int>(val.size());++i)
    {
      g[0] += deriv[i]*coord(0,i);
      g[1] += deriv[i]*coord(1,i);
      g[2] += deriv[i]*coord(2,i);
    }
    jac = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
  }
  
  // unknown case
  else
    dserror("ERROR: Jacobian1D called for unknown element type!");
    
  return jac;
}

/*----------------------------------------------------------------------*
 |  Compute length (in 3D area) of the element                popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::ComputeArea()
{
  double area = 0.0;
  
  // 2D linear case (2noded line element)
  if (Shape()==line2)
  {
    // no integration necessary (constant Jacobian)
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    coord = GetNodalCoords();
    
    // build vector between the two nodes
    double tang[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      tang[k]=coord(k,1)-coord(k,0);
    }
    area=sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
  }
  
  // 2D quadratic case (3noded line element)
  else if (Shape()==line3)
  {
    // Gauss quadrature
    CONTACT::Integrator integrator(CONTACTNGP,true);
    int nnodes = NumNode();
    double detg = 0.0;
    vector<double> val(nnodes);
    vector<double> deriv(nnodes);
      
    LINALG::SerialDenseMatrix coord(3,nnodes);
    coord = GetNodalCoords();
    
      // loop over all Gauss points, build Jacobian and compute area
    for (int j=0;j<integrator.nGP();++j)
    {
      double gpc[2] = {integrator.Coordinate(j), 0.0};
      EvaluateShape1D(gpc, val, deriv, nnodes);
      detg = Jacobian1D(val,deriv,coord);      
      area+= integrator.Weight(j)*detg;
    }  
  }
  
  // other cases (3D) not implemented yet
  else
    dserror("ERROR: Area computation not implemented for this type of CElement");
  
  return area;
}

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords                  popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::LocalToGlobal(const double* xi, double* globcoord,
                                      bool inttype)
{
  // check input
  if (!xi)
    dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord)
    dserror("ERROR: LocalToGlobal called with globcoord=NULL");
  if (Shape()!=line3 && Shape()!=line2)
    dserror("ERROR: LocalToGlobal called for CEl type != line3");
  
  // collect fundamental data
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");
  LINALG::SerialDenseMatrix coord(3,nnodes);
  vector<double> val(nnodes);
  vector<double> deriv(nnodes);
  
  // Evaluate shape, get nodal coords  and interpolate global coords
  EvaluateShape1D(xi, val, deriv, nnodes);
  for (int i=0;i<3;++i)
    globcoord[i]=0.0;
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: LocalToGlobal: Null pointer!");
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
    
    if (inttype)
    {
      // use shape function values for interpolation
      globcoord[0]+=val[i]*coord(0,i);
      globcoord[1]+=val[i]*coord(1,i);
      globcoord[2]+=val[i]*coord(2,i);
    }
    else
    {
      // use shape function derivatives for interpolation
      globcoord[0]+=deriv[i]*coord(0,i);
      globcoord[1]+=deriv[i]*coord(1,i);
      globcoord[2]+=deriv[i]*coord(2,i);
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Add CElements to potential contact partners               popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::AddSearchElements(const vector<int>& gids)
{
  // check input data and calling element type
  if ((int)gids.size()==0)
    dserror("ERROR: AddSearchElements called with vec of length zero!");
  if (!IsSlave())
    dserror("ERROR: AddSearchElements called for non-slave CElement!");
  
  // loop over all input gids
  for (int i=0;i<(int)gids.size();++i)
  {
    // loop over all search candidates already known
    bool found = false;
    for (int j=0;j<NumSearchElements();++j)
      if (gids[i]==searchelements_[j])
        found = true;
    
    // add new gid to vector of search candidates
    if (!found)
      searchelements_.push_back(gids[i]);
  }
  
  return true;
}

#endif  // #ifdef CCADISCRET
