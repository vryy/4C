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
  Area()=0.0;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement::CElement(const CONTACT::CElement& old) :
DRT::Element(old),
shape_(old.shape_),
isslave_(old.isslave_),
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
  if (isslave_) os << " Slave  ";
  else          os << " Master ";

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
  dserror("CONTACT::CElement::Evaluate not implemented!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  Get local coordinates for local node id                   popp 12/07|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::LocalCoordinatesOfNode(int lid, double* xi)
{
  // 2D linear case (2noded line element)
  // 2D quadratic case (3noded line element)
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
  
  // 3D linear case (2noded triangular element)
  // 3D quadratic case (3noded triangular element)
  else if (Shape()==tri3 || Shape()==tri6)
  {
    if (lid==0)
    {
      xi[0]=0.0; xi[1]=0.0;
    }
    else if (lid==1)
    {
      xi[0]=1.0; xi[1]=0.0;
    }
    else if (lid==2)
    {
      xi[0]=0.0; xi[1]=1.0;
    }
    else if (lid==3)
    {
      xi[0]=0.5; xi[1]=0.0;
    }
    else if (lid==4)
    {
      xi[0]=0.5; xi[1]=0.5;
    }
    else if (lid==5)
    {
      xi[0]=0.0; xi[1]=0.5;
    }
    else
      dserror("ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range",lid,Id());
  }
  
  // 3D bilinear case (4noded quadrilateral element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (Shape()==quad4 || Shape()==quad8 || Shape()==quad9)
  {
    if (lid==0)
    {
      xi[0]=-1.0; xi[1]=-1.0;
    }
    else if (lid==1)
    {
      xi[0]=1.0; xi[1]=-1.0;
    }
    else if (lid==2)
    {
      xi[0]=1.0; xi[1]=1.0;
    }
    else if (lid==3)
    {
      xi[0]=-1.0; xi[1]=1.0;
    }
    else if (lid==4)
    {
      xi[0]=0.0; xi[1]=-1.0;
    }
    else if (lid==5)
    {
      xi[0]=1.0; xi[1]=0.0;
    }
    else if (lid==6)
    {
      xi[0]=0.0; xi[1]=1.0;
    }
    else if (lid==7)
    {
      xi[0]=-1.0; xi[1]=0.0;
    }
    else if (lid==8)
    {
      xi[0]=0.0; xi[1]=0.0;
    }
    else
      dserror("ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range",lid,Id());
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
  
  // look for global ID nid in element's nodes
  for (int i=0;i<NumNode();++i)
    if (NodeIds()[i] == nid)
    {
      lid=i;
      break;
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
  if (n.size()!=3) dserror("ERROR: Normal vector must be of length 3!");
  
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
  // empty local basis vectors
  vector<double> gxi(3);
  vector<double> geta(3);
  DRT::Element::DiscretizationType dt = Shape();
  
  if (dt==line2 || dt==line3)
  {
    // metrics routine gives local basis vectors
    Metrics(xi,gxi,geta);
    
    // normal easy to compute from gxi in 2D
    n[0] =  gxi[1];
    n[1] = -gxi[0];
    n[2] =  0.0;
  }
  
  else if (dt==tri3 || dt==quad4 || dt==tri6 || dt==quad8 || dt==quad9)
  {
    // metrics routine gives local basis vectors
    Metrics(xi,gxi,geta);
    
    // n is cross product of gxi and geta
    n[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
    n[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
    n[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
  }
  
  else
    dserror("ERROR: ComPuteNormalAtXi not implemented for this element type");
  
  // compute length of element normal
  length = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (length==0.0) dserror("ERROR: ComputeNormalAtXi gives normal of length 0!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal derivative at loc. coord. xi       popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DerivNormalAtXi(double* xi, Epetra_SerialDenseMatrix& elens,
                                        vector<map<int,double> >& derivn)
{
  DRT::Element::DiscretizationType dt = Shape();
  
  if (dt==line2 || dt==line3)
  {
    int nnodes = NumNode();
    DRT::Node** mynodes = Nodes();
    if (!mynodes) dserror("ERROR: ComputeNormalAtXi: Null pointer!");
    
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    
    // get shape function values and derivatives at xi
    EvaluateShape(xi, val, deriv, nnodes);
  
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
    map<int,double>& derivnx = derivn[0];
    map<int,double>& derivny = derivn[1];
    
    //**********************************************************************
    // For the weighted normal case, the element lengths enter the nodal
    // normal formulation. They have to be linearized as well, which is an
    // element operation only. Thus, we can compute this part of th nodal
    // normal derivative before looping over the element nodes!
    //**********************************************************************
  #ifdef CONTACTWNORMAL
    // add directional derivative of element area
    typedef map<int,double>::const_iterator CI;
    map<int,double> derivarea;
    DerivArea(derivarea);
    
    // case1: weighted nodal normal (1 adjacent element)
    if (elens.N()==1)
    {
      for (CI p=derivarea.begin();p!=derivarea.end();++p)
      {
        derivnx[p->first] += elens(0,col)*(p->second);
        derivny[p->first] += elens(1,col)*(p->second);
      }
    }
    // case2: weighted nodal normal (2 adjacent elements)
    else if (elens.N()==2)
    {
      for (CI p=derivarea.begin();p!=derivarea.end();++p)
      {
        derivnx[p->first] += elens(4,ncol)*elens(0,col)*(p->second);
        derivny[p->first] += elens(4,ncol)*elens(1,col)*(p->second);
      }
    } 
  #endif // #ifdef CONTACTWNORMAL
    
    // loop over all nodes for directional derivative
    // (this is the part we in both the weighted and unweighted normal case)
    for (int i=0;i<nnodes;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      if (!mycnode) dserror("ERROR: ComputeNormalAtXi: Null pointer!");
      
  #ifdef CONTACTWNORMAL
      // case1: weighted nodal normal (1 adjacent element)
      if (elens.N()==1)
      {
        derivnx[mycnode->Dofs()[0]] += 0;
        derivnx[mycnode->Dofs()[1]] += elens(5,col)*deriv(i,0);
        derivny[mycnode->Dofs()[0]] -= elens(5,col)*deriv(i,0);
        derivny[mycnode->Dofs()[1]] += 0;
      }
      // case2: weighted nodal normal (2 adjacent elements)
      else if (elens.N()==2)
      {
        derivnx[mycnode->Dofs()[0]] -= (elens(0,ncol)*elens(1,col)*elens(5,ncol)/elens(4,col))*deriv(i,0);
        derivnx[mycnode->Dofs()[1]] += (elens(0,col)*elens(0,ncol)*elens(5,ncol)/elens(4,col)+elens(4,ncol)*elens(5,col))*deriv(i,0);
        derivny[mycnode->Dofs()[0]] -= (elens(1,col)*elens(1,ncol)*elens(5,ncol)/elens(4,col)+elens(4,ncol)*elens(5,col))*deriv(i,0);
        derivny[mycnode->Dofs()[1]] += (elens(0,col)*elens(1,ncol)*elens(5,ncol)/elens(4,col))*deriv(i,0);
      }
  #else
      // case3: unweighted nodal normal (1 adjacent element)
      if (elens.N()==1)
      {
        derivnx[mycnode->Dofs()[0]] += 0;
        derivnx[mycnode->Dofs()[1]] += deriv(i,0);
        derivny[mycnode->Dofs()[0]] -= deriv(i,0);
        derivny[mycnode->Dofs()[1]] += 0;
      }
      // case4: unweighted nodal normal (2 adjacent elements)
      else if (elens.N()==2)
      {
        derivnx[mycnode->Dofs()[0]] -= (elens(0,ncol)*elens(1,col)/elens(4,col))*deriv(i,0);
        derivnx[mycnode->Dofs()[1]] += (elens(0,col)*elens(0,ncol)/elens(4,col)+elens(4,ncol))*deriv(i,0);
        derivny[mycnode->Dofs()[0]] -= (elens(1,col)*elens(1,ncol)/elens(4,col)+elens(4,ncol))*deriv(i,0);
        derivny[mycnode->Dofs()[1]] += (elens(0,col)*elens(1,ncol)/elens(4,col))*deriv(i,0);
      }
      else
        dserror("ERROR: ComputeNormalAtXi: A 2D CNode can only have 1 or 2 adjacent CElements");
      
  #endif // #ifdef CONTACTWNORMAL
    }
  }
  
  else if (dt==tri3 || dt==quad4 || dt==tri6 || dt==quad8 || dt==quad9)
  {
    // not yet implemented
  }
  
  else
    dserror("ERROR: DerivNormalAtXi not implemented for this element type");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Get nodal coordinates of the element                      popp 01/08|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix CONTACT::CElement::GetNodalCoords()
{
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: GetNodalCoords: Null pointer!");
  LINALG::SerialDenseMatrix coord(3,nnodes);
  
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: GetNodalCoords: Null pointer!");
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
  }
  
  return coord;
}

/*----------------------------------------------------------------------*
 |  Evaluate element metrics (local basis vectors)            popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Metrics(double* xi, vector<double>& gxi,
                                vector<double>& geta)
{
  int nnodes = NumNode();
  int dim = 0;
  DRT::Element::DiscretizationType dt = Shape();
  if (dt==line2 || dt==line3) dim = 2;
  else if (dt==tri3 || dt==quad4 || dt==tri6 || dt==quad8 || dt==quad9) dim = 3;
  else dserror("ERROR: Metrics called for unknown element t<pe");
  
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,dim-1);
  
  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes);

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord = GetNodalCoords();
  
  // build basis vectors gxi and geta
  for (int i=0;i<nnodes;++i)
  {
    // first local basis vector
    gxi[0] += deriv(i,0)*coord(0,i);
    gxi[1] += deriv(i,0)*coord(1,i);
    gxi[2] += deriv(i,0)*coord(2,i);

    // second local basis vector
    if (dim==3)
    {
      geta[0] += deriv(i,1)*coord(0,i);
      geta[1] += deriv(i,1)*coord(1,i);
      geta[2] += deriv(i,1)*coord(2,i);
    }
  }
  
  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant                             popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::Jacobian(double* xi)
{
  double jac = 0.0;
  vector<double> gxi(3);
  vector<double> geta(3);
  DRT::Element::DiscretizationType dt = Shape();
    
  // 2D linear case (2noded line element)
  if (dt==line2)
    jac = Area()/2;
  
  // 2D quadratic case (3noded line element)
  else if (dt==line3)
  {
    // metrics routine gives local basis vectors
    Metrics(xi,gxi,geta);
    
    // length of gxi
    jac = sqrt(gxi[0]*gxi[0]+gxi[1]*gxi[1]+gxi[2]*gxi[2]);
  }
  
  // 3D linear case (3noded triangular element)
  else if (dt==tri3)
    jac = Area()*2;
  
  // 3D bilinear case (4noded quadrilateral element)
  // 3D quadratic case (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt==quad4 || dt==tri6 || dt==quad8 || dt==quad9)
  {
    // metrics routine gives local basis vectors
    Metrics(xi,gxi,geta);
    
    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
    cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
    cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
    jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  }
    
  // unknown case
  else
    dserror("ERROR: Jacobian called for unknown element type!");
    
  return jac;
}


/*----------------------------------------------------------------------*
 |  Evaluate derivative J,xi of Jacobian determinant          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DJacDXi(double* djacdxi,
                                const LINALG::SerialDenseVector& val,
                                const LINALG::SerialDenseMatrix& deriv,
                                const LINALG::SerialDenseMatrix& secderiv,
                                const LINALG::SerialDenseMatrix& coord)
{
  // the derivative dJacdXi
  djacdxi[0] = 0.0;
  djacdxi[1] = 0.0;
  
  // 2D linear case (2noded line element)
  if (Shape()==line2)
    djacdxi[0] = 0.0;
  
  // 2D quadratic case (3noded line element)
  else if (Shape()==line3)
  {
    double g[3] = {0.0, 0.0, 0.0};
    double gsec[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<val.Length();++i)
    {
      g[0] += deriv(i,0)*coord(0,i);
      g[1] += deriv(i,0)*coord(1,i);
      g[2] += deriv(i,0)*coord(2,i);
      
      gsec[0] += secderiv(i,0)*coord(0,i);
      gsec[1] += secderiv(i,0)*coord(1,i);
      gsec[2] += secderiv(i,0)*coord(2,i);
    }
    
    // the Jacobian itself
    double jac = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
    
    // compute dJacdXi
    for (int dim=0;dim<3;++dim)
      djacdxi[0] += g[dim]*gsec[dim]/jac;
  }
  
  // unknown case
  else
    dserror("ERROR: dJacdXi called for unknown element type!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det.              popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DerivJacobian(const LINALG::SerialDenseVector& val,
                                      const LINALG::SerialDenseMatrix& deriv,
                                      const LINALG::SerialDenseMatrix& coord,
                                      map<int,double>& derivjac)
{
  // get element nodes
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: DerivJacobian: Null pointer!");
  
  // loop over all nodes
  double g[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<nnodes;++i)
  {
    g[0] += deriv(i,0)*coord(0,i);
    g[1] += deriv(i,0)*coord(1,i);
    g[2] += deriv(i,0)*coord(2,i);
  }
  
  // the Jacobian itself
  double jac = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
  
  // *********************************************************************
  // compute Jacobian derivative
  // *********************************************************************
  // (loop over all nodes and over all nodal dofs to capture all
  // potential dependencies of the Jacobian. Note that here we only
  // need to compute the DIRECT derivative of Lin(J), as the current
  // GP coordinate does not change! The derivative DJacDXi is done in
  // a special function (see above)!
  // *********************************************************************
  for (int i=0;i<nnodes;++i)
  {
    CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[i]);
    if (!mycnode) dserror("ERROR: DerivJacobian: Null pointer!");
    
    for (int j=0;j<mycnode->NumDof();++j)
    {
      int col = mycnode->Dofs()[j];
      double tmp = (1/jac)*deriv(i,0)*g[j];
      derivjac[col] += tmp;
    } 
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute length / area of the element                      popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::ComputeArea()
{
  double area = 0.0;
  DRT::Element::DiscretizationType dt = Shape();
  
  // 2D linear case (2noded line element)
  if (dt==line2)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord = GetNodalCoords();
    
    // build vector between the two nodes
    double tang[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      tang[k]=coord(k,1)-coord(k,0);
    }
    area=sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
  }
  
  // 2D quadratic case (3noded line element)
  else if (dt==line3)
  {
    // Gauss quadrature with correct NumGP and Dim
    CONTACT::Integrator integrator(dt);
    double detg = 0.0;
      
    // loop over all Gauss points, build Jacobian and compute area
    for (int j=0;j<integrator.nGP();++j)
    {
      double gpc[2] = {integrator.Coordinate(j,0), 0.0};
      detg = Jacobian(gpc);      
      area+= integrator.Weight(j)*detg;
    }  
  }
  
  // 3D linear case (3noded triangular element)
  else if (dt==tri3)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord = GetNodalCoords();
    
    // build vectors between the three nodes
    double t1[3] = {0.0, 0.0, 0.0};
    double t2[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      t1[k]=coord(k,1)-coord(k,0);
      t2[k]=coord(k,2)-coord(k,0);
    }
    
    // cross product of t1 and t2
    double t1xt2[3] = {0.0, 0.0, 0.0};
    t1xt2[0] = t1[1]*t2[2]-t1[2]*t2[1];
    t1xt2[1] = t1[2]*t2[0]-t1[0]*t2[2];
    t1xt2[2] = t1[0]*t2[1]-t1[1]*t2[0];
    area=0.5*sqrt(t1xt2[0]*t1xt2[0]+t1xt2[1]*t1xt2[1]+t1xt2[2]*t1xt2[2]);
  }
  
  // 3D bilinear case    (4noded quadrilateral element)
  // 3D quadratic case   (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt==quad4 || dt==tri6 || dt==quad8 || dt==quad9)
  {
    // Gauss quadrature with correct NumGP and Dim
    CONTACT::Integrator integrator(dt);
    double detg = 0.0;
   
    // loop over all Gauss points, build Jacobian and compute area
    for (int j=0;j<integrator.nGP();++j)
    {
      double gpc[2] = {integrator.Coordinate(j,0), integrator.Coordinate(j,1)};
      detg = Jacobian(gpc);      
      area+= integrator.Weight(j)*detg;
    }
  }
  
  // other cases not implemented yet
  else
    dserror("ERROR: Area computation not implemented for this type of CElement");
  
  return area;
}

/*----------------------------------------------------------------------*
 |  Compute length / area linearization of the element        popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::DerivArea(map<int,double>& derivarea)
{
  // 2D linear case (2noded line element)
  if (Shape()==line2)
  {
    // no integration necessary (constant Jacobian)
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord = GetNodalCoords();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
        
    // build val, deriv etc. at xi=0.0 (arbitrary)
    double xi[2] = {0.0, 0.0};
    EvaluateShape(xi,val,deriv,nnodes);
    
    // get Jacobian derivative
    DerivJacobian(val,deriv,coord,derivarea);
    
    // multiply all entries with factor 2
    typedef map<int,double>::const_iterator CI;
    for (CI p=derivarea.begin();p!=derivarea.end();++p)
      derivarea[p->first] = 2*(p->second);
  }
  
  // 2D quadratic case (3noded line element)
  else if (Shape()==line3)
  {
    // Gauss quadrature with correct NumGP and Dim
    CONTACT::Integrator integrator(Shape());
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord = GetNodalCoords();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes);
    
    // loop over all Gauss points, build Jacobian derivative
    for (int j=0;j<integrator.nGP();++j)
    {
      double wgt = integrator.Weight(j);
      double gpc[2] = {integrator.Coordinate(j,0), 0.0};
      EvaluateShape(gpc, val, deriv, nnodes);
      
      // get Jacobian derivative
      map<int,double> derivjac;
      DerivJacobian(val,deriv,coord,derivjac);
      
      // add current GP to Area derivative map
      typedef map<int,double>::const_iterator CI;
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        derivarea[p->first] += wgt*(p->second);
    }  
  }
  
  // other cases not implemented yet
  else
    dserror("ERROR: Area derivative not implemented for this type of CElement");
    
  return;
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
    dserror("ERROR: LocalToGlobal called for CEl type != line2/3");
  
  // collect fundamental data
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");
  LINALG::SerialDenseMatrix coord(3,nnodes);
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,1);
  
  // Evaluate shape, get nodal coords  and interpolate global coords
  EvaluateShape(xi, val, deriv, nnodes);
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
      globcoord[0]+=deriv(i,0)*coord(0,i);
      globcoord[1]+=deriv(i,0)*coord(1,i);
      globcoord[2]+=deriv(i,0)*coord(2,i);
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
