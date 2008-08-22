/*!----------------------------------------------------------------------
\file drt_cnode.cpp
\brief A class for a contact node

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_cnode.H"
#include "../drt_lib/drt_dserror.H"
#include "drt_celement.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(int id, const double* coords, const int owner, 
                      const int numdof, const vector<int>& dofs, const bool isslave,
                      const bool initactive) :
DRT::Node(id,coords,owner),
isslave_(isslave),
initactive_(initactive),
isonbound_(false),
numdof_(numdof),
dofs_(dofs),
closestnode_(-1),
hasproj_(false),
active_(false),
slip_(false),
grow_(1.0e12)
{
  for (int i=0;i<3;++i)
  {
    dbc()[i]=false;
    n()[i]=0.0;
    u()[i]=0.0;
    xspatial()[i]=X()[i];
    lm()[i]=0.0;
    lmold()[i]=0.0;
    jump()[i]=0.0;
  }
   
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(const CONTACT::CNode& old) :
DRT::Node(old),
isslave_(old.isslave_),
initactive_(old.initactive_),
isonbound_(old.isonbound_),
numdof_(old.numdof_),
dofs_(old.dofs_),
closestnode_(old.closestnode_),
hasproj_(old.hasproj_),
active_(old.active_),
slip_(old.slip_),
drows_(old.drows_),
mrows_(old.mrows_),
mmodrows_(old.mmodrows_),
grow_(old.grow_)
{
  for (int i=0;i<3;++i)
  {
    dbc()[i]=old.dbc_[i];
    n()[i]=old.n_[i];
    u()[i]=old.u_[i];
    xspatial()[i]=old.xspatial_[i];
    lm()[i]=old.lm_[i];
    lmold()[i]=old.lmold_[i];
    jump()[i]=old.jump_[i];
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of CNode and return pointer to it (public)  |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode* CONTACT::CNode::Clone() const
{
  CONTACT::CNode* newnode = new CONTACT::CNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CNode& cnode)
{
  cnode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  DRT::Node::Print(os);
  if (IsSlave())
  {
    os << " Slave  ";
    if (IsInitActive()) os << " InitActive  ";
  }
  else           os << " Master ";
  if (IsOnBound()) os << " Boundary ";
  else             os << " Interior ";
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Node
  vector<char> basedata(0);
  DRT::Node::Pack(basedata);
  AddtoPack(data,basedata);
  // add isslave_
  AddtoPack(data,isslave_);
  // add initactive_
  AddtoPack(data,initactive_);
  // add isonbound_
  AddtoPack(data,isonbound_);
  // add dbc_
  AddtoPack(data,dbc_,3);
  // add numdof_
  AddtoPack(data,numdof_);
  // add dofs_
  AddtoPack(data,dofs_);
  // add xspatial_
  AddtoPack(data,xspatial_,3);
  // add n_
  AddtoPack(data,n_,3);
  // add u_
  AddtoPack(data,u_,3);
  // add lm_
  AddtoPack(data,lm_,3);
  // add lmold_
  AddtoPack(data,lmold_,3);
  // add jump_
  AddtoPack(data,jump_,3);
  // add closestnode_
  AddtoPack(data,closestnode_);
  // add hasproj_
  AddtoPack(data,hasproj_);
  // add active_
  AddtoPack(data,active_);
  // add slip_
  AddtoPack(data,slip_);
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Node
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Node::Unpack(basedata);
  // isslave_
  ExtractfromPack(position,data,isslave_);
  // isslave_
  ExtractfromPack(position,data,initactive_);
  // isonbound_
  ExtractfromPack(position,data,isonbound_);
  // dbc_
  ExtractfromPack(position,data,dbc_,3);
  // numdof_
  ExtractfromPack(position,data,numdof_);
  // dofs_
  ExtractfromPack(position,data,dofs_);
  // xspatial_
  ExtractfromPack(position,data,xspatial_,3);
  // n_
  ExtractfromPack(position,data,n_,3);
  // u_
  ExtractfromPack(position,data,u_,3);
  // lm_
  ExtractfromPack(position,data,lm_,3);
  // lmold_
  ExtractfromPack(position,data,lmold_,3);
  // jump_
  ExtractfromPack(position,data,jump_,3);
  // closestnode_
  ExtractfromPack(position,data,closestnode_);
  // hasproj_
  ExtractfromPack(position,data,hasproj_);
  // active_
  ExtractfromPack(position,data,active_);
  // slip_
  ExtractfromPack(position,data,slip_);

  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D' map                                popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddDValue(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());
  
  // check if this has been called before
  if ((int)drows_.size()==0)
    drows_.resize(NumDof());
  
  // check row index input
  if ((int)drows_.size()<=row)
    dserror("ERROR: AddDValue: tried to access invalid row index!");
  
  // add the pair (col,val) to the given row
  map<int,double>& dmap = drows_[row];
  dmap[col] += val;
    
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'M' map                                popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddMValue(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());
    
  // check if this has been called before
  if ((int)mrows_.size()==0)
    mrows_.resize(NumDof());
    
  // check row index input
  if ((int)mrows_.size()<=row)
    dserror("ERROR: AddMValue: tried to access invalid row index!");
    
  // add the pair (col,val) to the given row
  map<int,double>& mmap = mrows_[row];
  mmap[col] += val;
      
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'Mmod' map                             popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddMmodValue(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());
    
  // check if this has been called before
  if ((int)mmodrows_.size()==0)
    mmodrows_.resize(NumDof());
    
  // check row index input
  if ((int)mmodrows_.size()<=row)
    dserror("ERROR: AddMmodValue: tried to access invalid row index!");
    
  // add the pair (col,val) to the given row
  map<int,double>& mmodmap = mmodrows_[row];
  mmodmap[col] += val;
      
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the weighted gap                           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddgValue(double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());
  
  // initialize if called for the first time
  if (grow_==1.0e12) grow_=0;
  
  // add given value to grow_
  grow_+=val;
  return;
}

/*----------------------------------------------------------------------*
 |  Build averaged nodal normal                               popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::BuildAveragedNormal()
{
  // reset normal when this method is called
  for (int j=0;j<3;++j) n()[j]=0.0;
    
  int nseg = NumElement();
  DRT::Element** adjeles = Elements();
  
  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length of element itself
  //**********************************************************************
  Epetra_SerialDenseMatrix elens(6,nseg);
  
  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    CElement* adjcele = static_cast<CElement*> (adjeles[i]);

    // build element normal at current node
    vector<double> elen(3);
    double elenlength;
    adjcele->BuildNormalAtNode(Id(),elen,elenlength);
    double wgt = adjcele->Area();
    
    // add (weighted) element normal to nodal normal n
    for (int j=0;j<3;++j)
    {
      elens(j,i) = elen[j];
      
#ifdef CONTACTWNORMAL
      n()[j]+=wgt*elen[j]/elenlength;
#else
      n()[j]+=elen[j]/elenlength;
#endif // #ifdef CONTACTWNORMAL
    }
    
    // store some element info for normal linearization
    elens(3,i) = adjcele->Id();
    elens(4,i) = elenlength;
    elens(5,i) = wgt;
  }
  
  // create unit normal vector
  double length = sqrt(n()[0]*n()[0]+n()[1]*n()[1]+n()[2]*n()[2]);
  
  if (length==0.0)
    dserror("ERROR: Nodal normal of length zero, node ID %i",Id());
  else
    for (int j=0;j<3;++j) n()[j]/=length;
  
  // computation of nodal normal is finished here...!!!
  
  //**********************************************************************
  // Redefine length for directional derivative
  // (this step is necessary due to the fact that the unnormalized
  // nodal normal is scaled for making its linearization easier.
  // In both weighted and unweighted case this is done by muliplying
  // with the lengths of all adjacent elements!
  //**********************************************************************
  for (int i=0;i<nseg;++i)
    length *= elens(4,i);
  
  // build directional derivative of averaged nodal normal
  DerivAveragedNormal(elens,length);
 
  return;
}

/*----------------------------------------------------------------------*
 |  Build directional derivative of nodal normal              popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::DerivAveragedNormal(Epetra_SerialDenseMatrix& elens,
                                         double length)
{
  // prepare nodal storage maps for derivative
  if ((int)GetDerivN().size()==0) GetDerivN().resize(NumDof());
  if ((int)GetDerivT().size()==0) GetDerivT().resize(NumDof());
  
  int nseg = NumElement();
  DRT::Element** adjeles = Elements();
  
  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    CElement* adjcele = static_cast<CElement*> (adjeles[i]);
    
    // build element normal derivative at current node
    adjcele->DerivNormalAtNode(Id(),elens,GetDerivN());
  }
  
  // computation of directional derivative of unnormalized nodal
  // normal is finished here...!!!
  
  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with refernce / copy of derivative maps!)
  typedef map<int,double>::const_iterator CI;
  map<int,double>& derivnx = GetDerivN()[0];
  map<int,double>& derivny = GetDerivN()[1];
  map<int,double> copyderivnx = GetDerivN()[0];
  map<int,double> copyderivny = GetDerivN()[1];
  double nxnx = n()[0] * n()[0];
  double nxny = n()[0] * n()[1];
  double nyny = n()[1] * n()[1];
  
  // normalize x-components
  for (CI p=derivnx.begin();p!=derivnx.end();++p)
  {
    int col = p->first;
    double val = p->second;
    derivnx[col] = (val-nxnx*val-nxny*copyderivny[col])/length;
  }
  
  // normalize y-components
  for (CI p=derivny.begin();p!=derivny.end();++p)
  {
    int col = p->first;
    double val = p->second;
    derivny[col] = (val-nxny*copyderivnx[col]-nyny*val)/length;
  }
  
  // get directional derivative of nodal tangent "for free"
  // (we just have to use the orthogonality of n and t)
  if (NumDof()==3) dserror("ERROR: Not yet implemented for 3D case");
  map<int,double>& derivtx = GetDerivT()[0];
  map<int,double>& derivty = GetDerivT()[1];
  
  for (CI p=derivny.begin();p!=derivny.end();++p)
    derivtx[p->first] = -(p->second);
  for (CI p=derivnx.begin();p!=derivnx.end();++p)
    derivty[p->first] = (p->second);
  
  return;
}

/*----------------------------------------------------------------------*
 |  Find closest node from given node set                     popp 01/08|
 *----------------------------------------------------------------------*/
CONTACT::CNode* CONTACT::CNode::FindClosestNode(const RCP<DRT::Discretization> intdis,
                                                const RCP<Epetra_Map> nodesearchmap,
                                                double& mindist)
{
  CNode* closestnode = NULL;
  
  // loop over all nodes of the DRT::Discretization that are
  // included in the given Epetra_Map
  for(int i=0; i<nodesearchmap->NumMyElements();++i)
  {
    int gid = nodesearchmap->GID(i);
    DRT::Node* node = intdis->gNode(gid);
    if (!node) dserror("ERROR: FindClosestNode: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    
    // build distance between the two nodes
    double dist = 0.0;
    const double* p1 = xspatial();
    const double* p2 = cnode->xspatial();
    
    for (int j=0;j<3;++j)
      dist+=(p1[j]-p2[j])*(p1[j]-p2[j]);
    dist=sqrt(dist);
    
    // new closest node found, update
    if (dist <= mindist)
    {
      mindist=dist;
      closestnode=cnode;
    }
  }
  
  if (!closestnode)
    dserror("ERROR: FindClosestNode: No closest node found at all!");
  
  return closestnode;
}

#endif  // #ifdef CCADISCRET
