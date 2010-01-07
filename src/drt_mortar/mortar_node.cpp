/*!----------------------------------------------------------------------
\file mortar_node.cpp
\brief A class for a mortar coupling node

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

#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_defines.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarNode::MortarNode(int id, const double* coords, const int owner,
                const int numdof, const vector<int>& dofs, const bool isslave) :
DRT::Node(id,coords,owner),
isslave_(isslave),
isonbound_(false),
isdbc_(false),
numdof_(numdof),
dofs_(dofs),
closestnode_(-1),
hasproj_(false),
kappa_(1.0)
{
  for (int i=0;i<3;++i)
  {
    n()[i]=0.0;
    uold()[i]=0.0;
    xspatial()[i]=X()[i];
    lm()[i]=0.0;
    lmold()[i]=0.0;
    lmuzawa()[i]=0.0;
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarNode::MortarNode(const MORTAR::MortarNode& old) :
DRT::Node(old),
isslave_(old.isslave_),
isonbound_(old.isonbound_),
isdbc_(old.isdbc_),
numdof_(old.numdof_),
dofs_(old.dofs_),
closestnode_(old.closestnode_),
hasproj_(old.hasproj_),
drows_(old.drows_),
mrows_(old.mrows_),
mmodrows_(old.mmodrows_),
kappa_(old.kappa_)
{
  for (int i=0;i<3;++i)
  {
    n()[i]=old.n_[i];
    uold()[i]=old.uold_[i];
    xspatial()[i]=old.xspatial_[i];
    lm()[i]=old.lm_[i];
    lmold()[i]=old.lmold_[i];
    lmuzawa()[i]=old.lmuzawa_[i];
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarNode* MORTAR::MortarNode::Clone() const
{
  MORTAR::MortarNode* newnode = new MORTAR::MortarNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MORTAR::MortarNode& mrtrnode)
{
  mrtrnode.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this MortarNode (public)                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Mortar ";
  DRT::Node::Print(os);
  
  if (IsSlave()) os << " Slave  ";
  else           os << " Master ";
  
  if (IsOnBound()) os << " Boundary ";
  else             os << " Interior ";
  
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Pack(vector<char>& data) const
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
  // add isonbound_
  AddtoPack(data,isonbound_);
  // add isdbc_
  AddtoPack(data,isdbc_);
  // add numdof_
  AddtoPack(data,numdof_);
  // add dofs_
  AddtoPack(data,dofs_);
  // add xspatial_
  AddtoPack(data,xspatial_,3);
  // add n_
  AddtoPack(data,n_,3);
  // add uold_
  AddtoPack(data,uold_,3);
  // add lm_
  AddtoPack(data,lm_,3);
  // add lmold_
  AddtoPack(data,lmold_,3);
  // add lmuzawa_
  AddtoPack(data,lmuzawa_,3);
  // add closestnode_
  AddtoPack(data,closestnode_);
  // add hasproj_
  AddtoPack(data,hasproj_);
  // add kappa_
  AddtoPack(data,kappa_);
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Unpack(const vector<char>& data)
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
  // isonbound_
  ExtractfromPack(position,data,isonbound_);
  // isdbc_
  ExtractfromPack(position,data,isdbc_);
  // numdof_
  ExtractfromPack(position,data,numdof_);
  // dofs_
  ExtractfromPack(position,data,dofs_);
  // xspatial_
  ExtractfromPack(position,data,xspatial_,3);
  // n_
  ExtractfromPack(position,data,n_,3);
  // uold_
  ExtractfromPack(position,data,uold_,3);
  // lm_
  ExtractfromPack(position,data,lm_,3);
  // lmold_
  ExtractfromPack(position,data,lmold_,3);
  // lmuzawa_
  ExtractfromPack(position,data,lmuzawa_,3);
  // closestnode_
  ExtractfromPack(position,data,closestnode_);
  // hasproj_
  ExtractfromPack(position,data,hasproj_);
  // kappa_
  ExtractfromPack(position,data,kappa_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D' map                                popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddDValue(int& row, int& col, double& val)
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
void MORTAR::MortarNode::AddMValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddMValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddMValue: function called for boundary node %i", Id());

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
void MORTAR::MortarNode::AddMmodValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddMmodValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddMmodValue: function called for boundary node %i", Id());

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
 |  Build averaged nodal normal                               popp 12/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::BuildAveragedNormal()
{
  // reset normal and tangents when this method is called
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
  // elens(5,i): length/area of element itself
  //**********************************************************************
  Epetra_SerialDenseMatrix elens(6,nseg);

  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    MortarElement* adjmrtrele = static_cast<MortarElement*> (adjeles[i]);

    // build element normal at current node
    // (we have to pass in the index i to be able to store the
    // normal and other information at the right place in elens)
    adjmrtrele->BuildNormalAtNode(Id(),i,elens);

    // add (weighted) element normal to nodal normal n
    for (int j=0;j<3;++j)
      n()[j]+=elens(j,i)/elens(4,i);
  }

  // create unit normal vector
  double length = sqrt(n()[0]*n()[0]+n()[1]*n()[1]+n()[2]*n()[2]);
  if (length==0.0) dserror("ERROR: Nodal normal length 0, node ID %i",Id());
  else             for (int j=0;j<3;++j) n()[j]/=length;

  return;
}

/*----------------------------------------------------------------------*
 |  Find closest node from given node set                     popp 01/08|
 *----------------------------------------------------------------------*/
MORTAR::MortarNode* MORTAR::MortarNode::FindClosestNode(const RCP<DRT::Discretization> intdis,
                                                        const RCP<Epetra_Map> nodesearchmap, double& mindist)
{
  MortarNode* closestnode = NULL;

  // loop over all nodes of the DRT::Discretization that are
  // included in the given Epetra_Map ("brute force" search)
  for(int i=0; i<nodesearchmap->NumMyElements();++i)
  {
    int gid = nodesearchmap->GID(i);
    DRT::Node* node = intdis->gNode(gid);
    if (!node) dserror("ERROR: FindClosestNode: Cannot find node with gid %",gid);
    MortarNode* mrtrnode = static_cast<MortarNode*>(node);

    // build distance between the two nodes
    double dist = 0.0;
    const double* p1 = xspatial();
    const double* p2 = mrtrnode->xspatial();

    for (int j=0;j<3;++j)
      dist+=(p1[j]-p2[j])*(p1[j]-p2[j]);
    dist=sqrt(dist);

    // new closest node found, update
    if (dist <= mindist)
    {
      mindist=dist;
      closestnode=mrtrnode;
    }
  }

  if (!closestnode)
    dserror("ERROR: FindClosestNode: No closest node found at all!");

  return closestnode;
}

/*----------------------------------------------------------------------*
 | Check mesh re-initialization for this node                 popp 12/09|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarNode::CheckMeshDistortion(double& relocation, double& limit)
{
  // initialize return parameter
  bool ok = true;
  
  // loop over all adjacent elements of this node
  for (int i=0;i<NumElement();++i)
  {
    // get the current element
    DRT::Element* ele = Elements()[i];
    if (!ele) dserror("ERROR: Cannot find element with lid %",i);
    MortarElement* mrtrele = static_cast<MortarElement*>(ele);
    
    // minimal edge size of the current element
    double minedgesize = mrtrele->MinEdgeSize();
    
    // check whether relocation is not too large
    if (relocation > limit * minedgesize)
    {
      // print information to screen
      cout << "\n*****************WARNING***********************" << endl;
      cout << "Checking distortion for CNode:     " << Id() << endl;
      cout << "Relocation distance:               " << relocation << endl;
      cout << "AdjEle: " << mrtrele->Id() << "\tLimit*MinEdgeSize: " << limit*minedgesize << endl;
      cout << "*****************WARNING***********************" << endl;
      
      // set return parameter and stop
      ok = false;
      break;
    }
  }
  
  return ok; 
}

#endif  // #ifdef CCADISCRET
