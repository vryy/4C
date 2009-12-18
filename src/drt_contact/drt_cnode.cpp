/*!----------------------------------------------------------------------
\file drt_cnode.cpp
\brief A class for a contact node

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
isdbc_(false),
numdof_(numdof),
dofs_(dofs),
closestnode_(-1),
hasproj_(false),
active_(false),
activeold_(false),
slip_(false),
grow_(1.0e12)
{
  for (int i=0;i<3;++i)
  {
    n()[i]=0.0;
    txi()[i]=0.0;
    teta()[i]=0.0;
    uold()[i]=0.0;
    traction()[i]=0.0;
    tractionold()[i]=0.0;
    xspatial()[i]=X()[i];
    lm()[i]=0.0;
    lmold()[i]=0.0;
    lmuzawa()[i]=0.0;
    jump()[i]=0.0;
  }
  
  Kappa() = 1.0;

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
isdbc_(old.isdbc_),
numdof_(old.numdof_),
dofs_(old.dofs_),
closestnode_(old.closestnode_),
hasproj_(old.hasproj_),
active_(old.active_),
activeold_(old.activeold_),
slip_(old.slip_),
drows_(old.drows_),
mrows_(old.mrows_),
drowsold_(old.drowsold_),
mrowsold_(old.mrowsold_),
mmodrows_(old.mmodrows_),
mnodes_(old.mnodes_),
mnodesold_(old.mnodesold_),
grow_(old.grow_)
{
  for (int i=0;i<3;++i)
  {
    n()[i]=old.n_[i];
    txi()[i]=old.txi_[i];
    teta()[i]=old.teta_[i];
    uold()[i]=old.uold_[i];
    traction()[i]=old.traction_[i];
    tractionold()[i]=old.tractionold_[i];
    xspatial()[i]=old.xspatial_[i];
    lm()[i]=old.lm_[i];
    lmold()[i]=old.lmold_[i];
    lmuzawa()[i]=old.lmuzawa_[i];
    jump()[i]=old.jump_[i];
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of CNode and return pointer to it (public)  |
 |                                                           mwgee 10/07|
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
  // add txi_
  AddtoPack(data,txi_,3);
  // add teta_
  AddtoPack(data,teta_,3);
  // add uold_
  AddtoPack(data,uold_,3);
  // add traction_
  AddtoPack(data,traction_,3);
  // add tractionold_
  AddtoPack(data,tractionold_,3);
  // add lm_
  AddtoPack(data,lm_,3);
  // add lmold_
  AddtoPack(data,lmold_,3);
  // add lmuzawa_
  AddtoPack(data,lmuzawa_,3);
  // add jump_
  AddtoPack(data,jump_,3);
  // add closestnode_
  AddtoPack(data,closestnode_);
  // add hasproj_
  AddtoPack(data,hasproj_);
  // add active_
  AddtoPack(data,active_);
  // add activeold_
  AddtoPack(data,activeold_);
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
  // txi_
  ExtractfromPack(position,data,txi_,3);
  // teta_
  ExtractfromPack(position,data,teta_,3);
  // uold_
  ExtractfromPack(position,data,uold_,3);
  // traction_
  ExtractfromPack(position,data,traction_,3);
  // tractionold_
  ExtractfromPack(position,data,tractionold_,3);
  // lm_
  ExtractfromPack(position,data,lm_,3);
  // lmold_
  ExtractfromPack(position,data,lmold_,3);
  // lmuzawa_
  ExtractfromPack(position,data,lmuzawa_,3);
  // jump_
  ExtractfromPack(position,data,jump_,3);
  // closestnode_
  ExtractfromPack(position,data,closestnode_);
  // hasproj_
  ExtractfromPack(position,data,hasproj_);
  // active_
  ExtractfromPack(position,data,active_);
  // activeold_
  ExtractfromPack(position,data,activeold_);
  // slip_
  ExtractfromPack(position,data,slip_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D' map                                popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddDValue(int& row, int& col, double& val)
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
void CONTACT::CNode::AddMValue(int& row, int& col, double& val)
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
void CONTACT::CNode::AddMmodValue(int& row, int& col, double& val)
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
 |  Add a value to the 'SNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddSNode(int node)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddSnode: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddSNode: function called for boundary node %i", Id());

  snodes_.insert(node);
	
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'MNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddMNode(int node)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddMNode: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddMNode: function called for boundary node %i", Id());

  mnodes_.insert(node);
	
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the weighted gap                           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddgValue(double& val)
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
 |  Add a value to the 'DerivZ' map                           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddDerivZValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddZValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddZValue: function called for boundary node %i", Id());
    
  // check if this has been called before
  if ((int)derivz_.size()==0)
    derivz_.resize(NumDof());
    
  // check row index input
  if ((int)derivz_.size() <= row)
    dserror("ERROR: AddDerivZValue: tried to access invalid row index!");
    
  // add the pair (col,val) to the given row
  map<int,double>& zmap = derivz_[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivJump' map                     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::AddDerivJumpValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddJumpValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddJumpValue: function called for boundary node %i", Id());
    
  // check if this has been called before
  if ((int)derivjump_.size()==0)
    derivjump_.resize(NumDof());
    
  // check row index input
  if ((int)derivjump_.size() <= row)
    dserror("ERROR: AddDerivJumpValue: tried to access invalid row index!");
    
  // add the pair (col,val) to the given row
  map<int,double>& zmap = derivjump_[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries of D and M to old ones             gitterle 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::StoreDMOld()
{
  // copy drows_ to drowsold_

  // reset old nodal Mortar maps
  for (int j=0;j<(int)(GetDOld().size());++j)
  (GetDOld())[j].clear();
  for (int j=0;j<(int)((GetMOld()).size());++j)
  (GetMOld())[j].clear();

  // clear and zero nodal vectors
  drowsold_.clear();
  mrowsold_.clear();
  drowsold_.resize(0);
	mrowsold_.resize(0);
	
	// write drows_ to drowsold_
	drowsold_ = drows_;
  mrowsold_ = mrows_;
  
  // also vectors containing the according master nodes
  mnodesold_.clear();
  mnodesold_ = mnodes_;
  
  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries penalty tractions to old ones      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::StoreTracOld()
{
  // write entries to old ones
  for (int j=0;j<3;++j)
    tractionold_[j]=traction_[j];

  return;
}

/*----------------------------------------------------------------------*
 |  Build averaged nodal normal + tangents                    popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::BuildAveragedNormal()
{
  // reset normal and tangents when this method is called
  for (int j=0;j<3;++j)
  {
    n()[j]=0.0;
    txi()[j]=0.0;
    teta()[j]=0.0;
  }

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
    CElement* adjcele = static_cast<CElement*> (adjeles[i]);

    // build element normal at current node
    // (we have to pass in the index i to be able to store the
    // normal and other information at the right place in elens)
    adjcele->BuildNormalAtNode(Id(),i,elens);

    // add (weighted) element normal to nodal normal n
    for (int j=0;j<3;++j)
      n()[j]+=elens(j,i)/elens(4,i);
  }

  // create unit normal vector
  double length = sqrt(n()[0]*n()[0]+n()[1]*n()[1]+n()[2]*n()[2]);
  if (length==0.0) dserror("ERROR: Nodal normal length 0, node ID %i",Id());
  else             for (int j=0;j<3;++j) n()[j]/=length;

  // create unit tangent vectors
  // (note that this definition is not unique in 3D!)
  double ltxi = 1.0;

  if (NumDof()==2)
  {
    // simple definition for txi
    txi()[0] = -n()[1];
    txi()[1] =  n()[0];
    txi()[2] =  0.0;

    // teta is z-axis
    teta()[0] = 0.0;
    teta()[1] = 0.0;
    teta()[2] = 1.0;
  }
  else if (NumDof()==3)
  {
#ifdef CONTACTPSEUDO2D
    // we want to treat a 3D mesh as pseudo 2D contact problem
    // with all nodes fixed in z-direction
    // thus, the second tangent is fixed to (0,0,1)
    teta()[0] = 0.0;
    teta()[1] = 0.0;
    teta()[2] = 1.0;

    // txi follows from corkscrew rule (txi = teta x n)
    txi()[0] = teta()[1]*n()[2]-teta()[2]*n()[1];
    txi()[1] = teta()[2]*n()[0]-teta()[0]*n()[2];
    txi()[2] = teta()[0]*n()[1]-teta()[1]*n()[0];
#else
    // arbitrary definition for txi
    if (abs(n()[2])>1.0e-6)
    {
      txi()[0]=1.0;
      txi()[1]=1.0;
      txi()[2]=(-n()[0]-n()[1])/n()[2];
    }
    else if (abs(n()[1])>1.0e-6)
    {
      txi()[0]=1.0;
      txi()[2]=1.0;
      txi()[1]=(-n()[0]-n()[2])/n()[1];
    }
    else if (abs(n()[0])>1.0e-6)
    {
      txi()[1]=1.0;
      txi()[2]=1.0;
      txi()[0]=(-n()[1]-n()[2])/n()[0];
    }
    else
      dserror("ERROR: Something wrong with nodal normal");

    ltxi = sqrt(txi()[0]*txi()[0]+txi()[1]*txi()[1]+txi()[2]*txi()[2]);
    for (int j=0;j<3;++j) txi()[j]/=ltxi;

    // teta follows from corkscrew rule (teta = n x txi)
    teta()[0] = n()[1]*txi()[2]-n()[2]*txi()[1];
    teta()[1] = n()[2]*txi()[0]-n()[0]*txi()[2];
    teta()[2] = n()[0]*txi()[1]-n()[1]*txi()[0];
#endif // #ifdef CONTACTPSEUDO2D
  }
  else
    dserror("ERROR: Contact problems must be either 2D or 3D");

  // build linearization of averaged nodal normal and tangents
  DerivAveragedNormal(elens,length,ltxi);

  return;
}

/*----------------------------------------------------------------------*
 |  Build directional deriv. of nodal normal + tangents       popp 09/08|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::DerivAveragedNormal(Epetra_SerialDenseMatrix& elens,
                                         double length, double ltxi)
{
  // prepare nodal storage maps for derivative
  if ((int)GetDerivN().size()==0) GetDerivN().resize(3);
  if ((int)GetDerivTxi().size()==0) GetDerivTxi().resize(3);
  if ((int)GetDerivTeta().size()==0) GetDerivTeta().resize(3);

  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    CElement* adjcele = static_cast<CElement*> (adjeles[i]);

    // build element normal derivative at current node
    adjcele->DerivNormalAtNode(Id(),i,elens,GetDerivN());
  }

  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef map<int,double>::const_iterator CI;
  map<int,double>& derivnx = GetDerivN()[0];
  map<int,double>& derivny = GetDerivN()[1];
  map<int,double>& derivnz = GetDerivN()[2];
  map<int,double> cderivnx = GetDerivN()[0];
  map<int,double> cderivny = GetDerivN()[1];
  map<int,double> cderivnz = GetDerivN()[2];
  double nxnx = n()[0] * n()[0];
  double nxny = n()[0] * n()[1];
  double nxnz = n()[0] * n()[2];
  double nyny = n()[1] * n()[1];
  double nynz = n()[1] * n()[2];
  double nznz = n()[2] * n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  vector<int> allkeysn;
  for (CI p=derivnx.begin();p!=derivnx.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);

  }
  for (CI p=derivny.begin();p!=derivny.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);

  }
  for (CI p=derivnz.begin();p!=derivnz.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }

  // normalize x-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivnx[allkeysn[j]];
    derivnx[allkeysn[j]] = (val-nxnx*val-nxny*cderivny[allkeysn[j]]-nxnz*cderivnz[allkeysn[j]])/length;
  }

  // normalize y-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivny[allkeysn[j]];
    derivny[allkeysn[j]] = (val-nxny*cderivnx[allkeysn[j]]-nyny*val-nynz*cderivnz[allkeysn[j]])/length;
  }

  // normalize z-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivnz[allkeysn[j]];
    derivnz[allkeysn[j]] = (val-nxnz*cderivnx[allkeysn[j]]-nynz*cderivny[allkeysn[j]]-nznz*val)/length;
  }

  //**********************************************************************
  // tangent derivatives 2D
  //**********************************************************************
  if (NumDof()==2)
  {
    // get directional derivative of nodal tangent txi "for free"
    // (we just have to use the orthogonality of n and t)
    // the directional derivative of nodal tangent teta is 0
    map<int,double>& derivtxix = GetDerivTxi()[0];
    map<int,double>& derivtxiy = GetDerivTxi()[1];

    for (CI p=derivny.begin();p!=derivny.end();++p)
      derivtxix[p->first] = -(p->second);
    for (CI p=derivnx.begin();p!=derivnx.end();++p)
      derivtxiy[p->first] = (p->second);
  }

  //**********************************************************************
  // tangent derivatives 3D
  //**********************************************************************
  else
  {
#ifdef CONTACTPSEUDO2D
    // trivial tangent derivative teta
    // this is 0 as teta is fixed to (0,0,1)

    // get normalized tangent derivative txi
    // use corkscrew rule from BuildAveragedNormal()
    map<int,double>& derivtxix = GetDerivTxi()[0];
    map<int,double>& derivtxiy = GetDerivTxi()[1];
    map<int,double>& derivtxiz = GetDerivTxi()[2];

    for (CI p=derivnx.begin();p!=derivnx.end();++p)
    {
      derivtxiy[p->first] += teta()[2]*(p->second);
      derivtxiz[p->first] -= teta()[1]*(p->second);
    }
    for (CI p=derivny.begin();p!=derivny.end();++p)
    {
      derivtxix[p->first] -= teta()[2]*(p->second);
      derivtxiz[p->first] += teta()[0]*(p->second);
    }
    for (CI p=derivnz.begin();p!=derivnz.end();++p)
    {
      derivtxix[p->first] += teta()[1]*(p->second);
      derivtxiy[p->first] -= teta()[0]*(p->second);
    }
  }
#else
    // unnormalized tangent derivative txi
    // use definitions for txi from BuildAveragedNormal()
    if (abs(n()[2])>1.0e-6)
    {
      map<int,double>& derivtxiz = GetDerivTxi()[2];
      for (CI p=derivnx.begin();p!=derivnx.end();++p)
        derivtxiz[p->first] -= 1/n()[2]*(p->second);
      for (CI p=derivny.begin();p!=derivny.end();++p)
        derivtxiz[p->first] -= 1/n()[2]*(p->second);
      for (CI p=derivnz.begin();p!=derivnz.end();++p)
        derivtxiz[p->first] += (n()[0]+n()[1])/(n()[2]*n()[2])*(p->second);

    }
    else if (abs(n()[1])>1.0e-6)
    {
      map<int,double>& derivtxiy = GetDerivTxi()[1];
      for (CI p=derivnx.begin();p!=derivnx.end();++p)
        derivtxiy[p->first] -= 1/n()[1]*(p->second);
      for (CI p=derivny.begin();p!=derivny.end();++p)
        derivtxiy[p->first] += (n()[0]+n()[2])/(n()[1]*n()[1])*(p->second);
      for (CI p=derivnz.begin();p!=derivnz.end();++p)
        derivtxiy[p->first] -= 1/n()[1]*(p->second);
    }
    else if (abs(n()[0])>1.0e-6)
    {
      map<int,double>& derivtxix = GetDerivTxi()[0];
      for (CI p=derivnx.begin();p!=derivnx.end();++p)
        derivtxix[p->first] += (n()[1]+n()[2])/(n()[0]*n()[0])*(p->second);
      for (CI p=derivny.begin();p!=derivny.end();++p)
        derivtxix[p->first] -= 1/n()[0]*(p->second);
      for (CI p=derivnz.begin();p!=derivnz.end();++p)
        derivtxix[p->first] -= 1/n()[0]*(p->second);
    }
    else
      dserror("ERROR: Something wrong with nodal normal");

    // normalize txi directional derivative
    // (identical to normalization of normal derivative)
    typedef map<int,double>::const_iterator CI;
    map<int,double>& derivtxix = GetDerivTxi()[0];
    map<int,double>& derivtxiy = GetDerivTxi()[1];
    map<int,double>& derivtxiz = GetDerivTxi()[2];
    map<int,double> cderivtxix = GetDerivTxi()[0];
    map<int,double> cderivtxiy = GetDerivTxi()[1];
    map<int,double> cderivtxiz = GetDerivTxi()[2];
    double txtx = txi()[0] * txi()[0];
    double txty = txi()[0] * txi()[1];
    double txtz = txi()[0] * txi()[2];
    double tyty = txi()[1] * txi()[1];
    double tytz = txi()[1] * txi()[2];
    double tztz = txi()[2] * txi()[2];

    // build a vector with all keys from x,y,z maps
    // (we need this in order not to miss any entry!)
    vector<int> allkeyst;
    for (CI p=derivtxix.begin();p!=derivtxix.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);

    }
    for (CI p=derivtxiy.begin();p!=derivtxiy.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);

    }
    for (CI p=derivtxiz.begin();p!=derivtxiz.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }

    // normalize x-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val = cderivtxix[allkeyst[j]];
      derivtxix[allkeyst[j]] = (val-txtx*val-txty*cderivtxiy[allkeyst[j]]-txtz*cderivtxiz[allkeyst[j]])/ltxi;
    }

    // normalize y-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val =cderivtxiy[allkeyst[j]];
      derivtxiy[allkeyst[j]] = (val-txty*cderivtxix[allkeyst[j]]-tyty*val-tytz*cderivtxiz[allkeyst[j]])/ltxi;
    }

    // normalize z-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val = cderivtxiz[allkeyst[j]];
      derivtxiz[allkeyst[j]] = (val-txtz*cderivtxix[allkeyst[j]]-tytz*cderivtxiy[allkeyst[j]]-tztz*val)/ltxi;
    }

    // get normalized tangent derivative teta
    // use corkscrew rule from BuildAveragedNormal()
    map<int,double>& derivtetax = GetDerivTeta()[0];
    map<int,double>& derivtetay = GetDerivTeta()[1];
    map<int,double>& derivtetaz = GetDerivTeta()[2];

    for (CI p=derivnx.begin();p!=derivnx.end();++p)
    {
      derivtetay[p->first] -= txi()[2]*(p->second);
      derivtetaz[p->first] += txi()[1]*(p->second);
    }
    for (CI p=derivny.begin();p!=derivny.end();++p)
    {
      derivtetax[p->first] += txi()[2]*(p->second);
      derivtetaz[p->first] -= txi()[0]*(p->second);
    }
    for (CI p=derivnz.begin();p!=derivnz.end();++p)
    {
      derivtetax[p->first] -= txi()[1]*(p->second);
      derivtetay[p->first] += txi()[0]*(p->second);
    }
    for (CI p=derivtxix.begin();p!=derivtxix.end();++p)
    {
      derivtetay[p->first] += n()[2]*(p->second);
      derivtetaz[p->first] -= n()[1]*(p->second);
    }
    for (CI p=derivtxiy.begin();p!=derivtxiy.end();++p)
    {
      derivtetax[p->first] -= n()[2]*(p->second);
      derivtetaz[p->first] += n()[0]*(p->second);
    }
    for (CI p=derivtxiz.begin();p!=derivtxiz.end();++p)
    {
      derivtetax[p->first] += n()[1]*(p->second);
      derivtetay[p->first] -= n()[0]*(p->second);
    }
  }
#endif // #ifdef CONTACTPSEUDO2D

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
  // included in the given Epetra_Map ("brute force" search)
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
