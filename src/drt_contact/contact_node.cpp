/*!----------------------------------------------------------------------
\file contact_node.cpp
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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "contact_node.H"
#include "../drt_lib/drt_dserror.H"
#include "contact_element.H"
#include "contact_defines.H"

#include "../drt_nurbs_discret/drt_control_point.H"

CONTACT::CoNodeType CONTACT::CoNodeType::instance_;


DRT::ParObject* CONTACT::CoNodeType::Create( const std::vector<char> & data )
{
  double x[3];
  std::vector<int> dofs(0);
  CONTACT::CoNode* node = new CONTACT::CoNode(0,x,0,0,dofs,false,false);
  node->Unpack(data);
  return node;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 02/10 |
 *----------------------------------------------------------------------*/
CONTACT::CoNodeDataContainer::CoNodeDataContainer():
grow_(1.0e12),
activeold_(false),
kappa_(1.0),
derivn_(0,0),    //init deriv normal to length 0 with 0 entries per direction
derivtxi_(0,0),  //init deriv txi    to length 0 with 0 entries per direction
derivteta_(0,0)  //init deriv teta   to length 0 with 0 entries per direction
{
  // set all tangent entries to 0.0
  for (int i=0;i<3;++i)
  {
    txi()[i]=0.0;
    teta()[i]=0.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoNodeDataContainer::Pack(DRT::PackBuffer& data) const
{
  // add txi_
  DRT::ParObject::AddtoPack(data,txi_,3*sizeof(double));
  // add teta_
  DRT::ParObject::AddtoPack(data,teta_,3*sizeof(double));
  // add grow_
  DRT::ParObject::AddtoPack(data,grow_);
  // add kappa_
  DRT::ParObject::AddtoPack(data,kappa_);
  // add activeold_
  DRT::ParObject::AddtoPack(data,activeold_);

  // no need to pack derivs_
  // (these will evaluated anew anyway)

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoNodeDataContainer::Unpack(std::vector<char>::size_type& position,
                                          const std::vector<char>& data)
{
  // txi_
  DRT::ParObject::ExtractfromPack(position,data,txi_,3*sizeof(double));
  // teta_
  DRT::ParObject::ExtractfromPack(position,data,teta_,3*sizeof(double));
  // grow_
  DRT::ParObject::ExtractfromPack(position,data,grow_);
  // kappa_
  DRT::ParObject::ExtractfromPack(position,data,kappa_);
  // activeold_
  activeold_ = DRT::ParObject::ExtractInt(position,data);

  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoNode::CoNode(int id, const double* coords, const int owner,
                        const int numdof, const std::vector<int>& dofs, const bool isslave,
                        const bool initactive) :
MORTAR::MortarNode(id,coords,owner,numdof,dofs,isslave),
active_(false),
initactive_(initactive),
involvedm_(false),
linsize_(0)   // length of linearization
{

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoNode::CoNode(const CONTACT::CoNode& old) :
MORTAR::MortarNode(old),
active_(old.active_),
initactive_(old.initactive_),
involvedm_(false)
{
  // not yet used and thus not necessarily consistent
  dserror("ERROR: CoNode copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of CoNode and return pointer to it (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoNode* CONTACT::CoNode::Clone() const
{
  CONTACT::CoNode* newnode = new CONTACT::CoNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const CONTACT::CoNode& cnode)
{
  cnode.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  MORTAR::MortarNode::Print(os);
  if (IsSlave())
    if (IsInitActive()) os << " InitActive ";
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class MORTAR::MortarNode
  MORTAR::MortarNode::Pack(data);

  // add active_
  AddtoPack(data,active_);
  // add initactive_
  AddtoPack(data,initactive_);
  // add involved
  AddtoPack(data,involvedm_);
  // add linsize_
  AddtoPack(data,linsize_);
  // add data_
  bool hasdata = (codata_!=Teuchos::null);
  AddtoPack(data,hasdata);
  if (hasdata) codata_->Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class MORTAR::MortarNode
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  MORTAR::MortarNode::Unpack(basedata);

  // active_
  active_ = ExtractInt(position,data);
  // isslave_
  initactive_ = ExtractInt(position,data);
  // isslave_
  involvedm_ = ExtractInt(position,data);
  // isslave_
  linsize_ = ExtractInt(position,data);

  // data_
  bool hasdata = ExtractInt(position,data);
  if (hasdata)
  {
    codata_ = Teuchos::rcp(new CONTACT::CoNodeDataContainer());
    codata_->Unpack(position,data);
  }
  else
  {
    codata_ = Teuchos::null;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the weighted gap                           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::AddgValue(double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddgValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddgValue: function called for boundary node %i", Id());

  // initialize if called for the first time
  if (CoData().Getg()==1.0e12) CoData().Getg()=0.0;

  // add given value to grow_
  CoData().Getg()+=val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivZ' map                           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::AddDerivZValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddZValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddZValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)CoData().GetDerivZ().size()==0)
    CoData().GetDerivZ().resize(NumDof());

  // check row index input
  if ((int)CoData().GetDerivZ().size() <= row)
    dserror("ERROR: AddDerivZValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& zmap = CoData().GetDerivZ()[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::InitializeDataContainer()
{
  // get maximum size of lin vectors
  linsize_ = 0;
  for(int i=0;i<NumElement();++i)
    for(int j=0;j<Elements()[i]->NumNode();++j)
      linsize_ += Elements()[i]->NumDofPerNode(*(Elements()[i]->Nodes()[j]));

  // get maximum size of lin vectors
  dentries_ = 0;
  for(int i=0;i<NumElement();++i)
    for(int j=0;j<Elements()[i]->NumNode();++j)
      dentries_ += Elements()[i]->NumDofPerNode(*(Elements()[i]->Nodes()[j]));

  // only initialize if not yet done
  if (modata_==Teuchos::null && codata_==Teuchos::null)
  {
    codata_=Teuchos::rcp(new CONTACT::CoNodeDataContainer());
    modata_=Teuchos::rcp(new MORTAR::MortarNodeDataContainer());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::ResetDataContainer()
{
  // reset to Teuchos::null
  codata_  = Teuchos::null;
  modata_  = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  Build averaged nodal normal + tangents                    popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoNode::BuildAveragedNormal()
{
  // reset normal and tangents when this method is called
  for (int j=0;j<3;++j)
  {
    MoData().n()[j]=0.0;
    CoData().txi()[j]=0.0;
    CoData().teta()[j]=0.0;
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

  // temporary vector to store nodal normal
  double n_tmp[3]={0.,0.,0.};

  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    CoElement* adjcele = static_cast<CoElement*> (adjeles[i]);

    // build element normal at current node
    // (we have to pass in the index i to be able to store the
    // normal and other information at the right place in elens)
    adjcele->BuildNormalAtNode(Id(),i,elens);

    // add (weighted) element normal to nodal normal n
    for (int j=0;j<3;++j)
      n_tmp[j]+=elens(j,i)/elens(4,i);
  }

  // modify normal in case of symmetry condition
  for (int i=0; i<3; i++)
    if (DbcDofs()[i])
      n_tmp[i]=0.;

  // create unit normal vector
  double length = sqrt(n_tmp[0]*n_tmp[0]+n_tmp[1]*n_tmp[1]+n_tmp[2]*n_tmp[2]);
  if (length==0.0) dserror("ERROR: Nodal normal length 0, node ID %i",Id());
  else             for (int j=0;j<3;++j) MoData().n()[j]=n_tmp[j]/length;

  // create unit tangent vectors
  // (note that this definition is not unique in 3D!)
  double ltxi = 1.0;

  if (NumDof()==2)
  {
    // simple definition for txi
    CoData().txi()[0] = -MoData().n()[1];
    CoData().txi()[1] =  MoData().n()[0];
    CoData().txi()[2] =  0.0;

    // teta is z-axis
    CoData().teta()[0] = 0.0;
    CoData().teta()[1] = 0.0;
    CoData().teta()[2] = 1.0;
  }
  else if (NumDof()==3)
  {
#ifdef CONTACTPSEUDO2D
    // we want to treat a 3D mesh as pseudo 2D contact problem
    // with all nodes fixed in z-direction
    // thus, the second tangent is fixed to (0,0,1)
    CoData().teta()[0] = 0.0;
    CoData().teta()[1] = 0.0;
    CoData().teta()[2] = 1.0;

    // txi follows from corkscrew rule (txi = teta x n)
    CoData().txi()[0] = CoData().teta()[1]*MoData().n()[2]-CoData().teta()[2]*MoData().n()[1];
    CoData().txi()[1] = CoData().teta()[2]*MoData().n()[0]-CoData().teta()[0]*MoData().n()[2];
    CoData().txi()[2] = CoData().teta()[0]*MoData().n()[1]-CoData().teta()[1]*MoData().n()[0];
#else

    if (abs(MoData().n()[0])>1.0e-6 || abs(MoData().n()[1])>1.0e-6 )
    {
      CoData().txi()[0]=-MoData().n()[1];
      CoData().txi()[1]=MoData().n()[0];
      CoData().txi()[2]=0.0;
    }
    else
    {
      CoData().txi()[0]=0.0;
      CoData().txi()[1]=-MoData().n()[2];
      CoData().txi()[2]=MoData().n()[1];
    }

    ltxi = sqrt(CoData().txi()[0]*CoData().txi()[0]+CoData().txi()[1]*CoData().txi()[1]+CoData().txi()[2]*CoData().txi()[2]);
    for (int j=0;j<3;++j) CoData().txi()[j]/=ltxi;

    // teta follows from corkscrew rule (teta = n x txi)
    CoData().teta()[0] = MoData().n()[1]*CoData().txi()[2]-MoData().n()[2]*CoData().txi()[1];
    CoData().teta()[1] = MoData().n()[2]*CoData().txi()[0]-MoData().n()[0]*CoData().txi()[2];
    CoData().teta()[2] = MoData().n()[0]*CoData().txi()[1]-MoData().n()[1]*CoData().txi()[0];

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
void CONTACT::CoNode::DerivAveragedNormal(Epetra_SerialDenseMatrix& elens,
                                          double length, double ltxi)
{
  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  // prepare nodal storage maps for derivative
  if ((int)CoData().GetDerivN().size()==0) CoData().GetDerivN().resize(3,linsize_);
  if ((int)CoData().GetDerivTxi().size()==0) CoData().GetDerivTxi().resize(3,linsize_);
  if ((int)CoData().GetDerivTeta().size()==0) CoData().GetDerivTeta().resize(3,linsize_);

  // loop over all adjacent elements
  for (int i=0;i<nseg;++i)
  {
    CoElement* adjcele = static_cast<CoElement*> (adjeles[i]);

    // build element normal derivative at current node
    adjcele->DerivNormalAtNode(Id(),i,elens,CoData().GetDerivN());
  }

  // modify normal in case of symmetry condition
  for (int i=0; i<3; i++)
    if (DbcDofs()[i])
      CoData().GetDerivN()[i].clear();

  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef GEN::pairedvector<int,double>::const_iterator CI;
  GEN::pairedvector<int,double>& derivnx = CoData().GetDerivN()[0];
  GEN::pairedvector<int,double>& derivny = CoData().GetDerivN()[1];
  GEN::pairedvector<int,double>& derivnz = CoData().GetDerivN()[2];
  GEN::pairedvector<int,double> cderivnx = CoData().GetDerivN()[0];
  GEN::pairedvector<int,double> cderivny = CoData().GetDerivN()[1];
  GEN::pairedvector<int,double> cderivnz = CoData().GetDerivN()[2];
  const double nxnx = MoData().n()[0] * MoData().n()[0];
  const double nxny = MoData().n()[0] * MoData().n()[1];
  const double nxnz = MoData().n()[0] * MoData().n()[2];
  const double nyny = MoData().n()[1] * MoData().n()[1];
  const double nynz = MoData().n()[1] * MoData().n()[2];
  const double nznz = MoData().n()[2] * MoData().n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeysn;
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
    GEN::pairedvector<int,double>& derivtxix = CoData().GetDerivTxi()[0];
    GEN::pairedvector<int,double>& derivtxiy = CoData().GetDerivTxi()[1];

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
    std::map<int,double>& derivtxix = CoData().GetDerivTxi()[0];
    std::map<int,double>& derivtxiy = CoData().GetDerivTxi()[1];
    std::map<int,double>& derivtxiz = CoData().GetDerivTxi()[2];

    for (CI p=derivnx.begin();p!=derivnx.end();++p)
    {
      derivtxiy[p->first] += CoData().teta()[2]*(p->second);
      derivtxiz[p->first] -= CoData().teta()[1]*(p->second);
    }
    for (CI p=derivny.begin();p!=derivny.end();++p)
    {
      derivtxix[p->first] -= CoData().teta()[2]*(p->second);
      derivtxiz[p->first] += CoData().teta()[0]*(p->second);
    }
    for (CI p=derivnz.begin();p!=derivnz.end();++p)
    {
      derivtxix[p->first] += CoData().teta()[1]*(p->second);
      derivtxiy[p->first] -= CoData().teta()[0]*(p->second);
    }
  }
#else
  // unnormalized tangent derivative txi
  // use definitions for txi from BuildAveragedNormal()
  if (abs(MoData().n()[0])>1.0e-6 || abs(MoData().n()[1])>1.0e-6)
  {
    GEN::pairedvector<int,double>& derivtxix = CoData().GetDerivTxi()[0];
    GEN::pairedvector<int,double>& derivtxiy = CoData().GetDerivTxi()[1];

    for (CI p=derivny.begin();p!=derivny.end();++p)
      derivtxix[p->first] -= (p->second);

    for (CI p=derivnx.begin();p!=derivnx.end();++p)
      derivtxiy[p->first] += (p->second);

  }
  else
  {
    GEN::pairedvector<int,double>& derivtxiy = CoData().GetDerivTxi()[1];
    GEN::pairedvector<int,double>& derivtxiz = CoData().GetDerivTxi()[2];

    for (CI p=derivnz.begin();p!=derivnz.end();++p)
      derivtxiy[p->first] -= (p->second);

    for (CI p=derivny.begin();p!=derivny.end();++p)
      derivtxiz[p->first] += (p->second);
  }

  // normalize txi directional derivative
  // (identical to normalization of normal derivative)
  typedef GEN::pairedvector<int,double>::const_iterator CI;
  GEN::pairedvector<int,double>& derivtxix = CoData().GetDerivTxi()[0];
  GEN::pairedvector<int,double>& derivtxiy = CoData().GetDerivTxi()[1];
  GEN::pairedvector<int,double>& derivtxiz = CoData().GetDerivTxi()[2];
  GEN::pairedvector<int,double> cderivtxix = CoData().GetDerivTxi()[0];
  GEN::pairedvector<int,double> cderivtxiy = CoData().GetDerivTxi()[1];
  GEN::pairedvector<int,double> cderivtxiz = CoData().GetDerivTxi()[2];
  const double txtx = CoData().txi()[0] * CoData().txi()[0];
  const double txty = CoData().txi()[0] * CoData().txi()[1];
  const double txtz = CoData().txi()[0] * CoData().txi()[2];
  const double tyty = CoData().txi()[1] * CoData().txi()[1];
  const double tytz = CoData().txi()[1] * CoData().txi()[2];
  const double tztz = CoData().txi()[2] * CoData().txi()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeyst;
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
  GEN::pairedvector<int,double>& derivtetax = CoData().GetDerivTeta()[0];
  GEN::pairedvector<int,double>& derivtetay = CoData().GetDerivTeta()[1];
  GEN::pairedvector<int,double>& derivtetaz = CoData().GetDerivTeta()[2];

  for (CI p=derivnx.begin();p!=derivnx.end();++p)
  {
    derivtetay[p->first] -= CoData().txi()[2]*(p->second);
    derivtetaz[p->first] += CoData().txi()[1]*(p->second);
  }
  for (CI p=derivny.begin();p!=derivny.end();++p)
  {
    derivtetax[p->first] += CoData().txi()[2]*(p->second);
    derivtetaz[p->first] -= CoData().txi()[0]*(p->second);
  }
  for (CI p=derivnz.begin();p!=derivnz.end();++p)
  {
    derivtetax[p->first] -= CoData().txi()[1]*(p->second);
    derivtetay[p->first] += CoData().txi()[0]*(p->second);
  }
  for (CI p=derivtxix.begin();p!=derivtxix.end();++p)
  {
    derivtetay[p->first] += MoData().n()[2]*(p->second);
    derivtetaz[p->first] -= MoData().n()[1]*(p->second);
  }
  for (CI p=derivtxiy.begin();p!=derivtxiy.end();++p)
  {
    derivtetax[p->first] -= MoData().n()[2]*(p->second);
    derivtetaz[p->first] += MoData().n()[0]*(p->second);
  }
  for (CI p=derivtxiz.begin();p!=derivtxiz.end();++p)
  {
    derivtetax[p->first] += MoData().n()[1]*(p->second);
    derivtetay[p->first] -= MoData().n()[0]*(p->second);
  }
  }
#endif // #ifdef CONTACTPSEUDO2D

  return;
}

