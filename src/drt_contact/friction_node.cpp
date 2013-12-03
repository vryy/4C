/*!----------------------------------------------------------------------
\file friction_node.cpp
\brief A class for a frictional contact node

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

#include "friction_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "../drt_lib/drt_dserror.H"

CONTACT::FriNodeType CONTACT::FriNodeType::instance_;


DRT::ParObject* CONTACT::FriNodeType::Create( const std::vector<char> & data )
{
  double x[3];
  std::vector<int> dofs(0);

  // TODO: friplus = true for all nodes!!! change this with pack/unpack
  CONTACT::FriNode* node = new CONTACT::FriNode(0,x,0,0,dofs,false,false,true);
  node->Unpack(data);


  return node;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer::FriNodeDataContainer():
slip_(false),
slipold_(false)
{
  for (int i=0;i<3;++i)
  {
    jump()[i]=0.0;
    traction()[i]=0.0;
    tractionold()[i]=0.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::Pack(DRT::PackBuffer& data) const
{
  // add jump_
  DRT::ParObject::AddtoPack(data,jump_,3*sizeof(double));
  // add slip_
  DRT::ParObject::AddtoPack(data,slip_);
  // add slip_
  DRT::ParObject::AddtoPack(data,slipold_);
  // add traction_
  DRT::ParObject::AddtoPack(data,traction_,3*sizeof(double));
  // add tractionold_
  DRT::ParObject::AddtoPack(data,tractionold_,3*sizeof(double));

  // add drowsold_,mrowsold_,mnodesold_
  int hasdata = drowsold_.size();
  // check the sizes of vector/sets
  if(hasdata != (int) mrowsold_.size() or (hasdata == 0 and mnodesold_.size()!=0))
    dserror("Something wrong with sizes of vector/sets!");

  DRT::ParObject::AddtoPack(data,hasdata);

  if(hasdata!=0)
  {
    for (int i=0;i<hasdata;i++)
    {
      DRT::ParObject::AddtoPack(data,(drowsold_[i]));
      DRT::ParObject::AddtoPack(data,(mrowsold_[i]));
    }
    DRT::ParObject::AddtoPack(data,mnodesold_);
  }

  // add derivjump
  int hasdataderivjump = derivjump_.size();
  DRT::ParObject::AddtoPack(data,hasdataderivjump);

  if (hasdataderivjump!=0)
  {
    for (int i=0;i<hasdataderivjump;i++)
      DRT::ParObject::AddtoPack(data,(derivjump_[i]));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::Unpack(std::vector<char>::size_type& position,
                                           const std::vector<char>& data)
{
  // jump_
  DRT::ParObject::ExtractfromPack(position,data,jump_,3*sizeof(double));
  // slip_
  slip_ = DRT::ParObject::ExtractInt(position,data);
  // slipold_
  slipold_ = DRT::ParObject::ExtractInt(position,data);
  // traction_
  DRT::ParObject::ExtractfromPack(position,data,traction_,3*sizeof(double));
  // tractionold_
  DRT::ParObject::ExtractfromPack(position,data,tractionold_,3*sizeof(double));

  //drowsold_,mrowsold_,mnodesold_
  int hasdata;
  DRT::ParObject::ExtractfromPack(position,data,hasdata);

  if(hasdata!=0)
  {
    drowsold_.resize(hasdata);
    mrowsold_.resize(hasdata);
    for (int i=0;i<hasdata;i++)
    {
      DRT::ParObject::ExtractfromPack(position,data,drowsold_[i]);
      DRT::ParObject::ExtractfromPack(position,data,mrowsold_[i]);
    }
    DRT::ParObject::ExtractfromPack(position,data,mnodesold_);
  }

  //and derivjump_
  int hasdataderivjump;
  DRT::ParObject::ExtractfromPack(position,data,hasdataderivjump);

  if(hasdataderivjump!=0)
  {
    derivjump_.resize(hasdataderivjump);
    for (int i=0;i<hasdataderivjump;i++)
    {
      DRT::ParObject::ExtractfromPack(position,data,derivjump_[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 07/11|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainerPlus::FriNodeDataContainerPlus():
wear_(0.0),
deltawear_(0.0)
{
  wcurr_[0] = 0.0;
  wold_[0]  = 0.0;

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainerPlus::Pack(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data,wear_);
  DRT::ParObject::AddtoPack(data,deltawear_);

  // add d2row
  int hasdata = d2rows_.size();
  DRT::ParObject::AddtoPack(data,hasdata);

  if(hasdata!=0)
  {
    for (int i=0;i<hasdata;i++)
    {
      DRT::ParObject::AddtoPack(data,(d2rows_[i]));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainerPlus::Unpack(std::vector<char>::size_type& position,
                                               const std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(position,data,wear_);
  DRT::ParObject::ExtractfromPack(position,data,deltawear_);

  //d2rows_
  int hasdata;
  DRT::ParObject::ExtractfromPack(position,data,hasdata);

  if(hasdata!=0)
  {
    d2rows_.resize(hasdata);
    for (int i=0;i<hasdata;i++)
    {
      DRT::ParObject::ExtractfromPack(position,data,d2rows_[i]);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(int id, const double* coords, const int owner,
                          const int numdof, const std::vector<int>& dofs,
                          const bool isslave, const bool initactive,
                          const bool friplus) :
CONTACT::CoNode(id,coords,owner,numdof,dofs,isslave,initactive),
mechdiss_(0.0),
friplus_(friplus)
{

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(const CONTACT::FriNode& old) :
CONTACT::CoNode(old),
mechdiss_(old.mechdiss_)
{
  // not yet used and thus not necessarily consistent
  dserror("ERROR: FriNode copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of FriNode and return pointer to it (public)|
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode* CONTACT::FriNode::Clone() const
{
  CONTACT::FriNode* newnode = new CONTACT::FriNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                               mgit 02/10|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const CONTACT::FriNode& frinode)
{
  frinode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  CONTACT::CoNode::Print(os);
  if (IsSlave())
    if (IsInitActive()) os << " InitActive ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class MORTAR::MortarNode
  CONTACT::CoNode::Pack(data);

  // add data_
  bool hasdata = (fridata_!=Teuchos::null);
  AddtoPack(data,hasdata);
  if (hasdata) fridata_->Pack(data);

  bool hasdataplus = (fridataplus_!=Teuchos::null);
  AddtoPack(data,hasdataplus);
  if (hasdataplus) fridataplus_->Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class MORTAR::MortarNode
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  CONTACT::CoNode::Unpack(basedata);

  // **************************
  // FriData
  bool hasdata = ExtractInt(position,data);
  if (hasdata)
  {
    fridata_ = Teuchos::rcp(new CONTACT::FriNodeDataContainer());
    fridata_->Unpack(position,data);
  }
  else
    fridata_ = Teuchos::null;

  // **************************
  // FriDataPlus
  bool hasdataplus = ExtractInt(position,data);
  if (hasdataplus)
  {
    fridataplus_ = Teuchos::rcp(new CONTACT::FriNodeDataContainerPlus());
    fridataplus_->Unpack(position,data);
  }
  else
    fridataplus_ = Teuchos::null;

  // Check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'SNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddSNode(int node)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddSnode: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddSNode: function called for boundary node %i", Id());

  FriData().GetSNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'MNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddMNode(int node)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddMNode: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddMNode: function called for boundary node %i", Id());

  FriData().GetMNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D2' map                              farah 06/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddD2Value(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==true)
    dserror("ERROR: AddD2Value: function called for slave node %i", Id());

  // check if this has been called before
  if ((int)FriDataPlus().GetD2().size()==0)
    FriDataPlus().GetD2().resize(NumDof());

  // check row index input
  if ((int)FriDataPlus().GetD2().size()<=row)
    dserror("ERROR: AddD2Value: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& d2map = FriDataPlus().GetD2()[row];
  d2map[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'ANodes' set                       gitterle 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddANode(int node)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddANode: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddANode: function called for boundary node %i", Id());

  FriDataPlus().GetANodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'BNodes' set                        gitterle 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddBNode(int node)
{
  // check if this is a slave node
  if (IsSlave()==true)
    dserror("ERROR: AddBNode: function called for slave node %i", Id());

  GetBNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivJump' map                     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddDerivJumpValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddJumpValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddJumpValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)FriData().GetDerivJump().size()==0)
    FriData().GetDerivJump().resize(NumDof());

  // check row index input
  if ((int)FriData().GetDerivJump().size() <= row)
    dserror("ERROR: AddDerivJumpValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& zmap = FriData().GetDerivJump()[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add jump value from gp-wise integration of slip          farah 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddJumpValue(double val, int k)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddJumpValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddJumpValue: function called for boundary node %i", Id());

  FriData().jump_var()[k] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'A' map                             gitterle 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddAValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddAValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddAValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)FriDataPlus().GetA().size()==0)
    FriDataPlus().GetA().resize(NumDof());

  // check row index input
  if ((int)FriDataPlus().GetA().size()<=row)
    dserror("ERROR: AddAValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& amap = FriDataPlus().GetA()[row];
  amap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'B' map                             gitterle 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddBValue(int& row, int& col, double& val)
{

  // check if this has been called before
  if ((int) GetB().size()==0)
    GetB().resize(NumDof());

  // check row index input
  if ((int) GetB().size()<=row)
    dserror("ERROR: AddBValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& bmap = GetB()[row];
  bmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to mechanical dissipation                  gitterle 08/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddMechDissValue(double& val)
{
  // add given value to mechdiss_
  MechDiss()+=val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'T' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddTValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
//  if (IsSlave()==false)
//    dserror("ERROR: AddTValue: function called for master node %i", Id());

  // check if this has been called before
  if ((int)FriDataPlus().GetT().size()==0)
    FriDataPlus().GetT().resize(NumDof());

  // check row index input
  if ((int)FriDataPlus().GetT().size()<=row)
    dserror("ERROR: AddTValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& tmap = FriDataPlus().GetT()[row];
  tmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'E' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddEValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
//  if (IsSlave()==false)
//    dserror("ERROR: AddEValue: function called for master node %i", Id());

  // check if this has been called before
  if ((int)FriDataPlus().GetE().size()==0)
    FriDataPlus().GetE().resize(NumDof());

  // check row index input
  if ((int)FriDataPlus().GetE().size()<=row)
    dserror("ERROR: AddEValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int,double>& emap = FriDataPlus().GetE()[row];
  emap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries of D and M to old ones             gitterle 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreDMOld()
{
  // copy drows_ to drowsold_

  // reset old nodal Mortar maps
  for (int j=0;j<(int)(FriData().GetDOld().size());++j)
  (FriData().GetDOld())[j].clear();
  for (int j=0;j<(int)((FriData().GetMOld()).size());++j)
  (FriData().GetMOld())[j].clear();

  // clear and zero nodal vectors
  FriData().GetDOld().clear();
  FriData().GetMOld().clear();
  FriData().GetDOld().resize(0);
  FriData().GetMOld().resize(0);

  // write drows_ to drowsold_
  FriData().GetDOld() = MoData().GetD();
  FriData().GetMOld() = MoData().GetM();

  // also vectors containing the according master nodes
  FriData().GetMNodesOld().clear();
  FriData().GetMNodesOld() = FriData().GetMNodes();

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries penalty tractions to old ones      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreTracOld()
{
  // write entries to old ones
  for (int j=0;j<3;++j)
    FriData().tractionold()[j]=FriData().traction()[j];

  return;
}

/*-----------------------------------------------------------------------*
 |  Set the value of deltawear                             gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddDeltaWearValue(double& val)
{
  // add given value to deltawear_
  FriDataPlus().DeltaWear()+=val;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::InitializeDataContainer()
{
  // only initialize if not yet done
  if (modata_==Teuchos::null && codata_==Teuchos::null && fridata_==Teuchos::null)
  {
    modata_  = Teuchos::rcp(new MORTAR::MortarNodeDataContainer());
    codata_  = Teuchos::rcp(new CONTACT::CoNodeDataContainer());
    fridata_ = Teuchos::rcp(new CONTACT::FriNodeDataContainer());
  }
  
  // initialize data container for wear and tsi problems 
  if (friplus_==true)
  {
     if (fridataplus_==Teuchos::null)
      fridataplus_=Teuchos::rcp(new CONTACT::FriNodeDataContainerPlus());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::ResetDataContainer()
{
  // reset to Teuchos::null
  fridata_     = Teuchos::null;
  fridataplus_ = Teuchos::null; 
  codata_      = Teuchos::null;
  modata_      = Teuchos::null;

  return;
}

