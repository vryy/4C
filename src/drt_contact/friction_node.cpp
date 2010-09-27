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
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "contact_node.H"
#include "friction_node.H"
#include "../drt_lib/drt_dserror.H"
#include "contact_element.H"
#include "contact_defines.H"

CONTACT::FriNodeType CONTACT::FriNodeType::instance_;


DRT::ParObject* CONTACT::FriNodeType::Create( const std::vector<char> & data )
{
  double x[3];
  std::vector<int> dofs(0);
  CONTACT::FriNode* node = new CONTACT::FriNode(0,x,0,0,dofs,false,false);
  node->Unpack(data);
  return node;
}

/*----------------------------------------------------------------------*/
// METHODS RELATED TO FRINODEDATACONTAINER
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer::FriNodeDataContainer():
activeold_(false),
slip_(false)
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
 |  copy-ctor (public)                                        mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer::FriNodeDataContainer(const CONTACT::FriNodeDataContainer& old):
activeold_(old.activeold_),
slip_(old.slip_),
drowsold_(old.drowsold_),
mrowsold_(old.mrowsold_),
snodes_(old.snodes_),
mnodes_(old.mnodes_),
mnodesold_(old.mnodesold_)
{
  for (int i=0;i<3;++i)
  {
    jump()[i]=old.jump_[i];
    traction()[i]=old.traction_[i];
    tractionold()[i]=old.tractionold_[i];
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of FriNodeDataContainer and
    return pointer to it (public)                              mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer* CONTACT::FriNodeDataContainer::Clone() const
{
  CONTACT::FriNodeDataContainer* newnodedc = new CONTACT::FriNodeDataContainer(*this);
  return newnodedc;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::Pack(vector<char>& data) const
{
  // add jump_
  DRT::ParObject::AddtoPack(data,jump_,3*sizeof(double));
  // add activeold_
  DRT::ParObject::AddtoPack(data,activeold_);
  // add slip_
  DRT::ParObject::AddtoPack(data,slip_);
  // add traction_
  DRT::ParObject::AddtoPack(data,traction_,3*sizeof(double));
  // add tractionold_
  DRT::ParObject::AddtoPack(data,tractionold_,3*sizeof(double));
  // add drowsold_,mrowsold_,mnodesold_ and derivjump_
  int hasdata = drowsold_.size()!=0; 
  
  // check if the sizes of all vectors/sets to pack are nonzero
  if(hasdata && (mrowsold_.size()==0 or mnodesold_.size()==0 or derivjump_.size()==0))
    dserror("The sizes of all vectors/sets have to be nonzero"); 
  
  DRT::ParObject::AddtoPack(data,hasdata);
  if(hasdata)
  {    
    for (int i=0;i<3;i++)
    {
      DRT::ParObject::AddtoPack(data,(drowsold_[i]));
      DRT::ParObject::AddtoPack(data,(mrowsold_[i]));
      DRT::ParObject::AddtoPack(data,(derivjump_[i]));
    }
    DRT::ParObject::AddtoPack(data,mnodesold_);
  }  
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::Unpack(vector<char>::size_type& position, const vector<char>& data)
{
  // jump_
  DRT::ParObject::ExtractfromPack(position,data,jump_,3*sizeof(double));
  // activeold_
  DRT::ParObject::ExtractfromPack(position,data,activeold_);
  // slip_
  DRT::ParObject::ExtractfromPack(position,data,slip_);
  // traction_
  DRT::ParObject::ExtractfromPack(position,data,traction_,3*sizeof(double));
  // tractionold_
  DRT::ParObject::ExtractfromPack(position,data,tractionold_,3*sizeof(double));
  //drowsold_,mrowsold_,mnodesold_ and derivjump_
  int hasdata;
  DRT::ParObject::ExtractfromPack(position,data,hasdata);
  
  if(hasdata)
  {
    drowsold_.resize(3);
    mrowsold_.resize(3);
    derivjump_.resize(3);
    for (int i=0;i<3;i++)
    {
      DRT::ParObject::ExtractfromPack(position,data,drowsold_[i]);
      DRT::ParObject::ExtractfromPack(position,data,mrowsold_[i]);
      DRT::ParObject::ExtractfromPack(position,data,derivjump_[i]);
    }
    DRT::ParObject::ExtractfromPack(position,data,mnodesold_);
  }

  return;
}

/*----------------------------------------------------------------------*/
// METHODS RELATED TO FRINODE
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  ctor (public)                                              mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(int id, const double* coords, const int owner,
                          const int numdof, const vector<int>& dofs, const bool isslave,
                          const bool initactive) :
CONTACT::CoNode(id,coords,owner,numdof,dofs,isslave,initactive),
mechdiss_(0.0)
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
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of FriNode and return pointer to it (public) |
 |                                                           mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode* CONTACT::FriNode::Clone() const
{
  CONTACT::FriNode* newnode = new CONTACT::FriNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mgit 02/10|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::FriNode& frinode)
{
  frinode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  CONTACT::CoNode::Print(os);
  if (IsInitActive()) os << " InitActive ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class MORTAR::MortarNode
  vector<char> basedata(0);
  CONTACT::CoNode::Pack(basedata);
  AddtoPack(data,basedata);

  // data_
  int hasdata = fridata_!=Teuchos::null;
  AddtoPack(data,hasdata);
  if (hasdata)
    fridata_->Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class MORTAR::MortarNode
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  CONTACT::CoNode::Unpack(basedata);

  // data_
  int hasdata;
  ExtractfromPack(position,data,hasdata);
  if (hasdata)
  {
    fridata_ = Teuchos::rcp(new CONTACT::FriNodeDataContainer());
    fridata_->Unpack(position,data);
  }
  else
  {
    fridata_ = Teuchos::null;
  }

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
  map<int,double>& zmap = FriData().GetDerivJump()[row];
  zmap[col] += val;

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

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::InitializeDataContainer()
{
	// only initialize if not yet done
	if (modata_==Teuchos::null && codata_==Teuchos::null && fridata_==Teuchos::null)
	{
		modata_=rcp(new MORTAR::MortarNodeDataContainer());
		codata_ =rcp(new CONTACT::CoNodeDataContainer());
		fridata_=rcp(new CONTACT::FriNodeDataContainer());
	}

  return;
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::ResetDataContainer()
{
	// reset to null
  fridata_ = Teuchos::null;
  codata_  = Teuchos::null;
  modata_  = Teuchos::null;

  return;
}

#endif  // #ifdef CCADISCRET
