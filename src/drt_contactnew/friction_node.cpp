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


/*----------------------------------------------------------------------*
 |  ctor (public)                                              mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(int id, const double* coords, const int owner,
                          const int numdof, const vector<int>& dofs, const bool isslave,
                          const bool initactive) :
CONTACT::CoNode(id,coords,owner,numdof,dofs,isslave,initactive),
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
 |  copy-ctor (public)                                        mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(const CONTACT::FriNode& old) :
CONTACT::CoNode(old),
activeold_(old.activeold_),
slip_(old.slip_),
drowsold_(old.drowsold_),
mrowsold_(old.mrowsold_),
drowsPG_(old.drowsPG_),
mrowsPG_(old.mrowsPG_),
drowsoldPG_(old.drowsoldPG_),
mrowsoldPG_(old.mrowsoldPG_),
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
  
  // add jump_
  AddtoPack(data,jump_,3);
  // add activeold_
  AddtoPack(data,activeold_);
  // add slip_
  AddtoPack(data,slip_);
  // add traction_
  AddtoPack(data,traction_,3);
  // add tractionold_
  AddtoPack(data,tractionold_,3);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::Unpack(const vector<char>& data)
{
  int position = 0;
  
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  
  // extract base class MORTAR::MortarNode
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  CONTACT::CoNode::Unpack(basedata);
  
  // jump_
  ExtractfromPack(position,data,jump_,3);
  // activeold_
  ExtractfromPack(position,data,activeold_);
  // slip_
  ExtractfromPack(position,data,slip_);
  // traction_
  ExtractfromPack(position,data,traction_,3);
  // tractionold_
  ExtractfromPack(position,data,tractionold_,3);

  if (position != (int)data.size())
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

  GetSNodes().insert(node);
  
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

  GetMNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D' map (Petrov-Galerkin approach) gitterle 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddDValuePG(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)GetDPG().size()==0)
    GetDPG().resize(NumDof());

  // check row index input
  if ((int)GetDPG().size()<=row)
    dserror("ERROR: AddDValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  map<int,double>& dmap = GetDPG()[row];
  dmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'M' map (Petrov-Galerkin approach) gitterle 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddMValuePG(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave()==false)
    dserror("ERROR: AddDValue: function called for master node %i", Id());
  if (IsOnBound()==true)
    dserror("ERROR: AddDValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)GetMPG().size()==0)
    GetMPG().resize(NumDof());

  // check row index input
  if ((int)GetMPG().size()<=row)
    dserror("ERROR: AddMValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  map<int,double>& mmap = GetMPG()[row];
  mmap[col] += val;

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
  if ((int)GetDerivJump().size()==0)
    GetDerivJump().resize(NumDof());

  // check row index input
  if ((int)GetDerivJump().size() <= row)
    dserror("ERROR: AddDerivJumpValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  map<int,double>& zmap = GetDerivJump()[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries of D and M to old ones             gitterle 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreDMOld()
{
  // copy drows_ to drowsold_

  // reset old nodal Mortar maps
  for (int j=0;j<(int)(GetDOld().size());++j)
  (GetDOld())[j].clear();
  for (int j=0;j<(int)((GetMOld()).size());++j)
  (GetMOld())[j].clear();

  // clear and zero nodal vectors
  GetDOld().clear();
  GetMOld().clear();
  GetDOld().resize(0);
  GetMOld().resize(0);

  // write drows_ to drowsold_
  GetDOld() = GetD();
  GetMOld() = GetM();

  // also vectors containing the according master nodes
  GetMNodesOld().clear();
  GetMNodesOld() = GetMNodes();

  return;
}

/*----------------------------------------------------------------------*
 |  Store entries of D and M to old ones (PG-approach)     gitterle 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreDMOldPG()
{
  // copy drows_ to drowsold_

  // reset old nodal Mortar maps
  for (int j=0;j<(int)(GetDOldPG().size());++j)
  (GetDOldPG())[j].clear();
  for (int j=0;j<(int)((GetMOldPG()).size());++j)
  (GetMOldPG())[j].clear();

  // clear and zero nodal vectors
  GetDOldPG().clear();
  GetMOldPG().clear();
  GetDOldPG().resize(0);
  GetMOldPG().resize(0);

  // write drows_ to drowsold_
  GetDOldPG() = GetDPG();
  GetMOldPG() = GetMPG();

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries penalty tractions to old ones      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreTracOld()
{
  // write entries to old ones
  for (int j=0;j<3;++j)
    tractionold()[j]=traction()[j];

  return;
}

#endif  // #ifdef CCADISCRET
