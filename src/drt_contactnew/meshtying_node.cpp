/*!----------------------------------------------------------------------
\file meshtying_node.cpp
\brief A class for a meshtying node

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

#include "meshtying_node.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtNode::MtNode(int id, const double* coords, const int owner,
                        const int numdof, const vector<int>& dofs,
                        const bool isslave) :
MORTAR::MortarNode(id,coords,owner,numdof,dofs,isslave)
{
  // empty constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtNode::MtNode(const CONTACT::MtNode& old) :
MORTAR::MortarNode(old)
{
  // empty copy-constructor

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of MtNode and return pointer to it (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtNode* CONTACT::MtNode::Clone() const
{
  CONTACT::MtNode* newnode = new CONTACT::MtNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::MtNode& mtnode)
{
  mtnode.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this node (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Meshtying ";
  MORTAR::MortarNode::Print(os);
  
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtNode::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  
  // add base class MORTAR::MortarNode
  vector<char> basedata(0);
  MORTAR::MortarNode::Pack(basedata);
  AddtoPack(data,basedata);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtNode::Unpack(const vector<char>& data)
{
  int position = 0;
  
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  
  // extract base class MORTAR::MortarNode
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  MORTAR::MortarNode::Unpack(basedata);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

#endif  // #ifdef CCADISCRET
