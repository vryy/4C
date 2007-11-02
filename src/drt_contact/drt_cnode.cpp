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



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(int id, const double* coords, const int owner, 
                      const int numdof, const vector<int>& dofs, const bool isslave) :
DRT::Node(id,coords,owner),
isslave_(isslave),
numdof_(numdof),
dofs_(dofs)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(const CONTACT::CNode& old) :
DRT::Node(old),
isslave_(old.isslave_),
numdof_(old.numdof_),
dofs_(old.dofs_)
{
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
  os << "CNode " << setw(12) << Id()
     << " Owner " << setw(4) << Owner()
     << " Coords "
     << setw(12) << X()[0] << " "
     << setw(12) << X()[1] << " "
     << setw(12) << X()[2] << " "
     << " Dofs "; 
  for (int i=0; i<(int)dofs_.size(); ++i)
    os << dofs_[i] << " ";
  if (IsSlave()) os << " Slave Side  ";
  else           os << " Master Side ";

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
  // add numdof_
  AddtoPack(data,numdof_);
  // add dofs_
  AddtoPack(data,dofs_);

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
  // numdof_
  ExtractfromPack(position,data,numdof_);
  // dofs_
  ExtractfromPack(position,data,dofs_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


#endif  // #ifdef CCADISCRET
