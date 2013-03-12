/*!----------------------------------------------------------------------
\file drt_node.cpp
\brief A virtual class for a node

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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

#include "drt_node.H"
#include "drt_dserror.H"


DRT::NodeType DRT::NodeType::instance_;


DRT::ParObject* DRT::NodeType::Create( const std::vector<char> & data )
{
  double dummycoord[3] = {999.,999.,999.};
  DRT::Node* object = new DRT::Node(-1,dummycoord,-1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(int id, const double* coords, const int owner) :
ParObject(),
id_(id),
lid_(-1),
owner_(owner)
{
  for (int i=0; i<3; ++i) x_[i] = coords[i];
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(const DRT::Node& old) :
ParObject(old),
id_(old.id_),
lid_(old.lid_),
owner_(old.owner_),
element_(old.element_)
{
  for (int i=0; i<3; ++i) x_[i] = old.x_[i];

  // we do NOT want a deep copy of the condition_ a condition is
  // only a reference in the node anyway
  std::map<std::string,Teuchos::RCP<Condition> >::const_iterator fool;
  for (fool=old.condition_.begin(); fool!=old.condition_.end(); ++fool)
    SetCondition(fool->first,fool->second);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::~Node()
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Node::Clone() const
{
  DRT::Node* newnode = new DRT::Node(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Node& node)
{
  node.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Node::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Node " << std::setw(12) << Id()
     << " Owner " << std::setw(4) << Owner()
     << " Coords "
     << std::setw(12) << X()[0] << " "
     << std::setw(12) << X()[1] << " "
     << std::setw(12) << X()[2] << " ";

#if 0
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << endl << numcond << " Conditions:\n";
    std::map<std::string,RCP<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << curr->first << " ";
      os << *(curr->second) << endl;
    }
  }
#endif
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Node::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add id
  int id = Id();
  AddtoPack(data,id);
  // add owner
  int owner = Owner();
  AddtoPack(data,owner);
  // x_
  AddtoPack(data,x_,3*sizeof(double));

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Node::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // id_
  ExtractfromPack(position,data,id_);
  // owner_
  ExtractfromPack(position,data,owner_);
  // x_
  ExtractfromPack(position,data,x_,3*sizeof(double));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Node::GetCondition(const string& name,std::vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  std::multimap<std::string,Teuchos::RCP<Condition> >::const_iterator startit =
                                         condition_.lower_bound(name);
  std::multimap<std::string,Teuchos::RCP<Condition> >::const_iterator endit =
                                         condition_.upper_bound(name);
  int count=0;
  std::multimap<std::string,Teuchos::RCP<Condition> >::const_iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Node::GetCondition(const string& name) const
{
  std::multimap<std::string,Teuchos::RCP<Condition> >::const_iterator curr =
                                         condition_.find(name);
  if (curr==condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}


/*----------------------------------------------------------------------*
 |  Change position reference                                  (public) |
 |                                                            mc	06/11 |
 *----------------------------------------------------------------------*/
void DRT::Node::ChangePos(std::vector<double> nvector)
{

	for (int i=0; i<3; ++i) x_[i] = x_[i]+ nvector[i];

	return;
}


