/*---------------------------------------------------------------------*/
/*! \file

\brief A virtual class for a node

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_node.H"
#include "drt_dserror.H"


DRT::NodeType DRT::NodeType::instance_;


/*----------------------------------------------------------------------*
 |  kind of ctor (public)                                    mwgee 06/13|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::NodeType::Create(const std::vector<char>& data)
{
  double dummycoord[6] = {999., 999., 999., 999., 999., 999.};
  DRT::Node* object = new DRT::Node(-1, dummycoord, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(int id, const double* coords, const int owner, const bool iscosserat)
    : ParObject(), id_(id), lid_(-1), owner_(owner)
{
  if (!iscosserat)
    x_.resize(3);
  else
    x_.resize(6);
  for (unsigned i = 0; i < x_.size(); ++i) x_[i] = coords[i];
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(const DRT::Node& old)
    : ParObject(old),
      id_(old.id_),
      lid_(old.lid_),
      owner_(old.owner_),
      x_(old.x_),
      element_(old.element_)
{
  // we do NOT want a deep copy of the condition_ a condition is
  // only a reference in the node anyway
  std::map<std::string, Teuchos::RCP<Condition>>::const_iterator fool;
  for (fool = old.condition_.begin(); fool != old.condition_.end(); ++fool)
    SetCondition(fool->first, fool->second);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::~Node() { return; }


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
std::ostream& operator<<(std::ostream& os, const DRT::Node& node)
{
  node.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Node::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Node " << std::setw(12) << Id() << " Owner " << std::setw(4) << Owner() << " Coords "
     << std::setw(12) << X()[0] << " " << std::setw(12) << X()[1] << " " << std::setw(12) << X()[2]
     << " ";

#if 0
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << std::endl << numcond << " Conditions:\n";
    std::map<std::string,Teuchos::RCP<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << curr->first << " ";
      os << *(curr->second) << std::endl;
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
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add id
  int id = Id();
  AddtoPack(data, id);
  // add owner
  int owner = Owner();
  AddtoPack(data, owner);
  // x_
  //  AddtoPack(data,x_,3*sizeof(double));
  AddtoPack(data, x_);

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
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // id_
  ExtractfromPack(position, data, id_);
  // owner_
  ExtractfromPack(position, data, owner_);
  // x_
  //  ExtractfromPack(position,data,x_,3*sizeof(double));
  ExtractfromPack(position, data, x_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Node::GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  std::multimap<std::string, Teuchos::RCP<Condition>>::const_iterator startit =
      condition_.lower_bound(name);
  std::multimap<std::string, Teuchos::RCP<Condition>>::const_iterator endit =
      condition_.upper_bound(name);
  int count = 0;
  std::multimap<std::string, Teuchos::RCP<Condition>>::const_iterator curr;
  for (curr = startit; curr != endit; ++curr) out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Node::GetCondition(const std::string& name) const
{
  std::multimap<std::string, Teuchos::RCP<Condition>>::const_iterator curr = condition_.find(name);
  if (curr == condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}

/*----------------------------------------------------------------------*
 |  Change position reference                                  (public) |
 |                                                            mc  06/11 |
 *----------------------------------------------------------------------*/
void DRT::Node::ChangePos(std::vector<double> nvector)
{
  for (int i = 0; i < 3; ++i) x_[i] = x_[i] + nvector[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Set reference position                                     (public) |
 |                                                            jb  07/14 |
 *----------------------------------------------------------------------*/
void DRT::Node::SetPos(std::vector<double> nvector)
{
  for (int i = 0; i < 3; ++i) x_[i] = nvector[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Query data to be visualized by BINIO                       (public) |
 |                                                        mayr.mt 11/15 |
 *----------------------------------------------------------------------*/
bool DRT::Node::VisData(const std::string& name, std::vector<double>& data)
{
  if (name == "Nodeowner")
  {
    if ((int)data.size() < 1) dserror("Size mismatch");
    data[0] = Owner();
    return true;
  }
  return false;
}
