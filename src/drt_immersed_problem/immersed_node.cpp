/*----------------------------------------------------------------------*/
/*! \file

\brief specialized node for immersed problems.

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/

#include "immersed_node.H"


IMMERSED::ImmersedNodeType IMMERSED::ImmersedNodeType::instance_;


/*----------------------------------------------------------------------*
 |  kind of ctor (public)                                   rauch 11/14 |
 *----------------------------------------------------------------------*/
DRT::ParObject* IMMERSED::ImmersedNodeType::Create(const std::vector<char>& data)
{
  double dummycoord[6] = {999., 999., 999., 999., 999., 999.};
  DRT::Node* object = new IMMERSED::ImmersedNode(-1, dummycoord, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           rauch 11/14 |
 *----------------------------------------------------------------------*/
IMMERSED::ImmersedNode::ImmersedNode(int id, const double* coords, const int owner)
    : DRT::Node(id, coords, owner),
      ismatched_(false),
      IsBoundaryImmersed_(false),
      IsPseudoBoundary_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      rauch 11/14 |
 *----------------------------------------------------------------------*/
IMMERSED::ImmersedNode::ImmersedNode(const IMMERSED::ImmersedNode& old)
    : DRT::Node(old),
      ismatched_(old.ismatched_),
      IsBoundaryImmersed_(old.IsBoundaryImmersed_),
      IsPseudoBoundary_(old.IsPseudoBoundary_)
{
  dserror("ERROR: ImmersedNode copy-ctor has not been used before. Be careful when using it.");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                           rauch 11/14 |
 *----------------------------------------------------------------------*/
IMMERSED::ImmersedNode::~ImmersedNode() { return; }


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
IMMERSED::ImmersedNode* IMMERSED::ImmersedNode::Clone() const
{
  IMMERSED::ImmersedNode* newnode = new IMMERSED::ImmersedNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                             rauch 11/14 |
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const IMMERSED::ImmersedNode& immersednode)
{
  immersednode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             rauch 11/14 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedNode::Print(std::ostream& os) const
{
  os << "Immersed ";
  DRT::Node::Print(os);

  if (IsBoundaryImmersed())
    os << " Immersed Boundary  ";
  else
    os << " NOT Immersed Boundary ";

  if (IsMatched())
    os << " Matched ";
  else
    os << " NOT Matched ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class DRT::Node
  DRT::Node::Pack(data);

  // add IsBoundaryImmersed_
  AddtoPack(data, IsBoundaryImmersed_);
  // add ismatched_
  AddtoPack(data, ismatched_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data type %d != ObjId %d ", type, UniqueParObjectId());

  // extract base class DRT::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Node::Unpack(basedata);

  // isimersedboundary_
  IsBoundaryImmersed_ = ExtractInt(position, data);
  // ismatched_
  ismatched_ = ExtractInt(position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  Visualization Data                                         (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedNode::VisNames(std::map<std::string, int>& names)
{
  names.insert(std::pair<std::string, int>("IsBoundaryImmersedNode", 1));
  return;
}


/*----------------------------------------------------------------------*
 |  Query data to be visualized by BINIO                       (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
bool IMMERSED::ImmersedNode::VisData(const std::string& name, std::vector<double>& data)
{
  if (name == "IsBoundaryImmersedNode")
  {
    if ((int)data.size() < 1) dserror("Size mismatch");
    data[0] = IsBoundaryImmersed();
    return true;
  }
  return false;
}
