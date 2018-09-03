/*!----------------------------------------------------------------------
\file drt_fiber_node.cpp

\brief   This is basically a (3d-) node with an additional fiber direction.


\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264

\level 2


*----------------------------------------------------------------------*/

#include "drt_fiber_node.H"

DRT::FIBER::FiberNodeType DRT::FIBER::FiberNodeType::instance_;


DRT::ParObject* DRT::FIBER::FiberNodeType::Create(const std::vector<char>& data)
{
  double dummycoord[3] = {999., 999., 999.};
  double dummyfiber[3] = {1., 0., 0.};
  double dummycosy[6] = {1., 0., 0., 0., 1., 0.};
  double dummyangle[2] = {1., 1.};
  DRT::FIBER::FiberNode* object =
      new DRT::FIBER::FiberNode(-1, dummycoord, dummyfiber, dummycosy, dummyangle, -1);
  object->Unpack(data);
  return object;
}

/*
  Standard ctor
 */
DRT::FIBER::FiberNode::FiberNode(int id, const double* coords, const double* fiber,
    const double* cosy, const double* angle, const int owner)
    : DRT::Node(id, coords, owner), fiber_(3), cir_(3), tan_(3), helix_(0), transverse_(0)
{
  // store fiber information
  double norm = sqrt(fiber[0] * fiber[0] + fiber[1] * fiber[1] + fiber[2] * fiber[2]);
  if (norm < 1e-13) norm = 1.0;
  for (unsigned i = 0; i < 3; ++i) fiber_[i] = fiber[i] / norm;

  // store circumferential direction
  norm = sqrt(cosy[0] * cosy[0] + cosy[1] * cosy[1] + cosy[2] * cosy[2]);
  if (norm < 1e-13) norm = 1.0;
  for (unsigned i = 0; i < 3; ++i) cir_[i] = cosy[i] / norm;

  // store tangential direction
  norm = sqrt(cosy[3] * cosy[3] + cosy[4] * cosy[4] + cosy[5] * cosy[5]);
  if (norm < 1e-13) norm = 1.0;
  for (unsigned i = 0; i < 3; ++i) tan_[i] = cosy[i + 3] / norm;

  // store helix and transverse angle
  helix_ = angle[0];
  transverse_ = angle[1];

  return;
}

/*
  Copy Constructor

  Makes a deep copy of a fiber node

*/
DRT::FIBER::FiberNode::FiberNode(const DRT::FIBER::FiberNode& old)
    : DRT::Node(old),
      fiber_(old.fiber_),
      cir_(old.cir_),
      tan_(old.tan_),
      helix_(old.helix_),
      transverse_(old.transverse_)
{
  return;
}

/*
  Deep copy the derived class and return pointer to it

*/
DRT::FIBER::FiberNode* DRT::FIBER::FiberNode::Clone() const
{
  DRT::FIBER::FiberNode* newfn = new DRT::FIBER::FiberNode(*this);

  return newfn;
}

/*
  Destructor
*/
DRT::FIBER::FiberNode::~FiberNode() { return; }

/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this fiber node

*/
void DRT::FIBER::FiberNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  DRT::Node::AddtoPack(data, type);
  // add base class of fiber node
  DRT::Node::Pack(data);
  // add fiber
  DRT::Node::AddtoPack(data, fiber_);
  // add circumferential direction
  DRT::Node::AddtoPack(data, cir_);
  // add tangential direction
  DRT::Node::AddtoPack(data, tan_);
  // add helix angle
  DRT::Node::AddtoPack(data, helix_);
  // add transverse angle
  DRT::Node::AddtoPack(data, transverse_);

  return;
}

/*
  Unpack data from a char vector into this class

  Pack and Unpack are used to communicate this fiber node
*/
void DRT::FIBER::FiberNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Node::Unpack(basedata);
  // extract fiber
  DRT::Node::ExtractfromPack(position, data, fiber_);
  // extract circumferential direction
  DRT::Node::ExtractfromPack(position, data, cir_);
  // extract tangential direction
  DRT::Node::ExtractfromPack(position, data, tan_);
  // extract helix angle
  DRT::Node::ExtractfromPack(position, data, helix_);
  // extract transversse angle
  DRT::Node::ExtractfromPack(position, data, transverse_);

  return;
}

/*
  Print this fiber node
*/
void DRT::FIBER::FiberNode::Print(std::ostream& os) const
{
  os << "Fiber Node :";
  DRT::Node::Print(os);
  os << "\n+ additional fiber information" << std::setw(12) << Fiber()[0] << " " << std::setw(12)
     << Fiber()[1] << " " << std::setw(12) << Fiber()[2] << " ";

  return;
}
