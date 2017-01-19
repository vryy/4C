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


DRT::ParObject* DRT::FIBER::FiberNodeType::Create( const std::vector<char> & data )
{
  double dummycoord[3] = {999.,999.,999.};
  double dummyfiber[3]   = {1.,0.,0.};
  DRT::FIBER::FiberNode* object = new DRT::FIBER::FiberNode(-1,dummycoord,dummyfiber,-1);
  object->Unpack(data);
  return object;
}

/*
  Standard ctor
 */
DRT::FIBER::FiberNode::FiberNode(int           id    ,
               const double* coords,
               const double* fiber,
               const int     owner)
:
  DRT::Node(id,coords,owner),
  fiber_(3)
{
  double fibernorm = sqrt(fiber[0]*fiber[0]+fiber[1]*fiber[1]+fiber[2]*fiber[2]);
  if (fibernorm < 1e-13)
    fibernorm = 1.0;
  for (unsigned i=0; i<3; ++i)
    fiber_[i] = fiber[i]/fibernorm;
  return;
}

/*
  Copy Constructor

  Makes a deep copy of a fiber node

*/
DRT::FIBER::FiberNode::FiberNode(const DRT::FIBER::FiberNode& old)
  :
  DRT::Node(old),
  fiber_(old.fiber_)
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
DRT::FIBER::FiberNode::~FiberNode()
{
  return;
}

/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this fiber node

*/
void DRT::FIBER::FiberNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  DRT::Node::AddtoPack(data,type);
  // add base class of fiber node
  DRT::Node::Pack(data);
  // add fiber
  DRT::Node::AddtoPack(data,fiber_);

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
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Node
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Node::Unpack(basedata);
  // extract fiber
  DRT::Node::ExtractfromPack(position,data,fiber_);

  return;
}

/*
  Print this fiber node
*/
void DRT::FIBER::FiberNode::Print(std::ostream& os) const
{
  os << "Fiber Node :";
  DRT::Node::Print(os);
  os << "\n+ additional fiber information"
  << std::setw(12) << Fiber()[0] << " "
  << std::setw(12) << Fiber()[1] << " "
  << std::setw(12) << Fiber()[2] << " ";

  return;

}

