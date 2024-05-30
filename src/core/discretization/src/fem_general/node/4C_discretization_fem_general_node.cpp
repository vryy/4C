/*---------------------------------------------------------------------*/
/*! \file

\brief A virtual class for a node

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_fem_general_node.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::Nodes::NodeType CORE::Nodes::NodeType::instance_;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* CORE::Nodes::NodeType::Create(const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* object = new CORE::Nodes::Node(-1, dummycoord, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Nodes::Node::Node(const int id, const std::vector<double>& coords, const int owner)
    : ParObject(), id_(id), lid_(-1), owner_(owner), x_(coords)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Nodes::Node::Node(const CORE::Nodes::Node& old)
    : ParObject(old),
      id_(old.id_),
      lid_(old.lid_),
      owner_(old.owner_),
      x_(old.x_),
      element_(old.element_)
{
  // we do NOT want a deep copy of the condition_ a condition is
  // only a reference in the node anyway
  std::map<std::string, Teuchos::RCP<CORE::Conditions::Condition>>::const_iterator fool;
  for (fool = old.condition_.begin(); fool != old.condition_.end(); ++fool)
    SetCondition(fool->first, fool->second);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Nodes::Node* CORE::Nodes::Node::Clone() const
{
  auto* newnode = new CORE::Nodes::Node(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CORE::Nodes::Node& node)
{
  node.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Node " << std::setw(12) << Id() << " Owner " << std::setw(4) << Owner() << " Coords "
     << std::setw(12) << X()[0] << " " << std::setw(12) << X()[1] << " " << std::setw(12) << X()[2]
     << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
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
  AddtoPack(data, x_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // id_
  ExtractfromPack(position, data, id_);
  // owner_
  ExtractfromPack(position, data, owner_);
  // x_
  ExtractfromPack(position, data, x_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::GetCondition(
    const std::string& name, std::vector<CORE::Conditions::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  auto startit = condition_.lower_bound(name);
  auto endit = condition_.upper_bound(name);
  int count = 0;
  std::multimap<std::string, Teuchos::RCP<CORE::Conditions::Condition>>::const_iterator curr;
  for (curr = startit; curr != endit; ++curr) out[count++] = curr->second.get();
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Conditions::Condition* CORE::Nodes::Node::GetCondition(const std::string& name) const
{
  auto curr = condition_.find(name);
  if (curr == condition_.end()) return nullptr;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::ChangePos(std::vector<double> nvector)
{
  FOUR_C_ASSERT(x_.size() == nvector.size(),
      "Mismatch in size of the nodal coordinates vector and the vector to change the nodal "
      "position");
  for (std::size_t i = 0; i < x_.size(); ++i) x_[i] = x_[i] + nvector[i];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Nodes::Node::SetPos(std::vector<double> nvector)
{
  FOUR_C_ASSERT(x_.size() == nvector.size(),
      "Mismatch in size of the nodal coordinates vector and the vector to set the new nodal "
      "position");
  for (std::size_t i = 0; i < x_.size(); ++i) x_[i] = nvector[i];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::Nodes::Node::VisData(const std::string& name, std::vector<double>& data)
{
  if (name == "Nodeowner")
  {
    if (static_cast<int>(data.size()) < 1) FOUR_C_THROW("Size mismatch");
    data[0] = Owner();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
