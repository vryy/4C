/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_condition.hpp"

#include "4C_discretization_fem_general_element.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


Core::Conditions::Condition::Condition(const int id, const Core::Conditions::ConditionType type,
    const bool buildgeometry, const Core::Conditions::GeometryType gtype)
    : id_(id), buildgeometry_(buildgeometry), type_(type), gtype_(gtype)
{
}

std::ostream& operator<<(std::ostream& os, const Core::Conditions::Condition& cond)
{
  cond.Print(os);
  return os;
}


void Core::Conditions::Condition::Print(std::ostream& os) const
{
  os << "Condition " << id_ << " " << to_string(type_) << ": ";
  container_.Print(os);
  os << std::endl;
  if (nodes_.size() != 0)
  {
    os << "Nodes of this condition:";
    for (const auto& node_gid : nodes_) os << " " << node_gid;
    os << std::endl;
  }
  if (geometry_ != Teuchos::null and (int) geometry_->size())
  {
    os << "Elements of this condition:";
    for (const auto& [ele_id, ele] : *geometry_) os << " " << ele_id;
    os << std::endl;
  }
}

void Core::Conditions::Condition::AdjustId(const int shift)
{
  std::map<int, Teuchos::RCP<Core::Elements::Element>> geometry;
  std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator iter;

  for (const auto& [ele_id, ele] : *geometry_)
  {
    ele->SetId(ele_id + shift);
    geometry[ele_id + shift] = (*geometry_)[ele_id];
  }

  swap(*geometry_, geometry);
}

Teuchos::RCP<Core::Conditions::Condition> Core::Conditions::Condition::copy_without_geometry() const
{
  Teuchos::RCP<Core::Conditions::Condition> copy(new Condition(*this));
  copy->clear_geometry();
  return copy;
}


FOUR_C_NAMESPACE_CLOSE
