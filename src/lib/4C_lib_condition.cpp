/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#include "4C_lib_condition.hpp"

#include "4C_lib_element.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


DRT::Condition::Condition(const int id, const CORE::Conditions::ConditionType type,
    const std::string name, const bool buildgeometry, const CORE::Conditions::GeometryType gtype)
    : id_(id), buildgeometry_(buildgeometry), type_(type), name_(std::move(name)), gtype_(gtype)
{
}

std::ostream& operator<<(std::ostream& os, const DRT::Condition& cond)
{
  cond.Print(os);
  return os;
}

std::string DRT::Condition::Name() const
{
  std::stringstream condition_name;
  condition_name << "Condition " << id_ << " " << name_ << ": ";

  return condition_name.str();
}

void DRT::Condition::Print(std::ostream& os) const
{
  os << DRT::Condition::Name();
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

void DRT::Condition::AdjustId(const int shift)
{
  std::map<int, Teuchos::RCP<DRT::Element>> geometry;
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter;

  for (const auto& [ele_id, ele] : *geometry_)
  {
    ele->SetId(ele_id + shift);
    geometry[ele_id + shift] = (*geometry_)[ele_id];
  }

  swap(*geometry_, geometry);
}

Teuchos::RCP<DRT::Condition> DRT::Condition::copy_without_geometry() const
{
  Teuchos::RCP<Condition> copy(new Condition(*this));
  copy->ClearGeometry();
  return copy;
}


FOUR_C_NAMESPACE_CLOSE
