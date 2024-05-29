/*! \file
 *
\brief Line element associated to the shell 7-Parameter element

\level 3
*/

#include "4C_shell7p_line.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::Shell7pLineType DRT::ELEMENTS::Shell7pLineType::instance_;

DRT::ELEMENTS::Shell7pLineType& DRT::ELEMENTS::Shell7pLineType::Instance() { return instance_; }

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Shell7pLineType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

DRT::ELEMENTS::Shell7pLine::Shell7pLine(int id, int owner, int nnode, const int* nodeids,
    CORE::Nodes::Node** nodes, CORE::Elements::Element* parent, const int lline)
    : CORE::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lline);
  // type of gaussian integration
  switch (Shape())
  {
    case CORE::FE::CellType::line2:
      gaussrule_ = CORE::FE::GaussRule1D::line_2point;
      break;
    case CORE::FE::CellType::line3:
      gaussrule_ = CORE::FE::GaussRule1D::line_3point;
      break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
  }
}

DRT::ELEMENTS::Shell7pLine::Shell7pLine(const DRT::ELEMENTS::Shell7pLine& old)
    : CORE::Elements::FaceElement(old), gaussrule_(old.gaussrule_)
{
}

CORE::Elements::Element* DRT::ELEMENTS::Shell7pLine::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Shell7pLine(*this);
  return newelement;
}

CORE::FE::CellType DRT::ELEMENTS::Shell7pLine::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

void DRT::ELEMENTS::Shell7pLine::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("this Shell7line element does not support communication");
}


void DRT::ELEMENTS::Shell7pLine::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Shell line element does not support communication");
}


void DRT::ELEMENTS::Shell7pLine::Print(std::ostream& os) const
{
  os << "Shell7pLine ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
