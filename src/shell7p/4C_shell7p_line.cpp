/*! \file
 *
\brief Line element associated to the shell 7-Parameter element

\level 3
*/

#include "4C_shell7p_line.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Shell7pLineType Discret::ELEMENTS::Shell7pLineType::instance_;

Discret::ELEMENTS::Shell7pLineType& Discret::ELEMENTS::Shell7pLineType::instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Shell7pLineType::create(
    const int id, const int owner)
{
  return Teuchos::null;
}

Discret::ELEMENTS::Shell7pLine::Shell7pLine(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  // type of gaussian integration
  switch (shape())
  {
    case Core::FE::CellType::line2:
      gaussrule_ = Core::FE::GaussRule1D::line_2point;
      break;
    case Core::FE::CellType::line3:
      gaussrule_ = Core::FE::GaussRule1D::line_3point;
      break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
  }
}

Discret::ELEMENTS::Shell7pLine::Shell7pLine(const Discret::ELEMENTS::Shell7pLine& old)
    : Core::Elements::FaceElement(old), gaussrule_(old.gaussrule_)
{
}

Core::Elements::Element* Discret::ELEMENTS::Shell7pLine::clone() const
{
  auto* newelement = new Discret::ELEMENTS::Shell7pLine(*this);
  return newelement;
}

Core::FE::CellType Discret::ELEMENTS::Shell7pLine::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

void Discret::ELEMENTS::Shell7pLine::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Shell7line element does not support communication");
}


void Discret::ELEMENTS::Shell7pLine::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Shell line element does not support communication");
}


void Discret::ELEMENTS::Shell7pLine::print(std::ostream& os) const
{
  os << "Shell7pLine ";
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
