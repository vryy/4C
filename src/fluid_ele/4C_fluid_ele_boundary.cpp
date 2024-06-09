/*----------------------------------------------------------------------*/
/*! \file

\brief fluid boundary element

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidBoundaryType Discret::ELEMENTS::FluidBoundaryType::instance_;

Discret::ELEMENTS::FluidBoundaryType& Discret::ELEMENTS::FluidBoundaryType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::FluidBoundaryType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::FluidBoundary* object = new Discret::ELEMENTS::FluidBoundary(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidBoundary::FluidBoundary(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Fluid* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1)
{
  set_parent_master_element(parent, lsurface);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  distype_ = Core::FE::getShapeOfBoundaryElement(num_node(), ParentMasterElement()->Shape());

  numdofpernode_ = ParentMasterElement()->NumDofPerNode(*Nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < num_node(); ++nlid)
  {
    if (numdofpernode_ != ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]))
      FOUR_C_THROW(
          "You need different NumDofPerNode for each node on this fluid boundary? (%d != %d)",
          numdofpernode_, ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]));
  }
  return;
}

/*------------------------------------------------------------------------*
 |  ctor (private) - used by FluidBoundaryType                  ager 12/16|
 *-----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidBoundary::FluidBoundary(int id, int owner)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidBoundary::FluidBoundary(const Discret::ELEMENTS::FluidBoundary& old)
    : Core::Elements::FaceElement(old), distype_(old.distype_), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::FluidBoundary::Clone() const
{
  Discret::ELEMENTS::FluidBoundary* newelement = new Discret::ELEMENTS::FluidBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           ager 12/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidBoundary::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  FaceElement::Pack(data);
  // Discretisation type
  add_to_pack(data, distype_);
  // add numdofpernode_
  add_to_pack(data, numdofpernode_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           ager 12/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidBoundary::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  FaceElement::Unpack(basedata);
  // distype
  distype_ = static_cast<Core::FE::CellType>(ExtractInt(position, data));
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidBoundary::Print(std::ostream& os) const
{
  os << "FluidBoundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidBoundary::Lines()
{
  FOUR_C_THROW("Lines of FluidBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          ager 12/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidBoundary::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
