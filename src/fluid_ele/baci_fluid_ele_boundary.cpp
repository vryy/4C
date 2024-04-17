/*----------------------------------------------------------------------*/
/*! \file

\brief fluid boundary element

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::FluidBoundaryType DRT::ELEMENTS::FluidBoundaryType::instance_;

DRT::ELEMENTS::FluidBoundaryType& DRT::ELEMENTS::FluidBoundaryType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::FluidBoundaryType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::FluidBoundary* object = new DRT::ELEMENTS::FluidBoundary(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidBoundaryType::Create(const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundary::FluidBoundary(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Fluid* parent, const int lsurface)
    : DRT::FaceElement(id, owner), distype_(CORE::FE::CellType::dis_none), numdofpernode_(-1)
{
  SetParentMasterElement(parent, lsurface);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  distype_ = CORE::FE::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());

  numdofpernode_ = ParentMasterElement()->NumDofPerNode(*Nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < NumNode(); ++nlid)
  {
    if (numdofpernode_ != ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]))
      dserror("You need different NumDofPerNode for each node on this fluid boundary? (%d != %d)",
          numdofpernode_, ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]));
  }
  return;
}

/*------------------------------------------------------------------------*
 |  ctor (private) - used by FluidBoundaryType                  ager 12/16|
 *-----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundary::FluidBoundary(int id, int owner)
    : DRT::FaceElement(id, owner), distype_(CORE::FE::CellType::dis_none), numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundary::FluidBoundary(const DRT::ELEMENTS::FluidBoundary& old)
    : DRT::FaceElement(old), distype_(old.distype_), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidBoundary::Clone() const
{
  DRT::ELEMENTS::FluidBoundary* newelement = new DRT::ELEMENTS::FluidBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           ager 12/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidBoundary::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  FaceElement::Pack(data);
  // Discretisation type
  AddtoPack(data, distype_);
  // add numdofpernode_
  AddtoPack(data, numdofpernode_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           ager 12/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidBoundary::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  FaceElement::Unpack(basedata);
  // distype
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidBoundary::Print(std::ostream& os) const
{
  os << "FluidBoundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidBoundary::Lines()
{
  dserror("Lines of FluidBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          ager 12/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidBoundary::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
