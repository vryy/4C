/*----------------------------------------------------------------------*/
/*! \file
\brief line element

\level 1

*----------------------------------------------------------------------*/

#include "baci_so3_line.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::StructuralLineType DRT::ELEMENTS::StructuralLineType::instance_;

DRT::ELEMENTS::StructuralLineType& DRT::ELEMENTS::StructuralLineType::Instance()
{
  return instance_;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::StructuralLineType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new StructuralLine( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralLine::StructuralLine(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::Element* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
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
      dserror("shape type unknown!\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralLine::StructuralLine(const DRT::ELEMENTS::StructuralLine& old)
    : DRT::FaceElement(old), gaussrule_(old.gaussrule_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::StructuralLine::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::StructuralLine(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::StructuralLine::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Pack(CORE::COMM::PackBuffer& data) const
{
  dserror("StructuralLine element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Unpack(const std::vector<char>& data)
{
  dserror("StructuralLine element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Print(std::ostream& os) const
{
  os << "StructuralLine ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
