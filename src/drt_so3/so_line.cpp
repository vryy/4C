/*----------------------------------------------------------------------*/
/*! \file
\brief line element

\level 1

*----------------------------------------------------------------------*/

#include "so_line.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_dserror.H"


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
    case line2:
      gaussrule_ = DRT::UTILS::intrule_line_2point;
      break;
    case line3:
      gaussrule_ = DRT::UTILS::intrule_line_3point;
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
  DRT::ELEMENTS::StructuralLine* newelement = new DRT::ELEMENTS::StructuralLine(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::StructuralLine::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    case 3:
      return line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Pack(DRT::PackBuffer& data) const
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
