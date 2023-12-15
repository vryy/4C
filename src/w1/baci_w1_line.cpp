/*----------------------------------------------------------------------------*/
/*! \file
\brief Line definitions for wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_utils_exceptions.H"
#include "baci_w1.H"

BACI_NAMESPACE_OPEN


DRT::ELEMENTS::Wall1LineType DRT::ELEMENTS::Wall1LineType::instance_;

DRT::ELEMENTS::Wall1LineType& DRT::ELEMENTS::Wall1LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
  *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Line::Wall1Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Wall1* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Line::Wall1Line(const DRT::ELEMENTS::Wall1Line& old) : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1Line::Clone() const
{
  DRT::ELEMENTS::Wall1Line* newelement = new DRT::ELEMENTS::Wall1Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          farah 02/14 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Wall1Line::Shape() const
{
  CORE::FE::CellType distype_line = CORE::FE::CellType::dis_none;

  switch (ParentMasterElement()->Shape())
  {
    case CORE::FE::CellType::tri3:
    {
      distype_line = CORE::FE::CellType::line2;
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      distype_line = CORE::FE::CellType::line3;
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      distype_line = CORE::FE::CellType::line2;
      break;
    }
    case CORE::FE::CellType::quad8:
    {
      distype_line = CORE::FE::CellType::line3;
      break;
    }
    case CORE::FE::CellType::quad9:
    {
      distype_line = CORE::FE::CellType::line3;
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      distype_line = CORE::FE::CellType::nurbs2;
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      distype_line = CORE::FE::CellType::nurbs3;
      break;
    }
    default:
      dserror("DRT::ELEMENTS::Wall1Line::Wall1Line: Unknown parent shape!");
  }

  return distype_line;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::Pack(CORE::COMM::PackBuffer& data) const
{
  dserror("this Wall1Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::Unpack(const std::vector<char>& data)
{
  dserror("this line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::Print(std::ostream& os) const
{
  os << "Wall1Line ";
  Element::Print(os);
  return;
}

BACI_NAMESPACE_CLOSE
