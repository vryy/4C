/*----------------------------------------------------------------------------*/
/*! \file

\brief 2D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_ale2.H"
#include "baci_lib_discret.H"
#include "baci_utils_exceptions.H"


DRT::ELEMENTS::Ale2LineType DRT::ELEMENTS::Ale2LineType::instance_;

DRT::ELEMENTS::Ale2LineType& DRT::ELEMENTS::Ale2LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale2Line::Ale2Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Ale2* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale2Line::Ale2Line(const DRT::ELEMENTS::Ale2Line& old) : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Ale2Line::Clone() const
{
  DRT::ELEMENTS::Ale2Line* newelement = new DRT::ELEMENTS::Ale2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Ale2Line::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return DiscretizationType::line2;
    case 3:
      return DiscretizationType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Pack(DRT::PackBuffer& data) const
{
  dserror("this Ale2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Unpack(const std::vector<char>& data)
{
  dserror("this Ale2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale2Line::~Ale2Line() { return; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Print(std::ostream& os) const
{
  os << "Ale2Line ";
  Element::Print(os);
  return;
}
