/*! \file
 *
\brief Line element associated to the shell 7-Parameter element

\level 3
*/

#include "baci_shell7p_line.H"

DRT::ELEMENTS::Shell7pLineType DRT::ELEMENTS::Shell7pLineType::instance_;

DRT::ELEMENTS::Shell7pLineType& DRT::ELEMENTS::Shell7pLineType::Instance() { return instance_; }

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell7pLineType::Create(const int id, const int owner)
{
  return Teuchos::null;
}

DRT::ELEMENTS::Shell7pLine::Shell7pLine(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::Element* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  // type of gaussian integration
  switch (Shape())
  {
    case DRT::Element::DiscretizationType::line2:
      gaussrule_ = CORE::DRT::UTILS::GaussRule1D::line_2point;
      break;
    case DRT::Element::DiscretizationType::line3:
      gaussrule_ = CORE::DRT::UTILS::GaussRule1D::line_3point;
      break;
    default:
      dserror("shape type unknown!\n");
  }
}

DRT::ELEMENTS::Shell7pLine::Shell7pLine(const DRT::ELEMENTS::Shell7pLine& old)
    : DRT::FaceElement(old), gaussrule_(old.gaussrule_)
{
}

DRT::Element* DRT::ELEMENTS::Shell7pLine::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Shell7pLine(*this);
  return newelement;
}

DRT::Element::DiscretizationType DRT::ELEMENTS::Shell7pLine::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return DRT::Element::DiscretizationType::line2;
    case 3:
      return DRT::Element::DiscretizationType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}

void DRT::ELEMENTS::Shell7pLine::Pack(DRT::PackBuffer& data) const
{
  dserror("this Shell7line element does not support communication");
}


void DRT::ELEMENTS::Shell7pLine::Unpack(const std::vector<char>& data)
{
  dserror("this Shell line element does not support communication");
}


void DRT::ELEMENTS::Shell7pLine::Print(std::ostream& os) const
{
  os << "Shell7pLine ";
  Element::Print(os);
}
