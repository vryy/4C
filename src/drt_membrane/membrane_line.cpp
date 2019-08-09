/*!----------------------------------------------------------------------
\brief

\level 3

\maintainer Sebastian Fuchs

\brief Nonlinear Membrane Finite Element line

*----------------------------------------------------------------------*/
#include "membrane.H"


/*----------------------------------------------------------------------*
 |  LINE 2 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_line2Type DRT::ELEMENTS::Membrane_line2Type::instance_;

DRT::ELEMENTS::Membrane_line2Type& DRT::ELEMENTS::Membrane_line2Type::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  LINE 3 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_line3Type DRT::ELEMENTS::Membrane_line3Type::instance_;

DRT::ELEMENTS::Membrane_line3Type& DRT::ELEMENTS::Membrane_line3Type::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  constructor (public)                                   fbraeu 06/16 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MembraneLine<distype>::MembraneLine(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Membrane<distype>* parent, const int lline)
    : DRT::FaceElement(id, owner), intpointsline_(DRT::UTILS::intrule_line_2point)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  switch (DRT::UTILS::DisTypeToFaceShapeType<distype>::shape)
  {
    case line2:
    {
      DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule_line_2point;
      // get gauss integration points
      intpointsline_ = DRT::UTILS::IntegrationPoints1D(gaussrule);
      break;
    }
    case line3:
    {
      DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule_line_3point;
      // get gauss integration points
      intpointsline_ = DRT::UTILS::IntegrationPoints1D(gaussrule);
      break;
    }
    break;
    default:
      dserror("shape type unknown!\n");
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                              fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MembraneLine<distype>::MembraneLine(const DRT::ELEMENTS::MembraneLine<distype>& old)
    : DRT::FaceElement(old), intpointsline_(old.intpointsline_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::MembraneLine<distype>::Clone() const
{
  DRT::ELEMENTS::MembraneLine<distype>* newelement =
      new DRT::ELEMENTS::MembraneLine<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::Element::DiscretizationType DRT::ELEMENTS::MembraneLine<distype>::Shape() const
{
  return DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             fb 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Pack(DRT::PackBuffer& data) const
{
  dserror("this membrane line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        fbraeu 06/165 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Unpack(const std::vector<char>& data)
{
  dserror("this membrane line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  destructor (public)                                     fbraeu 06/16|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MembraneLine<distype>::~MembraneLine()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             fbraeu 06/16|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Print(std::ostream& os) const
{
  os << "MembraneLine ";
  os << " Discretization type: "
     << DRT::DistypeToString(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape).c_str();
  Element::Print(os);
  return;
}

template class DRT::ELEMENTS::MembraneLine<DRT::Element::tri3>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::tri6>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::quad4>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::quad9>;
