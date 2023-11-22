/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element line

*----------------------------------------------------------------------*/
#include "baci_membrane.H"


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
template <CORE::FE::CellType distype>
DRT::ELEMENTS::MembraneLine<distype>::MembraneLine(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Membrane<distype>* parent, const int lline)
    : DRT::FaceElement(id, owner), intpointsline_(CORE::DRT::UTILS::GaussRule1D::line_2point)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  switch (CORE::DRT::UTILS::DisTypeToFaceShapeType<distype>::shape)
  {
    case CORE::FE::CellType::line2:
    {
      CORE::DRT::UTILS::GaussRule1D gaussrule = CORE::DRT::UTILS::GaussRule1D::line_2point;
      // get gauss integration points
      intpointsline_ = CORE::DRT::UTILS::IntegrationPoints1D(gaussrule);
      break;
    }
    case CORE::FE::CellType::line3:
    {
      CORE::DRT::UTILS::GaussRule1D gaussrule = CORE::DRT::UTILS::GaussRule1D::line_3point;
      // get gauss integration points
      intpointsline_ = CORE::DRT::UTILS::IntegrationPoints1D(gaussrule);
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
template <CORE::FE::CellType distype>
DRT::ELEMENTS::MembraneLine<distype>::MembraneLine(const DRT::ELEMENTS::MembraneLine<distype>& old)
    : DRT::FaceElement(old), intpointsline_(old.intpointsline_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
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
template <CORE::FE::CellType distype>
CORE::FE::CellType DRT::ELEMENTS::MembraneLine<distype>::Shape() const
{
  return CORE::DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             fb 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  dserror("this membrane line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        fbraeu 06/165 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Unpack(const std::vector<char>& data)
{
  dserror("this membrane line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             fbraeu 06/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneLine<distype>::Print(std::ostream& os) const
{
  os << "MembraneLine ";
  os << " Discretization type: "
     << CORE::FE::CellTypeToString(CORE::DRT::UTILS::DisTypeToFaceShapeType<distype>::shape)
            .c_str();
  Element::Print(os);
  return;
}

template class DRT::ELEMENTS::MembraneLine<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::MembraneLine<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::MembraneLine<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::MembraneLine<CORE::FE::CellType::quad9>;
