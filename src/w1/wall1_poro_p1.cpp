/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "wall1_poro_p1.H"
#include "wall1_poro_p1_eletypes.H"

#include "utils_factory.H"

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1<distype>::Wall1_PoroP1(int id, int owner)
    : DRT::ELEMENTS::Wall1_Poro<distype>(id, owner)
{
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1<distype>::Wall1_PoroP1(const DRT::ELEMENTS::Wall1_PoroP1<distype>& old)
    : DRT::ELEMENTS::Wall1_Poro<distype>(old)
{
}

template <DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::Wall1_PoroP1<distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Wall1_PoroP1<distype>(*this);
  return newelement;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  Base::AddtoPack(data, type);

  // add base class Element
  Base::Pack(data);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  Base::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  Base::ExtractfromPack(position, data, basedata);
  Base::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1_PoroP1<distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line, Wall1_PoroP1>(DRT::UTILS::buildLines, this);
}

template <DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1_PoroP1<distype>::Surfaces()
{
  std::vector<Teuchos::RCP<Element>> surfaces(1);
  surfaces[0] = Teuchos::rcp(this, false);
  return surfaces;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Print(std::ostream& os) const
{
  os << "Wall1_PoroP1 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << Base::data_;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case DRT::Element::tri3:
      return DRT::ELEMENTS::WallTri3PoroP1Type::Instance().UniqueParObjectId();
    case DRT::Element::quad4:
      return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance().UniqueParObjectId();
    case DRT::Element::quad9:
      return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance().UniqueParObjectId();
    default:
      dserror("unknown element type");
  }
  return -1;
}

template <DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::Wall1_PoroP1<distype>::ElementType() const
{
  {
    switch (distype)
    {
      case DRT::Element::tri3:
        return DRT::ELEMENTS::WallTri3PoroP1Type::Instance();
      case DRT::Element::quad4:
        return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
      case DRT::Element::quad9:
        return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance();
      default:
        dserror("unknown element type");
    }
    return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
  }
}

template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>;
